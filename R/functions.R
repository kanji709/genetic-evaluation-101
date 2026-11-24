# =============================================================================
# functions.R  --  helper functions used by every module of the course
#
# Plain base R (no packages), so you can read exactly what each step does.
# Real evaluation software (BLUPF90, MiX99, ASReml, ...) does the same maths
# with sparse matrices and iterative solvers, so it can handle millions of
# animals. The logic is identical.
#
# The same functions are written out inside tutorial.qmd, so that document
# renders on its own. If you change a function here, change it there too.
# =============================================================================


## ---- prepare_pedigree ----
# Input : data.frame with columns id, sire, dam (unknown parent = "0" or NA),
#         plus any other columns (sex, generation, ...) which are kept.
# Output: the pedigree sorted so parents come before their offspring, with
#         integer codes n (row number), s (sire row), d (dam row); 0 = unknown.
prepare_pedigree <- function(ped) {
  ped <- as.data.frame(ped, stringsAsFactors = FALSE)
  for (v in c("id", "sire", "dam")) ped[[v]] <- as.character(ped[[v]])
  ped$sire[is.na(ped$sire) | ped$sire == ""] <- "0"
  ped$dam[is.na(ped$dam)   | ped$dam  == ""] <- "0"
  if (any(duplicated(ped$id))) stop("Duplicated IDs in pedigree: ",
                                    paste(head(ped$id[duplicated(ped$id)]), collapse = ", "))

  # Parents listed as sire/dam but without a row of their own become founders
  missing_parents <- setdiff(unique(c(ped$sire, ped$dam)), c("0", ped$id))
  if (length(missing_parents) > 0) {
    add <- ped[rep(NA_integer_, length(missing_parents)), , drop = FALSE]
    add$id <- missing_parents; add$sire <- "0"; add$dam <- "0"
    ped <- rbind(add, ped)
  }

  # Generation number: founders = 0, otherwise 1 + the younger parent's number
  n  <- nrow(ped)
  is <- match(ped$sire, ped$id)
  id <- match(ped$dam,  ped$id)
  gen <- integer(n)
  for (iter in seq_len(n + 1)) {
    new_gen <- pmax(ifelse(is.na(is), 0L, gen[is] + 1L),
                    ifelse(is.na(id), 0L, gen[id] + 1L))
    if (all(new_gen == gen)) break
    if (iter > n) stop("Pedigree loop detected: an animal is its own ancestor.")
    gen <- new_gen
  }

  ped <- ped[order(gen, seq_len(n)), , drop = FALSE]
  ped$n <- seq_len(nrow(ped))
  ped$s <- match(ped$sire, ped$id, nomatch = 0L)
  ped$d <- match(ped$dam,  ped$id, nomatch = 0L)
  rownames(ped) <- NULL
  ped
}

## ---- make_A ----
# Numerator relationship matrix A (the "tabular method").
# A[i,j] = additive genetic relationship: the expected proportion of genes
#          that i and j share identical by descent (inherited from the same
#          ancestor). Parent-offspring 0.5, full sibs 0.5, half sibs 0.25.
#          Formally, twice the coancestry (kinship) coefficient.
# A[i,i] = 1 + F_i, where F_i is the inbreeding coefficient of i.
make_A <- function(s, d) {
  n <- length(s)
  A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    if (i > 1) {
      j <- seq_len(i - 1)
      a <- numeric(i - 1)
      if (s[i] > 0) a <- a + 0.5 * A[j, s[i]]
      if (d[i] > 0) a <- a + 0.5 * A[j, d[i]]
      A[j, i] <- a
      A[i, j] <- a
    }
    A[i, i] <- 1 + if (s[i] > 0 && d[i] > 0) 0.5 * A[s[i], d[i]] else 0
  }
  A
}

## ---- make_Ainv ----
# Inverse of A written directly from the pedigree: Henderson (1976) rules,
# extended for inbreeding by Quaas (1976). F = inbreeding coefficients.
# Here we take F from diag(A) - 1 for simplicity; large programs compute F
# without building A (Meuwissen & Luo 1992).
make_Ainv <- function(s, d, F = rep(0, length(s))) {
  n <- length(s)
  Ainv <- matrix(0, n, n)
  for (i in seq_len(n)) {
    p <- c(s[i], d[i])
    p <- p[p > 0]                         # known parents
    # variance of the Mendelian sampling term (the "luck of the draw")
    v <- if (length(p) == 2) 0.5 - 0.25 * (F[p[1]] + F[p[2]])
         else if (length(p) == 1) 0.75 - 0.25 * F[p]
         else 1
    k <- c(i, p)
    w <- c(1, rep(-0.5, length(p)))
    Ainv[k, k] <- Ainv[k, k] + (1 / v) * outer(w, w)
  }
  Ainv
}

## ---- read_genotypes ----
# BLUPF90-style genotype file: one animal per line, "<id> <0/1/2 string>".
# 0 and 2 = homozygous (two identical copies), 1 = heterozygous, 5 = missing.
read_genotypes <- function(file) {
  lines <- readLines(file)
  ids   <- sub("^\\s*(\\S+)\\s+.*$", "\\1", lines)
  calls <- sub("^\\s*\\S+\\s+(\\S+)\\s*$", "\\1", lines)
  M <- do.call(rbind, lapply(strsplit(calls, ""), as.integer))
  M[M == 5] <- NA
  dimnames(M) <- list(ids, paste0("snp", seq_len(ncol(M))))
  M
}

## ---- qc_genotypes ----
# Quality control: drop SNPs with low minor-allele frequency (MAF) or too
# many missing calls, then fill remaining missing calls with the SNP mean.
qc_genotypes <- function(M, min_maf = 0.05, max_missing = 0.10, verbose = TRUE) {
  p    <- colMeans(M, na.rm = TRUE) / 2
  maf  <- pmin(p, 1 - p)
  miss <- colMeans(is.na(M))
  keep <- maf >= min_maf & miss <= max_missing
  if (verbose) cat(sprintf("QC: kept %d of %d SNPs (%d low MAF, %d high missing)\n",
              sum(keep), ncol(M), sum(maf < min_maf), sum(miss > max_missing)))
  M <- M[, keep, drop = FALSE]
  for (j in which(colSums(is.na(M)) > 0)) M[is.na(M[, j]), j] <- mean(M[, j], na.rm = TRUE)
  M
}

## ---- make_G ----
# Genomic relationship matrix G (VanRaden 2008, method 1). Same meaning as A,
# but measured from DNA. p = allele frequencies used to centre the genotypes;
# by default the frequencies in the genotyped animals themselves, which makes
# G relative to the CURRENT population (many off-diagonals are negative).
make_G <- function(M, p = colMeans(M) / 2) {
  W <- sweep(M, 2, 2 * p)                 # centred genotypes: 0/1/2 minus 2p
  G <- tcrossprod(W) / (2 * sum(p * (1 - p)))
  dimnames(G) <- list(rownames(M), rownames(M))
  G
}

## ---- tune_G ----
# Put G on the same scale as A22 (pedigree relationships among genotyped
# animals): G* = a + b*G with the same average diagonal and the same average
# off-diagonal as A22 (Vitezica et al. 2011; BLUPF90 "TunedG 2"). Then blend
# in a little A22 so G* can be inverted (BLUPF90 default: 95% G, 5% A22).
tune_G <- function(G, A22, blend = 0.05, verbose = TRUE) {
  m  <- rbind(c(1, mean(diag(G))), c(1, mean(G)))
  ab <- solve(m, c(mean(diag(A22)), mean(A22)))
  Gt <- ab[1] + ab[2] * G
  if (verbose) cat(sprintf("Tuning G: G* = %.3f + %.3f * G, then %.0f%% A22 blended in\n",
                           ab[1], ab[2], 100 * blend))
  (1 - blend) * Gt + blend * A22
}

## ---- solve_mme ----
# Mixed model equations (Henderson):
#   y = Xb + Zu + e,  var(u) = K * sigma2_a,  var(e) = I * sigma2_e
#   [ X'X   X'Z              ] [b]   [X'y]
#   [ Z'X   Z'Z + Kinv*alpha ] [u] = [Z'y]      alpha = sigma2_e / sigma2_a
# K is A (BLUP), G (GBLUP) or H (ssGBLUP); only Kinv changes between methods.
# Reliability = 1 - PEV / (K_ii * sigma2_a), PEV = prediction error variance.
# Kdiag = diag(K): 1 + F for A, diag(G) for GBLUP, diag(H) for ssGBLUP.
solve_mme <- function(y, X, Z, Kinv, sigma2_a, sigma2_e, Kdiag = 1, reliability = TRUE) {
  alpha <- sigma2_e / sigma2_a
  XtX <- crossprod(X); XtZ <- crossprod(X, Z); ZtZ <- crossprod(Z)
  LHS <- rbind(cbind(XtX,    XtZ),
               cbind(t(XtZ), ZtZ + Kinv * alpha))
  RHS <- c(crossprod(X, y), crossprod(Z, y))
  nb  <- ncol(X)
  iu  <- nb + seq_len(ncol(Z))
  rel <- NULL
  if (reliability) {
    C   <- solve(LHS)                     # inverse of LHS (fine for small data)
    sol <- drop(C %*% RHS)
    pev <- diag(C)[iu] * sigma2_e
    rel <- 1 - pev / (Kdiag * sigma2_a)
    if (any(rel < -1e-6)) warning("Negative reliabilities: check Kdiag, variances and the model.")
  } else {
    sol <- solve(LHS, RHS)
  }
  list(fixed = setNames(sol[seq_len(nb)], colnames(X)),
       ebv = sol[iu], rel = rel)
}

## ---- make_Z ----
# Z: one row per record, one column per animal, a 1 where they match.
make_Z <- function(record_ids, animal_ids) {
  col <- match(record_ids, animal_ids)
  if (anyNA(col)) stop("Records for animals not in the relationship matrix: ",
                       paste(head(record_ids[is.na(col)]), collapse = ", "))
  Z <- matrix(0, length(record_ids), length(animal_ids))
  Z[cbind(seq_along(col), col)] <- 1
  Z
}

## ---- reml_animal ----
# REML estimates of sigma2_a and sigma2_e for the animal model
# y = Xb + Za + e, var(a) = K sigma2_a. Maximises the restricted likelihood
# numerically. Fine for a few thousand records; real data use AIREMLF90,
# ASReml, WOMBAT, etc.
reml_animal <- function(y, X, Z, K) {
  ev <- eigen(Z %*% K %*% t(Z), symmetric = TRUE)
  d  <- pmax(ev$values, 0)
  Uy <- drop(crossprod(ev$vectors, y))
  UX <- crossprod(ev$vectors, X)
  neg_loglik <- function(log_var) {
    s2a <- exp(log_var[1]); s2e <- exp(log_var[2])
    w   <- 1 / (s2a * d + s2e)                       # V = U diag(1/w) U'
    XVX <- crossprod(UX, UX * w)
    XVy <- crossprod(UX, Uy * w)
    yPy <- sum(Uy^2 * w) - sum(XVy * solve(XVX, XVy))
    0.5 * (-sum(log(w)) + determinant(XVX)$modulus[1] + yPy)
  }
  start <- log(rep(var(y) / 2, 2))
  fit <- optim(start, neg_loglik, method = "BFGS")
  v <- exp(fit$par)
  c(sigma2_a = v[1], sigma2_e = v[2], h2 = v[1] / sum(v),
    converged = fit$convergence == 0)
}

## ---- is_half_sib ----
# TRUE for pairs that share exactly one known parent (half sibs).
is_half_sib <- function(sire, dam) {
  ss <- outer(sire, sire, "==") & sire != "0"
  dd <- outer(dam,  dam,  "==") & dam  != "0"
  (ss & !dd) | (dd & !ss)
}

## ---- evaluate_prediction ----
# Validation statistics for one method:
#  accuracy   = correlation between prediction and true BV (young animals)
#  slope      = regression of true BV on prediction (1 = right spread)
#  gain_error = estimated minus true genetic gain from the base generation
#               to the young generation (0 = no level bias), in kg
evaluate_prediction <- function(pred, true_bv, generation,
                                young = max(generation), base = 0) {
  y <- generation == young & !is.na(pred)
  b <- generation == base  & !is.na(pred)
  c(accuracy   = cor(pred[y], true_bv[y]),
    slope      = unname(coef(lm(true_bv[y] ~ pred[y]))[2]),
    gain_error = (mean(pred[y]) - mean(pred[b])) - (mean(true_bv[y]) - mean(true_bv[b])))
}

## ---- simulate_herd ----
# Simulated beef herd recorded for weaning weight (kg). Returns a list with
# pedigree, phenotypes, genotypes (matrix), true breeding values, and the
# founders' allele frequencies (unknowable in real life; used in Section 7).
#  * 5 chromosomes x 600 SNPs, 100 of them affect the trait (QTL)
#  * h2 = 0.30 (sigma2_a = 90, sigma2_e = 210)
#  * generation 0: 10 bulls + 100 cows; generations 1-3: 300 calves each
#  * each generation: the 10 heaviest bull calves become sires, 100 random heifers dams
#  * contemporary groups = herd x generation, with random environmental effects
#  * generation 3 calves have no weights: they are the selection candidates
#  * genotyped: all sires used, all generation 3 calves, 40% of generations 1-2
simulate_herd <- function(seed = 2026, n_chr = 5, snp_per_chr = 600, n_qtl = 100,
                          sigma2_a = 90, sigma2_e = 210, sd_cg = 8) {
  set.seed(seed)
  n_snp <- n_chr * snp_per_chr
  chr <- rep(seq_len(n_chr), each = snp_per_chr)
  pos <- unlist(lapply(seq_len(n_chr), function(i) sort(runif(snp_per_chr))))  # Morgans
  make_gamete <- function(h1, h2) {          # one recombined copy of a parent's genome
    out <- integer(n_snp)
    for (c in seq_len(n_chr)) {
      k  <- which(chr == c)
      xo <- sort(runif(rpois(1, 1)))         # about 1 crossover per chromosome
      strand <- (sample(0:1, 1) + findInterval(pos[k], xo)) %% 2
      out[k] <- ifelse(strand == 0, h1[k], h2[k])
    }
    out
  }

  mu <- 240
  herd_eff <- c(H1 = 0, H2 = 12, H3 = -8, H4 = 5)
  sex_eff  <- c(F = 0, M = 15)
  n_sires <- 10; n_dams <- 100; calves_per_dam <- 3; n_gen <- 3
  n_per_gen <- n_dams * calves_per_dam
  nf <- n_sires + n_dams
  n_total <- nf + n_gen * n_per_gen

  id   <- sprintf("A%04d", seq_len(n_total))
  sire <- dam <- rep("0", n_total)
  sex  <- c(rep("M", n_sires), rep("F", n_dams), rep(NA, n_gen * n_per_gen))
  generation <- c(rep(0L, nf), rep(NA_integer_, n_gen * n_per_gen))

  H1 <- H2 <- matrix(0L, n_total, n_snp)
  p0 <- runif(n_snp, 0.05, 0.95)              # founder allele frequencies
  H1[1:nf, ] <- matrix(rbinom(nf * n_snp, 1, rep(p0, each = nf)), nf, n_snp)
  H2[1:nf, ] <- matrix(rbinom(nf * n_snp, 1, rep(p0, each = nf)), nf, n_snp)

  qtl <- sort(sample(n_snp, n_qtl))
  eff <- rnorm(n_qtl)
  eff <- eff * sqrt(sigma2_a / sum(2 * p0[qtl] * (1 - p0[qtl]) * eff^2))
  tbv_of <- function(rows) drop((H1[rows, qtl, drop = FALSE] + H2[rows, qtl, drop = FALSE] -
                                 rep(2 * p0[qtl], each = length(rows))) %*% eff)

  herd <- sample(names(herd_eff), n_total, replace = TRUE)
  cg_eff <- herd_eff + matrix(rnorm(4 * n_gen, 0, sd_cg), 4, n_gen,
                              dimnames = list(names(herd_eff), 0:(n_gen - 1)))  # herd x generation
  tbv <- numeric(n_total)
  y   <- rep(NA_real_, n_total)
  tbv[1:nf] <- tbv_of(1:nf)
  phen_of <- function(rows) mu + cg_eff[cbind(herd[rows], as.character(generation[rows]))] +
    sex_eff[sex[rows]] + tbv[rows] + rnorm(length(rows), 0, sqrt(sigma2_e))
  founder_cows <- (n_sires + 1):nf
  y[founder_cows] <- phen_of(founder_cows)

  sires <- 1:n_sires; dams <- founder_cows; last <- nf
  for (g in 1:n_gen) {
    rows <- last + seq_len(n_per_gen)
    dam_of  <- rep(dams, each = calves_per_dam)
    sire_of <- sample(sires, n_per_gen, replace = TRUE)
    for (k in seq_along(rows)) {
      H1[rows[k], ] <- make_gamete(H1[sire_of[k], ], H2[sire_of[k], ])
      H2[rows[k], ] <- make_gamete(H1[dam_of[k], ],  H2[dam_of[k], ])
    }
    sire[rows] <- id[sire_of]; dam[rows] <- id[dam_of]
    sex[rows]  <- sample(c("M", "F"), n_per_gen, replace = TRUE)
    generation[rows] <- g
    tbv[rows] <- tbv_of(rows)
    if (g < n_gen) y[rows] <- phen_of(rows)
    males <- rows[sex[rows] == "M"]
    sires <- males[order(-y[males])][1:n_sires]
    dams  <- sample(rows[sex[rows] == "F"], n_dams)
    last  <- max(rows)
  }

  used_sires <- unique(match(sire[sire != "0"], id))
  young <- which(generation == n_gen)
  mid   <- which(generation %in% 1:(n_gen - 1))
  genotyped <- sort(unique(c(used_sires, young, mid[runif(length(mid)) < 0.40])))
  M <- H1[genotyped, ] + H2[genotyped, ]
  dimnames(M) <- list(id[genotyped], paste0("snp", seq_len(n_snp)))

  has_y <- which(!is.na(y))
  list(pedigree   = data.frame(id, sire, dam, sex, generation, stringsAsFactors = FALSE),
       phenotypes = data.frame(id = id[has_y], herd = herd[has_y],
                               generation = generation[has_y], sex = sex[has_y],
                               cg = paste0(herd[has_y], "_G", generation[has_y]),
                               weaning_wt = round(y[has_y], 1), stringsAsFactors = FALSE),
       genotypes  = M,
       true_bv    = data.frame(id, generation, true_bv = round(tbv, 3), stringsAsFactors = FALSE),
       founder_freq = setNames(p0, colnames(M)))
}

## ---- evaluate_herd ----
# Runs pedigree BLUP, GBLUP, ssGBLUP, and ssGBLUP with G built from the
# founders' allele frequencies, on one simulated herd. Returns the
# validation statistics for the youngest generation (used in Section 7).
evaluate_herd <- function(seed, sigma2_a = 90, sigma2_e = 210) {
  h <- simulate_herd(seed)
  ped  <- prepare_pedigree(h$pedigree)
  phen <- h$phenotypes
  tbv  <- h$true_bv$true_bv[match(ped$id, h$true_bv$id)]
  A    <- make_A(ped$s, ped$d)
  Ainv <- make_Ainv(ped$s, ped$d, diag(A) - 1)
  X <- model.matrix(~ 0 + cg + sex, data = phen)
  Z <- make_Z(phen$id, ped$id)
  y <- phen$weaning_wt
  out <- list()
  out$BLUP <- solve_mme(y, X, Z, Ainv, sigma2_a, sigma2_e, reliability = FALSE)$ebv

  M <- qc_genotypes(h$genotypes, verbose = FALSE)
  geno <- match(rownames(M), ped$id)
  G <- make_G(M)
  pg <- phen[phen$id %in% rownames(M), ]
  fg <- solve_mme(pg$weaning_wt, model.matrix(~ 0 + cg + sex, data = pg),
                  make_Z(pg$id, rownames(M)), solve(G + diag(0.01, nrow(G))),
                  sigma2_a, sigma2_e, reliability = FALSE)
  out$GBLUP <- rep(NA_real_, nrow(ped)); out$GBLUP[geno] <- fg$ebv

  A22 <- A[geno, geno]; iA22 <- solve(A22)
  ss <- function(Gx) {
    Hinv <- Ainv; Hinv[geno, geno] <- Hinv[geno, geno] + solve(Gx) - iA22
    solve_mme(y, X, Z, Hinv, sigma2_a, sigma2_e, reliability = FALSE)$ebv
  }
  out$ssGBLUP <- ss(tune_G(G, A22, verbose = FALSE))
  G0 <- make_G(M, p = h$founder_freq[colnames(M)])
  out$`ssGBLUP, founder frequencies` <- ss(0.95 * G0 + 0.05 * A22)

  res <- t(sapply(out, evaluate_prediction, true_bv = tbv, generation = ped$generation))
  data.frame(seed = seed, method = rownames(res), res, row.names = NULL)
}
