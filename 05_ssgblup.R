# =============================================================================
# LESSON 5  --  Single-step GBLUP (ssGBLUP): pedigree + DNA together
#
# Goal: combine A (everyone) and G (genotyped animals) into one matrix H, so
#       ALL records and ALL animals are used in one evaluation. This is how
#       most national evaluations run today.
# =============================================================================
source("R/functions.R")
sigma2_a <- 90; sigma2_e <- 210

ped  <- prepare_pedigree(read.csv("data/pedigree.csv", colClasses = c(id = "character", sire = "character", dam = "character")))
phen <- read.csv("data/phenotypes.csv", colClasses = c(id = "character"))
M    <- qc_genotypes(read_genotypes("data/genotypes.txt"))

# -- 1. Pedigree part --------------------------------------------------------------------
A    <- make_A(ped$s, ped$d); dimnames(A) <- list(ped$id, ped$id)
Ainv <- make_Ainv(ped$s, ped$d, diag(A) - 1)

# -- 2. Genomic part ------------------------------------------------------------------------
# A22 = the part of A for genotyped animals only (rows and columns 'geno').
geno <- match(rownames(M), ped$id)
stopifnot(!anyNA(geno))
A22  <- A[geno, geno]
G    <- tune_G(make_G(M), A22, blend = 0.05)   # put G on the same scale as A22

# -- 3. H-inverse ---------------------------------------------------------------------------
#   H-inverse = A-inverse, plus (G-inverse - A22-inverse) in the rows and
#               columns of the genotyped animals.
# Start from the pedigree; for genotyped animals ADD what the DNA says and
# SUBTRACT what the pedigree already said, so nothing is counted twice.
Hinv <- Ainv
Hinv[geno, geno] <- Hinv[geno, geno] + solve(G) - solve(A22)

# -- 4. Solve: same code as Lesson 3, only Ainv -> Hinv ------------------------------------
X <- model.matrix(~ 0 + cg + sex, data = phen)[, , drop = FALSE]
Z <- make_Z(phen$id, ped$id)
fit <- solve_mme(phen$weaning_wt, X, Z, Hinv, sigma2_a, sigma2_e,
                 Kdiag = diag(solve(Hinv)))      # diag(H) for reliabilities
print(round(fit$fixed, 2))

out <- data.frame(id = ped$id, generation = ped$generation,
                  genotyped = seq_along(ped$id) %in% geno,
                  has_record = ped$id %in% phen$id,
                  gebv = round(fit$ebv, 2), reliability = round(fit$rel, 3))
dir.create("output", showWarnings = FALSE)
write.csv(out, "output/gebv_ssgblup.csv", row.names = FALSE)

# Reliability compared with pedigree BLUP (run Lesson 3 first)
if (file.exists("output/ebv_blup.csv")) {
  blup  <- read.csv("output/ebv_blup.csv", colClasses = c(id = "character"))
  group <- ifelse(out$generation == 3, "generation-3 calves (genotyped, no record)",
           ifelse(!out$genotyped, "not genotyped",
           ifelse(out$has_record, "genotyped, with record", "genotyped founder bulls")))
  cmp <- data.frame(BLUP    = tapply(blup$reliability[match(out$id, blup$id)], group, mean),
                    ssGBLUP = tapply(out$reliability, group, mean))
  cat("\nAverage reliability by group:\n"); print(round(cmp, 3))
}

# -- EXERCISES -------------------------------------------------------------------------
# 1. Founder cows (generation 0, sex F) are not genotyped. Compare their EBV
#    from Lesson 3 with their GEBV here. Do they change? Compare the cg
#    estimates of the two runs as well.
# 2. Try blend = 0.20 in tune_G(), then re-run Lesson 6. What changes?
#    (Set it back to 0.05 and re-run afterwards.)
# 3. GBLUP (Lesson 4) used far fewer records than ssGBLUP. Whose records
#    did GBLUP throw away, and why does that matter for the young calves?
