# Module 03 worked example: BLUP for a simple animal model (base R only)
#
# Model:  y = Xb + Za + e
#   y : weaning weight records
#   b : fixed effects (intercept, herd, sex)
#   a : random additive genetic effects (breeding values), a ~ N(0, A * sigma2_a)
#   e : residuals, e ~ N(0, I * sigma2_e)
#
# Run from the repository root:  source("modules/03-blup-animal-model.R")

ped   <- read.csv("data/toy_pedigree.csv")
pheno <- read.csv("data/toy_phenotypes.csv")

# Assumed variance components (in practice estimated by REML - see Module 04)
sigma2_a <- 90    # additive genetic variance (kg^2)
sigma2_e <- 210   # residual variance (kg^2)
h2       <- sigma2_a / (sigma2_a + sigma2_e)
alpha    <- sigma2_e / sigma2_a
cat(sprintf("Heritability = %.2f, alpha = %.3f\n\n", h2, alpha))

# ---- Step 1: numerator relationship matrix A (tabular method) ---------------
# Requires parents to be listed before offspring.
n <- nrow(ped)
A <- matrix(0, n, n, dimnames = list(ped$id, ped$id))
for (i in seq_len(n)) {
  s <- match(ped$sire[i], ped$id)   # NA if unknown
  d <- match(ped$dam[i],  ped$id)
  # Diagonal: 1 + F_i, where F_i = 0.5 * a(sire, dam)
  A[i, i] <- 1 + if (!is.na(s) && !is.na(d)) 0.5 * A[s, d] else 0
  if (i > 1) {
    for (j in seq_len(i - 1)) {
      a_js <- if (!is.na(s)) A[j, s] else 0
      a_jd <- if (!is.na(d)) A[j, d] else 0
      A[i, j] <- A[j, i] <- 0.5 * (a_js + a_jd)
    }
  }
}
cat("Relationship matrix A:\n")
print(round(A, 4))
cat("\nInbreeding coefficients (F = diag(A) - 1):\n")
print(round(diag(A) - 1, 4))

# ---- Step 2: design matrices ------------------------------------------------
y <- pheno$weaning_weight_kg
X <- model.matrix(~ herd + sex, data = pheno)
Z <- matrix(0, nrow(pheno), n, dimnames = list(NULL, ped$id))
Z[cbind(seq_len(nrow(pheno)), match(pheno$id, ped$id))] <- 1

# ---- Step 3: Henderson's mixed model equations ------------------------------
#  [ X'X   X'Z             ] [b]   [X'y]
#  [ Z'X   Z'Z + A^-1*alpha ] [a] = [Z'y]
Ainv <- solve(A)
LHS <- rbind(
  cbind(crossprod(X),    crossprod(X, Z)),
  cbind(crossprod(Z, X), crossprod(Z) + Ainv * alpha)
)
RHS <- c(crossprod(X, y), crossprod(Z, y))
sol <- solve(LHS, RHS)

p    <- ncol(X)
bhat <- sol[1:p]
ebv  <- sol[-(1:p)]
names(bhat) <- colnames(X)

cat("\nFixed effect solutions (BLUE):\n")
print(round(bhat, 3))

# ---- Step 4: accuracy of EBVs -----------------------------------------------
# Prediction error variance PEV = diag(C22) * sigma2_e, where C22 is the
# random-effect block of the inverted LHS. Accuracy r = sqrt(1 - PEV / var(a_i)),
# where var(a_i) = (1 + F_i) * sigma2_a, so inbred animals are handled correctly.
C   <- solve(LHS)
PEV <- diag(C)[-(1:p)] * sigma2_e
acc <- sqrt(pmax(0, 1 - PEV / (diag(A) * sigma2_a)))

result <- data.frame(
  id       = ped$id,
  has_record = ped$id %in% pheno$id,
  F        = round(diag(A) - 1, 4),
  EBV_kg   = round(ebv, 2),
  accuracy = round(acc, 3)
)
cat("\nEstimated breeding values:\n")
print(result[order(-result$EBV_kg), ], row.names = FALSE)

# Exercises:
# 1. Check A[6, 6] by hand. Why is animal 6 inbred?
# 2. Change sigma2_a to 30 (lower h2). What happens to the spread of EBVs and to accuracy?
# 3. Animals 1-3 have no records. Explain where their EBVs come from.
# 4. Remove sex from the model. Which EBVs change most, and why?
