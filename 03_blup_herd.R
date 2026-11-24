# =============================================================================
# LESSON 3  --  Pedigree BLUP on a herd of 1,010 animals
#
# Goal: run a complete genetic evaluation the way you would at work:
#       read files -> check data -> estimate variances -> solve -> check EBVs.
# Trait: weaning weight (kg).
# =============================================================================
source("R/functions.R")

# -- 1. Read and CHECK the data (never skip this) ------------------------------------
ped_raw <- read.csv("data/pedigree.csv", colClasses = c(id = "character", sire = "character", dam = "character"))
phen    <- read.csv("data/phenotypes.csv", colClasses = c(id = "character"))

cat("Pedigree rows:", nrow(ped_raw), "  Records:", nrow(phen), "\n")
print(summary(phen$weaning_wt))
print(table(phen$cg, phen$sex))            # every contemporary group has data?
hist(phen$weaning_wt, main = "Weaning weight", xlab = "kg")   # odd values?
stopifnot(all(phen$id %in% ped_raw$id))    # every recorded animal is in the pedigree

# A CONTEMPORARY GROUP (cg) = animals raised together under the same conditions
# (here: same herd, same generation). Animals are only compared within a group.

# -- 2. Pedigree -----------------------------------------------------------------------
ped  <- prepare_pedigree(ped_raw)          # keeps sex and generation columns
A    <- make_A(ped$s, ped$d)
F    <- diag(A) - 1
cat(sprintf("Average inbreeding: %.4f  (max %.4f)\n", mean(F), max(F)))
Ainv <- make_Ainv(ped$s, ped$d, F)

# -- 3. Model and variance components ------------------------------------------------
# weaning_wt = contemporary group + sex + animal + residual
X <- model.matrix(~ 0 + cg + sex, data = phen)[, , drop = FALSE]
Z <- make_Z(phen$id, ped$id)
y <- phen$weaning_wt

# sigma2_a and sigma2_e are estimated by REML before any BLUP run.
vc <- reml_animal(y, X, Z, A)
print(round(vc, 3))
# The simulation used sigma2_a = 90, sigma2_e = 210 (h2 = 0.30). With only 700
# records the estimate is rough. We use the true values below so that every
# lesson gives the same numbers; with real data you would use the estimates.
sigma2_a <- 90; sigma2_e <- 210

# -- 4. Solve --------------------------------------------------------------------------
fit <- solve_mme(y, X, Z, Ainv, sigma2_a, sigma2_e, Kdiag = 1 + F)
print(round(fit$fixed, 2))
# cg... = expected weight of a female calf of BASE genetic merit (the average
#         of the generation-0 founders) in that contemporary group
# sexM  = how much heavier males are than females

# -- 5. Check the results ----------------------------------------------------------------
ebv <- data.frame(id = ped$id, sire = ped$sire, dam = ped$dam, sex = ped$sex,
                  generation = ped$generation,
                  has_record = ped$id %in% phen$id,
                  n_offspring = tabulate(c(ped$s, ped$d), nbins = nrow(ped)),
                  ebv = round(fit$ebv, 2), reliability = round(fit$rel, 3))
dir.create("output", showWarnings = FALSE)
write.csv(ebv, "output/ebv_blup.csv", row.names = FALSE)

cat("\nTop 10 animals by EBV:\n")
print(head(ebv[order(-ebv$ebv), ], 10), row.names = FALSE)

cat("\nMean EBV and reliability by generation:\n")
print(aggregate(cbind(ebv, reliability) ~ generation, data = ebv, FUN = mean))
boxplot(ebv ~ generation, data = ebv, xlab = "Generation", ylab = "EBV (kg)")
# Generation-3 calves have no records: their EBV is just the average of their
# parents' EBVs ("parent average"), so full sibs get IDENTICAL EBVs.

# -- EXERCISES -------------------------------------------------------------------------------
# 1. Pick two full sibs in generation 3 (same sire AND dam). Compare their
#    EBVs. Why are they identical? What would you need to tell them apart?
# 2. Is there genetic progress? How many kg per generation?
#    Hint: coef(lm(ebv ~ generation, data = ebv))
# 3. Is a high reliability the same as a high EBV? Find a bull (sex "M")
#    with high reliability but an average EBV. What gives him the high
#    reliability? (Look at n_offspring.)
