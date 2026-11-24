# =============================================================================
# LESSON 2  --  BLUP by hand (8 animals)
#
# BLUP = Best Linear Unbiased Prediction.
# Goal: build the mixed model equations yourself, step by step.
#
# Data: Mrode (2014) Example 3.1, pre-weaning gain (kg) of 5 beef calves.
#       sigma2_a = 20, sigma2_e = 40  ->  heritability = 20/60 = 0.33
#
# R you need for this lesson:
#   t(X)        transpose: rows become columns
#   X %*% Y     matrix multiplication
#   solve(L, R) solve the system of equations L x = R for x
#   solve(L)    the inverse of L
#   model.matrix(~ 0 + sex, data)  one 0/1 column per level of sex
# =============================================================================
source("R/functions.R")

ped <- prepare_pedigree(data.frame(
  id   = c("1", "2", "3", "4", "5", "6", "7", "8"),
  sire = c("0", "0", "0", "1", "3", "1", "4", "3"),
  dam  = c("0", "0", "0", "0", "2", "2", "5", "6")))

phen <- data.frame(id  = c("4", "5", "6", "7", "8"),
                   sex = c("M", "F", "F", "M", "M"),
                   wwg = c(4.5, 2.9, 3.9, 3.5, 5.0))
sigma2_a <- 20; sigma2_e <- 40

# -- The model ------------------------------------------------------------------
#   y  =  X b  +  Z a  +  e
#   gain = sex effect + breeding value + everything else
# b: "fixed effects" = groups we correct for (sex, herd, year ...)
# a: "random animal effects" = the breeding values we want to predict

X <- model.matrix(~ 0 + sex, data = phen)[, , drop = FALSE]  # [, ] drops clutter
print(X)

Z <- make_Z(phen$id, ped$id)
dimnames(Z) <- list(phen$id, ped$id)
print(Z)
# Animals 1-3 have no record (all-zero columns) but still get a breeding
# value, predicted through their offspring.

A    <- make_A(ped$s, ped$d)
Ainv <- make_Ainv(ped$s, ped$d, F = diag(A) - 1)

# -- Build the mixed model equations piece by piece --------------------------------
alpha <- sigma2_e / sigma2_a          # = 2. Bigger alpha = less trust in records
y <- phen$wwg
LHS <- rbind(cbind(t(X) %*% X, t(X) %*% Z),
             cbind(t(Z) %*% X, t(Z) %*% Z + Ainv * alpha))
RHS <- c(t(X) %*% y, t(Z) %*% y)
print(round(LHS, 2)); print(RHS)

sol <- solve(LHS, RHS)
names(sol) <- c(colnames(X), paste0("animal_", ped$id))
print(round(sol, 3))
# Textbook answer (Mrode Example 3.1), which we match within 0.001:
#   sexF 3.404, sexM 4.358
#   animals 1..8:  0.098 -0.019 -0.041 -0.009 -0.186  0.177 -0.249  0.183

# -- Reliability -----------------------------------------------------------------
# PEV (prediction error variance) = how uncertain an EBV still is.
# Reliability = 1 - PEV / ((1 + F) * sigma2_a): 0 = no information, 1 = certain.
fit <- solve_mme(y, X, Z, Ainv, sigma2_a, sigma2_e, Kdiag = diag(A))
print(data.frame(id = ped$id, has_record = ped$id %in% phen$id,
                 ebv = round(fit$ebv, 3), reliability = round(fit$rel, 3)))

# Things to notice:
#  * Calf 4 gained 4.5 kg and heifer 6 only 3.9 kg, yet calf 6 has the higher
#    EBV (0.177 vs -0.009). Males gain about 1 kg more on average; the sex
#    effect is removed before judging genetics.
#  * Founders 1, 2 and 3 have no record but still get an EBV, learned from
#    their offspring, with the lowest reliabilities.

# -- EXERCISES ------------------------------------------------------------------
# 1. Set sigma2_a = 40, sigma2_e = 20 (h2 = 0.67) and re-run. Do EBVs get
#    bigger or smaller? Why?
# 2. Calf 7's gain was mis-typed; it should be 5.5. Change it and re-run.
#    Which EBVs change? Why does calf 4 (its sire) go DOWN slightly?
#    Hint: look at the sexM estimate.
# 3. Delete the sex effect (X <- matrix(1, 5, 1), an overall mean only).
#    How does the ranking of calves 4 and 6 change?
