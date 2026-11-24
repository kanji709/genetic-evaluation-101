# =============================================================================
# LESSON 1  --  Pedigrees and the relationship matrix A
#
# Goal: turn a pedigree into numbers that say how closely two animals are
#       related.
#
# Example pedigree from Mrode (2014) "Linear Models for the Prediction of
# Animal Breeding Values", 3rd ed., Example 3.1 -- so you can check the book.
# =============================================================================
source("R/functions.R")

# -- 1. A tiny pedigree -------------------------------------------------------
# Animals 1, 2, 3 are founders (parents unknown = "0").
ped <- data.frame(
  id   = c("1", "2", "3", "4", "5", "6", "7", "8"),
  sire = c("0", "0", "0", "1", "3", "1", "4", "3"),
  dam  = c("0", "0", "0", "0", "2", "2", "5", "6")
)
ped <- prepare_pedigree(ped)
print(ped)
# Columns s and d are the ROW NUMBERS of the parents. Software works with
# row numbers, not names; this step is called "renumbering".

# -- 2. Build A -----------------------------------------------------------------
# A[i, j] = the expected proportion of genes that i and j share because they
# inherited them from the same ancestors ("identical by descent").
A <- make_A(ped$s, ped$d)
dimnames(A) <- list(ped$id, ped$id)
print(A)

# How to read it:
#   A["1","4"] = 0.5     1 is the sire of 4 (parent-offspring)
#   A["5","6"] = 0.25    half sibs: same dam (2), different sires
#   A["4","5"] = 0       no known common ancestor
#   diag = 1 + F         F = inbreeding coefficient. Nobody here is inbred
#                        (all diagonals are 1) because no animal's parents are
#                        related to each other. Exercise 2 creates one.

F <- diag(A) - 1
cat("\nInbreeding coefficients:\n"); print(F)

# -- 3. A-inverse, two ways ------------------------------------------------------
# Henderson's rules (with Quaas' extension for inbreeding) write A-inverse
# straight from the pedigree. Here F comes from A; big programs get F
# another way and never build A.
Ainv_rules <- make_Ainv(ped$s, ped$d, F)
Ainv_solve <- solve(A)
cat("\nMax difference between the two inverses:",
    signif(max(abs(Ainv_rules - Ainv_solve)), 3), "\n")
print(round(Ainv_rules, 3))
# A-inverse is mostly ZEROS ("sparse"): each animal only links to its own
# parents and mates. That is why software can handle millions of animals.
cat(sprintf("\nShare of zeros: A = %.0f%%, A-inverse = %.0f%%\n",
            100 * mean(A == 0), 100 * mean(abs(Ainv_rules) < 1e-12)))

# -- EXERCISES -----------------------------------------------------------------
# 1. Which pairs have the largest relationship (ignore the diagonal)?
#    What kind of relatives are they all?
#    Hint: A2 <- A; diag(A2) <- 0; which(A2 == max(A2), arr.ind = TRUE)
# 2. Add a 9th animal whose sire is 7 and dam is 8. Rebuild A. What is
#    animal 9's inbreeding coefficient, and why is it above zero?
# 3. Change animal 5's sire from "3" to "0" (unknown). Which relationships
#    in A change? What does that tell you about missing pedigree data?
