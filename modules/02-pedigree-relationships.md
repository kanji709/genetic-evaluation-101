# Module 02 — Pedigrees and relationships

## Learning objectives

- Clean and validate a pedigree (ordering, missing parents, sex and date checks)
- Compute additive relationships and inbreeding coefficients
- Build the numerator relationship matrix **A** with the tabular method
- Understand why A⁻¹ (not A) is used in practice and how Henderson's rules build it directly

## Core concepts

**Pedigree checks.** Parents must appear before offspring; an animal cannot be its own parent; sires should be male and dams female; parents must be born before offspring.

**Additive relationship** a(i, j): twice the coefficient of coancestry — the expected proportion of alleles shared identical by descent. Parent–offspring and full sibs ≈ 0.5, half sibs ≈ 0.25.

**Tabular method** (animals ordered oldest first):

- Diagonal: a(i, i) = 1 + F(i), where F(i) = ½ a(sire, dam)
- Off-diagonal: a(i, j) = ½ [a(j, sire) + a(j, dam)]

**Inbreeding.** F(i) is the probability that the two alleles at a locus are identical by descent. It rises when relatives are mated.

**A⁻¹.** A is dense, but its inverse is sparse and can be written down directly from the pedigree (Henderson 1976; Quaas 1976). That is what makes evaluations of millions of animals feasible.

## Worked example

Use `data/toy_pedigree.csv`. Animals 4 and 5 are half sibs (common sire 1), so a(4, 5) = 0.25 and their offspring 6 has F = ½ × 0.25 = 0.125.

The function in `modules/03-blup-animal-model.R` (Step 1) builds A with the tabular method. Compare it with:

```r
library(nadiv)
ped <- read.csv("data/toy_pedigree.csv")
ped[ped == 0] <- NA
A <- as.matrix(makeA(ped[, c("id", "sire", "dam")]))
```

## Exercises

1. By hand, compute a(6, 7) and the inbreeding coefficient of animal 8.
2. Add a new animal 9 with sire 8 and dam 5. Compute F(9).
3. Count the non-zero elements of A and of A⁻¹. Why does the difference matter for large populations?

## Further reading

- Mrode, *Linear Models for the Prediction of Animal Breeding Values*, chapter 2
