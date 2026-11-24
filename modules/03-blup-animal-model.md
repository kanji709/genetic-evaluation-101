# Module 03 — BLUP and the animal model

## Learning objectives

- Write an animal model in matrix form and identify fixed vs random effects
- Set up Henderson's mixed model equations (MME)
- Solve the MME for BLUE of fixed effects and BLUP of breeding values
- Explain how animals without records receive EBVs

## Core concepts

**Animal model.**

  y = Xb + Za + e,  a ~ N(0, Aσ²a),  e ~ N(0, Iσ²e)

- **Fixed effects (b):** systematic environmental effects estimated from the data — herd, year, sex, contemporary group
- **Random effects (a):** breeding values, assumed to come from a distribution with known covariance (A)

**Mixed model equations.** With α = σ²e / σ²a:

```
[ X'X    X'Z          ] [ b̂ ]   [ X'y ]
[ Z'X    Z'Z + A⁻¹α   ] [ â ] = [ Z'y ]
```

**Why BLUP?** Best (minimum prediction error variance), Linear (in y), Unbiased, Prediction. It adjusts for fixed effects and uses information from all relatives at once through A⁻¹.

**Shrinkage.** The A⁻¹α term pulls EBVs towards zero; the lower the heritability (higher α), the stronger the shrinkage.

## Worked example

Run the script:

```r
source("modules/03-blup-animal-model.R")
```

It reads the toy data, builds A, sets up X and Z, solves the MME and prints EBVs and their accuracies. Animals 1–3 have no records but still get EBVs through their offspring.

## Exercises

See the exercise list at the end of `03-blup-animal-model.R`, plus:

1. Write out X and Z for the toy data by hand and check them against the script.
2. Why would ranking animals on raw weaning weight be misleading here?

## Further reading

- Henderson, C.R. (1975) Best linear unbiased estimation and prediction under a selection model. *Biometrics* 31: 423–447
- Mrode, chapter 3
