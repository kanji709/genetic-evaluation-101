# Module 04 — Estimating variance components

## Learning objectives

- Explain why BLUP needs variance components and where they come from
- Contrast ML and REML and explain why REML is preferred
- Fit an animal model with REML in R and report heritability with a standard error
- Recognise common problems: small data, confounding, convergence failures

## Core concepts

**Why it matters.** α = σ²e / σ²a controls shrinkage in the MME. Wrong variance components give biased rankings and misleading accuracies.

**Maximum likelihood (ML)** estimates fixed effects and variances together, but ignores the degrees of freedom used by fixed effects, so variances are biased downwards.

**Restricted/residual maximum likelihood (REML)** maximises the likelihood of error contrasts (data with fixed effects removed), giving less biased estimates. It is the industry standard.

**Algorithms and software.** AI-REML, EM-REML; common tools include ASReml, BLUPF90 (AIREMLF90), DMU, WOMBAT and R packages such as `sommer`.

## Worked example

The toy data are too small for reliable REML. Simulate a larger dataset (e.g. 500 animals), then fit:

```r
library(sommer)
# ped: id, sire, dam; dat: id, herd, sex, y
A <- nadiv::makeA(ped)
fit <- mmer(y ~ herd + sex,
            random = ~ vsr(id, Gu = as.matrix(A)),
            rcov   = ~ units,
            data   = dat)
summary(fit)$varcomp
vpredict(fit, h2 ~ V1 / (V1 + V2))
```

*(Template — adapt to the version of `sommer` installed; newer versions also provide `mmes()`.)*

## Exercises

1. Simulate data with a known h² of 0.3. How close is the REML estimate with 100, 500 and 2000 animals?
2. Why is the standard error of h² so large with small datasets?
3. What happens if herd is confounded with sire (each sire used in only one herd)?

## Further reading

- Patterson & Thompson (1971) Recovery of inter-block information when block sizes are unequal. *Biometrika* 58: 545–554
- Mrode, chapter 11
