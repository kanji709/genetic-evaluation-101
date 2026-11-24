# Module 05 — Genomic evaluation

## Learning objectives

- Explain how SNP markers capture realised (rather than expected) relationships
- Build a genomic relationship matrix **G** (VanRaden method 1)
- Fit GBLUP by replacing A with G
- Describe single-step GBLUP (ssGBLUP) and why it combines genotyped and non-genotyped animals

## Core concepts

**Markers.** SNP genotypes coded 0/1/2 (copies of one allele). Dense SNP panels are in linkage disequilibrium with the causal variants.

**Expected vs realised relationships.** A full-sib pair has expected a = 0.5, but the realised proportion shared varies (≈ 0.4–0.6) due to Mendelian sampling. Genotypes measure it.

**VanRaden G (method 1).**

  G = ZZ′ / [2 Σ pⱼ(1 − pⱼ)],  where Z = M − 2p (centred genotype matrix)

**GBLUP.** The animal model with G in place of A — only genotyped animals.

**Single-step (ssGBLUP).** Uses H, which blends A and G, so all animals (genotyped or not) are evaluated together. In practice H⁻¹ = A⁻¹ + [0 0; 0 G⁻¹ − A₂₂⁻¹].

**Reference population.** Accuracy of genomic prediction depends on the number of phenotyped and genotyped animals, h² and relatedness to the selection candidates.

## Worked example

```r
set.seed(1)
n <- 10; m <- 200
M <- matrix(rbinom(n * m, 2, 0.3), n, m)   # simulated genotypes (0/1/2)
p <- colMeans(M) / 2
Zc <- sweep(M, 2, 2 * p)
G <- tcrossprod(Zc) / (2 * sum(p * (1 - p)))
round(G[1:5, 1:5], 3)
```

Compare with `AGHmatrix::Gmatrix(M, method = "VanRaden")`.

## Exercises

1. Why do unrelated animals in the simulation still show non-zero off-diagonals in G?
2. Simulate two full sibs from known parental haplotypes and compare their G element with A = 0.5.
3. List two reasons genotyped young animals get more accurate EBVs than with pedigree BLUP.

## Further reading

- Meuwissen, Hayes & Goddard (2001) Prediction of total genetic value using genome-wide dense marker maps. *Genetics* 157: 1819–1829
- VanRaden (2008) Efficient methods to compute genomic predictions. *J. Dairy Sci.* 91: 4414–4423
- Legarra, Aguilar & Misztal (2009) A relationship matrix including full pedigree and genomic information. *J. Dairy Sci.* 92: 4656–4663
