# Breeding Values from Scratch

A hands-on course in genetic evaluation for interns and summer students: pedigree **BLUP**,
genomic **GBLUP** and **single-step GBLUP**, in plain R. Basic R is assumed; no genetics.
Every method is checked against a simulated beef herd whose true breeding values are known.

## Getting started

1. Clone or download this repository.
2. Open `blup_tutorial.Rproj` in RStudio (or open the folder in Positron), so the working
   directory is the repository root.
3. Render the tutorial and read it alongside the lesson scripts:

   ```bash
   quarto render tutorial.qmd    # -> tutorial.html
   ```

4. For each lesson, open the matching script, run it line by line (Ctrl/Cmd + Enter) and do the
   exercises. Each exercise set in the tutorial has a collapsed answer box.

**Requirements:** R ≥ 4.0. The lessons use base R only. Rendering needs
[Quarto](https://quarto.org) and the R packages `knitr` and `rmarkdown`. The first render takes
about a minute (Lesson 6 simulates 20 herds); later renders reuse `tutorial_cache/`.

## Lessons

| Lesson | Topic | Script |
|---|---|---|
| 1 | Breeding values and relatives | `01_relationships.R` |
| 2 | BLUP by hand (Mrode 2014, Example 3.1) | `02_blup_by_hand.R` |
| 3 | A full herd evaluation, with REML | `03_blup_herd.R` |
| 4 | GBLUP: using DNA | `04_gblup.R` |
| 5 | Single-step GBLUP | `05_ssgblup.R` |
| 6 | Which method wins? | `06_compare.R`, `06b_replicates.R` |
| 7 | Your own data and BLUPF90 (optional) | `07_export_blupf90.R` |

Lessons 3–5 write results to `output/`, which Lesson 6 reads, so run them in order. If an
exercise asks you to change a script, change it back afterwards.

## Repository layout

```
tutorial.qmd            the course; self-contained (simulates its data, defines all functions)
01_…07_*.R              lesson scripts (read data/, write output/)
00_simulate_data.R      regenerates data/ (seed 2026, same herd as the tutorial)
R/functions.R           all functions used by the scripts
data/                   practice herd: pedigree, phenotypes, genotypes, true breeding values
tools/sync_functions.R  copies R/functions.R into tutorial.qmd
```

Not tracked (see `.gitignore`): `tutorial.html` and Quarto build folders, `output/` and
`blupf90/` (written by the scripts).

## Data

Simulated beef herd, weaning weight (kg), h² = 0.30.

| File | Contents |
|---|---|
| `pedigree.csv` | 1,010 animals: id, sire, dam (0 = unknown), sex, generation |
| `phenotypes.csv` | 700 weaning weights with herd, generation, sex and contemporary group (`cg`) |
| `genotypes.txt` | 562 animals × 3,000 SNPs, BLUPF90 format (`id 0121...`) |
| `true_bv.csv` | Every animal's true breeding value. Real data never has this. |

## BLUPF90

Lesson 7 exports the herd and two parameter files for `renumf90` and `blupf90+`
(download from <https://nce.ads.uga.edu/html/projects/programs/>; free for research). Tested
with BLUPF90+ on macOS: pedigree BLUP and ssGBLUP solutions match the R lessons within 0.006 kg.

## Editing the course

- **Functions** live in `R/functions.R`. `tutorial.qmd` contains copies so it renders on its own.
  After changing a function, run `Rscript tools/sync_functions.R`, then re-render.
- **Text and exercises** are edited directly in `tutorial.qmd`. Numbers in the text are computed
  when it renders. If you change the simulation, delete `tutorial_cache/` and run
  `Rscript 00_simulate_data.R` so the scripts' data match.

## References

Mrode, R. A. (2014). *Linear Models for the Prediction of Animal Breeding Values*, 3rd ed. CABI.
The full reading list is at the end of the tutorial.
