# Genetic Evaluation 101

A hands-on introduction to genetic evaluation for students and interns: how we go from pedigrees, phenotypes and genotypes to **estimated breeding values (EBVs)** that breeders can act on.

The tutorial is built around short reading notes plus small, runnable R examples on toy data, so every equation can be checked by hand and then by code.

## Who this is for

- Undergraduate / postgraduate students new to quantitative genetics
- Interns joining an animal or plant breeding analytics team
- Data scientists who know statistics but not breeding

Assumed background: basic statistics (mean, variance, regression) and some exposure to R. Matrix algebra helps but is introduced as needed.

## Learning outcomes

By the end, you should be able to:

1. Explain phenotype = genotype + environment, heritability and breeding value
2. Build a numerator relationship matrix (**A**) from a pedigree
3. Set up and solve Henderson's mixed model equations for an animal model (BLUP)
4. Describe how variance components are estimated (REML) and why they matter
5. Explain genomic relationships (**G**), GBLUP and single-step evaluation
6. Interpret EBV accuracy and predict response to selection

## Modules

| # | Module | Key ideas |
|---|--------|-----------|
| 01 | [Foundations of quantitative genetics](modules/01-foundations.md) | P = G + E, additive effects, heritability, breeding values |
| 02 | [Pedigrees and relationships](modules/02-pedigree-relationships.md) | Pedigree checks, inbreeding, the A matrix and its inverse |
| 03 | [BLUP and the animal model](modules/03-blup-animal-model.md) | Mixed models, Henderson's MME, fixed vs random effects |
| 04 | [Estimating variance components](modules/04-variance-components.md) | ML vs REML, heritability estimates, software |
| 05 | [Genomic evaluation](modules/05-genomic-evaluation.md) | Markers, G matrix, GBLUP, single-step (ssGBLUP) |
| 06 | [Accuracy, selection and reporting](modules/06-accuracy-selection.md) | Reliability, breeder's equation, communicating EBVs |

Each module has: learning objectives, core concepts, a worked example, exercises and further reading.

## Getting started

1. Install [R](https://cran.r-project.org/) (≥ 4.2) and, optionally, [RStudio](https://posit.co/download/rstudio-desktop/) or [Quarto](https://quarto.org/).
2. Clone this repository:
   ```bash
   git clone https://github.com/kanji709/genetic-evaluation-101.git
   cd genetic-evaluation-101
   ```
3. Install the packages used in later modules:
   ```r
   source("install.R")
   ```
4. Run the first worked example (base R only, no packages needed):
   ```r
   source("modules/03-blup-animal-model.R")
   ```

## Repository structure

```
genetic-evaluation-101/
├── README.md            # You are here
├── install.R            # Installs R packages used in the tutorial
├── CONTRIBUTING.md      # How students/interns can suggest fixes and additions
├── LICENSE              # MIT (code); tutorial text and figures are CC BY 4.0
├── data/                # Small simulated datasets used in examples
│   ├── README.md
│   ├── toy_pedigree.csv
│   └── toy_phenotypes.csv
└── modules/             # One file of notes per module (+ worked R examples)
    ├── 01-foundations.md
    ├── ...
    └── 03-blup-animal-model.R
```

## Suggested pace

| Format | Plan |
|--------|------|
| Intensive (interns) | One module per day over a week, with day 6 for a mini project |
| Semester (students) | One module per week, alternating lecture and practical |

## Further reading

- Falconer, D.S. & Mackay, T.F.C. *Introduction to Quantitative Genetics* (4th ed.)
- Mrode, R.A. *Linear Models for the Prediction of Animal Breeding Values* (3rd ed.)
- Lynch, M. & Walsh, B. *Genetics and Analysis of Quantitative Traits*

## Licence

Code is released under the [MIT License](LICENSE). Tutorial text and figures are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); please credit "Genetic Evaluation 101, Jie Kang".

## Contact

Maintained by Jie Kang ([@kanji709](https://github.com/kanji709)). Questions and corrections are welcome as [GitHub issues](https://github.com/kanji709/genetic-evaluation-101/issues).
