# Data

All datasets here are **simulated toy data** for teaching. They are deliberately tiny so every calculation can be checked by hand.

## `toy_pedigree.csv`

8 animals across three generations. Unknown parents are coded `0`.

| Column | Description |
|--------|-------------|
| `id` | Animal ID (parents always appear before their offspring) |
| `sire` | Sire ID (`0` = unknown) |
| `dam` | Dam ID (`0` = unknown) |
| `sex` | `M` / `F` |
| `birth_year` | Year of birth |

Note: animal 6 is the offspring of half-sibs (4 and 5 share sire 1), so it is inbred (F = 0.125). Animal 8 is also inbred — a good exercise for Module 02.

## `toy_phenotypes.csv`

Weaning weight records for animals 4–8. Animals 1–3 have no records but still receive EBVs through their relatives.

| Column | Description |
|--------|-------------|
| `id` | Animal ID (matches the pedigree) |
| `herd` | Herd (fixed effect): `H1`, `H2` |
| `sex` | Sex (fixed effect): `M`, `F` |
| `weaning_weight_kg` | Phenotype (kg) |

## Rules for adding data

- Simulated or openly licensed public data only — never real farm, client or animal records
- Keep files small (< 1 MB); link to larger public datasets instead of committing them
- Document every new file in this README
