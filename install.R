# Installs the R packages used in Genetic Evaluation 101.
# Modules 01-03 use base R only; later modules use the packages below.

pkgs <- c(
  "nadiv",      # pedigree tools: A matrix, inverse, inbreeding (Module 02)
  "sommer",     # mixed models and REML (Modules 04-05)
  "AGHmatrix"   # genomic relationship matrices (Module 05)
)

missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All tutorial packages are already installed.")
}
