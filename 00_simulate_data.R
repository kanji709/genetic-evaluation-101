# =============================================================================
# 00_simulate_data.R  --  writes the practice data set to data/
#
# You don't need to run this: the files are already there. The simulation
# itself is simulate_herd() in R/functions.R (read it if you want to know how
# the "truth" was made).
#
# If you change the seed to get a different herd, first make a backup copy of
# data/, and change the seed back to 2026 afterwards: every lesson reads data/.
# =============================================================================
source("R/functions.R")
herd <- simulate_herd(seed = 2026)

dir.create("data", showWarnings = FALSE)
write.csv(herd$pedigree,   "data/pedigree.csv",   row.names = FALSE, quote = FALSE)
write.csv(herd$phenotypes, "data/phenotypes.csv", row.names = FALSE, quote = FALSE)
write.csv(herd$true_bv,    "data/true_bv.csv",    row.names = FALSE, quote = FALSE)
M <- herd$genotypes
writeLines(paste(rownames(M), apply(M, 1, paste, collapse = "")), "data/genotypes.txt")

cat(sprintf("Animals: %d   weights: %d   genotyped: %d   SNPs: %d\n",
            nrow(herd$pedigree), nrow(herd$phenotypes), nrow(M), ncol(M)))
