# =============================================================================
# LESSON 6b  --  One herd is not enough: repeat on 20 simulated herds
#
# Each run simulates a new herd and fits pedigree BLUP, GBLUP, ssGBLUP and
# ssGBLUP with G built from the founders' allele frequencies (possible only
# in a simulation). Takes about a minute. Does not touch data/ or output/.
# =============================================================================
source("R/functions.R")
seeds <- 1:20
res <- do.call(rbind, lapply(seeds, function(s) { cat("herd", s, "\n"); evaluate_herd(s) }))

summ <- function(x) sprintf("%.2f (%.2f)", mean(x), sd(x))
print(data.frame(accuracy   = tapply(res$accuracy,   res$method, summ),
                 slope      = tapply(res$slope,      res$method, summ),
                 gain_error = tapply(res$gain_error, res$method, summ)))
# Mean (SD) over herds.

acc  <- reshape(res[, c("seed", "method", "accuracy")], idvar = "seed",
                timevar = "method", direction = "wide")
cat("\nGBLUP beats BLUP in", sum(acc$accuracy.GBLUP > acc$accuracy.BLUP), "of", length(seeds), "herds\n")
cat("ssGBLUP beats BLUP in", sum(acc$accuracy.ssGBLUP > acc$accuracy.BLUP), "of", length(seeds), "herds\n")
