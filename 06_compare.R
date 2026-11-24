# =============================================================================
# LESSON 6  --  Which method predicted best?
#
# Goal: check the predictions for the generation-3 calves (no weights yet, all
#       genotyped) against their TRUE breeding values. In real life the truth
#       is never known; breeders hide the latest records, predict them, and
#       compare (validation, e.g. the LR method of Legarra & Reverter 2018).
#
# Run Lessons 3, 4 and 5 first: they write the files read here.
# =============================================================================
source("R/functions.R")
truth <- read.csv("data/true_bv.csv",        colClasses = c(id = "character"))
blup  <- read.csv("output/ebv_blup.csv",     colClasses = c(id = "character"))
gblup <- read.csv("output/gebv_gblup.csv",   colClasses = c(id = "character"))
ss    <- read.csv("output/gebv_ssgblup.csv", colClasses = c(id = "character"))

pred <- list(BLUP    = blup$ebv[match(truth$id, blup$id)],
             GBLUP   = gblup$gebv[match(truth$id, gblup$id)],   # NA if not genotyped
             ssGBLUP = ss$gebv[match(truth$id, ss$id)])
res <- t(sapply(pred, evaluate_prediction, true_bv = truth$true_bv, generation = truth$generation))
print(round(res, 3))

# How to read this:
#  accuracy   = correlation between predicted and true breeding value
#               for the young calves. Higher = better ranking.
#  slope      = regression of true on predicted value. 1 = right spread;
#               below 1 = predictions too spread out (over-promising).
#  gain_error = estimated minus true genetic gain from generation 0 to 3, kg.
#               0 = no level bias. Accuracy and slope cannot see this.
#
# One small herd is a noisy test: these numbers move a lot from herd to herd.
# Lesson 6b (06b_replicates.R) repeats everything on 20 herds.

young <- truth$generation == 3
op <- par(mfrow = c(1, 3))
for (m in names(pred)) {
  plot(pred[[m]][young], truth$true_bv[young], pch = 16, col = rgb(0.1, 0.3, 0.6, 0.5),
       xlab = paste("Predicted,", m), ylab = "True breeding value",
       main = sprintf("%s   r = %.2f", m, res[m, "accuracy"]))
  abline(lm(truth$true_bv[young] ~ pred[[m]][young]), col = "red")
}
par(op)

# -- EXERCISES -------------------------------------------------------------------------
# 1. Which method wins in this herd? Is the difference bigger than the
#    herd-to-herd spread you see in 06b_replicates.R?
# 2. In what situation would you expect GBLUP to beat pedigree BLUP?
# 3. Which method has the largest gain_error? What could cause it?
