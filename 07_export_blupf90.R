# =============================================================================
# LESSON 7 (optional)  --  The same evaluation in BLUPF90
#
# The R code in Lessons 3-5 stores full matrices, so it stops working above
# roughly 10,000 animals. Production evaluations use dedicated software.
# BLUPF90 (University of Georgia; free for research, check the licence for
# commercial use) is the most common for single-step. This script writes the
# practice data in BLUPF90 format plus two parameter files.
#
# Install: download renumf90 and blupf90+ for your system from
#   https://nce.ads.uga.edu/html/projects/programs/  (Linux, Windows, Mac_OSX/64bit)
# into the blupf90/ folder. macOS/Linux: chmod +x renumf90 blupf90+
# The Mac programs are Intel builds; Apple Silicon Macs run them through
# Rosetta 2 (softwareupdate --install-rosetta). Then, in a terminal inside blupf90/:
#
#     ./renumf90 renum_blup.par      # BLUP: renumbers IDs, writes renf90.par
#     ./blupf90+ renf90.par          # solves; writes solutions.orig
#
#     ./renumf90 renum_ssgblup.par   # ssGBLUP: same, plus genotypes
#     ./blupf90+ renf90.par
#
# Tested on the practice herd: BLUP and ssGBLUP solutions agree with Lessons
# 3 and 5 within 0.006 kg. Notes:
#  * solutions.orig columns: trait, effect, level, original_id, solution, s.e.
#    Effect 1 = contemporary group, 2 = sex, 3 = animal.
#  * BLUPF90 sets one fixed-effect level to zero (males), so cg and sex
#    solutions look different from R's; differences between levels are equal.
#  * Inbreeding is included in A-inverse by default (renf90.inb).
#  * G is tuned, then blended 95% G + 5% A22, as in Lesson 5.
#  * PED_DEPTH 0 keeps the whole pedigree (default: 3 generations of ancestors).
# =============================================================================
dir.create("blupf90", showWarnings = FALSE)
ped  <- read.csv("data/pedigree.csv",  colClasses = "character")
phen <- read.csv("data/phenotypes.csv", colClasses = "character")

# BLUPF90 wants plain space-separated text, no header, missing parent = 0
write.table(phen[, c("id", "cg", "sex", "weaning_wt")], "blupf90/phen.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ped[, c("id", "sire", "dam")], "blupf90/ped.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
invisible(file.copy("data/genotypes.txt", "blupf90/geno.txt", overwrite = TRUE))

par_text <- function(with_snp) c(
  "DATAFILE", "phen.txt",
  "TRAITS", "4",                       # column 4 = weaning_wt
  "FIELDS_PASSED TO OUTPUT", "",
  "WEIGHT(S)", "",
  "RESIDUAL_VARIANCE", "210",
  "EFFECT", "2 cross alpha",           # column 2 = contemporary group (fixed)
  "EFFECT", "3 cross alpha",           # column 3 = sex (fixed)
  "EFFECT", "1 cross alpha",           # column 1 = animal ID
  "RANDOM", "animal",
  "FILE", "ped.txt",
  "FILE_POS", "1 2 3 0 0",             # animal, sire, dam columns in ped.txt
  if (with_snp) c("SNP_FILE", "geno.txt"),
  "PED_DEPTH", "0",                    # keep all generations of the pedigree
  "(CO)VARIANCES", "90",
  "OPTION origID",                     # original IDs in the solutions output
  "OPTION sol se")                     # also write standard errors
writeLines(par_text(FALSE), "blupf90/renum_blup.par")
writeLines(par_text(TRUE),  "blupf90/renum_ssgblup.par")
cat("Wrote blupf90/phen.txt, ped.txt, geno.txt, renum_blup.par, renum_ssgblup.par\n")

# Afterwards, compare with the R results (run Lessons 3 and 5 first):
#   sol <- read.table("blupf90/solutions.orig", header = TRUE,
#                     colClasses = c(original_id = "character"))
#   ani <- sol[sol$effect == 3, ]
#   r   <- read.csv("output/gebv_ssgblup.csv", colClasses = c(id = "character"))
#   cor(ani$solution, r$gebv[match(ani$original_id, r$id)])   # about 1
