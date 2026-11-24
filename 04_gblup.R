# =============================================================================
# LESSON 4  --  GBLUP: replace the pedigree with DNA
#
# Goal: build the genomic relationship matrix G from SNP genotypes and use it
#       in place of A. Only GENOTYPED animals take part.
#
# Words: an ALLELE is one version of a gene or DNA position. Each animal has
# two copies of every position (one from each parent). The SNP code counts
# copies of one chosen allele: 0 or 2 = HOMOZYGOUS (both copies the same),
# 1 = HETEROZYGOUS (one of each).
# =============================================================================
source("R/functions.R")
sigma2_a <- 90; sigma2_e <- 210

# -- 1. Genotypes ----------------------------------------------------------------------
M <- read_genotypes("data/genotypes.txt")
cat("Genotype matrix:", nrow(M), "animals x", ncol(M), "SNPs\n")
print(M[1:5, 1:10])
M <- qc_genotypes(M, min_maf = 0.05)

# -- 2. G: relationships measured from DNA ---------------------------------------------
G <- make_G(M)
# G is centred on the genotyped animals' own allele frequencies, so it
# measures relationships relative to the CURRENT population: the average pair
# is about 0, and many pairs are negative ("less related than average").

ped <- prepare_pedigree(read.csv("data/pedigree.csv", colClasses = "character"))
A   <- make_A(ped$s, ped$d); dimnames(A) <- list(ped$id, ped$id)
g   <- match(rownames(G), ped$id)
A22 <- A[g, g]                              # pedigree relationships, same animals
up  <- upper.tri(G)
plot(A22[up], G[up], pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.2),
     xlab = "Pedigree relationship (A)", ylab = "Genomic relationship (G)")
abline(0, 1, col = "red")

hs <- is_half_sib(ped$sire[g], ped$dam[g])[up]   # true half-sib pairs
cat(sprintf("Half-sib pairs: A says 0.25 for all; G ranges %.2f to %.2f (sd %.3f)\n",
            min(G[up][hs]), max(G[up][hs]), sd(G[up][hs])))
# The spread comes from Mendelian sampling: which half of each parent's DNA
# each calf happened to receive. That detail is what genomics adds.

# -- 3. GBLUP --------------------------------------------------------------------------
# Only genotyped animals can be in the model. Records of non-genotyped
# animals are DROPPED -- the main weakness of GBLUP.
phen   <- read.csv("data/phenotypes.csv", colClasses = c(id = "character"))
phen_g <- phen[phen$id %in% rownames(G), ]
cat(sprintf("Records usable by GBLUP: %d of %d\n", nrow(phen_g), nrow(phen)))

Gstar <- G + diag(0.01, nrow(G))           # tiny ridge so G can be inverted
X <- model.matrix(~ 0 + cg + sex, data = phen_g)[, , drop = FALSE]
Z <- make_Z(phen_g$id, rownames(G))
fit <- solve_mme(phen_g$weaning_wt, X, Z, solve(Gstar), sigma2_a, sigma2_e, Kdiag = diag(Gstar))

gebv <- data.frame(id = rownames(G), gebv = round(fit$ebv, 2), reliability = round(fit$rel, 3))
dir.create("output", showWarnings = FALSE)
write.csv(gebv, "output/gebv_gblup.csv", row.names = FALSE)
print(head(gebv[order(-gebv$gebv), ], 10), row.names = FALSE)

# -- EXERCISES -------------------------------------------------------------------------
# 1. Take the two full sibs from Lesson 3 (all generation-3 calves are
#    genotyped). Are their GEBVs still identical? Why not?
# 2. Look at diag(G). Values are not exactly 1. What does a value above 1 mean?
# 3. Re-run with min_maf = 0.20. How many SNPs are left? Measure how far G
#    is from A with sd(G[up] - A22[up]) for both settings. Does G get noisier?
#    (Afterwards set min_maf back to 0.05 and re-run, because Lesson 6 reads
#    the output file.)
