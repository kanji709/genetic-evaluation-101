# =============================================================================
# sync_functions.R  --  copy every function from R/functions.R into
#                       tutorial.qmd, so the standalone tutorial and the lesson
#                       scripts always run the same code.
#
# Run from the repository root after editing R/functions.R:
#     Rscript tools/sync_functions.R
#
# How it works: each function in R/functions.R starts with a marker line
#     ## ---- name ----
# and tutorial.qmd shows it in a chunk labelled fn-name (or fn-name-<suffix>).
# The code of those chunks is replaced; chunk options (#| lines) are kept.
# =============================================================================
src    <- readLines("R/functions.R")
starts <- grep("^## ---- [A-Za-z_]+ ----$", src)
ends   <- c(starts[-1] - 1, length(src))
fns <- Map(function(a, b) {
  x <- src[(a + 1):b]
  while (length(x) > 0 && x[length(x)] == "") x <- x[-length(x)]
  x
}, starts, ends)
names(fns) <- sub("^## ---- ([A-Za-z_]+) ----$", "\\1", src[starts])

qmd <- readLines("tutorial.qmd")
out <- character(0)
used <- character(0)
i <- 1
while (i <= length(qmd)) {
  out <- c(out, qmd[i])
  if (grepl("^```\\{r\\}", qmd[i])) {
    j <- i + 1
    while (j <= length(qmd) && grepl("^#\\|", qmd[j])) j <- j + 1
    opts  <- if (j > i + 1) qmd[(i + 1):(j - 1)] else character(0)
    label <- sub("^#\\|\\s*label:\\s*", "", grep("^#\\|\\s*label:", opts, value = TRUE))
    name  <- sub("^fn-([A-Za-z_]+)(-.*)?$", "\\1", label)
    if (length(label) == 1 && grepl("^fn-", label) && name %in% names(fns)) {
      k <- j
      while (qmd[k] != "```") k <- k + 1          # end of the chunk
      out  <- c(out, opts, fns[[name]], "```")
      used <- c(used, name)
      i <- k + 1
      next
    }
  }
  i <- i + 1
}

missing <- setdiff(names(fns), used)
writeLines(out, "tutorial.qmd")
cat(sprintf("Updated %d chunks in tutorial.qmd.\n", length(used)))
if (length(missing) > 0)
  cat("Not shown in tutorial.qmd (scripts only):", paste(missing, collapse = ", "), "\n")
