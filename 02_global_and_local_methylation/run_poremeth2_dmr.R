#!/usr/bin/env Rscript
# PoreMeth2 segmentation of the methylation difference (Figure 11a).
# Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). 
suppressWarnings(suppressMessages({
  .libPaths(c(path.expand("~/R/x86_64-pc-linux-gnu-library/4.5"), .libPaths()))
  library(PoreMeth2)
}))

args  <- commandArgs(trailingOnly = TRUE)
omega <- if (length(args) >= 1) as.numeric(args[1]) else 0.1
eta   <- if (length(args) >= 2) as.numeric(args[2]) else 1e-5
FW    <- if (length(args) >= 3) as.integer(args[3])  else 3L
tag   <- if (length(args) >= 4) args[4] else sprintf("omega%s_eta%s_FW%d", omega, format(eta, scientific = TRUE), FW)

P   <- "/data/8oxo_project/AdMSC_DMR"
IN  <- Sys.getenv("POREMETH_IN", unset = file.path(P, "readlevel_calls", "poremeth"))
OUT <- file.path(P, "results", paste0("poremeth2_dmr_", tag, ".tsv"))
cat(sprintf("input: %s\n", IN))

cn <- c("chrom", "pos", "entropy", "entropy_cov", "beta", "beta_cov")
rd <- function(arm) {
  f <- file.path(IN, paste0(arm, ".modkit_adapted_sorted.entropy.file.tsv"))
  d <- read.table(f, sep = "\t", header = FALSE, col.names = cn,
                  colClasses = c("character", "integer", "numeric", "integer", "numeric", "integer"))
  cat(sprintf("  %s: %d positions, %d contigs\n", arm, nrow(d), length(unique(d$chrom))))
  d
}

cat(sprintf("PoreMeth2DMR  omega=%g  eta=%g  FW=%d\n", omega, eta, FW))
cat("reading input\n")
TableControl <- rd("NT")          # untreated
TableTest    <- rd("PA")          # palmitate

cat("calling DMRs (test = PA, control = NT)\n")
t0  <- Sys.time()
dmr <- PoreMeth2DMR(TableTest, TableControl, omega = omega, eta = eta, FW = FW)
cat(sprintf("elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

dmr <- as.data.frame(dmr, stringsAsFactors = FALSE)
cat("\ncolumns returned by PoreMeth2DMR:\n"); print(colnames(dmr)); print(utils::head(dmr, 3))

write.table(dmr, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nwrote %d segments -> %s\n", nrow(dmr), OUT))

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
db  <- NULL
for (nm in c("DeltaBeta", "deltaBeta", "Delta_Beta", "dBeta")) if (nm %in% colnames(dmr)) db <- num(dmr[[nm]])
if (is.null(db)) { i <- grep("beta", colnames(dmr), ignore.case = TRUE); if (length(i)) db <- num(dmr[[i[1]]]) }

if (!is.null(db)) {
  keep <- !is.na(db) & abs(db) >= 0.2
  cat(sprintf("\nsegments total                      %d\n", nrow(dmr)))
  cat(sprintf("segments with |delta beta| >= 0.2   %d (%.1f%%)\n", sum(keep), 100 * sum(keep) / nrow(dmr)))
  cat(sprintf("  of those, hypermethylated in PA   %d\n", sum(keep & db > 0)))
  cat(sprintf("  of those, hypomethylated  in PA   %d\n", sum(keep & db < 0)))
} else {
  cat("\n(could not identify the delta-beta column automatically; inspect the table)\n")
}
