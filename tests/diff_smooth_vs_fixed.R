# ============================================================
# Cell-diff: smooth manifest (T_run_smooth) vs greedy manifest (T_run_fixed)
# ============================================================
# For every sweep summary present in BOTH dirs, aligns cells on their grid
# identifiers and compares each <metric>_mean using the combined SE:
#     z = (mean_smooth - mean_fixed) / sqrt(se_smooth^2 + se_fixed^2)
# Reports cells that moved beyond a z threshold, with special attention to the
# forward-vs-myopic contrasts (fwd_vs_myo_P / _PG / _PS) that hosted the greedy
# resonance: those SHOULD shrink at the resonance cells (artifact removed) while
# the genuine forward advantage elsewhere is preserved. A move >Z in any OTHER
# headline metric is a flag to investigate.
#
# Usage:  Rscript tests/diff_smooth_vs_fixed.R [Zthresh=3] [dirA=T_run_smooth] [dirB=T_run_fixed]
# Non-destructive: reads only. Writes a summary to stdout + tests/diff_smooth_vs_fixed.log.

setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
args  <- commandArgs(trailingOnly = TRUE)
Z     <- if (length(args) >= 1) as.numeric(args[1]) else 3
dirS  <- file.path("sweep_results", if (length(args) >= 2) args[2] else "T_run_smooth")
dirF  <- file.path("sweep_results", if (length(args) >= 3) args[3] else "T_run_fixed")

CONTRASTS <- c("fwd_vs_myo_P", "fwd_vs_myo_PG", "fwd_vs_myo_PS")

files <- sort(intersect(
  sub("_summary\\.rds$", "", basename(Sys.glob(file.path(dirS, "*_summary.rds")))),
  sub("_summary\\.rds$", "", basename(Sys.glob(file.path(dirF, "*_summary.rds"))))
))
if (length(files) == 0) stop("No overlapping *_summary.rds in\n  ", dirS, "\n  ", dirF,
                             "\n(has the smooth re-run finished writing summaries?)")

cat(sprintf("Cell-diff  smooth[%s]  vs  fixed[%s]   z-threshold=%.1f\n", dirS, dirF, Z))
cat(sprintf("%d overlapping sweeps: %s\n\n", length(files), paste(files, collapse = ", ")))

grand_flags <- 0; grand_cmp <- 0
worst <- list()   # per non-contrast metric worst mover

for (nm in files) {
  S <- as.data.frame(readRDS(file.path(dirS, paste0(nm, "_summary.rds"))))
  F <- as.data.frame(readRDS(file.path(dirF, paste0(nm, "_summary.rds"))))
  # id cols = everything that is not <x>_mean / <x>_se / n_trials
  is_ms <- grepl("_(mean|se)$", names(S)) | names(S) == "n_trials"
  idcols <- names(S)[!is_ms]
  # align rows by grid identity
  key <- function(df) do.call(paste, c(df[idcols], sep = "|"))
  F <- F[match(key(S), key(F)), , drop = FALSE]
  if (anyNA(F[[1]])) { cat(sprintf("  [%s] SKIP — grid mismatch\n", nm)); next }

  bases <- sub("_mean$", "", grep("_mean$", names(S), value = TRUE))
  bases <- bases[paste0(bases, "_se") %in% names(S)]

  ncells <- nrow(S); nflag <- 0; con_lines <- c(); oth_worst_z <- 0; oth_worst <- ""
  for (b in bases) {
    ms <- S[[paste0(b, "_mean")]]; mf <- F[[paste0(b, "_mean")]]
    ss <- S[[paste0(b, "_se")]];   sf <- F[[paste0(b, "_se")]]
    d  <- ms - mf
    se <- sqrt(pmax(ss, 0)^2 + pmax(sf, 0)^2)
    z  <- ifelse(se > 0, d / se, ifelse(abs(d) > 1e-9, sign(d) * Inf, 0))
    z[is.na(z)] <- 0
    grand_cmp <- grand_cmp + sum(!is.na(ms) & !is.na(mf))
    fl <- which(abs(z) > Z)
    nflag <- nflag + length(fl)
    if (b %in% CONTRASTS) {
      # show the biggest-magnitude fixed cells for this contrast (resonance/real signal)
      ord <- order(-abs(mf))[seq_len(min(3, length(mf)))]
      for (i in ord) con_lines <- c(con_lines, sprintf(
        "      %-14s cell%2d: fixed %+.4f (se %.4f) -> smooth %+.4f (se %.4f)  z=%+.1f",
        b, i, mf[i], sf[i], ms[i], ss[i], z[i]))
    } else {
      wi <- which.max(abs(z))
      if (length(wi) && abs(z[wi]) > oth_worst_z) { oth_worst_z <- abs(z[wi]); oth_worst <- sprintf("%s@cell%d (z=%+.1f, %+.4f->%+.4f)", b, wi, z[wi], mf[wi], ms[wi]) }
    }
  }
  grand_flags <- grand_flags + nflag
  worst[[nm]] <- list(z = oth_worst_z, desc = oth_worst)
  cat(sprintf("[%s] %d cells x %d metrics | flags(|z|>%.0f)=%d | worst non-contrast: %s\n",
              nm, ncells, length(bases), Z, nflag, if (nzchar(oth_worst)) oth_worst else "-"))
  if (length(con_lines)) cat("    contrasts (artifact should shrink; real signal persist):\n",
                             paste(con_lines, collapse = "\n"), "\n", sep = "")
}

cat(sprintf("\n=== SUMMARY ===\n%d total mean-cell comparisons, %d moved |z|>%.0f (%.1f%%; ~%.1f%% expected by chance).\n",
            grand_cmp, grand_flags, Z, 100 * grand_flags / max(grand_cmp, 1),
            100 * 2 * pnorm(-Z)))
oz <- sapply(worst, `[[`, "z")
cat("Largest NON-contrast movers across sweeps (these are the ones to scrutinize):\n")
for (nm in names(sort(oz, decreasing = TRUE))[seq_len(min(5, length(oz)))])
  cat(sprintf("  %-20s %s\n", nm, worst[[nm]]$desc))
cat("\nInterpretation: contrast metrics dropping at high-|fixed| cells = artifact removed (good).\n",
    "Non-contrast movers >~4 SE that are NOT the resonance cells warrant a look.\n", sep = "")
