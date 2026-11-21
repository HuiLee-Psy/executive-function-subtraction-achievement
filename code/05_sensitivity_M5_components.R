# ============================================================================
# code/05_sensitivity_M5_components.R
# Reproduces Section 3.6 & Supplementary Table S1:
#   Exploratory parallel mediation model (M5) decomposing EF into three manifest
#   component scores (WM_pro, FLEX, INH), with 5,000-resample percentile bootstrap.
# ============================================================================

suppressPackageStartupMessages({
  library(lavaan)
  library(readxl)
})
source("code/00_load_and_prep.R")

out_dir <- "results/05_sensitivity_M5"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "M5_console.log"), split = TRUE)

cat("05_sensitivity_M5_components.R | ", as.character(Sys.time()),
    " | lavaan ", as.character(packageVersion("lavaan")), "\n", sep = "")
df <- load_analysis_data()

m5 <- "
  Fluency  =~ FL_SIMP + FL_COMP
  Academic =~ AC_MATH + AC_CHIN + AC_ENG
  Fluency  ~ a1*WM_pro + a2*FLEX + a3*INH
  Academic ~ b*Fluency + c1*WM_pro + c2*FLEX + c3*INH
  ind_WM    := a1 * b
  ind_FLEX  := a2 * b
  ind_INH   := a3 * b
  ind_total := a1 * b + a2 * b + a3 * b
"

fit5 <- lavaan::sem(m5, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
fm5 <- lavaan::fitMeasures(fit5, c("chisq", "df", "pvalue", "cfi", "tli",
                                   "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"))
cat(sprintf(
  "\n[M5 Global Fit] chi2(%.0f) = %.2f, p = %.4f, CFI = %.4f, TLI = %.4f, RMSEA = %.4f [%.4f, %.4f], SRMR = %.4f\n",
  fm5["df"], fm5["chisq"], fm5["pvalue"], fm5["cfi"], fm5["tli"],
  fm5["rmsea"], fm5["rmsea.ci.lower"], fm5["rmsea.ci.upper"], fm5["srmr"]
))

pe5 <- lavaan::parameterEstimates(fit5, standardized = TRUE, ci = TRUE, level = 0.95)
paths5 <- pe5[pe5$op %in% c("~", ":=") &
              pe5$lhs %in% c("Fluency", "Academic", "ind_WM", "ind_FLEX", "ind_INH", "ind_total"),
              c("lhs", "op", "rhs", "est", "se", "ci.lower", "ci.upper", "pvalue", "std.all")]
cat("\n[M5 structural paths (delta method)]\n")
print(round_df(paths5, 4), row.names = FALSE)
write.csv(round_df(paths5, 4), file.path(out_dir, "M5_paths_delta.csv"), row.names = FALSE)

# 5,000-resample nonparametric percentile bootstrap (seed = 20260923)
set.seed(20260923)
n <- nrow(df); idx <- seq_len(n); R <- 5000L
labs <- c("ind_WM", "ind_FLEX", "ind_INH", "ind_total")
boot_beta <- matrix(NA_real_, R, length(labs), dimnames = list(NULL, paste0(labs, "_beta")))

cat(sprintf("\nRunning %d bootstrap resamples for M5...\n", R))
for (i in seq_len(R)) {
  d_b <- df[sample(idx, n, replace = TRUE), , drop = FALSE]
  f_b <- tryCatch(
    suppressWarnings(lavaan::sem(m5, data = d_b, std.lv = TRUE, estimator = "ML", missing = "fiml")),
    error = function(e) NULL
  )
  if (!is.null(f_b) && isTRUE(lavaan::lavInspect(f_b, "converged"))) {
    pe_b <- lavaan::parameterEstimates(f_b, standardized = TRUE)
    boot_beta[i, ] <- sapply(labs, function(l) pe_b$std.all[pe_b$label == l][1])
  }
}

pt_beta <- sapply(labs, function(l) pe5$std.all[pe5$op == ":=" & pe5$lhs == l][1])
boot_summ <- data.frame(
  Effect     = colnames(boot_beta),
  Point_beta = as.numeric(pt_beta),
  Boot_Mean  = colMeans(boot_beta, na.rm = TRUE),
  Boot_SD    = apply(boot_beta, 2, sd, na.rm = TRUE),
  CI_LL      = apply(boot_beta, 2, function(x) unname(quantile(x, 0.025, na.rm = TRUE))),
  CI_UL      = apply(boot_beta, 2, function(x) unname(quantile(x, 0.975, na.rm = TRUE))),
  NonConv    = apply(boot_beta, 2, function(x) sum(is.na(x))),
  stringsAsFactors = FALSE
)
cat("\n[M5 5,000-resample percentile bootstrap standardized indirect associations (Table S1)]\n")
print(round_df(boot_summ, 4), row.names = FALSE)
write.csv(round_df(boot_summ, 4), file.path(out_dir, "TableS1_M5_bootstrap_beta.csv"), row.names = FALSE)

sink()
cat("Completed 05_sensitivity_M5_components.R\n")
