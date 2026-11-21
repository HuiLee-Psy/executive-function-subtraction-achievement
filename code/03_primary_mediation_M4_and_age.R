# ============================================================================
# code/03_primary_mediation_M4_and_age.R
# Reproduces:
#   - Section 3.3 & Table 4: Primary hybrid SEM (M4) global fit, R², and
#     5,000-resample nonparametric percentile bootstrap decomposition
#   - Section 3.6 & Supplementary Table S1: Age-adjusted M4 sensitivity model
#     (standardized student age predicting both Fluency and Academic Achievement)
# ============================================================================

suppressPackageStartupMessages({
  library(lavaan)
  library(readxl)
})
source("code/00_load_and_prep.R")

out_dir <- "results/03_primary_mediation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "PRIMARY_M4_console.log"), split = TRUE)

cat("03_primary_mediation_M4_and_age.R | ", as.character(Sys.time()),
    " | lavaan ", as.character(packageVersion("lavaan")), "\n", sep = "")
df <- load_analysis_data()

# ---- 1. Primary Model M4 ----
m4 <- "
  Fluency  =~ FL_SIMP + FL_COMP
  Academic =~ AC_MATH + AC_CHIN + AC_ENG
  Fluency  ~ a*EF_COMP
  Academic ~ b*Fluency + c_prime*EF_COMP
  ind   := a * b
  total := c_prime + a * b
  prop  := ind / total
"

f4 <- lavaan::sem(m4, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
fm4 <- lavaan::fitMeasures(f4, c("chisq", "df", "pvalue", "cfi", "tli",
                                 "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr"))
cat(sprintf(
  "\n[M4 Global Fit] chi2(%.0f) = %.2f, p = %.4f, CFI = %.4f, TLI = %.4f, RMSEA = %.4f [%.4f, %.4f], SRMR = %.4f\n",
  fm4["df"], fm4["chisq"], fm4["pvalue"], fm4["cfi"], fm4["tli"],
  fm4["rmsea"], fm4["rmsea.ci.lower"], fm4["rmsea.ci.upper"], fm4["srmr"]
))

r2_m4 <- lavaan::lavInspect(f4, "rsquare")[c("Fluency", "Academic")]
cat(sprintf("[M4 R2] Fluency = %.4f | Academic Achievement = %.4f\n", r2_m4["Fluency"], r2_m4["Academic"]))

pe0 <- lavaan::parameterEstimates(f4, standardized = TRUE)
g0  <- function(l) pe0$est[pe0$label == l][1]
gb0 <- function(l) pe0$std.all[pe0$label == l][1]

# 5,000-resample participant-level nonparametric percentile bootstrap (seed = 20260925)
set.seed(20260925)
n <- nrow(df); idx <- seq_len(n); R <- 5000L
cols <- c("a_b", "b_b", "cp_b", "ind_b", "tot_b",
          "a_beta", "b_beta", "cp_beta", "ind_beta", "tot_beta", "prop_beta")
boot_mat <- matrix(NA_real_, R, length(cols), dimnames = list(NULL, cols))

cat(sprintf("\nRunning %d bootstrap resamples for M4...\n", R))
for (i in seq_len(R)) {
  d_boot <- df[sample(idx, n, replace = TRUE), , drop = FALSE]
  f_boot <- tryCatch(
    suppressWarnings(lavaan::sem(m4, data = d_boot, std.lv = TRUE, estimator = "ML", missing = "fiml")),
    error = function(e) NULL
  )
  if (!is.null(f_boot) && isTRUE(lavaan::lavInspect(f_boot, "converged"))) {
    pe <- lavaan::parameterEstimates(f_boot, standardized = TRUE)
    g  <- function(l) pe$est[pe$label == l][1]
    gb <- function(l) pe$std.all[pe$label == l][1]
    boot_mat[i, ] <- c(
      g("a"), g("b"), g("c_prime"), g("ind"), g("total"),
      gb("a"), gb("b"), gb("c_prime"), gb("ind"), gb("total"),
      gb("ind") / gb("total")
    )
  }
}

ci_pct <- function(x) unname(quantile(x, c(0.025, 0.975), na.rm = TRUE))
point <- c(
  g0("a"), g0("b"), g0("c_prime"), g0("ind"), g0("total"),
  gb0("a"), gb0("b"), gb0("c_prime"), gb0("ind"), gb0("total"),
  gb0("ind") / gb0("total")
)

t4_summary <- data.frame(
  Effect    = cols,
  Point     = point,
  Boot_Mean = colMeans(boot_mat, na.rm = TRUE),
  Boot_SD   = apply(boot_mat, 2, sd, na.rm = TRUE),
  CI_LL     = apply(boot_mat, 2, function(x) ci_pct(x)[1]),
  CI_UL     = apply(boot_mat, 2, function(x) ci_pct(x)[2]),
  NonConv   = apply(boot_mat, 2, function(x) sum(is.na(x))),
  stringsAsFactors = FALSE
)
t4_summary <- round_df(t4_summary, 4)

cat("\n[Table 4: Primary M4 structural estimates & 5,000-resample percentile bootstrap]\n")
print(t4_summary, row.names = FALSE)
write.csv(t4_summary, file.path(out_dir, "Table4_M4_bootstrap_summary.csv"), row.names = FALSE)

# ---- 2. Age-Adjusted M4 Sensitivity Model (Section 3.6 & Table S1) ----
m4_age <- "
  Fluency  =~ FL_SIMP + FL_COMP
  Academic =~ AC_MATH + AC_CHIN + AC_ENG
  Fluency  ~ a*EF_COMP + d1*AGE_Z
  Academic ~ b*Fluency + c_prime*EF_COMP + d2*AGE_Z
  ind   := a * b
  total := c_prime + a * b
"
f4_age <- lavaan::sem(m4_age, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
ss_age <- lavaan::standardizedSolution(f4_age, ci = TRUE, level = 0.95)
age_keep <- ss_age[ss_age$label %in% c("a", "b", "c_prime", "d1", "d2", "ind", "total"),
                   c("label", "est.std", "se", "ci.lower", "ci.upper", "pvalue")]
cat("\n[Sensitivity: Age-adjusted M4 standardized delta-method estimates (Table S1)]\n")
print(round_df(age_keep, 4), row.names = FALSE)
write.csv(round_df(age_keep, 4), file.path(out_dir, "TableS1_M4_age_adjusted.csv"), row.names = FALSE)

sink()
cat("Completed 03_primary_mediation_M4_and_age.R\n")
