# ============================================================================
# code/06_sex_invariance.R
# Reproduces Section 3.5, Supplementary Table S2 (multigroup measurement and
# structural invariance sequence), and Supplementary Table S3 (sex-stratified
# parameter estimates and full-covariance Wald contrasts under partial scalar
# invariance).
#
# CRITICAL FIXES vs. v1.0.0:
#   1. Explicitly sets group.label = c("Female", "Male") in lavaan::sem() so
#      Group 1 = Girls (n = 238) and Group 2 = Boys (n = 280), preventing
#      lavaan's default ordering by first row appearance (Male in row 1).
#   2. Extracts group-specific structural slopes, indirect associations, R²,
#      and full-covariance Wald contrasts (Girls - Boys) directly from the
#      partial scalar model (Model 3b) rather than the configural model.
# ============================================================================

suppressPackageStartupMessages({
  library(lavaan)
  library(readxl)
})
source("code/00_load_and_prep.R")

out_dir <- "results/06_sex_invariance"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "SEX_INVARIANCE_console.log"), split = TRUE)

cat("06_sex_invariance.R | ", as.character(Sys.time()),
    " | lavaan ", as.character(packageVersion("lavaan")), "\n", sep = "")
df <- load_analysis_data()
cat("Group counts:\n")
print(table(df$sex_group))

GROUP_ORDER <- c("Female", "Male")

fit_row <- function(f, step_label) {
  fm <- lavaan::fitMeasures(f, c("chisq", "df", "cfi", "tli", "rmsea", "srmr"))
  data.frame(
    Model = step_label,
    chisq = unname(fm["chisq"]),
    df    = unname(fm["df"]),
    CFI   = unname(fm["cfi"]),
    TLI   = unname(fm["tli"]),
    RMSEA = unname(fm["rmsea"]),
    SRMR  = unname(fm["srmr"]),
    stringsAsFactors = FALSE
  )
}

# ---- 1. Hierarchical Multigroup Invariance Sequence (Supplementary Table S2) ----
m_base <- "
  Fluency  =~ FL_SIMP + FL_COMP
  Academic =~ AC_MATH + AC_CHIN + AC_ENG
  Fluency  ~ EF_COMP
  Academic ~ Fluency + EF_COMP
"

f_config <- lavaan::sem(
  m_base, data = df, group = "sex_group", group.label = GROUP_ORDER,
  std.lv = TRUE, estimator = "ML", missing = "fiml"
)
f_metric <- lavaan::sem(
  m_base, data = df, group = "sex_group", group.label = GROUP_ORDER,
  std.lv = TRUE, estimator = "ML", missing = "fiml",
  group.equal = c("loadings")
)
f_scalar <- lavaan::sem(
  m_base, data = df, group = "sex_group", group.label = GROUP_ORDER,
  std.lv = TRUE, estimator = "ML", missing = "fiml",
  group.equal = c("loadings", "intercepts")
)
f_pscalar <- lavaan::sem(
  m_base, data = df, group = "sex_group", group.label = GROUP_ORDER,
  std.lv = TRUE, estimator = "ML", missing = "fiml",
  group.equal = c("loadings", "intercepts"),
  group.partial = c("AC_CHIN ~ 1", "AC_ENG ~ 1")
)
f_struct <- lavaan::sem(
  m_base, data = df, group = "sex_group", group.label = GROUP_ORDER,
  std.lv = TRUE, estimator = "ML", missing = "fiml",
  group.equal = c("loadings", "intercepts", "regressions"),
  group.partial = c("AC_CHIN ~ 1", "AC_ENG ~ 1")
)

# Verify group order explicitly
stopifnot(identical(lavaan::lavInspect(f_pscalar, "group.label"), GROUP_ORDER))
cat("\nVerified lavaan group order:", paste(lavaan::lavInspect(f_pscalar, "group.label"), collapse = " -> "), "\n")

# Build Table S2 with exact baseline comparisons:
#   Metric vs. Configural; Scalar vs. Metric; Partial scalar vs. Metric; Structural vs. Partial scalar
t_s2 <- rbind(
  fit_row(f_config,  "1. Configural"),
  fit_row(f_metric,  "2. Metric (equal loadings)"),
  fit_row(f_scalar,  "3a. Full scalar"),
  fit_row(f_pscalar, "3b. Partial scalar"),
  fit_row(f_struct,  "4. Equal structural slopes")
)

ref_idx <- c(NA, 1L, 2L, 2L, 4L)
t_s2$dCFI   <- ifelse(is.na(ref_idx), NA_real_, t_s2$CFI - t_s2$CFI[ref_idx])
t_s2$dchisq <- ifelse(is.na(ref_idx), NA_real_, t_s2$chisq - t_s2$chisq[ref_idx])
t_s2$ddf    <- ifelse(is.na(ref_idx), NA_real_, t_s2$df - t_s2$df[ref_idx])
t_s2$p_diff <- ifelse(is.na(ref_idx), NA_real_, pchisq(t_s2$dchisq, t_s2$ddf, lower.tail = FALSE))

cat("\n[Supplementary Table S2: Multigroup invariance sequence]\n")
print(round_df(t_s2, 4), row.names = FALSE)
write.csv(round_df(t_s2, 4), file.path(out_dir, "TableS2_sex_invariance_fit.csv"), row.names = FALSE)

# ---- 2. Partial Scalar Group Estimates & Full-Covariance Wald Contrasts (Table S3) ----
m_ps_labeled <- "
  Fluency  =~ FL_SIMP + FL_COMP
  Academic =~ AC_MATH + AC_CHIN + AC_ENG
  Fluency  ~ c(a_g, a_b)*EF_COMP
  Academic ~ c(b_g, b_b)*Fluency + c(c_g, c_b)*EF_COMP

  ind_g    := a_g * b_g
  ind_b    := a_b * b_b
  diff_a   := a_g - a_b
  diff_b   := b_g - b_b
  diff_c   := c_g - c_b
  diff_ind := ind_g - ind_b
"

f_ps_lab <- lavaan::sem(
  m_ps_labeled, data = df, group = "sex_group", group.label = GROUP_ORDER,
  std.lv = TRUE, estimator = "ML", missing = "fiml",
  group.equal = c("loadings", "intercepts"),
  group.partial = c("AC_CHIN ~ 1", "AC_ENG ~ 1")
)

# Confirm identical log-likelihood / chi-square to unlabeled partial scalar model
stopifnot(abs(lavaan::fitMeasures(f_ps_lab, "chisq") - lavaan::fitMeasures(f_pscalar, "chisq")) < 1e-6)

pe_unstd <- lavaan::parameterEstimates(f_ps_lab, ci = TRUE, level = 0.95)
pe_std   <- lavaan::standardizedSolution(f_ps_lab, ci = TRUE, level = 0.95)
r2_grp   <- lavaan::lavInspect(f_ps_lab, "rsquare") # named list: Female, Male

pick_u <- function(lbl) pe_unstd[pe_unstd$label == lbl, ][1, ]
pick_s <- function(lbl) pe_std[pe_std$label == lbl, ][1, ]

t_s3 <- data.frame(
  Parameter = c(
    "Path a (EF -> fluency), b",
    "Path b (Fluency -> achievement), b",
    "Path c' (Direct EF -> achievement), b",
    "Indirect association (a x b), b",
    "Standardized indirect association, beta",
    "Academic achievement R2"
  ),
  Girls_est = c(
    pick_u("a_g")$est, pick_u("b_g")$est, pick_u("c_g")$est, pick_u("ind_g")$est,
    pick_s("ind_g")$est.std, unname(r2_grp[["Female"]]["Academic"])
  ),
  Girls_CI_LL = c(NA, NA, NA, NA, pick_s("ind_g")$ci.lower, NA),
  Girls_CI_UL = c(NA, NA, NA, NA, pick_s("ind_g")$ci.upper, NA),
  Boys_est = c(
    pick_u("a_b")$est, pick_u("b_b")$est, pick_u("c_b")$est, pick_u("ind_b")$est,
    pick_s("ind_b")$est.std, unname(r2_grp[["Male"]]["Academic"])
  ),
  Boys_CI_LL = c(NA, NA, NA, NA, pick_s("ind_b")$ci.lower, NA),
  Boys_CI_UL = c(NA, NA, NA, NA, pick_s("ind_b")$ci.upper, NA),
  Diff_Girls_minus_Boys = c(
    pick_u("diff_a")$est, pick_u("diff_b")$est, pick_u("diff_c")$est, pick_u("diff_ind")$est,
    NA_real_, NA_real_
  ),
  Diff_CI_LL = c(
    pick_u("diff_a")$ci.lower, pick_u("diff_b")$ci.lower, pick_u("diff_c")$ci.lower, pick_u("diff_ind")$ci.lower,
    NA_real_, NA_real_
  ),
  Diff_CI_UL = c(
    pick_u("diff_a")$ci.upper, pick_u("diff_b")$ci.upper, pick_u("diff_c")$ci.upper, pick_u("diff_ind")$ci.upper,
    NA_real_, NA_real_
  ),
  p = c(
    pick_u("diff_a")$pvalue, pick_u("diff_b")$pvalue, pick_u("diff_c")$pvalue, pick_u("diff_ind")$pvalue,
    NA_real_, NA_real_
  ),
  stringsAsFactors = FALSE
)

cat("\n[Supplementary Table S3: Sex-stratified estimates under partial scalar invariance]\n")
print(round_df(t_s3, 4), row.names = FALSE)
write.csv(round_df(t_s3, 4), file.path(out_dir, "TableS3_sex_partial_scalar_paths.csv"), row.names = FALSE)

sink()
cat("Completed 06_sex_invariance.R\n")
