# ============================================================================
# code/02_measurement_and_latent_sensitivity.R
# Reproduces:
#   - Section 3.2 & Table 3: Single-factor EF CFA (M0) and 3-factor measurement
#     model loadings, model-based composite reliability (omega), and AVE
#   - Section 3.2.1: Multidimensional latent EF admissibility checks
#     (M2 orthogonal bifactor, M3 two-factor oblique ESEM)
#   - Section 3.6 & Supplementary Table S1: Latent-factor mediation sensitivity
#     models M0 (5-indicator EF) and M1 (4-indicator EF, SSRT omitted) with
#     standardized delta-method 95% CIs
# ============================================================================

suppressPackageStartupMessages({
  library(lavaan)
  library(semTools)
  library(readxl)
})
source("code/00_load_and_prep.R")

out_dir <- "results/02_measurement"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "MEASUREMENT_console.log"), split = TRUE)

cat("02_measurement_and_latent_sensitivity.R | ", as.character(Sys.time()),
    " | lavaan ", as.character(packageVersion("lavaan")), "\n", sep = "")
df <- load_analysis_data()

safe_fit <- function(model, data, ...) {
  withCallingHandlers(
    tryCatch(lavaan::sem(model, data = data, ...),
             error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }),
    warning = function(w) { cat("WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") }
  )
}

# Model-based McDonald's omega and AVE from standardized factor loadings
calc_omega_ave <- function(fit, factor_name) {
  ss <- lavaan::standardizedSolution(fit)
  lam <- ss$est.std[ss$op == "=~" & ss$lhs == factor_name]
  omega <- (sum(lam)^2) / (sum(lam)^2 + sum(1 - lam^2))
  ave   <- mean(lam^2)
  data.frame(
    Factor = factor_name, Indicators = length(lam),
    Omega = omega, AVE = ave,
    Min_Loading = min(lam), Max_Loading = max(lam),
    stringsAsFactors = FALSE
  )
}

fit_summary_row <- function(fit, name) {
  fm <- lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli",
                                   "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr"))
  data.frame(
    Model = name,
    chisq = unname(fm["chisq"]), df = unname(fm["df"]), p = unname(fm["pvalue"]),
    CFI = unname(fm["cfi"]), TLI = unname(fm["tli"]),
    RMSEA = unname(fm["rmsea"]), RMSEA_LL = unname(fm["rmsea.ci.lower"]),
    RMSEA_UL = unname(fm["rmsea.ci.upper"]), SRMR = unname(fm["srmr"]),
    stringsAsFactors = FALSE
  )
}

# ---- 1. M0: Single-factor EF CFA (5 indicators) ----
m0_cfa <- "EF =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST"
f_m0_cfa <- safe_fit(m0_cfa, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")

# ---- 2. Three-factor measurement model (Table 3) ----
m_3f <- "
  EF       =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST
  Fluency  =~ FL_SIMP + FL_COMP
  Academic =~ AC_MATH + AC_CHIN + AC_ENG
"
f_3f <- safe_fit(m_3f, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")

ss_m0 <- lavaan::standardizedSolution(f_m0_cfa)
ss_m0 <- ss_m0[ss_m0$op == "=~", c("lhs", "rhs", "est.std", "se", "z", "pvalue", "ci.lower", "ci.upper")]
ss_3f <- lavaan::standardizedSolution(f_3f)
ss_3f <- ss_3f[ss_3f$op == "=~", c("lhs", "rhs", "est.std", "se", "z", "pvalue", "ci.lower", "ci.upper")]

cat("\n[Table 3: Single-factor EF CFA (M0) standardized loadings]\n")
print(round_df(ss_m0, 4), row.names = FALSE)
cat("\n[Table 3: Three-factor measurement model standardized loadings]\n")
print(round_df(ss_3f, 4), row.names = FALSE)

rel_tbl <- rbind(
  cbind(Model = "M0_CFA", calc_omega_ave(f_m0_cfa, "EF")),
  cbind(Model = "3Factor", calc_omega_ave(f_3f, "EF")),
  cbind(Model = "3Factor", calc_omega_ave(f_3f, "Fluency")),
  cbind(Model = "3Factor", calc_omega_ave(f_3f, "Academic"))
)
cat("\n[Reliability (omega) and AVE]\n")
print(round_df(rel_tbl, 4), row.names = FALSE)

write.csv(round_df(ss_m0, 4), file.path(out_dir, "Table3_M0_loadings.csv"), row.names = FALSE)
write.csv(round_df(ss_3f, 4), file.path(out_dir, "Table3_3factor_loadings.csv"), row.names = FALSE)
write.csv(round_df(rel_tbl, 4), file.path(out_dir, "Reliability_AVE.csv"), row.names = FALSE)

# ---- 3. M2 (Orthogonal Bifactor) & M3 (2-Factor Oblique ESEM) Admissibility Checks ----
cat("\n[M2: Orthogonal bifactor EF model admissibility check]\n")
m2_bifactor <- "
  G   =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST
  Inh =~ EF_FLANK + EF_SST
  G ~~ 0*Inh
"
f_m2 <- safe_fit(m2_bifactor, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
if (!is.null(f_m2)) {
  pe_m2 <- lavaan::standardizedSolution(f_m2)
  print(round_df(pe_m2[pe_m2$op == "=~", c("lhs", "rhs", "est.std", "se", "pvalue")], 4), row.names = FALSE)
}

cat("\n[M3: Two-factor oblique ESEM admissibility check]\n")
m3_esem <- '
  efa("efa1")*F1 + efa("efa1")*F2 =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST
'
f_m3 <- safe_fit(m3_esem, data = df, estimator = "ML", missing = "fiml", rotation = "geomin")
if (!is.null(f_m3)) {
  pe_m3 <- lavaan::parameterEstimates(f_m3, standardized = TRUE)
  res_vars <- pe_m3[pe_m3$op == "~~" & pe_m3$lhs == pe_m3$rhs & pe_m3$lhs %in% c("EF_AX","EF_BACKWARD","EF_WCST","EF_FLANK","EF_SST"),
                    c("lhs", "est", "std.all")]
  cat("Residual variances (negative value indicates Heywood case):\n")
  print(round_df(res_vars, 4), row.names = FALSE)
}

# ---- 4. Latent-EF Mediation Sensitivity Models: M0 (5-ind) and M1 (4-ind, omit SSRT) ----
m0_med <- "
  EF       =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST
  Fluency  =~ FL_SIMP + FL_COMP
  Academic =~ AC_MATH + AC_CHIN + AC_ENG
  Fluency  ~ a*EF
  Academic ~ b*Fluency + c_prime*EF
  indirect := a * b
  total    := c_prime + a * b
"
m1_med <- "
  EF       =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK
  Fluency  =~ FL_SIMP + FL_COMP
  Academic =~ AC_MATH + AC_CHIN + AC_ENG
  Fluency  ~ a*EF
  Academic ~ b*Fluency + c_prime*EF
  indirect := a * b
  total    := c_prime + a * b
"

f_m0_med <- safe_fit(m0_med, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
f_m1_med <- safe_fit(m1_med, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")

extract_std_effects <- function(fit, tag) {
  ss <- lavaan::standardizedSolution(fit, ci = TRUE, level = 0.95)
  keep <- ss[ss$label %in% c("a", "b", "c_prime", "indirect", "total"),
             c("label", "est.std", "se", "ci.lower", "ci.upper", "pvalue")]
  cbind(Model = tag, keep)
}

latent_sens <- rbind(
  extract_std_effects(f_m0_med, "M0_latent_5ind"),
  extract_std_effects(f_m1_med, "M1_latent_4ind_noSSRT")
)
cat("\n[Latent EF sensitivity models: standardized delta-method estimates (Table S1)]\n")
print(round_df(latent_sens, 4), row.names = FALSE)
write.csv(round_df(latent_sens, 4), file.path(out_dir, "TableS1_latent_EF_sensitivity.csv"), row.names = FALSE)

cat("\n[M0 latent mediation R2]\n")
print(round(lavaan::lavInspect(f_m0_med, "rsquare")[c("Fluency", "Academic")], 4))
cat("\n[M1 latent mediation R2]\n")
print(round(lavaan::lavInspect(f_m1_med, "rsquare")[c("Fluency", "Academic")], 4))

fits_all <- rbind(
  fit_summary_row(f_m0_cfa, "M0_5ind_CFA"),
  fit_summary_row(f_3f,     "Three_factor_CFA"),
  fit_summary_row(f_m0_med, "M0_latent_mediation"),
  fit_summary_row(f_m1_med, "M1_latent_mediation_noSSRT")
)
cat("\n[Global fit summary]\n")
print(round_df(fits_all, 4), row.names = FALSE)
write.csv(round_df(fits_all, 4), file.path(out_dir, "Measurement_and_Latent_Fit.csv"), row.names = FALSE)

sink()
cat("Completed 02_measurement_and_latent_sensitivity.R\n")
