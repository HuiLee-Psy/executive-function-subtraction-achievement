# ============================================================================
# SEM-recalc-v3.R
#
# Purpose: Re-run the EF -> Fluency -> Academic latent mediation model on the
# UPDATED data (data_analysis.xlsx) and COMPARE EF indicator configurations to
# decide whether to drop Flanker, SST, and/or AX-CPT.
#
# Key adaptations vs SEM-final-v2.R (old data):
#   1. Input = data_analysis.xlsx, sheet "Sheet1" (old was data-0501.xlsx sheet "data").
#   2. Column rename: Backward_max_span -> Backward_max_span (values identical).
#   3. Sex coding: 1 = Male, 0 = Female (OLD script assumed 0 = Male, 1 = Female).
#   4. SST: stopsignal_ssrt_ms updated for 384/518 rows.
#   5. N = 518 (was 523).
#
# Output: ./results/
# ============================================================================

# 相对路径已配置：请从 V1_result/ 根目录运行本脚本（Rscript code/SEM-recalc-v3.R）

suppressPackageStartupMessages({
  req <- c("lavaan", "semTools", "readxl", "dplyr", "tibble", "tidyr", "psych", "moments")
  for (p in req) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
})

set.seed(20260921)
options(stringsAsFactors = FALSE)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
log_path <- file.path(results_dir, "SEM-recalc_console.log")
sink(log_path, split = TRUE)

cat(strrep("=", 88), "\n", sep = "")
cat("SEM-recalc-v3.R | run at ", as.character(Sys.time()), "\n", sep = "")
cat("R ", R.version.string, " | lavaan ", as.character(packageVersion("lavaan")),
    " | semTools ", as.character(packageVersion("semTools")), "\n", sep = "")
cat(strrep("=", 88), "\n", sep = "")

round_num <- function(x, k = 3) {
  if (is.data.frame(x)) {
    for (nm in names(x)) if (is.numeric(x[[nm]])) x[[nm]] <- round(x[[nm]], k)
    return(x)
  }
  if (is.numeric(x)) return(round(x, k))
  x
}

save_csv <- function(x, filename) {
  fp <- file.path(results_dir, filename)
  write.csv(x, fp, row.names = FALSE, na = "")
  cat("  -> wrote ", fp, "\n", sep = "")
}

safe_fit <- function(model, data, ...) {
  out <- withCallingHandlers(
    tryCatch(
      lavaan::sem(model, data = data, ...),
      error = function(e) { cat("MODEL ERROR: ", conditionMessage(e), "\n", sep = ""); NULL }
    ),
    warning = function(w) { cat("MODEL WARNING: ", conditionMessage(w), "\n", sep = ""); invokeRestart("muffleWarning") }
  )
  out
}

is_ok_fit <- function(fit) !is.null(fit) && isTRUE(lavaan::lavInspect(fit, "converged"))

fit_table_one <- function(fit, model_id, estimator = "ML") {
  if (!is_ok_fit(fit)) {
    return(data.frame(Model = model_id, Estimator = estimator, N = NA, chisq = NA, df = NA,
                      pvalue = NA, CFI = NA, TLI = NA, RMSEA = NA, SRMR = NA, AIC = NA, BIC = NA))
  }
  fm <- lavaan::fitMeasures(fit, c("npar", "chisq", "df", "pvalue", "cfi", "tli",
                                   "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"))
  data.frame(Model = model_id, Estimator = estimator, N = lavaan::nobs(fit),
             chisq = unname(fm["chisq"]), df = unname(fm["df"]), pvalue = unname(fm["pvalue"]),
             CFI = unname(fm["cfi"]), TLI = unname(fm["tli"]), RMSEA = unname(fm["rmsea"]),
             RMSEA_LL = unname(fm["rmsea.ci.lower"]), RMSEA_UL = unname(fm["rmsea.ci.upper"]),
             SRMR = unname(fm["srmr"]), AIC = unname(fm["aic"]), BIC = unname(fm["bic"]))
}

extract_loadings <- function(fit, model_id) {
  if (!is_ok_fit(fit)) return(data.frame())
  ss <- lavaan::standardizedSolution(fit, type = "std.all", ci = TRUE)
  ss <- ss[ss$op == "=~", , drop = FALSE]
  if (!nrow(ss)) return(data.frame())
  data.frame(Model = model_id, Factor = ss$lhs, Indicator = ss$rhs,
             Loading = ss$est.std, SE = ss$se, CI_LL = ss$ci.lower, CI_UL = ss$ci.upper, p = ss$pvalue)
}

calc_rel <- function(fit, model_id) {
  if (!is_ok_fit(fit)) return(data.frame())
  ss <- lavaan::standardizedSolution(fit, type = "std.all")
  ss <- ss[ss$op == "=~", c("lhs", "rhs", "est.std"), drop = FALSE]
  factors <- unique(ss$lhs)
  do.call(rbind, lapply(factors, function(f) {
    lam <- ss$est.std[ss$lhs == f]; lam <- lam[is.finite(lam)]
    k <- length(lam); theta <- 1 - lam^2
    omega <- if (sum(lam) == 0) NA_real_ else (sum(lam)^2) / ((sum(lam)^2) + sum(theta))
    ave <- mean(lam^2)
    data.frame(Model = model_id, Factor = f, Indicators = k, Omega = omega, AVE = ave,
               Min_Loading = min(lam), Max_Loading = max(lam))
  }))
}

extract_structural <- function(fit, model_id, ci_level = .95) {
  if (!is_ok_fit(fit)) return(data.frame())
  pe <- lavaan::parameterEstimates(fit, standardized = TRUE, ci = TRUE, level = ci_level)
  pe <- pe[pe$op %in% c("~", ":="), , drop = FALSE]
  if (!nrow(pe)) return(data.frame())
  data.frame(Model = model_id,
             Effect = ifelse(is.na(pe$label) | pe$label == "", pe$lhs, pe$label),
             lhs = pe$lhs, op = pe$op, rhs = pe$rhs,
             b = pe$est, SE = pe$se, CI_LL = pe$ci.lower, CI_UL = pe$ci.upper,
             p = pe$pvalue, beta = pe$std.all)
}

# ----------------------------------------------------------------------------
# 1. DATA LOAD
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 88), "\n1. DATA LOAD\n", strrep("=", 88), "\n", sep = "")

df <- readxl::read_excel("data/data_analysis.xlsx", sheet = "Sheet1")
cat("N rows = ", nrow(df), "; columns = ", ncol(df), "\n", sep = "")

required <- c("sex", "age", "math_score", "chinese_score", "english_score",
              "axcpt_d_prime_context", "Backward_max_span", "wcst_wcst_completed_categories",
              "flanker_flanker_effect_ms", "stopsignal_ssrt_ms",
              "simple_subtraction_adjusted_score", "complex_subtraction_adjusted_score")
missing_req <- setdiff(required, names(df))
if (length(missing_req)) stop("Missing required variables: ", paste(missing_req, collapse = ", "))

num_vars <- setdiff(required, "sex")
for (v in num_vars) df[[v]] <- suppressWarnings(as.numeric(df[[v]]))

# Sex coding: 1 = Male, 0 = Female (corrected vs old script)
df$sex_group <- factor(ifelse(df$sex %in% c(0, 1), df$sex, NA),
                       levels = c(0, 1), labels = c("Female", "Male"))
cat("Sex coding (raw): 1 = Male, 0 = Female\n")
print(table(df$sex, useNA = "ifany"))
cat("Sex-group N:\n")
print(table(df$sex_group, useNA = "ifany"))

# ----------------------------------------------------------------------------
# 2. VARIABLE PREPARATION
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 88), "\n2. VARIABLE PREPARATION\n", strrep("=", 88), "\n", sep = "")

df$Flanker_Rev <- -df$flanker_flanker_effect_ms
df$SST_Rev     <- -df$stopsignal_ssrt_ms

df$EF_AX       <- as.numeric(scale(df$axcpt_d_prime_context))
df$EF_BACKWARD <- as.numeric(scale(df$Backward_max_span))          # renamed from Backward_max_span
df$EF_WCST     <- as.numeric(scale(df$wcst_wcst_completed_categories))
df$EF_FLANK    <- as.numeric(scale(df$Flanker_Rev))
df$EF_SST      <- as.numeric(scale(df$SST_Rev))
df$FL_SIMP     <- as.numeric(scale(df$simple_subtraction_adjusted_score))
df$FL_COMP     <- as.numeric(scale(df$complex_subtraction_adjusted_score))
df$AC_MATH     <- as.numeric(scale(df$math_score))
df$AC_CHIN     <- as.numeric(scale(df$chinese_score))
df$AC_ENG      <- as.numeric(scale(df$english_score))

all_z <- c("EF_AX", "EF_BACKWARD", "EF_WCST", "EF_FLANK", "EF_SST",
           "FL_SIMP", "FL_COMP", "AC_MATH", "AC_CHIN", "AC_ENG")

# Missingness on the analysis variables
miss <- data.frame(Variable = all_z,
                   N = sapply(df[all_z], function(x) sum(!is.na(x))),
                   Missing = sapply(df[all_z], function(x) sum(is.na(x))),
                   Missing_pct = sapply(df[all_z], function(x) mean(is.na(x)) * 100))
cat("Missingness:\n"); print(round_num(miss, 3), row.names = FALSE)

# ----------------------------------------------------------------------------
# 3. MODEL DEFINITIONS
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 88), "\n3. MODEL DEFINITIONS\n", strrep("=", 88), "\n", sep = "")

base_measurement <- function(ef_indicators) {
  paste0(
    "EF =~ ", paste(ef_indicators, collapse = " + "), "\n",
    "Fluency =~ FL_SIMP + FL_COMP\n",
    "Academic =~ AC_MATH + AC_CHIN + AC_ENG\n"
  )
}
structural_part <- function() {
  paste0(
    "Fluency ~ a*EF\n",
    "Academic ~ b*Fluency + c_prime*EF\n",
    "indirect := a*b\n",
    "direct := c_prime\n",
    "total := c_prime + a*b\n"
  )
}

models <- list(
  M_full       = paste0(base_measurement(c("EF_AX", "EF_BACKWARD", "EF_WCST", "EF_FLANK", "EF_SST")), structural_part()),
  M_noFlanker  = paste0(base_measurement(c("EF_AX", "EF_BACKWARD", "EF_WCST", "EF_SST")),       structural_part()),
  M_noSST      = paste0(base_measurement(c("EF_AX", "EF_BACKWARD", "EF_WCST", "EF_FLANK")),     structural_part()),
  M_noAX       = paste0(base_measurement(c("EF_BACKWARD", "EF_WCST", "EF_FLANK", "EF_SST")),    structural_part()),
  M_noFlankSST = paste0(base_measurement(c("EF_AX", "EF_BACKWARD", "EF_WCST")),                 structural_part()),
  M_noAXFlank  = paste0(base_measurement(c("EF_BACKWARD", "EF_WCST", "EF_SST")),                structural_part())
)

model_desc <- data.frame(
  Model = names(models),
  EF_indicators = c(
    "AX + Corsi + WCST + Flanker + SST",
    "AX + Corsi + WCST + SST",
    "AX + Corsi + WCST + Flanker",
    "Corsi + WCST + Flanker + SST",
    "AX + Corsi + WCST",
    "Corsi + WCST + SST"),
  Change = c("Baseline (all 5)", "Drop Flanker", "Drop SST", "Drop AX-CPT",
             "Drop Flanker + SST", "Drop AX-CPT + Flanker"))
print(model_desc, row.names = FALSE)
save_csv(model_desc, "04_model_definitions.csv")

# ----------------------------------------------------------------------------
# 4. ESTIMATE ALL MODELS
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 88), "\n4. ESTIMATE ALL MODELS (ML + FIML)\n", strrep("=", 88), "\n", sep = "")

fits <- lapply(models, function(m) safe_fit(m, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml"))
names(fits) <- names(models)

fit_tbl <- do.call(rbind, Map(function(f, id) fit_table_one(f, id, "ML+FIML"), fits, names(fits)))
fit_tbl <- round_num(fit_tbl, 4)
cat("Fit comparison:\n"); print(fit_tbl, row.names = FALSE)
save_csv(fit_tbl, "05_fit_comparison.csv")

load_tbl_all <- do.call(rbind, Map(function(f, id) extract_loadings(f, id), fits, names(fits)))
load_tbl_all <- round_num(load_tbl_all, 4)
cat("\nLoadings:\n"); print(load_tbl_all, row.names = FALSE)
save_csv(load_tbl_all, "06_loadings.csv")

rel_tbl_all <- do.call(rbind, Map(function(f, id) calc_rel(f, id), fits, names(fits)))
rel_tbl_all <- round_num(rel_tbl_all, 4)
cat("\nReliability / AVE:\n"); print(rel_tbl_all, row.names = FALSE)
save_csv(rel_tbl_all, "07_reliability_AVE.csv")

struct_tbl_all <- do.call(rbind, Map(function(f, id) extract_structural(f, id), fits, names(fits)))
struct_tbl_all <- round_num(struct_tbl_all, 4)
cat("\nStructural effects:\n"); print(struct_tbl_all, row.names = FALSE)
save_csv(struct_tbl_all, "08_structural_effects.csv")

# ----------------------------------------------------------------------------
# 5. SIDE-BY-SIDE SUMMARY (key rows only)
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 88), "\n5. SIDE-BY-SIDE SUMMARY\n", strrep("=", 88), "\n", sep = "")

summary_rows <- lapply(names(fits), function(id) {
  ft <- fit_tbl[fit_tbl$Model == id, ]
  rel_ef <- rel_tbl_all[rel_tbl_all$Model == id & rel_tbl_all$Factor == "EF", , drop = FALSE]
  st <- struct_tbl_all[struct_tbl_all$Model == id, , drop = FALSE]
  ind_b <- st$b[st$Effect == "indirect"]; ind_beta <- st$beta[st$Effect == "indirect"]
  dir_b <- st$b[st$Effect == "direct"]; dir_beta <- st$beta[st$Effect == "direct"]
  data.frame(Model = id, EF_k = if (nrow(rel_ef)) rel_ef$Indicators[1] else NA,
             chisq = ft$chisq, df = ft$df, p = ft$pvalue,
             CFI = ft$CFI, TLI = ft$TLI, RMSEA = ft$RMSEA, SRMR = ft$SRMR,
             EF_Omega = if (nrow(rel_ef)) rel_ef$Omega[1] else NA,
             EF_AVE = if (nrow(rel_ef)) rel_ef$AVE[1] else NA,
             a = st$b[st$Effect == "a"][1], b = st$b[st$Effect == "b"][1],
             c_prime = st$b[st$Effect == "c_prime"][1],
             indirect_b = ind_b[1], indirect_beta = ind_beta[1],
             direct_beta = dir_beta[1])
})
summary_tbl <- do.call(rbind, summary_rows)
summary_tbl <- round_num(summary_tbl, 4)
cat("Side-by-side summary:\n"); print(summary_tbl, row.names = FALSE)
save_csv(summary_tbl, "09_side_by_side_summary.csv")

# ----------------------------------------------------------------------------
# 6. DECISION SCREEN: which indicator(s) to drop?
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 88), "\n6. DROP-DECISION SCREEN\n", strrep("=", 88), "\n", sep = "")

# Full-model standardized loadings of the three candidate indicators
full_ss <- lavaan::standardizedSolution(fits$M_full, type = "std.all")
full_ef_load <- full_ss[full_ss$op == "=~" & full_ss$lhs == "EF", c("rhs", "est.std")]
rownames(full_ef_load) <- NULL
cat("Full-model (M_full) EF loadings:\n")
print(round_num(full_ef_load, 4), row.names = FALSE)

cand <- c("Flanker" = "M_noFlanker", "SST" = "M_noSST", "AXCPT" = "M_noAX")
decision <- do.call(rbind, lapply(names(cand), function(nm) {
  full_id <- "M_full"; drop_id <- cand[[nm]]
  f0 <- fit_tbl[fit_tbl$Model == full_id, ]; f1 <- fit_tbl[fit_tbl$Model == drop_id, ]
  st0 <- struct_tbl_all[struct_tbl_all$Model == full_id, ]; st1 <- struct_tbl_all[struct_tbl_all$Model == drop_id, ]
  load_in_full <- full_ef_load$est.std[full_ef_load$rhs == switch(nm,
    Flanker = "EF_FLANK", SST = "EF_SST", AXCPT = "EF_AX")]
  data.frame(Candidate = nm,
             Loading_in_full = load_in_full,
             dCFI = f1$CFI - f0$CFI,
             dRMSEA = f1$RMSEA - f0$RMSEA,
             dSRMR = f1$SRMR - f0$SRMR,
             dEF_Omega = rel_tbl_all$Omega[rel_tbl_all$Model == drop_id & rel_tbl_all$Factor == "EF"][1] -
                         rel_tbl_all$Omega[rel_tbl_all$Model == full_id & rel_tbl_all$Factor == "EF"][1],
             d_indirect_beta = st1$beta[st1$Effect == "indirect"][1] - st0$beta[st0$Effect == "indirect"][1])
}))
decision <- round_num(decision, 4)
cat("Drop-decision screen (relative to M_full):\n")
print(decision, row.names = FALSE)
save_csv(decision, "10_drop_decision_screen.csv")

cat("\nDecision screen interpretation:\n")
cat(" - A candidate indicator is 'weak' if its full-model loading < .30 (or .40).\n")
cat(" - Dropping a weak indicator is acceptable only if fit does not meaningfully\n")
cat("   deteriorate (dCFI >= -.01, dRMSEA <= +.015) and the focal indirect effect\n")
cat("   is stable (|d_indirect_beta| small).\n")

sink()
cat("\nDONE. Results in ", results_dir, "/\n", sep = "")
