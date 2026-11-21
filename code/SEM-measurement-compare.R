# ============================================================================
# SEM-measurement-compare.R
# EF 测量结构比较：M0 五指标单因子 / M1 删SST / M2 General+Inhibition 双因子 /
# M3 ESEM / M4 composite。含中介模型（delta 法 CI）。
# 输出：results/
# ============================================================================

# 相对路径已配置：请从 V1_result/ 根目录运行本脚本（Rscript code/SEM-measurement-compare.R）

suppressPackageStartupMessages({
  req <- c("lavaan","semTools","readxl","dplyr","tibble","tidyr","psych","moments","GPArotation")
  for (p in req) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
})

set.seed(20260922)
options(stringsAsFactors = FALSE)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
log_path <- file.path(results_dir, "SEM-measurement-compare_console.log")
sink(log_path, split = TRUE)

cat(strrep("=", 90), "\n", sep = "")
cat("SEM-measurement-compare.R | ", as.character(Sys.time()), "\n", sep = "")
cat("R ", R.version.string, " | lavaan ", as.character(packageVersion("lavaan")), "\n", sep = "")
cat(strrep("=", 90), "\n", sep = "")

round_num <- function(x, k = 3) {
  if (is.data.frame(x)) { for (nm in names(x)) if (is.numeric(x[[nm]])) x[[nm]] <- round(x[[nm]], k); return(x) }
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
    tryCatch(lavaan::sem(model, data = data, ...),
             error = function(e) { cat("MODEL ERROR: ", conditionMessage(e), "\n", sep = ""); NULL }),
    warning = function(w) { cat("MODEL WARNING: ", conditionMessage(w), "\n", sep = ""); invokeRestart("muffleWarning") }
  )
  out
}
is_ok_fit <- function(fit) !is.null(fit) && isTRUE(lavaan::lavInspect(fit, "converged"))

fit_table_one <- function(fit, model_id, estimator = "ML") {
  if (!is_ok_fit(fit)) return(data.frame(Model = model_id, Estimator = estimator, N = NA,
    chisq = NA, df = NA, pvalue = NA, CFI = NA, TLI = NA, RMSEA = NA,
    RMSEA_LL = NA, RMSEA_UL = NA, SRMR = NA))
  fm <- lavaan::fitMeasures(fit, c("chisq","df","pvalue","cfi","tli","rmsea",
    "rmsea.ci.lower","rmsea.ci.upper","srmr"))
  data.frame(Model = model_id, Estimator = estimator, N = lavaan::nobs(fit),
    chisq = unname(fm["chisq"]), df = unname(fm["df"]), pvalue = unname(fm["pvalue"]),
    CFI = unname(fm["cfi"]), TLI = unname(fm["tli"]), RMSEA = unname(fm["rmsea"]),
    RMSEA_LL = unname(fm["rmsea.ci.lower"]), RMSEA_UL = unname(fm["rmsea.ci.upper"]),
    SRMR = unname(fm["srmr"]))
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
  ss <- ss[ss$op == "=~", c("lhs","rhs","est.std"), drop = FALSE]
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
# 1. DATA
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 90), "\n1. DATA\n", strrep("=", 90), "\n", sep = "")
df <- readxl::read_excel("data/data_analysis.xlsx", sheet = "Sheet1")
cat("N rows = ", nrow(df), "\n", sep = "")

df$Flanker_Rev <- -df$flanker_flanker_effect_ms
df$SST_Rev     <- -df$stopsignal_ssrt_ms
df$EF_AX       <- as.numeric(scale(df$axcpt_d_prime_context))
df$EF_BACKWARD <- as.numeric(scale(df$Backward_max_span))
df$EF_WCST     <- as.numeric(scale(df$wcst_wcst_completed_categories))
df$EF_FLANK    <- as.numeric(scale(df$Flanker_Rev))
df$EF_SST      <- as.numeric(scale(df$SST_Rev))
df$FL_SIMP     <- as.numeric(scale(df$simple_subtraction_adjusted_score))
df$FL_COMP     <- as.numeric(scale(df$complex_subtraction_adjusted_score))
df$AC_MATH     <- as.numeric(scale(df$math_score))
df$AC_CHIN     <- as.numeric(scale(df$chinese_score))
df$AC_ENG      <- as.numeric(scale(df$english_score))

# --- 流畅性测量说明（第4节要求）---
cat("\n[Fluency measurement note]\n")
cat("Using ONLY simple/complex ADJUSTED scores as the 2 fluency indicators.\n")
cat("correct/incorrect/total are derived (adjusted = correct - incorrect; total = correct + incorrect) and NOT used as separate indicators.\n")
cat("FL_SIMP ~ FL_COMP correlation = ",
    round(cor(df$FL_SIMP, df$FL_COMP, use = "pairwise.complete.obs"), 3),
    " (Spearman-Brown = ",
    round(2*cor(df$FL_SIMP, df$FL_COMP, use = "pairwise.complete.obs")/
          (1+cor(df$FL_SIMP, df$FL_COMP, use = "pairwise.complete.obs")), 3), ")\n", sep = "")

# --- EF 指标相关矩阵 + 特征值（回答"为何 EF 任务相关低"）---
cat("\n[EF indicator correlation matrix (raw, higher=better for all)]\n")
ef_raw <- c("axcpt_d_prime_context","Backward_max_span","wcst_wcst_completed_categories","Flanker_Rev","SST_Rev")
names(ef_raw) <- c("AX-CPT d'","Corsi span","WCST cat","Flanker(rev)","SST SSRT(rev)")
cor_ef <- cor(df[ef_raw], use = "pairwise.complete.obs")
colnames(cor_ef) <- rownames(cor_ef) <- names(ef_raw)
print(round(cor_ef, 3))
save_csv(round_num(as.data.frame(cor_ef) |> tibble::rownames_to_column("Indicator"), 3), "E0_EF_correlations.csv")

eig <- eigen(cor_ef)$values
cat("\nEigenvalues of EF correlation matrix:\n")
print(round(eig, 3))
cat("(Eigenvalues > 1 suggests more than one dimension; note first eigenvalue vs rest.)\n")

# --- psych::fa 交叉验证（2因子 ML oblimin）---
fa2 <- psych::fa(cor_ef, nfactors = 2, n.obs = nrow(df), fm = "ml", rotate = "oblimin")
cat("\n[psych::fa cross-check: 2-factor ML oblimin]\n")
print(round(loadings(fa2)[, 1:2], 3))
cat("Factor correlation:\n"); print(round(fa2$Phi, 3))
cat("Communalities:\n"); print(round(fa2$communality, 3))

# ----------------------------------------------------------------------------
# 2. M0: 五指标单因子 CFA
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 90), "\n2. M0: 5-INDICATOR SINGLE-FACTOR EF (baseline)\n", strrep("=", 90), "\n", sep = "")

m0 <- 'EF =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST'
fit0 <- safe_fit(m0, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
cat("M0 fit:\n"); print(fit_table_one(fit0, "M0"), row.names = FALSE)
load0 <- extract_loadings(fit0, "M0")
rel0  <- calc_rel(fit0, "M0")
cat("M0 loadings (std.all + 95% CI):\n"); print(load0, row.names = FALSE)
cat("M0 reliability/AVE:\n"); print(rel0, row.names = FALSE)
# 指标级 R² = lambda^2
load0$R2_indicator <- load0$Loading^2
cat("M0 indicator R² (= loading^2):\n")
print(round_num(load0[, c("Factor","Indicator","Loading","R2_indicator")], 3), row.names = FALSE)
save_csv(round_num(load0, 4), "M0_loadings.csv")
save_csv(round_num(rel0, 4), "M0_reliability.csv")

# ----------------------------------------------------------------------------
# 3. M1: 删除 SST 的四指标 CFA
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 90), "\n3. M1: 4-INDICATOR EF (drop SST)\n", strrep("=", 90), "\n", sep = "")
m1 <- 'EF =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK'
fit1 <- safe_fit(m1, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
cat("M1 fit:\n"); print(fit_table_one(fit1, "M1"), row.names = FALSE)
load1 <- extract_loadings(fit1, "M1")
rel1  <- calc_rel(fit1, "M1")
cat("M1 loadings:\n"); print(load1, row.names = FALSE)
cat("M1 reliability/AVE:\n"); print(rel1, row.names = FALSE)
save_csv(round_num(load1, 4), "M1_loadings.csv")
save_csv(round_num(rel1, 4), "M1_reliability.csv")

# ----------------------------------------------------------------------------
# 4. M2: General EF + Specific Inhibition (bifactor)
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 90), "\n4. M2: GENERAL EF + INHIBITION (bifactor)\n", strrep("=", 90), "\n", sep = "")
m2 <- paste0(
  "G =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST\n",
  "Inh =~ EF_FLANK + EF_SST\n",
  "G ~~ 0*Inh\n"
)
fit2 <- safe_fit(m2, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")

m2_status <- data.frame(Model = "M2", Converged = is_ok_fit(fit2),
                        Note = if (is_ok_fit(fit2)) "converged" else "NOT converged")
if (is_ok_fit(fit2)) {
  # 检查负残差方差 / Heywood / 非正定
  est <- lavaan::parameterEstimates(fit2)
  neg_var <- sum(est$op == "~~" & est$lhs == est$rhs & est$est < 0)
  vcov_pd <- tryCatch({
    ev <- eigen(lavaan::vcov(fit2), symmetric = TRUE, only.values = TRUE)$values
    all(is.finite(ev)) && min(ev) > 1e-8
  }, error = function(e) FALSE)
  cat("M2 fit:\n"); print(fit_table_one(fit2, "M2"), row.names = FALSE)
  cat("M2 negative residual variances (Heywood): ", neg_var, "\n", sep = "")
  cat("M2 parameter vcov positive-definite: ", vcov_pd, "\n", sep = "")
  load2 <- extract_loadings(fit2, "M2")
  cat("M2 loadings (G = general, Inh = inhibition-specific):\n"); print(load2, row.names = FALSE)
  cat("M2 factor correlation: fixed to 0 (orthogonal by design)\n")
  save_csv(round_num(load2, 4), "M2_loadings.csv")
  m2_status$Note <- paste0("converged; neg_var=", neg_var, "; vcov_pd=", vcov_pd)
} else {
  cat("M2 did NOT converge / not identified. Reporting honestly.\n")
}
m2_status$neg_var <- if (is_ok_fit(fit2)) neg_var else NA
m2_status$vcov_pd <- if (is_ok_fit(fit2)) vcov_pd else NA
save_csv(m2_status, "M2_identification.csv")

# ----------------------------------------------------------------------------
# 5. M3: ESEM (2-factor geomin) + 3-factor attempt
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 90), "\n5. M3: ESEM\n", strrep("=", 90), "\n", sep = "")

run_esem <- function(nf, label) {
  inds <- c("EF_AX", "EF_BACKWARD", "EF_WCST", "EF_FLANK", "EF_SST")
  facs <- paste0('efa("efa")*F', seq_len(nf), collapse = " + ")
  model <- paste0(facs, " =~ ", paste(inds, collapse = " + "))
  fit <- safe_fit(model, data = df, rotation = "geomin", estimator = "ML",
                  missing = "fiml", std.lv = TRUE, orthogonal = FALSE)
  list(model = model, fit = fit)
}

esem2 <- run_esem(2, "2-factor")
if (is_ok_fit(esem2$fit)) {
  cat("M3 (2-factor ESEM, geomin) fit:\n")
  print(fit_table_one(esem2$fit, "M3_ESEM2"), row.names = FALSE)
  pe <- lavaan::parameterEstimates(esem2$fit, standardized = TRUE, ci = TRUE)
  ld <- pe[pe$op == "=~", c("lhs","rhs","est","std.all","se","pvalue")]
  cat("M3 ESEM pattern loadings (main + cross):\n"); print(round_num(ld, 4), row.names = FALSE)
  # 因子相关
  corlv <- tryCatch(lavaan::lavInspect(esem2$fit, "cor.lv"), error = function(e) NULL)
  cat("M3 ESEM factor correlation:\n"); print(round(corlv, 4))
  save_csv(round_num(ld, 4), "M3_ESEM2_loadings.csv")
} else {
  cat("M3 2-factor ESEM failed; reporting honestly.\n")
}

esem3 <- run_esem(3, "3-factor")
if (is_ok_fit(esem3$fit)) {
  cat("\nM3b (3-factor ESEM, geomin) fit:\n")
  print(fit_table_one(esem3$fit, "M3_ESEM3"), row.names = FALSE)
  pe <- lavaan::parameterEstimates(esem3$fit, standardized = TRUE, ci = TRUE)
  ld <- pe[pe$op == "=~", c("lhs","rhs","est","std.all","se","pvalue")]
  cat("M3b ESEM pattern loadings:\n"); print(round_num(ld, 4), row.names = FALSE)
} else {
  cat("M3b 3-factor ESEM failed/unstable (expected with only 5 indicators).\n")
}

# ----------------------------------------------------------------------------
# 6. M4: EF composite (equal-weight z-mean) + mediation
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 90), "\n6. M4: EF COMPOSITE\n", strrep("=", 90), "\n", sep = "")
df$EF_COMP <- apply(df[, c("EF_AX","EF_BACKWARD","EF_WCST","EF_FLANK","EF_SST")], 1, function(x) {
  if (sum(!is.na(x)) < 4) return(NA_real_)
  mean(x, na.rm = TRUE)
})
cat("EF composite: N = ", sum(!is.na(df$EF_COMP)), ", M = ",
    round(mean(df$EF_COMP, na.rm=TRUE), 3), ", SD = ", round(sd(df$EF_COMP, na.rm=TRUE), 3), "\n", sep = "")

m4 <- paste0(
  "Fluency =~ FL_SIMP + FL_COMP\n",
  "Academic =~ AC_MATH + AC_CHIN + AC_ENG\n",
  "Fluency ~ a*EF_COMP\n",
  "Academic ~ b*Fluency + c_prime*EF_COMP\n",
  "indirect := a*b\n",
  "direct := c_prime\n",
  "total := c_prime + a*b\n"
)
fit4 <- safe_fit(m4, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
cat("M4 composite-mediation fit:\n"); print(fit_table_one(fit4, "M4"), row.names = FALSE)
if (is_ok_fit(fit4)) {
  st4 <- extract_structural(fit4, "M4")
  cat("M4 structural effects (delta CI):\n"); print(st4, row.names = FALSE)
  r2_4 <- lavaan::inspect(fit4, "r2")
  cat("M4 R²:\n"); print(round(r2_4, 3))
  save_csv(round_num(st4, 4), "M4_composite_effects.csv")
}

# ----------------------------------------------------------------------------
# 7. MEDIATION for M0 / M1 / M2(G)  (delta CI)
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 90), "\n7. MEDIATION MODELS (EF -> Fluency -> Achievement)\n", strrep("=", 90), "\n", sep = "")

mediation_model <- function(ef_block, predictor) {
  paste0(ef_block, "\n",
         "Fluency =~ FL_SIMP + FL_COMP\n",
         "Academic =~ AC_MATH + AC_CHIN + AC_ENG\n",
         "Fluency ~ a*", predictor, "\n",
         "Academic ~ b*Fluency + c_prime*", predictor, "\n",
         "indirect := a*b\n",
         "direct := c_prime\n",
         "total := c_prime + a*b\n")
}

med_specs <- list(
  M0 = mediation_model("EF =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST", "EF"),
  M1 = mediation_model("EF =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK", "EF"),
  M2 = mediation_model("G =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST\nInh =~ EF_FLANK + EF_SST\nG ~~ 0*Inh", "G")
)

med_fits <- lapply(med_specs, function(m) safe_fit(m, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml"))

med_struct <- do.call(rbind, lapply(names(med_fits), function(id) {
  f <- med_fits[[id]]
  if (!is_ok_fit(f)) return(data.frame(Model = id, Effect = "FAILED", b = NA, beta = NA, CI_LL = NA, CI_UL = NA, p = NA))
  pe <- lavaan::parameterEstimates(f, standardized = TRUE, ci = TRUE, level = .95)
  z <- pe[pe$op %in% c("~", ":=") & (pe$lhs %in% c("Fluency","Academic","indirect","direct","total")), , drop = FALSE]
  data.frame(Model = id,
             Effect = ifelse(is.na(z$label) | z$label == "", paste(z$lhs, z$op, z$rhs), z$label),
             b = z$est, beta = z$std.all, SE = z$se, CI_LL = z$ci.lower, CI_UL = z$ci.upper, p = z$pvalue)
}))
cat("Mediation structural effects (delta CI):\n")
print(round_num(med_struct, 4), row.names = FALSE)
save_csv(round_num(med_struct, 4), "MED_structural_delta_CI.csv")

# R² for mediation models
med_r2 <- do.call(rbind, lapply(names(med_fits), function(id) {
  f <- med_fits[[id]]
  if (!is_ok_fit(f)) return(data.frame(Model = id, Variable = NA, R2 = NA))
  r2 <- lavaan::inspect(f, "r2")
  data.frame(Model = id, Variable = names(r2), R2 = as.numeric(r2))
}))
cat("Mediation R²:\n"); print(round_num(med_r2, 3), row.names = FALSE)
save_csv(round_num(med_r2, 4), "MED_R2.csv")

# ----------------------------------------------------------------------------
# 8. CONSOLIDATED COMPARISON
# ----------------------------------------------------------------------------
cat("\n", strrep("=", 90), "\n8. CONSOLIDATED COMPARISON\n", strrep("=", 90), "\n", sep = "")

fit_all <- rbind(
  fit_table_one(fit0, "M0_5ind_CFA"),
  fit_table_one(fit1, "M1_dropSST"),
  fit_table_one(fit2, "M2_GenInh"),
  fit_table_one(esem2$fit, "M3_ESEM2"),
  fit_table_one(fit4, "M4_composite_med")
)
cat("Fit comparison:\n"); print(round_num(fit_all, 4), row.names = FALSE)
save_csv(round_num(fit_all, 4), "COMPARISON_fit.csv")

rel_all <- rbind(rel0, rel1, calc_rel(fit2, "M2_GenInh"))
cat("Reliability/AVE (M0/M1/M2):\n"); print(round_num(rel_all, 4), row.names = FALSE)
save_csv(round_num(rel_all, 4), "COMPARISON_reliability.csv")

cat("\nDONE.\n")
sink()
