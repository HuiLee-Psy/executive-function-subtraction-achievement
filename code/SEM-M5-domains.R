# ============================================================================
# SEM-M5-domains.R
# M5 理论竞争模型：三个理论 domain 的等权 composite（不建 latent）
#   WM_pro   = mean(EF_AX, EF_BACKWARD)    工作记忆/主动控制
#   FLEX     = EF_WCST                     认知灵活性（单指标）
#   INH      = mean(EF_FLANK, EF_SST)      抑制控制
# 结构：WM_pro/FLEX/INH -> Fluency -> Academic（+ 直接路径）
# 对比对象：M4（单一 general EF composite）
# 输出：results/
# ============================================================================

# 相对路径已配置：请从 V1_result/ 根目录运行本脚本（Rscript code/SEM-M5-domains.R）

suppressPackageStartupMessages({
  req <- c("lavaan", "readxl")
  for (p in req) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
})

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(results_dir, "M5_console.log"), split = TRUE)

round_df <- function(x, k = 4) {
  x[] <- lapply(x, function(c) if (is.numeric(c)) round(c, k) else c)
  x
}

cat("SEM-M5-domains.R | ", as.character(Sys.time()), " | lavaan ",
    as.character(packageVersion("lavaan")), "\n", sep = "")

df <- readxl::read_excel("data/data_analysis.xlsx", sheet = "Sheet1")
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

comp2 <- function(x, y, min_n = 1L) {
  m <- cbind(x, y)
  apply(m, 1, function(r) if (sum(!is.na(r)) < min_n) NA_real_ else mean(r, na.rm = TRUE))
}
df$WM_pro <- comp2(df$EF_AX, df$EF_BACKWARD)       # 无缺失
df$FLEX   <- df$EF_WCST                              # 单指标
df$INH    <- comp2(df$EF_FLANK, df$EF_SST)           # 2 行仅 SST

# 报告 composite 缺失与描述
cat("\n[Composite descriptive]\n")
for (v in c("WM_pro","FLEX","INH")) {
  cat(v, ": N = ", sum(!is.na(df[[v]])), ", M = ", round(mean(df[[v]], na.rm=TRUE),3),
      ", SD = ", round(sd(df[[v]], na.rm=TRUE),3), "\n", sep = "")
}

# domain 相关矩阵（含流畅性/学业）
cat("\n[Domain composite correlations]\n")
dom_cor <- cor(df[, c("WM_pro","FLEX","INH","FL_SIMP","FL_COMP","AC_MATH","AC_CHIN","AC_ENG")],
               use = "pairwise.complete.obs")
print(round(dom_cor, 3))
write.csv(round(dom_cor, 4), file.path(results_dir, "M5_domain_correlations.csv"))

# 内部一致性（domain 内两指标相关）
cat("\n[Within-domain indicator correlations]\n")
cat("WM_pro (AX ~ Corsi): ", round(cor(df$EF_AX, df$EF_BACKWARD, use="pairwise.complete.obs"), 3), "\n")
cat("INH (Flanker ~ SST): ", round(cor(df$EF_FLANK, df$EF_SST, use="pairwise.complete.obs"), 3), "\n")

# ----------------------------------------------------------------------------
# M5 结构模型
# ----------------------------------------------------------------------------
m5 <- paste0(
  "Fluency =~ FL_SIMP + FL_COMP\n",
  "Academic =~ AC_MATH + AC_CHIN + AC_ENG\n",
  "Fluency ~ a1*WM_pro + a2*FLEX + a3*INH\n",
  "Academic ~ b*Fluency + c1*WM_pro + c2*FLEX + c3*INH\n",
  "ind_WM := a1*b\n",
  "ind_FLEX := a2*b\n",
  "ind_INH := a3*b\n",
  "ind_total := a1*b + a2*b + a3*b\n"
)

fit5 <- lavaan::sem(m5, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
cat("\n[M5 fit]\n")
fm5 <- lavaan::fitMeasures(fit5, c("chisq","df","pvalue","cfi","tli","rmsea",
  "rmsea.ci.lower","rmsea.ci.upper","srmr","aic","bic"))
print(round(fm5, 4))

# M4 重新拟合（同条件下，取 AIC/BIC 对比）
m4 <- paste0(
  "Fluency =~ FL_SIMP + FL_COMP\n",
  "Academic =~ AC_MATH + AC_CHIN + AC_ENG\n",
  "Fluency ~ a*EF_COMP\n",
  "Academic ~ b*Fluency + c_prime*EF_COMP\n",
  "indirect := a*b\n"
)
df$EF_COMP <- apply(df[, c("EF_AX","EF_BACKWARD","EF_WCST","EF_FLANK","EF_SST")], 1, function(x) {
  if (sum(!is.na(x)) < 4) return(NA_real_); mean(x, na.rm = TRUE)
})
fit4 <- lavaan::sem(m4, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
fm4 <- lavaan::fitMeasures(fit4, c("chisq","df","pvalue","cfi","tli","rmsea","srmr","aic","bic"))

cat("\n[M4 vs M5 comparison]\n")
comp <- data.frame(
  Model = c("M4_unity(1 composite)", "M5_diversity(3 domains)"),
  chisq = c(unname(fm4["chisq"]), unname(fm5["chisq"])),
  df = c(unname(fm4["df"]), unname(fm5["df"])),
  p = c(unname(fm4["pvalue"]), unname(fm5["pvalue"])),
  CFI = c(unname(fm4["cfi"]), unname(fm5["cfi"])),
  TLI = c(unname(fm4["tli"]), unname(fm5["tli"])),
  RMSEA = c(unname(fm4["rmsea"]), unname(fm5["rmsea"])),
  SRMR = c(unname(fm4["srmr"]), unname(fm5["srmr"])),
  AIC = c(unname(fm4["aic"]), unname(fm5["aic"])),
  BIC = c(unname(fm4["bic"]), unname(fm5["bic"])),
  Fluency_R2 = c(unname(lavaan::inspect(fit4, "r2")["Fluency"]), unname(lavaan::inspect(fit5, "r2")["Fluency"])),
  Academic_R2 = c(unname(lavaan::inspect(fit4, "r2")["Academic"]), unname(lavaan::inspect(fit5, "r2")["Academic"]))
)
comp <- round_df(comp, 4)
print(comp, row.names = FALSE)
write.csv(comp, file.path(results_dir, "M5_vs_M4_comparison.csv"), row.names = FALSE)

# M5 路径 + delta CI
pe5 <- lavaan::parameterEstimates(fit5, standardized = TRUE, ci = TRUE, level = .95)
paths <- pe5[pe5$op %in% c("~", ":=") &
             pe5$lhs %in% c("Fluency","Academic","ind_WM","ind_FLEX","ind_INH","ind_total"),
             c("lhs","op","rhs","est","se","ci.lower","ci.upper","pvalue","std.all")]
paths <- round_df(data.frame(paths), 4)
cat("\n[M5 paths + delta CI]\n"); print(paths, row.names = FALSE)
write.csv(paths, file.path(results_dir, "M5_paths_delta_CI.csv"), row.names = FALSE)

# R²
cat("\n[M5 R²]\n"); print(round(lavaan::inspect(fit5, "r2"), 4))

# ----------------------------------------------------------------------------
# Bootstrap 5000（间接效应 ind_WM / ind_FLEX / ind_INH / ind_total）
# ----------------------------------------------------------------------------
cat("\n=== Bootstrap M5 (R = 5000) ===\n")
set.seed(20260923)
n <- nrow(df); idx <- seq_len(n)
R <- 5000
boot_mat <- matrix(NA_real_, R, 8)
colnames(boot_mat) <- c("a1","a2","a3","b","ind_WM","ind_FLEX","ind_INH","ind_total")
colnames_beta <- c("ind_WM_beta","ind_FLEX_beta","ind_INH_beta","ind_total_beta")
boot_beta <- matrix(NA_real_, R, 4)
for (i in seq_len(R)) {
  rows <- sample(idx, n, replace = TRUE)
  d <- df[rows, , drop = FALSE]
  f <- tryCatch(suppressWarnings(lavaan::sem(m5, data = d, std.lv = TRUE, estimator = "ML", missing = "fiml")),
                error = function(e) NULL)
  if (!is.null(f) && isTRUE(lavaan::lavInspect(f, "converged"))) {
    pe <- lavaan::parameterEstimates(f, standardized = TRUE)
    g  <- function(lab) pe$est[pe$label == lab][1]
    gb <- function(lab) pe$std.all[pe$label == lab][1]
    boot_mat[i, ] <- c(g("a1"), g("a2"), g("a3"), g("b"), g("ind_WM"), g("ind_FLEX"), g("ind_INH"), g("ind_total"))
    boot_beta[i, ] <- c(gb("ind_WM"), gb("ind_FLEX"), gb("ind_INH"), gb("ind_total"))
  }
}
write.csv(boot_mat, file.path(results_dir, "M5_BOOT_raw.csv"), row.names = FALSE, na = "")
write.csv(boot_beta, file.path(results_dir, "M5_BOOT_raw_beta.csv"), row.names = FALSE, na = "")

ci_p <- function(x) unname(quantile(x, c(.025, .975), na.rm = TRUE))
# point (beta) 与 bootstrap beta
pt_beta <- sapply(c("ind_WM","ind_FLEX","ind_INH","ind_total"), function(lab)
  pe5$std.all[pe5$op == ":=" & pe5$lhs == lab][1])
boot_summ <- data.frame(
  Effect = colnames_beta,
  Point_beta = as.numeric(pt_beta),
  Boot_Mean = colMeans(boot_beta, na.rm = TRUE),
  Boot_SD = apply(boot_beta, 2, sd, na.rm = TRUE),
  CI_LL = apply(boot_beta, 2, function(x) ci_p(x)[1]),
  CI_UL = apply(boot_beta, 2, function(x) ci_p(x)[2]),
  NonConv = apply(boot_beta, 2, function(x) sum(is.na(x)))
)
boot_summ[] <- lapply(boot_summ, function(c) if (is.numeric(c)) round(c, 4) else c)
cat("\n[M5 bootstrap indirect effects (beta)]\n"); print(boot_summ, row.names = FALSE)
write.csv(boot_summ, file.path(results_dir, "M5_BOOT_indirect_beta.csv"), row.names = FALSE)

sink()
cat("DONE\n")
