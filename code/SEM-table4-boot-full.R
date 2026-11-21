# ============================================================================
# SEM-table4-boot-full.R
# M4 主模型的完整 nonparametric bootstrap：覆盖 a/b/c_prime/indirect/total
# 的 unstd(b) 与 std(β) 全部参数，以及 prop(=ind/total)。
# 用于统一 Table 4 的 SE 与 95% CI 口径（全表 nonparametric percentile CI）。
# seed 与 SEM-bootstrap.R 的 M4 一致（20260925），保证与已有 BOOT_summary_M4
# 完全对齐，仅新增此前漏算的 a_beta/b_beta/cp_beta/tot_b/prop_beta。
# ============================================================================
# 相对路径已配置：请从 V1_result/ 根目录运行本脚本（Rscript code/SEM-table4-boot-full.R）

suppressPackageStartupMessages({
  for (p in c("lavaan", "readxl")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
})

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

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
df$EF_COMP <- apply(df[, c("EF_AX","EF_BACKWARD","EF_WCST","EF_FLANK","EF_SST")], 1, function(x) {
  if (sum(!is.na(x)) < 4) return(NA_real_)
  mean(x, na.rm = TRUE)
})

m4 <- "Fluency =~ FL_SIMP + FL_COMP
Academic =~ AC_MATH + AC_CHIN + AC_ENG
Fluency ~ a*EF_COMP
Academic ~ b*Fluency + c_prime*EF_COMP
ind := a*b
total := c_prime + a*b
prop := ind/total"

f4 <- lavaan::sem(m4, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")

pe0 <- lavaan::parameterEstimates(f4, standardized = TRUE)
g0  <- function(l) pe0$est[pe0$label == l][1]
gb0 <- function(l) pe0$std.all[pe0$label == l][1]

cat("Point estimates (delta):\n")
cat(sprintf("  a       b=%.4f beta=%.4f\n", g0("a"), gb0("a")))
cat(sprintf("  b       b=%.4f beta=%.4f\n", g0("b"), gb0("b")))
cat(sprintf("  c_prime b=%.4f beta=%.4f\n", g0("c_prime"), gb0("c_prime")))
cat(sprintf("  ind     b=%.4f beta=%.4f\n", g0("ind"), gb0("ind")))
cat(sprintf("  total   b=%.4f beta=%.4f\n", g0("total"), gb0("total")))
cat(sprintf("  prop    beta=%.4f\n", gb0("ind")/gb0("total")))

# 完整 nonparametric bootstrap（R=5000）
set.seed(20260925)
n <- nrow(df); idx <- seq_len(n); R <- 5000
cols <- c("a_b","b_b","cp_b","ind_b","tot_b","a_beta","b_beta","cp_beta","ind_beta","tot_beta","prop_beta")
out <- matrix(NA_real_, R, length(cols)); colnames(out) <- cols
for (i in seq_len(R)) {
  rows <- sample(idx, n, replace = TRUE)
  d <- df[rows, , drop = FALSE]
  f <- tryCatch(suppressWarnings(lavaan::sem(m4, data = d, std.lv = TRUE, estimator = "ML", missing = "fiml")),
                error = function(e) NULL)
  if (!is.null(f) && isTRUE(lavaan::lavInspect(f, "converged"))) {
    pe <- lavaan::parameterEstimates(f, standardized = TRUE)
    g   <- function(l) pe$est[pe$label == l][1]
    gb  <- function(l) pe$std.all[pe$label == l][1]
    ind_b   <- g("ind");   ind_beta <- gb("ind")
    tot_b   <- g("total"); tot_beta <- gb("total")
    out[i, ] <- c(g("a"), g("b"), g("c_prime"), ind_b, tot_b,
                  gb("a"), gb("b"), gb("c_prime"), ind_beta, tot_beta,
                  ind_beta / tot_beta)
  }
}

ci_percentile <- function(x, level = .95) unname(quantile(x, c((1-level)/2, 1-(1-level)/2), na.rm = TRUE))

point <- c(g0("a"), g0("b"), g0("c_prime"), g0("ind"), g0("total"),
           gb0("a"), gb0("b"), gb0("c_prime"), gb0("ind"), gb0("total"),
           gb0("ind")/gb0("total"))

sm <- data.frame(
  Effect = cols,
  Point = point,
  Boot_Mean = colMeans(out, na.rm = TRUE),
  Boot_SD = apply(out, 2, sd, na.rm = TRUE),
  CI_LL = apply(out, 2, function(x) ci_percentile(x)[1]),
  CI_UL = apply(out, 2, function(x) ci_percentile(x)[2]),
  NonConv = apply(out, 2, function(x) sum(is.na(x))),
  stringsAsFactors = FALSE
)
sm[] <- lapply(sm, function(c) if (is.numeric(c)) round(c, 4) else c)

write.csv(out, file.path(results_dir, "BOOT_full_M4_raw.csv"), row.names = FALSE, na = "")
write.csv(sm,  file.path(results_dir, "BOOT_full_M4_summary.csv"), row.names = FALSE, na = "")

cat("\n=== Full bootstrap summary (M4, R=5000, nonparametric percentile) ===\n")
print(sm, row.names = FALSE)
cat("\nWrote BOOT_full_M4_raw.csv / BOOT_full_M4_summary.csv\n")
