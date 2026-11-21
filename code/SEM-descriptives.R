# ============================================================================
# SEM-descriptives.R
# N=518 描述性统计：样本特征 + 各指标原始均值/SD + 性别比较 + 相关矩阵
# 用于替换论文 §3.1 的 Table 1 / Table 2
# sex: 1 = Male, 0 = Female
# 输出：results/
# ============================================================================

# 相对路径已配置：请从 V1_result/ 根目录运行本脚本（Rscript code/SEM-descriptives.R）

suppressPackageStartupMessages({
  req <- c("readxl")
  for (p in req) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
})

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(results_dir, "DESCRIPTIVES_console.log"), split = TRUE)

cat("SEM-descriptives.R | ", as.character(Sys.time()), "\n", sep = "")

df <- readxl::read_excel("data/data_analysis.xlsx", sheet = "Sheet1")
cat("nrow = ", nrow(df), "\n", sep = "")
cat("columns:\n"); print(names(df))

# ---- 性别 ----
cat("\n[sex table]\n")
print(table(df$sex, useNA = "ifany"))

# ---- 年龄 ----
cat("\n[age]\n")
cat("non-missing n = ", sum(!is.na(df$age)), ", missing = ", sum(is.na(df$age)), "\n", sep = "")
cat("M = ", round(mean(df$age, na.rm=TRUE), 3), ", SD = ", round(sd(df$age, na.rm=TRUE), 3),
    ", range = [", min(df$age, na.rm=TRUE), ", ", max(df$age, na.rm=TRUE), "]\n", sep = "")

# ---- 各指标缺失 ----
cat("\n[missing per indicator]\n")
ind_names <- c("axcpt_d_prime_context","Backward_max_span","wcst_wcst_completed_categories",
               "flanker_flanker_effect_ms","stopsignal_ssrt_ms",
               "simple_subtraction_adjusted_score","complex_subtraction_adjusted_score",
               "math_score","chinese_score","english_score")
print(sapply(ind_names, function(v) sum(is.na(df[[v]]))))

# ---- 反向计分（用于相关与描述的正确方向）----
df$Flanker_Rev <- -df$flanker_flanker_effect_ms
df$SST_Rev     <- -df$stopsignal_ssrt_ms

# 指标列：用于描述与相关（原始 + 反向）
vars <- c("axcpt_d_prime_context","Backward_max_span","wcst_wcst_completed_categories",
          "Flanker_Rev","SST_Rev",
          "simple_subtraction_adjusted_score","complex_subtraction_adjusted_score",
          "math_score","chinese_score","english_score")
labels <- c("AX-CPT d'-context","Backward digit span","WCST categories",
            "Flanker (reversed, ms)","SSRT (reversed, ms)",
            "Simple subtraction (adj.)","Complex subtraction (adj.)",
            "Math score","Chinese score","English score")

# ---- 原始均值/SD + 性别比较（Welch t）----
cat("\n[descriptive + sex comparison]\n")
desc <- data.frame(Indicator = labels,
                   M_total = NA, SD_total = NA,
                   M_boys = NA, SD_boys = NA,
                   M_girls = NA, SD_girls = NA,
                   t = NA, p = NA, d_girls_minus_boys = NA,
                   stringsAsFactors = FALSE)
for (i in seq_along(vars)) {
  v <- vars[i]
  x <- df[[v]]
  male <- df$sex == 1; female <- df$sex == 0
  xt <- tryCatch(t.test(x[male], x[female]), error = function(e) NULL)
  # Cohen d（用合并 SD）
  nm <- sum(!is.na(x[male])); nf <- sum(!is.na(x[female]))
  sm <- sd(x[male], na.rm=TRUE); sf <- sd(x[female], na.rm=TRUE)
  sp <- sqrt(((nm-1)*sm^2 + (nf-1)*sf^2) / (nm+nf-2))
  d <- (mean(x[female], na.rm=TRUE) - mean(x[male], na.rm=TRUE)) / sp
  desc[i, c("M_total","SD_total")] <- c(mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
  desc[i, c("M_boys","SD_boys")] <- c(mean(x[male], na.rm=TRUE), sd(x[male], na.rm=TRUE))
  desc[i, c("M_girls","SD_girls")] <- c(mean(x[female], na.rm=TRUE), sd(x[female], na.rm=TRUE))
  desc[i, c("t","p")] <- if (!is.null(xt)) c(xt$statistic, xt$p.value) else c(NA, NA)
  desc[i, "d_girls_minus_boys"] <- d
}
desc[] <- lapply(desc, function(c) if (is.numeric(c)) round(c, 3) else c)
print(desc, row.names = FALSE)
write.csv(desc, file.path(results_dir, "DESCRIPTIVES_table1.csv"), row.names = FALSE)

# ---- 相关矩阵（10 指标，z 标准化后）----
cat("\n[correlation matrix (z-scored indicators)]\n")
zv <- sapply(vars, function(v) as.numeric(scale(df[[v]])))
colnames(zv) <- labels
cormat <- cor(zv, use = "pairwise.complete.obs")
print(round(cormat, 3))
write.csv(round(cormat, 4), file.path(results_dir, "DESCRIPTIVES_correlations.csv"))

# ---- 关键内部相关 ----
cat("\n[key internal correlations]\n")
cat("Simple ~ Complex subtraction: ", round(cor(zv[,"Simple subtraction (adj.)"], zv[,"Complex subtraction (adj.)"]), 3), "\n")
cat("Flanker(rev) ~ SSRT(rev): ", round(cor(zv[,"Flanker (reversed, ms)"], zv[,"SSRT (reversed, ms)"]), 3), "\n")
cat("Achievement alpha (math/chin/eng): ",
    round((3/(3-1)) * (1 - (sum(apply(zv[, c("Math score","Chinese score","English score")], 2, var)) / var(rowSums(zv[, c("Math score","Chinese score","English score")])))), 3), "\n")

# ---- EF 指标相关范围 ----
ef_cols <- c("AX-CPT d'-context","Backward digit span","WCST categories","Flanker (reversed, ms)","SSRT (reversed, ms)")
ef_cor <- cormat[ef_cols, ef_cols]
diag(ef_cor) <- NA
cat("\n[EF indicator inter-correlation range]\n")
cat("min = ", round(min(ef_cor, na.rm=TRUE), 3), ", max = ", round(max(ef_cor, na.rm=TRUE), 3), "\n", sep = "")

# ---- 年龄按性别 ----
cat("\n[age by sex]\n")
cat("Boys: M=", round(mean(df$age[df$sex==1], na.rm=TRUE),3),
    " SD=", round(sd(df$age[df$sex==1], na.rm=TRUE),3), "\n", sep="")
cat("Girls: M=", round(mean(df$age[df$sex==0], na.rm=TRUE),3),
    " SD=", round(sd(df$age[df$sex==0], na.rm=TRUE),3), "\n", sep="")
cat("age sex t-test p = ", round(tryCatch(t.test(df$age[df$sex==1], df$age[df$sex==0])$p.value, error=function(e) NA), 4), "\n", sep="")

sink()
cat("DONE\n")
