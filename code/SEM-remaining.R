# 补全剩余数字：composite M4 拟合+标准化路径+R2；单结果 composite 拟合；composite 分性别 R2
# 相对路径已配置：请从 V1_result/ 根目录运行本脚本（Rscript code/SEM-remaining.R）
suppressPackageStartupMessages({library(lavaan); library(readxl)})
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
df$EF_COMP <- apply(df[, c("EF_AX","EF_BACKWARD","EF_WCST","EF_FLANK","EF_SST")], 1,
                    function(x){if(sum(!is.na(x))<4) NA_real_ else mean(x,na.rm=TRUE)})
df$sex_group <- factor(ifelse(df$sex %in% c(0,1), df$sex, NA), levels=c(0,1), labels=c("Female","Male"))

prfit <- function(tag, f) {
  fm <- fitMeasures(f, c("chisq","df","pvalue","cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","srmr"))
  cat(sprintf("%-20s chi2=%.3f df=%d p=%.3f CFI=%.4f TLI=%.4f RMSEA=%.3f [%.3f,%.3f] SRMR=%.3f\n",
      tag, fm["chisq"], fm["df"], fm["pvalue"], fm["cfi"], fm["tli"], fm["rmsea"],
      fm["rmsea.ci.lower"], fm["rmsea.ci.upper"], fm["srmr"]))
}

# 1) composite M4
m4 <- "Fluency =~ FL_SIMP + FL_COMP
Academic =~ AC_MATH + AC_CHIN + AC_ENG
Fluency ~ a*EF_COMP
Academic ~ b*Fluency + c_prime*EF_COMP
ind := a*b"
f4 <- sem(m4, data=df, std.lv=TRUE, estimator="ML", missing="fiml")
prfit("M4 composite", f4)
pe <- standardizedSolution(f4)
cat("[M4 standardized paths]\n")
for (lbl in c("a","b","c_prime","ind")) {
  r <- pe[pe$label==lbl, ]; cat(sprintf("%s: beta=%.4f se=%.4f p=%.4f\n", lbl, r$est.std, r$se, r$pvalue))
}
r2 <- inspect(f4, "rsquare")
cat(sprintf("M4 R2 Fluency=%.4f Academic=%.4f\n", r2["Fluency"], r2["Academic"]))

# 2) 单结果 composite fits
cat("\n[Single-outcome composite fits]\n")
for (out in c("AC_MATH","AC_CHIN","AC_ENG")) {
  m <- sprintf("Fluency =~ FL_SIMP + FL_COMP\nFluency ~ EF_COMP\n%s ~ Fluency + EF_COMP", out)
  f <- sem(m, data=df, std.lv=TRUE, estimator="ML", missing="fiml")
  prfit(out, f)
  cat(sprintf("  R2 %s = %.4f\n", out, inspect(f,"rsquare")[out]))
}

# 3) composite 分性别 R2（configural）
cat("\n[Composite configural R2 by sex]\n")
mcfg <- "Fluency =~ FL_SIMP + FL_COMP
Academic =~ AC_MATH + AC_CHIN + AC_ENG
Fluency ~ EF_COMP
Academic ~ Fluency + EF_COMP"
fc <- sem(mcfg, data=df, group="sex_group", std.lv=TRUE, estimator="ML", missing="fiml")
for (g in 1:2) {
  rr <- inspect(fc, "rsquare")[[g]]
  cat(sprintf("group %d (%s): R2 Fluency=%.4f Academic=%.4f\n",
      g, c("Female","Male")[g], rr["Fluency"], rr["Academic"]))
}
