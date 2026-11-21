# 补齐论文所需的两个拟合值（N=518）
# 相对路径已配置：请从 V1_result/ 根目录运行本脚本（Rscript code/SEM-paper-fits.R）
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

pr <- function(tag, f) {
  fm <- fitMeasures(f, c("chisq","df","pvalue","cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","srmr"))
  cat(sprintf("%-28s chi2=%.3f df=%d p=%.3f CFI=%.4f TLI=%.4f RMSEA=%.3f [%.3f, %.3f] SRMR=%.3f\n",
      tag, fm["chisq"], fm["df"], fm["pvalue"], fm["cfi"], fm["tli"], fm["rmsea"],
      fm["rmsea.ci.lower"], fm["rmsea.ci.upper"], fm["srmr"]))
}

# 1. 完整测量模型（三因子自由相关，无结构路径）
mm <- "EF =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST
Fluency =~ FL_SIMP + FL_COMP
Academic =~ AC_MATH + AC_CHIN + AC_ENG"
f_mm <- sem(mm, data=df, std.lv=TRUE, estimator="ML", missing="fiml")
pr("Measurement (3-factor)", f_mm)

# 2. 潜变量中介 M0（EF→Fluency→Academic + 直接）
m0 <- paste0(mm, "\nFluency ~ EF\nAcademic ~ Fluency + EF")
f_m0 <- sem(m0, data=df, std.lv=TRUE, estimator="ML", missing="fiml")
pr("Latent mediation M0", f_m0)

# 3. composite 中介 M4（EF_COMP 观测 → Fluency → Academic）
m4 <- "Fluency =~ FL_SIMP + FL_COMP
Academic =~ AC_MATH + AC_CHIN + AC_ENG
Fluency ~ EF_COMP
Academic ~ Fluency + EF_COMP"
f_m4 <- sem(m4, data=df, std.lv=TRUE, estimator="ML", missing="fiml")
pr("Composite mediation M4", f_m4)

# 3-factor measurement model 的因子间相关（供参考）
cat("\n[Measurement model factor correlations]\n")
pe <- standardizedSolution(f_mm)
print(pe[pe$op=="~~" & pe$lhs != pe$rhs & pe$lhs %in% c("EF","Fluency","Academic") & pe$rhs %in% c("EF","Fluency","Academic"),
         c("lhs","op","rhs","est.std")], row.names=FALSE)
