# ============================================================================
# SEM-crossdomain-tost-std.R
# 修正版跨域 TOST：标准化口径（est.std）。
# 原 SEM-crossdomain-tost.R 用 coef()/vcov()（未标准化）算 Δβ，并注释假设
# "未标准化≈标准化"——对 composite 不成立（EF_COMP 是 5 个 z 分数的均值，SD≈0.54≠1）。
# 本脚本改用 standardizedSolution() 的 est.std 与 lavInspect(.,"vcov.std.all")，
# 报告标准化的间接效应差 Δβ 及其 90% CI（Monte Carlo, 10,000 draws, |Δβ|=.10）。
# 不动 M0–M5 主模型与结构估计，仅重算这一处探索性等价检验。
# ============================================================================
# 相对路径已配置：请从 V1_result/ 根目录运行本脚本（Rscript code/SEM-crossdomain-tost-std.R）
suppressPackageStartupMessages({library(lavaan); library(readxl); library(MASS)})

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

# ---- joint 三结果模型（LR 检验 + TOST） ----
m_free <- "Fluency =~ FL_SIMP + FL_COMP
Fluency ~ a*EF_COMP
AC_MATH ~ b_Math*Fluency + c_Math*EF_COMP
AC_CHIN ~ b_Chin*Fluency + c_Chin*EF_COMP
AC_ENG  ~ b_Eng*Fluency  + c_Eng*EF_COMP
ind_Math := a*b_Math
ind_Chin := a*b_Chin
ind_Eng  := a*b_Eng"
f_free <- sem(m_free, data=df, std.lv=TRUE, estimator="ML", missing="fiml")

m_eq <- "Fluency =~ FL_SIMP + FL_COMP
Fluency ~ a*EF_COMP
AC_MATH ~ b*Fluency + c_Math*EF_COMP
AC_CHIN ~ b*Fluency + c_Chin*EF_COMP
AC_ENG  ~ b*Fluency + c_Eng*EF_COMP"
f_eq <- sem(m_eq, data=df, std.lv=TRUE, estimator="ML", missing="fiml")

fm <- function(f) fitMeasures(f, c("chisq","df","cfi","rmsea"))
cat("FREE chi2=", fm(f_free)["chisq"], " df=", fm(f_free)["df"], " CFI=", fm(f_free)["cfi"],
    " RMSEA=", fm(f_free)["rmsea"], "\n")
cat("EQUAL chi2=", fm(f_eq)["chisq"], " df=", fm(f_eq)["df"], " CFI=", fm(f_eq)["cfi"], "\n")
dch <- fm(f_eq)["chisq"] - fm(f_free)["chisq"]; ddf <- fm(f_eq)["df"] - fm(f_free)["df"]
cat(sprintf("LR test: dchi2=%.4f ddf=%d p=%.4f\n", dch, ddf, pchisq(dch, ddf, lower.tail=FALSE)))

# ---- 标准化点估计（est.std） ----
pe <- standardizedSolution(f_free, type="std.all")
get_std <- function(lbl) pe$est.std[pe$label==lbl][1]
cat("\n[Standardized estimates (est.std)]\n")
for (e in c("a","b_Math","b_Chin","b_Eng","ind_Math","ind_Chin","ind_Eng")) {
  r <- pe[pe$label==e, ]
  cat(sprintf("  %-10s est.std=%.4f se=%.4f\n", e, r$est.std, r$se))
}

indM <- get_std("ind_Math"); indC <- get_std("ind_Chin"); indE <- get_std("ind_Eng")
dMC <- indM-indC; dME <- indM-indE; dCE <- indC-indE
cat(sprintf("\n[Standardized pairwise differences (point)]\n  Math-Chin=%.4f  Math-Eng=%.4f  Chin-Eng=%.4f\n",
            dMC, dME, dCE))

# ---- 标准化参数的 vcov（std.all） ----
V <- lavInspect(f_free, "vcov.std.all")
dn <- dimnames(V)[[1]]
cat("\n[vcov.std.all] dims=", paste(dim(V), collapse="x"),
    " has_names=", !is.null(dn), "\n", sep="")
if (!is.null(dn)) cat("  labels:", paste(dn, collapse=", "), "\n")

need <- c("a","b_Math","b_Chin","b_Eng")
if (is.null(dn) || !all(need %in% dn)) stop("vcov.std.all 标签缺失，需手动映射")
Vs <- V[need, need, drop=FALSE]

# ---- 配对 Wald（delta 法，标准化） ----
wald_pair <- function(n1, n2) {
  ai <- get_std("a"); bi <- get_std(paste0("b_", n1)); bj <- get_std(paste0("b_", n2))
  g <- c(bi-bj, ai, -ai)
  nm <- c("a", paste0("b_", n1), paste0("b_", n2))
  Vw <- Vs[nm, nm, drop=FALSE]
  d <- ai*bi - ai*bj
  se <- sqrt(as.numeric(t(g) %*% Vw %*% g))
  c(d=d, se=se, z=d/se, p=2*(1-pnorm(abs(d/se))))
}
cat("\n[Pairwise Wald (standardized delta)]\n")
for (pair in list(c("Math","Chin"), c("Math","Eng"), c("Chin","Eng"))) {
  r <- wald_pair(pair[1], pair[2])
  cat(sprintf("  Wald %s vs %s: d=%.4f se=%.4f z=%.3f p=%.3f\n",
              pair[1], pair[2], r["d"], r["se"], r["z"], r["p"]))
}

# ---- TOST（Monte Carlo, 标准化, 等价界 |Δβ|=.10） ----
set.seed(20260531)
mu <- c(a=get_std("a"), b_Math=get_std("b_Math"), b_Chin=get_std("b_Chin"), b_Eng=get_std("b_Eng"))
Sig <- Vs[c("a","b_Math","b_Chin","b_Eng"), c("a","b_Math","b_Chin","b_Eng"), drop=FALSE]
R <- 10000
draws <- mvrnorm(R, mu, Sig)
ind <- cbind(draws[,"a"]*draws[,"b_Math"], draws[,"a"]*draws[,"b_Chin"], draws[,"a"]*draws[,"b_Eng"])
diffs <- cbind(ind[,1]-ind[,2], ind[,1]-ind[,3], ind[,2]-ind[,3])
colnames(diffs) <- c("Math-Chin","Math-Eng","Chin-Eng")
cat("\n[TOST standardized (Monte Carlo 10,000 draws, |Δβ|=.10)]\n")
for (j in 1:3) {
  ci <- quantile(diffs[,j], c(.05,.95))
  cat(sprintf("  %-10s: mean=%.4f 90%%CI=[%.4f, %.4f]\n",
              colnames(diffs)[j], mean(diffs[,j]), ci[1], ci[2]))
}

# 写入 CSV 供追溯
out <- data.frame(
  pair = colnames(diffs),
  dbeta = c(dMC, dME, dCE),
  mc_mean = apply(diffs, 2, mean),
  ci_ll = apply(diffs, 2, function(x) quantile(x, .05)),
  ci_ul = apply(diffs, 2, function(x) quantile(x, .95)),
  stringsAsFactors = FALSE
)
write.csv(out, "results/TOST_std_pairwise.csv", row.names = FALSE)
cat("\nWrote results/TOST_std_pairwise.csv\n")
