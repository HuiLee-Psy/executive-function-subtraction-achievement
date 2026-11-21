library(lavaan)
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

df$grp <- factor(df$sex,levels=c(0,1),labels=c("F","M"))
s <- "Fluency =~ FL_SIMP+FL_COMP
Academic =~ AC_MATH+AC_CHIN+AC_ENG
Fluency ~ c(af,am)*EF_COMP
Academic ~ c(bf,bm)*Fluency + c(cf,cm)*EF_COMP
da:=af-am
db:=bf-bm
dc:=cf-cm
di:=af*bf-am*bm"
ff <- sem(s,data=df,group="grp",std.lv=TRUE,missing="fiml",group.equal=c("loadings","intercepts"),group.partial=c("AC_CHIN~1","AC_ENG~1"))
print(subset(parameterEstimates(ff),op==":=")[,c("lhs","est","se","pvalue","ci.lower","ci.upper")])
print(lavTestWald(ff,"af==am\nbf==bm\ncf==cm"))
