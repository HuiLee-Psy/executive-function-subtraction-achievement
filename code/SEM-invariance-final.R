# ============================================================================
# SEM-invariance-final.R  (修正版：多组模型去掉回归标签，避免提前跨组约束)
# 完整性别不变性序列（N=518），主轴 composite + latent 对照：
#   configural → metric → scalar → partial scalar（释放中文/英语截距）→ structural
#   分组路径估计 + Wald 差异检验（Δa, Δb, Δc', Δindirect）
# sex: 1 = Male, 0 = Female（分组标签 Female / Male）
# 输出：results/
# ============================================================================

# 相对路径已配置：请从 V1_result/ 根目录运行本脚本（Rscript code/SEM-invariance-final.R）

suppressPackageStartupMessages({
  req <- c("lavaan", "readxl")
  for (p in req) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
})

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(results_dir, "INVARIANCE_FINAL_console.log"), split = TRUE)

cat("SEM-invariance-final.R | ", as.character(Sys.time()), " | lavaan ",
    as.character(packageVersion("lavaan")), "\n", sep = "")

round_df <- function(x, k = 4) {
  x[] <- lapply(x, function(c) if (is.numeric(c)) round(c, k) else c)
  x
}

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
  if (sum(!is.na(x)) < 4) return(NA_real_); mean(x, na.rm = TRUE)
})
df$sex_group <- factor(ifelse(df$sex %in% c(0, 1), df$sex, NA),
                       levels = c(0, 1), labels = c("Female", "Male"))

safe_fit <- function(model, data, ...) {
  withCallingHandlers(
    tryCatch(lavaan::sem(model, data = data, ...),
             error = function(e) { cat("MODEL ERROR: ", conditionMessage(e), "\n", sep = ""); NULL }),
    warning = function(w) invokeRestart("muffleWarning")
  )
}

fit_row <- function(f, tag, step) {
  if (is.null(f) || !lavaan::lavInspect(f, "converged"))
    return(data.frame(Tag = tag, Step = step, chisq = NA, df = NA, CFI = NA, TLI = NA, RMSEA = NA, SRMR = NA))
  fm <- lavaan::fitMeasures(f, c("chisq","df","cfi","tli","rmsea","srmr"))
  data.frame(Tag = tag, Step = step, chisq = unname(fm["chisq"]), df = unname(fm["df"]),
             CFI = unname(fm["cfi"]), TLI = unname(fm["tli"]),
             RMSEA = unname(fm["rmsea"]), SRMR = unname(fm["srmr"]))
}

# 从 configural 模型提取分组 a/b/c'（按 lhs~rhs 定位，不依赖标签）
group_paths <- function(f, predictor) {
  pe <- lavaan::parameterEstimates(f, standardized = TRUE)
  grp_lab <- c("Female", "Male")
  rows <- list()
  for (g in 1:2) {
    sub <- pe[pe$group == g, ]
    pick <- function(lhs, rhs) {
      sub[sub$op == "~" & sub$lhs == lhs & sub$rhs == rhs, ]
    }
    a <- pick("Fluency", predictor)
    b <- pick("Academic", "Fluency")
    c <- pick("Academic", predictor)
    rows[[g]] <- data.frame(Group = grp_lab[g],
      a_b = a$est, a_beta = a$std.all, a_se = a$se,
      b_b = b$est, b_beta = b$std.all, b_se = b$se,
      c_b = c$est, c_beta = c$std.all, c_se = c$se)
  }
  out <- do.call(rbind, rows)
  out$ind_b <- out$a_b * out$b_b
  out$ind_beta <- out$a_beta * out$b_beta
  # Sobel delta-method SE for indirect effect（假设 a、b 独立）
  out$ind_se <- sqrt(out$a_b^2 * out$b_se^2 + out$b_b^2 * out$a_se^2)
  out
}

wald_diff <- function(gp, lab, name) {
  est <- c(gp[gp$Group == "Female", paste0(lab, "_b")], gp[gp$Group == "Male", paste0(lab, "_b")])
  se  <- c(gp[gp$Group == "Female", paste0(lab, "_se")], gp[gp$Group == "Male", paste0(lab, "_se")])
  d <- est[1] - est[2]
  sed <- sqrt(se[1]^2 + se[2]^2)  # 独立分组
  z <- d / sed
  data.frame(Path = name, Delta = d, SE = sed, z = z, p = 2 * (1 - pnorm(abs(z))))
}

run_inv_full <- function(ef_block, predictor, tag) {
  cat("\n===== ", tag, " full invariance (N=518) =====\n", sep = "")
  base <- paste0(ef_block, "\n",
    "Fluency =~ FL_SIMP + FL_COMP\n",
    "Academic =~ AC_MATH + AC_CHIN + AC_ENG\n",
    "Fluency ~ ", predictor, "\n",
    "Academic ~ Fluency + ", predictor, "\n")
  seq_def <- list(
    configural = list(equal = c(), partial = c()),
    metric     = list(equal = c("loadings"), partial = c()),
    scalar     = list(equal = c("loadings", "intercepts"), partial = c()),
    partial_scalar = list(equal = c("loadings", "intercepts"), partial = c("AC_CHIN ~ 1", "AC_ENG ~ 1")),
    structural = list(equal = c("loadings", "intercepts", "regressions"), partial = c("AC_CHIN ~ 1", "AC_ENG ~ 1"))
  )
  fits <- list()
  for (nm in names(seq_def)) {
    ge <- seq_def[[nm]]$equal; gp <- seq_def[[nm]]$partial
    f <- safe_fit(base, data = df, group = "sex_group", std.lv = TRUE,
                  estimator = "ML", missing = "fiml",
                  group.equal = ge, group.partial = gp)
    fits[[nm]] <- f
  }
  tbl <- do.call(rbind, lapply(names(fits), function(nm) fit_row(fits[[nm]], tag, nm)))
  tbl$dCFI <- c(NA, diff(tbl$CFI))
  tbl$dchisq <- c(NA, diff(tbl$chisq))
  tbl$ddf <- c(NA, diff(tbl$df))
  # 新增：明确标注参照模型的 ΔCFI/Δχ²/Δdf 派生列（保留原有 dCFI/dchisq/ddf 逐步变化，不改变既有分析逻辑）
  tbl$dCFI_vs_previous <- c(NA, diff(tbl$CFI))          # 相对上一模型（与 dCFI 同值，仅显式命名）
  metric_idx <- which(tbl$Step == "metric")               # 按名称识别 metric 基准模型，避免假设固定行位置
  if (length(metric_idx) == 1) {
    tbl$dCFI_vs_metric   <- tbl$CFI   - tbl$CFI[metric_idx]    # 相对 metric 模型（正文 scalar/partial scalar 口径）
    tbl$dchisq_vs_metric <- tbl$chisq - tbl$chisq[metric_idx]  # 相对 metric 模型（正文 partial scalar Δχ²=2.15 口径）
    tbl$ddf_vs_metric    <- tbl$df    - tbl$df[metric_idx]     # 相对 metric 模型（正文 partial scalar Δdf=1 口径）
    if (metric_idx > 1) {
      tbl$dCFI_vs_metric[seq_len(metric_idx - 1)]   <- NA   # metric 之前的行不适用
      tbl$dchisq_vs_metric[seq_len(metric_idx - 1)] <- NA
      tbl$ddf_vs_metric[seq_len(metric_idx - 1)]    <- NA
    }
  } else {
    tbl$dCFI_vs_metric   <- NA_real_
    tbl$dchisq_vs_metric <- NA_real_
    tbl$ddf_vs_metric    <- NA_real_
  }
  cat("[delta references] dCFI/dchisq/ddf = vs previous step; *_vs_metric = vs metric model\n")
  gp <- group_paths(fits$configural, predictor)
  cat("[fit sequence]\n"); print(round_df(tbl, 4), row.names = FALSE)
  cat("[group paths (configural)]\n"); print(round_df(gp, 4), row.names = FALSE)
  w <- do.call(rbind, list(
    wald_diff(gp, "a", "a (EF->Fluency)"),
    wald_diff(gp, "b", "b (Fluency->Academic)"),
    wald_diff(gp, "c", "c' (EF->Academic)"),
    wald_diff(gp, "ind", "indirect (a*b)")
  ))
  cat("[Wald difference tests (Female - Male)]\n"); print(round_df(w, 4), row.names = FALSE)
  list(tbl = tbl, group = gp, wald = w)
}

inv_comp <- run_inv_full("", "EF_COMP", "Composite")
inv_lat  <- run_inv_full("EF =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST", "EF", "Latent")

cat("\n[Combined fit sequence]\n")
print(round_df(rbind(inv_comp$tbl, inv_lat$tbl), 4), row.names = FALSE)
write.csv(rbind(inv_comp$tbl, inv_lat$tbl), file.path(results_dir, "INVARIANCE_FINAL_fit.csv"), row.names = FALSE)

cat("\n[Combined group paths]\n")
print(round_df(rbind(cbind(Tag="Composite", inv_comp$group), cbind(Tag="Latent", inv_lat$group)), 4), row.names = FALSE)
write.csv(rbind(cbind(Tag="Composite", inv_comp$group), cbind(Tag="Latent", inv_lat$group)),
          file.path(results_dir, "INVARIANCE_FINAL_group_paths.csv"), row.names = FALSE)

cat("\n[Combined Wald]\n")
print(round_df(rbind(cbind(Tag="Composite", inv_comp$wald), cbind(Tag="Latent", inv_lat$wald)), 4), row.names = FALSE)
write.csv(rbind(cbind(Tag="Composite", inv_comp$wald), cbind(Tag="Latent", inv_lat$wald)),
          file.path(results_dir, "INVARIANCE_FINAL_wald.csv"), row.names = FALSE)

sink()
cat("DONE\n")
