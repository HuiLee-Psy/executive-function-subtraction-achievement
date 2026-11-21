# ============================================================================
# SEM-downstream.R
# 在 data_analysis.xlsx (N=518) 上重跑下游分析，主轴改用 composite (M4)：
#   1. composite 跨域等价性（3 个单结局 + 联合三结局 + LR 检验）
#   2. latent EF (M0) 跨域等价性（稳健性对照）
#   3. composite 性别不变性（多组序列）
#   4. latent EF (M0) 性别不变性（稳健性对照）
# sex 编码：1 = Male, 0 = Female
# 输出：results/
# ============================================================================

# 相对路径已配置：请从 V1_result/ 根目录运行本脚本（Rscript code/SEM-downstream.R）

suppressPackageStartupMessages({
  req <- c("lavaan", "readxl")
  for (p in req) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
})

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(results_dir, "DOWNSTREAM_console.log"), split = TRUE)

cat("SEM-downstream.R | ", as.character(Sys.time()), " | lavaan ",
    as.character(packageVersion("lavaan")), "\n", sep = "")

round_df <- function(x, k = 4) {
  x[] <- lapply(x, function(c) if (is.numeric(c)) round(c, k) else c)
  x
}

df <- readxl::read_excel("data/data_analysis.xlsx", sheet = "Sheet1")
df$Flanker_Rev <- -df$flanker_flanker_effect_ms
df$SST_Rev     <- -df$stopsignal_ssrt_ms
df$EF_AX       <- as.numeric(scale(df$axcpt_d_prime_context))
df$EF_BACKWARD <- as.numeric(scale(df$Backward_max_span))       # 实为 Backward digit span
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
df$sex_group <- factor(ifelse(df$sex %in% c(0, 1), df$sex, NA),
                       levels = c(0, 1), labels = c("Female", "Male"))
cat("N = ", nrow(df), "; sex: ", paste(names(table(df$sex_group)), table(df$sex_group), collapse = ", "), "\n", sep = "")

safe_fit <- function(model, data, ...) {
  withCallingHandlers(
    tryCatch(lavaan::sem(model, data = data, ...),
             error = function(e) { cat("MODEL ERROR: ", conditionMessage(e), "\n", sep = ""); NULL }),
    warning = function(w) invokeRestart("muffleWarning")
  )
}

# ----------------------------------------------------------------------------
# 1. 跨域等价性（单结局 + 联合三结局 + LR）
# ----------------------------------------------------------------------------
run_crossdomain <- function(ef_block, predictor, tag) {
  cat("\n=== ", tag, " cross-domain (N=518) ===\n", sep = "")
  outcomes <- c(AC_MATH = "AC_MATH", AC_CHIN = "AC_CHIN", AC_ENG = "AC_ENG")
  single <- lapply(names(outcomes), function(onm) {
    m <- paste0(ef_block, "\n",
      "Fluency =~ FL_SIMP + FL_COMP\n",
      "Fluency ~ a*", predictor, "\n",
      outcomes[[onm]], " ~ b*Fluency + c*", predictor, "\n",
      "indirect := a*b\n")
    safe_fit(m, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
  })
  names(single) <- names(outcomes)
  tbl <- do.call(rbind, lapply(names(single), function(onm) {
    f <- single[[onm]]
    if (is.null(f) || !lavaan::lavInspect(f, "converged"))
      return(data.frame(Tag = tag, Outcome = onm, a = NA, b = NA, c_prime = NA, indirect_b = NA, indirect_beta = NA, R2_outcome = NA))
    pe <- lavaan::parameterEstimates(f, standardized = TRUE)
    g <- function(l) pe$est[pe$label == l][1]; gb <- function(l) pe$std.all[pe$label == l][1]
    data.frame(Tag = tag, Outcome = onm, a = g("a"), b = g("b"), c_prime = g("c"),
               indirect_b = g("indirect"), indirect_beta = gb("indirect"),
               R2_outcome = unname(lavaan::inspect(f, "r2")[outcomes[[onm]]]))
  }))
  # 联合三结局：自由 vs 约束 b 路径
  joint_free <- paste0(ef_block, "\n",
    "Fluency =~ FL_SIMP + FL_COMP\n",
    "Fluency ~ a*", predictor, "\n",
    "AC_MATH ~ b1*Fluency + c1*", predictor, "\n",
    "AC_CHIN ~ b2*Fluency + c2*", predictor, "\n",
    "AC_ENG ~ b3*Fluency + c3*", predictor, "\n",
    "ind_M := a*b1\nind_C := a*b2\nind_E := a*b3\n")
  joint_eq <- paste0(ef_block, "\n",
    "Fluency =~ FL_SIMP + FL_COMP\n",
    "Fluency ~ a*", predictor, "\n",
    "AC_MATH ~ b*Fluency + c1*", predictor, "\n",
    "AC_CHIN ~ b*Fluency + c2*", predictor, "\n",
    "AC_ENG ~ b*Fluency + c3*", predictor, "\n")
  jf <- safe_fit(joint_free, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
  je <- safe_fit(joint_eq, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
  lr <- if (!is.null(jf) && !is.null(je)) lavaan::anova(jf, je) else NULL
  cat("Joint model LR test (free vs equal b):\n")
  if (!is.null(lr)) print(round_df(as.data.frame(lr), 4), row.names = FALSE)
  # 单结局间接效应 bootstrap
  cat("Bootstrap (R=2000) for single-outcome indirect effects...\n")
  boot_out <- lapply(names(single), function(onm) {
    m <- paste0(ef_block, "\n",
      "Fluency =~ FL_SIMP + FL_COMP\n",
      "Fluency ~ a*", predictor, "\n",
      outcomes[[onm]], " ~ b*Fluency + c*", predictor, "\n",
      "indirect := a*b\n")
    set.seed(20260923 + match(onm, names(outcomes)))
    n <- nrow(df); idx <- seq_len(n); R <- 2000
    bb <- numeric(R)
    for (i in seq_len(R)) {
      rows <- sample(idx, n, replace = TRUE)
      d <- df[rows, , drop = FALSE]
      f <- tryCatch(suppressWarnings(lavaan::sem(m, data = d, std.lv = TRUE, estimator = "ML", missing = "fiml")),
                    error = function(e) NULL)
      if (!is.null(f) && isTRUE(lavaan::lavInspect(f, "converged"))) {
        pe <- lavaan::parameterEstimates(f, standardized = TRUE)
        bb[i] <- pe$std.all[pe$label == "indirect"][1]
      } else bb[i] <- NA_real_
    }
    q <- quantile(bb, c(.025, .975), na.rm = TRUE)
    data.frame(Tag = tag, Outcome = onm, Boot_Mean = mean(bb, na.rm=TRUE),
               CI_LL = unname(q[1]), CI_UL = unname(q[2]), NonConv = sum(is.na(bb)))
  })
  boot_df <- do.call(rbind, boot_out)
  list(single = tbl, boot = boot_df)
}

cd_comp <- run_crossdomain("", "EF_COMP", "Composite")
cd_lat  <- run_crossdomain("EF =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST", "EF", "Latent")

cat("\n[Cross-domain summary: single-outcome indirect effects]\n")
print(round_df(rbind(cd_comp$single, cd_lat$single), 4), row.names = FALSE)
cat("\n[Cross-domain bootstrap (beta)]\n")
print(round_df(rbind(cd_comp$boot, cd_lat$boot), 4), row.names = FALSE)
write.csv(rbind(cd_comp$single, cd_lat$single), file.path(results_dir, "DOWNSTREAM_crossdomain_single.csv"), row.names = FALSE)
write.csv(rbind(cd_comp$boot, cd_lat$boot), file.path(results_dir, "DOWNSTREAM_crossdomain_boot.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 2. 性别不变性（多组序列）
# ----------------------------------------------------------------------------
run_invariance <- function(ef_block, predictor, tag) {
  cat("\n=== ", tag, " sex invariance (N=518) ===\n", sep = "")
  base <- paste0(ef_block, "\n",
    "Fluency =~ FL_SIMP + FL_COMP\n",
    "Academic =~ AC_MATH + AC_CHIN + AC_ENG\n",
    "Fluency ~ a*", predictor, "\n",
    "Academic ~ b*Fluency + c_prime*", predictor, "\n")
  steps <- list(
    configural = c(),
    metric = c("loadings"),
    scalar = c("loadings", "intercepts"),
    structural = c("loadings", "intercepts", "regressions")
  )
  fits <- list()
  for (nm in names(steps)) {
    ge <- steps[[nm]]
    f <- safe_fit(base, data = df, group = "sex_group", std.lv = TRUE,
                  estimator = "ML", missing = "fiml", group.equal = ge)
    fits[[nm]] <- f
  }
  tbl <- do.call(rbind, lapply(names(fits), function(nm) {
    f <- fits[[nm]]
    if (is.null(f) || !lavaan::lavInspect(f, "converged"))
      return(data.frame(Tag = tag, Step = nm, chisq = NA, df = NA, CFI = NA, TLI = NA, RMSEA = NA, SRMR = NA))
    fm <- lavaan::fitMeasures(f, c("chisq","df","cfi","tli","rmsea","srmr"))
    data.frame(Tag = tag, Step = nm, chisq = unname(fm["chisq"]), df = unname(fm["df"]),
               CFI = unname(fm["cfi"]), TLI = unname(fm["tli"]),
               RMSEA = unname(fm["rmsea"]), SRMR = unname(fm["srmr"]))
  }))
  # ΔCFI / Δχ² 相对上一级
  if (nrow(tbl) == 4 && all(!is.na(tbl$chisq))) {
    tbl$dCFI <- c(NA, diff(tbl$CFI)); tbl$dchisq <- c(NA, diff(tbl$chisq)); tbl$ddf <- c(NA, diff(tbl$df))
  }
  # 结构模型下 分性别 路径（从 configural 提取）
  fg <- fits$configural
  grp_eff <- NULL
  if (!is.null(fg) && lavaan::lavInspect(fg, "converged")) {
    pe <- lavaan::parameterEstimates(fg, standardized = TRUE)
    # 提取每个 group 的 a, b, c_prime, indirect（label 含 group 后缀由 lavaan 自动处理）
    # 这里直接取两个 group 的 std.all
    grp_eff <- data.frame(Tag = tag)
  }
  list(tbl = tbl)
}

inv_comp <- run_invariance("", "EF_COMP", "Composite")
inv_lat  <- run_invariance("EF =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST", "EF", "Latent")

cat("\n[Invariance sequence]\n")
print(round_df(rbind(inv_comp$tbl, inv_lat$tbl), 4), row.names = FALSE)
write.csv(rbind(inv_comp$tbl, inv_lat$tbl), file.path(results_dir, "DOWNSTREAM_invariance.csv"), row.names = FALSE)

sink()
cat("DONE\n")
