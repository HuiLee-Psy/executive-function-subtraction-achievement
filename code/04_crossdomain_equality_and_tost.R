# ============================================================================
# code/04_crossdomain_equality_and_tost.R
# Reproduces:
#   - Section 3.4.1 & Table 5: Single-outcome SEMs (Math, L1 Chinese, L2 English)
#     with 2,000-resample percentile bootstrap 95% CIs (Composite + Latent EF)
#   - Section 3.4.2: Joint 3-outcome SEM LR equality tests (Composite & Latent),
#     pairwise standardized Wald tests with nominal AND Holm-adjusted p values,
#     and exploratory Monte Carlo TOST equivalence analysis (10,000 draws, |Δβ| = 0.10)
# ============================================================================

suppressPackageStartupMessages({
  library(lavaan)
  library(MASS)
  library(readxl)
})
source("code/00_load_and_prep.R")

out_dir <- "results/04_crossdomain"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "CROSSDOMAIN_console.log"), split = TRUE)

cat("04_crossdomain_equality_and_tost.R | ", as.character(Sys.time()),
    " | lavaan ", as.character(packageVersion("lavaan")), "\n", sep = "")
df <- load_analysis_data()

# ---- 1. Single-Outcome Replications & Joint LR Equality Test ----
run_crossdomain <- function(ef_block, predictor, tag, R_boot = 2000L) {
  cat("\n=== Cross-domain analysis:", tag, "===\n")
  outcomes <- c("AC_MATH", "AC_CHIN", "AC_ENG")
  rows <- list()

  for (k in seq_along(outcomes)) {
    onm <- outcomes[k]
    m <- paste0(
      ef_block, "\n",
      "Fluency =~ FL_SIMP + FL_COMP\n",
      "Fluency ~ a*", predictor, "\n",
      onm, " ~ b*Fluency + c_prime*", predictor, "\n",
      "indirect := a * b\n"
    )
    f <- lavaan::sem(m, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
    fm <- lavaan::fitMeasures(f, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"))
    pe <- lavaan::parameterEstimates(f, standardized = TRUE)
    gb <- function(l) pe$std.all[pe$label == l][1]
    r2 <- unname(lavaan::lavInspect(f, "rsquare")[onm])

    # Nonparametric percentile bootstrap (R = 2,000)
    set.seed(20260923 + k)
    n <- nrow(df); idx <- seq_len(n)
    bb <- numeric(R_boot)
    for (i in seq_len(R_boot)) {
      d_b <- df[sample(idx, n, replace = TRUE), , drop = FALSE]
      f_b <- tryCatch(
        suppressWarnings(lavaan::sem(m, data = d_b, std.lv = TRUE, estimator = "ML", missing = "fiml")),
        error = function(e) NULL
      )
      if (!is.null(f_b) && isTRUE(lavaan::lavInspect(f_b, "converged"))) {
        pe_b <- lavaan::parameterEstimates(f_b, standardized = TRUE)
        bb[i] <- pe_b$std.all[pe_b$label == "indirect"][1]
      } else {
        bb[i] <- NA_real_
      }
    }
    q <- unname(quantile(bb, c(0.025, 0.975), na.rm = TRUE))

    rows[[k]] <- data.frame(
      Tag = tag, Outcome = onm,
      chisq = unname(fm["chisq"]), df = unname(fm["df"]), p_chisq = unname(fm["pvalue"]),
      CFI = unname(fm["cfi"]), RMSEA = unname(fm["rmsea"]), SRMR = unname(fm["srmr"]),
      a_beta = gb("a"), b_beta = gb("b"), c_prime_beta = gb("c_prime"),
      indirect_beta = gb("indirect"), boot_CI_LL = q[1], boot_CI_UL = q[2],
      R2_outcome = r2, NonConv = sum(is.na(bb)),
      stringsAsFactors = FALSE
    )
  }
  single_tbl <- do.call(rbind, rows)

  # Joint 3-outcome SEM (free vs. equal fluency-to-outcome slopes b)
  m_free <- paste0(
    ef_block, "\n",
    "Fluency =~ FL_SIMP + FL_COMP\n",
    "Fluency ~ a*", predictor, "\n",
    "AC_MATH ~ b_Math*Fluency + c_Math*", predictor, "\n",
    "AC_CHIN ~ b_Chin*Fluency + c_Chin*", predictor, "\n",
    "AC_ENG  ~ b_Eng*Fluency  + c_Eng*", predictor, "\n",
    "ind_Math := a * b_Math\n",
    "ind_Chin := a * b_Chin\n",
    "ind_Eng  := a * b_Eng\n"
  )
  m_eq <- paste0(
    ef_block, "\n",
    "Fluency =~ FL_SIMP + FL_COMP\n",
    "Fluency ~ a*", predictor, "\n",
    "AC_MATH ~ b*Fluency + c_Math*", predictor, "\n",
    "AC_CHIN ~ b*Fluency + c_Chin*", predictor, "\n",
    "AC_ENG  ~ b*Fluency + c_Eng*", predictor, "\n"
  )
  f_free <- lavaan::sem(m_free, data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
  f_eq   <- lavaan::sem(m_eq,   data = df, std.lv = TRUE, estimator = "ML", missing = "fiml")
  lr     <- lavaan::anova(f_free, f_eq)

  cat("Joint LR equality test (free vs. equal b):\n")
  print(round_df(as.data.frame(lr), 4), row.names = FALSE)

  list(single = single_tbl, f_free = f_free, f_eq = f_eq, lr = lr)
}

cd_comp <- run_crossdomain("", "EF_COMP", "Composite", R_boot = 2000L)
cd_lat  <- run_crossdomain("EF =~ EF_AX + EF_BACKWARD + EF_WCST + EF_FLANK + EF_SST", "EF", "Latent", R_boot = 2000L)

t5_all <- round_df(rbind(cd_comp$single, cd_lat$single), 4)
cat("\n[Table 5: Single-outcome domain replications (Composite & Latent)]\n")
print(t5_all, row.names = FALSE)
write.csv(t5_all, file.path(out_dir, "Table5_crossdomain_single.csv"), row.names = FALSE)

# ---- 2. Pairwise Standardized Wald Contrasts (Nominal + Holm) & Exploratory TOST ----
f_free <- cd_comp$f_free
pe_std <- lavaan::standardizedSolution(f_free, type = "std.all")
get_s  <- function(lbl) pe_std$est.std[pe_std$label == lbl][1]

V_std <- lavaan::lavInspect(f_free, "vcov.std.all")
need  <- c("a", "b_Math", "b_Chin", "b_Eng")
Vs    <- V_std[need, need, drop = FALSE]

wald_pair <- function(n1, n2) {
  ai <- get_s("a"); bi <- get_s(paste0("b_", n1)); bj <- get_s(paste0("b_", n2))
  g  <- c(bi - bj, ai, -ai)
  nm <- c("a", paste0("b_", n1), paste0("b_", n2))
  Vw <- Vs[nm, nm, drop = FALSE]
  d  <- ai * bi - ai * bj
  se <- sqrt(as.numeric(t(g) %*% Vw %*% g))
  z  <- d / se
  p  <- 2 * (1 - pnorm(abs(z)))
  data.frame(Pair = paste0(n1, "-", n2), dbeta = d, SE = se, z = z, p_nominal = p, stringsAsFactors = FALSE)
}

wald_tbl <- rbind(
  wald_pair("Math", "Chin"),
  wald_pair("Math", "Eng"),
  wald_pair("Chin", "Eng")
)
wald_tbl$p_holm <- p.adjust(wald_tbl$p_nominal, method = "holm")

cat("\n[Pairwise standardized Wald contrasts with Holm correction]\n")
print(round_df(wald_tbl, 4), row.names = FALSE)
write.csv(round_df(wald_tbl, 4), file.path(out_dir, "Crossdomain_Wald_Holm.csv"), row.names = FALSE)

# Exploratory TOST equivalence check (Monte Carlo 10,000 draws, |Δβ| = 0.10, 90% CI)
set.seed(20260531)
mu    <- c(a = get_s("a"), b_Math = get_s("b_Math"), b_Chin = get_s("b_Chin"), b_Eng = get_s("b_Eng"))
draws <- MASS::mvrnorm(10000L, mu = mu, Sigma = Vs)
ind_d <- cbind(draws[, "a"] * draws[, "b_Math"],
               draws[, "a"] * draws[, "b_Chin"],
               draws[, "a"] * draws[, "b_Eng"])
diffs <- cbind(ind_d[, 1] - ind_d[, 2], ind_d[, 1] - ind_d[, 3], ind_d[, 2] - ind_d[, 3])
colnames(diffs) <- c("Math-Chin", "Math-Eng", "Chin-Eng")

tost_tbl <- data.frame(
  Pair     = colnames(diffs),
  dbeta    = wald_tbl$dbeta,
  MC_Mean  = colMeans(diffs),
  CI90_LL  = apply(diffs, 2, function(x) unname(quantile(x, 0.05))),
  CI90_UL  = apply(diffs, 2, function(x) unname(quantile(x, 0.95))),
  Within_0.10 = apply(diffs, 2, function(x) all(abs(quantile(x, c(0.05, 0.95))) < 0.10)),
  stringsAsFactors = FALSE
)
cat("\n[Exploratory TOST standardized equivalence check (10,000 MC draws, bounds = ±0.10)]\n")
print(round_df(tost_tbl, 4), row.names = FALSE)
write.csv(round_df(tost_tbl, 4), file.path(out_dir, "Crossdomain_TOST_std.csv"), row.names = FALSE)

sink()
cat("Completed 04_crossdomain_equality_and_tost.R\n")
