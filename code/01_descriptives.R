# ============================================================================
# code/01_descriptives.R
# Reproduces Section 3.1, Table 1 (descriptive statistics and sex comparisons),
# and Table 2 (bivariate Pearson correlations).
# Note: Both Welch's t and Cohen's d are consistently oriented as Girls - Boys.
# ============================================================================

suppressPackageStartupMessages({
  library(readxl)
})
source("code/00_load_and_prep.R")

out_dir <- "results/01_descriptives"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "DESCRIPTIVES_console.log"), split = TRUE)

cat("01_descriptives.R | ", as.character(Sys.time()), "\n", sep = "")
df <- load_analysis_data()
cat("N =", nrow(df), "\n")
print(table(df$sex_group, useNA = "ifany"))

# Age summary (overall and by sex)
cat(sprintf(
  "\n[Age] Overall: M = %.3f, SD = %.3f, range = [%.2f, %.2f]\n",
  mean(df$age, na.rm = TRUE), sd(df$age, na.rm = TRUE),
  min(df$age, na.rm = TRUE), max(df$age, na.rm = TRUE)
))
age_t <- t.test(df$age[df$sex == 0], df$age[df$sex == 1])
cat(sprintf(
  "  Girls (n=%d): M = %.3f, SD = %.3f | Boys (n=%d): M = %.3f, SD = %.3f | Welch t = %.3f, p = %.4f\n",
  sum(df$sex == 0), mean(df$age[df$sex == 0], na.rm = TRUE), sd(df$age[df$sex == 0], na.rm = TRUE),
  sum(df$sex == 1), mean(df$age[df$sex == 1], na.rm = TRUE), sd(df$age[df$sex == 1], na.rm = TRUE),
  age_t$statistic, age_t$p.value
))

# Missing counts per indicator
raw_vars <- c(
  "axcpt_d_prime_context", "Backward_max_span", "wcst_wcst_completed_categories",
  "Flanker_Rev", "SST_Rev",
  "simple_subtraction_adjusted_score", "complex_subtraction_adjusted_score",
  "math_score", "chinese_score", "english_score"
)
labels <- c(
  "AX-CPT d'-context", "Backward digit span", "WCST categories",
  "Flanker (reversed, ms)", "SSRT (reversed, ms)",
  "Simple subtraction (adj.)", "Complex subtraction (adj.)",
  "Math score", "Chinese score", "English score"
)

cat("\n[Missing observations per indicator]\n")
print(sapply(raw_vars, function(v) sum(is.na(df[[v]]))))

# Table 1: Raw descriptive statistics & sex comparisons (Girls - Boys)
desc <- data.frame(
  Indicator = labels,
  M_total = NA_real_, SD_total = NA_real_,
  M_boys  = NA_real_, SD_boys  = NA_real_,
  M_girls = NA_real_, SD_girls = NA_real_,
  t_girls_minus_boys = NA_real_,
  p = NA_real_,
  d_girls_minus_boys = NA_real_,
  stringsAsFactors = FALSE
)

female <- df$sex == 0
male   <- df$sex == 1
for (i in seq_along(raw_vars)) {
  x <- df[[raw_vars[i]]]
  # Welch's t oriented as Girls minus Boys (matching Cohen's d sign)
  xt <- t.test(x[female], x[male])
  nf <- sum(!is.na(x[female])); nm <- sum(!is.na(x[male]))
  sf <- sd(x[female], na.rm = TRUE); sm <- sd(x[male], na.rm = TRUE)
  sp <- sqrt(((nf - 1) * sf^2 + (nm - 1) * sm^2) / (nf + nm - 2))
  d  <- (mean(x[female], na.rm = TRUE) - mean(x[male], na.rm = TRUE)) / sp

  desc$M_total[i]  <- mean(x, na.rm = TRUE)
  desc$SD_total[i] <- sd(x, na.rm = TRUE)
  desc$M_boys[i]   <- mean(x[male], na.rm = TRUE)
  desc$SD_boys[i]  <- sm
  desc$M_girls[i]  <- mean(x[female], na.rm = TRUE)
  desc$SD_girls[i] <- sf
  desc$t_girls_minus_boys[i] <- unname(xt$statistic)
  desc$p[i]                  <- xt$p.value
  desc$d_girls_minus_boys[i] <- d
}

desc_out <- round_df(desc, 3)
cat("\n[Table 1: Descriptive statistics and sex comparisons (Girls - Boys)]\n")
print(desc_out, row.names = FALSE)
write.csv(desc_out, file.path(out_dir, "Table1_descriptives.csv"), row.names = FALSE)

# Table 2: Bivariate Pearson correlation matrix (pairwise complete)
zv <- sapply(raw_vars, function(v) as.numeric(scale(df[[v]])))
colnames(zv) <- labels
cormat <- cor(zv, use = "pairwise.complete.obs")
cat("\n[Table 2: Bivariate Pearson correlations]\n")
print(round(cormat, 3))
write.csv(round(cormat, 4), file.path(out_dir, "Table2_correlations.csv"))

# Key internal correlations and two-task Spearman-Brown coefficient
r_sub <- cor(zv[, "Simple subtraction (adj.)"], zv[, "Complex subtraction (adj.)"], use = "pairwise.complete.obs")
sb_sub <- (2 * r_sub) / (1 + r_sub)
r_inh <- cor(zv[, "Flanker (reversed, ms)"], zv[, "SSRT (reversed, ms)"], use = "pairwise.complete.obs")
cat(sprintf("\nSimple ~ Complex subtraction: r = %.4f | Spearman-Brown = %.4f\n", r_sub, sb_sub))
cat(sprintf("Flanker (rev) ~ SSRT (rev):   r = %.4f\n", r_inh))

ef_cor <- cormat[1:5, 1:5]; diag(ef_cor) <- NA
cat(sprintf("EF intercorrelation range: [%.3f, %.3f]\n", min(ef_cor, na.rm = TRUE), max(ef_cor, na.rm = TRUE)))

sink()
cat("Completed 01_descriptives.R\n")
