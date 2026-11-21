# ============================================================================
# code/00_load_and_prep.R
# Shared data loading, reverse scoring, z-standardization, and EF composite
# construction for all statistical analyses (N = 518).
#
# Expected input: data/data_analysis.xlsx (sheet = "Sheet1")
# Sex coding in raw data: 0 = Female (Girls, n = 238), 1 = Male (Boys, n = 280)
# ============================================================================

load_analysis_data <- function(path = "data/data_analysis.xlsx", sheet = "Sheet1") {
  if (!file.exists(path)) {
    stop(sprintf(
      "Input data file not found at '%s'. Run scripts from the repository root and place the anonymized dataset at '%s'.",
      path, path
    ))
  }
  df <- readxl::read_excel(path, sheet = sheet)

  # Reverse-score RT indicators so higher values denote stronger executive control
  df$Flanker_Rev <- -df$flanker_flanker_effect_ms
  df$SST_Rev     <- -df$stopsignal_ssrt_ms

  # z-standardize all 10 indicators + age across the analytic sample (M = 0, SD = 1)
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
  df$AGE_Z       <- as.numeric(scale(df$age))

  # Primary broad EF composite: unweighted mean of available standardized EF tasks (min 4 of 5)
  ef_mat <- df[, c("EF_AX", "EF_BACKWARD", "EF_WCST", "EF_FLANK", "EF_SST")]
  df$EF_COMP <- apply(ef_mat, 1, function(x) {
    if (sum(!is.na(x)) < 4L) return(NA_real_)
    mean(x, na.rm = TRUE)
  })

  # Exploratory M5 component composites (Section 2.4.4 / Section 3.6)
  comp_mean <- function(mat, min_n = 1L) {
    apply(mat, 1, function(x) if (sum(!is.na(x)) < min_n) NA_real_ else mean(x, na.rm = TRUE))
  }
  df$WM_pro <- comp_mean(cbind(df$EF_AX, df$EF_BACKWARD), min_n = 1L)
  df$FLEX   <- df$EF_WCST
  df$INH    <- comp_mean(cbind(df$EF_FLANK, df$EF_SST), min_n = 1L)

  # Sex factor: explicitly ordered Female (0) then Male (1)
  df$sex_group <- factor(
    ifelse(df$sex %in% c(0, 1), df$sex, NA),
    levels = c(0, 1),
    labels = c("Female", "Male")
  )

  df
}

round_df <- function(x, k = 4) {
  x[] <- lapply(x, function(col) if (is.numeric(col)) round(col, k) else col)
  x
}
