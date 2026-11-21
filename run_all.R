# ============================================================================
# run_all.R
# Master reproduction script for:
#   "Computational Fluency in the Association Between Executive Functions and
#    Academic Achievement in Chinese School-Aged Children" (BMC Psychology)
#
# Run from the repository root:
#   Rscript run_all.R
# ============================================================================

scripts <- c(
  "code/01_descriptives.R",
  "code/02_measurement_and_latent_sensitivity.R",
  "code/03_primary_mediation_M4_and_age.R",
  "code/04_crossdomain_equality_and_tost.R",
  "code/05_sensitivity_M5_components.R",
  "code/06_sex_invariance.R"
)

for (s in scripts) {
  cat("\n", strrep("=", 70), "\nRunning: ", s, "\n", strrep("=", 70), "\n", sep = "")
  source(s, local = new.env(parent = globalenv()))
}

cat("\nAll statistical analyses completed successfully.\n")
