# Executive Functions, Computational Fluency, and Academic Achievement

Reproducible statistical analysis code for the manuscript:

> **Computational Fluency in the Association Between Executive Functions and Academic Achievement in Chinese School-Aged Children**

## Repository Structure

```text
.
├── README.md
├── run_all.R                                     # Master script to run all analyses in sequence
└── code/
    ├── 00_load_and_prep.R                        # Shared data loading, reverse scoring, z-standardization, EF composites
    ├── 01_descriptives.R                         # Section 3.1, Table 1 (Girls - Boys Welch t & Cohen's d), Table 2
    ├── 02_measurement_and_latent_sensitivity.R   # Section 3.2, Table 3, omega/AVE, M2/M3 checks, M0/M1 latent sensitivity (Table S1)
    ├── 03_primary_mediation_M4_and_age.R         # Section 3.3, Table 4 (5,000-resample bootstrap), age-adjusted M4 (Table S1)
    ├── 04_crossdomain_equality_and_tost.R        # Section 3.4, Table 5, joint LR test, Holm-adjusted Wald tests, Monte Carlo TOST
    ├── 05_sensitivity_M5_components.R            # Section 3.6, Table S1 (M5 parallel mediation, 5,000-resample bootstrap)
    └── 06_sex_invariance.R                       # Section 3.5, Tables S2 & S3 (multigroup invariance with explicit group order)
```

## Software Requirements

Analyses were conducted in **R 4.5.2** using the following CRAN packages:

- `lavaan` (0.6-21 or later)
- `semTools` (0.5-8)
- `MASS`
- `readxl`

## Data Availability

Individual-level data are not publicly deposited because the dataset links cognitive task performance to school examination grades of minor primary-school students. De-identified data may be made available from the corresponding author upon reasonable request and subject to institutional ethics approval. Once obtained, place the dataset at `data/data_analysis.xlsx` (sheet `"Sheet1"`) and run:

```bash
Rscript run_all.R
```

## Key Updates in This Release (v2.0.0 — Submission Version)

1. **Explicit sex-group ordering and partial-scalar contrasts (`code/06_sex_invariance.R`)**: Fixed `group.label = c("Female", "Male")` in `lavaan::sem()` to prevent default ordering by first row appearance, and extracted group-specific slopes, standardized indirect-effect delta-method 95% CIs, achievement $R^2$, and full-covariance Wald contrasts directly under the partial scalar model (Tables S2–S3).
2. **Age-adjusted sensitivity analysis (`code/03_primary_mediation_M4_and_age.R`)**: Added standardized student age as a covariate predicting both computational fluency and academic achievement in Model M4 (Table S1).
3. **Holm-adjusted cross-domain contrasts (`code/04_crossdomain_equality_and_tost.R`)**: Added post hoc Holm correction across the three standardized pairwise Wald contrasts alongside the joint LR equality test and exploratory TOST check.
4. **Consistent descriptive sign orientation (`code/01_descriptives.R`)**: Oriented both Welch's $t$ and Cohen's $d$ consistently as Girls minus Boys (Table 1) and fixed pairwise-complete handling for internal correlations.
5. **Removed unused packages and plotting scripts**: Removed unused package declarations (`simsem`) and consolidated all statistical analyses into six modular scripts (01–06).
