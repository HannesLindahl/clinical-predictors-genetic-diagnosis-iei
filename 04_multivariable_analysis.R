# ==============================================================================
# 04_multivariable_analysis.R
#
# Study:
#   Clinical predictors of genetic diagnoses in adults evaluated for
#   inborn errors of immunity
#
# Purpose:
#   Fit the primary multivariable logistic regression model for achieving a
#   molecular diagnosis and generate the numerical results reported in Table 2
#   and Figure 4.
#
# Input:
#   data/derived/analysis_data.rds
#
# Outputs:
#   data/derived/primary_model.rds
#   results/models/primary_model_results.csv
#
# Notes:
#   The preceding variable-prioritization analysis was exploratory and was not
#   intended to develop or validate a clinical prediction model.
#
#   The final primary model includes age at symptom onset, sex, serum IgG1
#   concentration, and eczema.
#
#   Age at symptom onset and IgG1 are modelled as continuous variables.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tibble)
library(readr)


# ------------------------------------------------------------------------------
# 2. Load analysis dataset
# ------------------------------------------------------------------------------

dat <- readRDS(
  "data/derived/analysis_data.rds"
)


# ==============================================================================
# PART A. PREPARE PRIMARY MODEL DATA
# ==============================================================================


# ------------------------------------------------------------------------------
# 3. Define variables
# ------------------------------------------------------------------------------

primary_model_vars <- c(
  "molecular_diagnosis",
  "age_onset",
  "male",
  "igg1",
  "eczema_binary"
)


stopifnot(
  all(primary_model_vars %in% names(dat))
)


# ------------------------------------------------------------------------------
# 4. Construct complete-case analysis dataset
# ------------------------------------------------------------------------------

primary_data <- dat %>%
  select(all_of(primary_model_vars)) %>%
  filter(complete.cases(.))


cat("\nPRIMARY MODEL ANALYSIS POPULATION\n")
cat("---------------------------------\n")
cat("Complete cases:", nrow(primary_data), "\n")
cat(
  "Molecular diagnoses:",
  sum(primary_data$molecular_diagnosis == 1),
  "\n"
)
cat(
  "No molecular diagnosis:",
  sum(primary_data$molecular_diagnosis == 0),
  "\n\n"
)


# Both outcome classes must be represented.

stopifnot(
  length(unique(primary_data$molecular_diagnosis)) == 2
)


# ==============================================================================
# PART B. PRIMARY MULTIVARIABLE LOGISTIC REGRESSION
# ==============================================================================


# ------------------------------------------------------------------------------
# 5. Fit primary model
# ------------------------------------------------------------------------------

primary_model <- glm(
  molecular_diagnosis ~
    age_onset +
    male +
    igg1 +
    eczema_binary,
  data = primary_data,
  family = binomial(link = "logit")
)


# Basic convergence check.

stopifnot(primary_model$converged)


cat("\nPRIMARY MODEL\n")
cat("-------------\n")

print(
  summary(primary_model)
)


# ==============================================================================
# PART C. EXTRACT MODEL RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 6. Obtain odds ratios, confidence intervals, and p-values
# ------------------------------------------------------------------------------

# Confidence intervals are calculated using profile likelihood.
#
# These may differ slightly from Wald confidence intervals based directly on
# the standard errors in summary(primary_model).

coef_table <- summary(primary_model)$coefficients


ci_logit <- suppressMessages(
  confint(primary_model)
)


primary_results <- tibble(
  term = rownames(coef_table),
  estimate_log_odds = coef_table[, "Estimate"],
  standard_error = coef_table[, "Std. Error"],
  p_value = coef_table[, "Pr(>|z|)"],
  conf_low_log_odds = ci_logit[, 1],
  conf_high_log_odds = ci_logit[, 2]
) %>%
  mutate(
    odds_ratio = exp(estimate_log_odds),
    conf_low = exp(conf_low_log_odds),
    conf_high = exp(conf_high_log_odds)
  )


# ------------------------------------------------------------------------------
# 7. Add manuscript-friendly variable labels
# ------------------------------------------------------------------------------

primary_results <- primary_results %>%
  mutate(
    variable = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "age_onset" ~ "Age at symptom onset (years)",
      term == "male" ~ "Male sex",
      term == "igg1" ~ "IgG1 (g/L)",
      term == "eczema_binary" ~ "Eczema",
      TRUE ~ term
    )
  ) %>%
  select(
    variable,
    term,
    estimate_log_odds,
    standard_error,
    odds_ratio,
    conf_low,
    conf_high,
    p_value
  )


cat("\nPRIMARY MODEL RESULTS\n")
cat("---------------------\n")

print(
  primary_results,
  n = Inf
)


# ==============================================================================
# PART D. CHECK AGAINST REPORTED RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 8. Compare estimates with manuscript Table 2
# ------------------------------------------------------------------------------

# The following values are the adjusted odds ratios reported in the manuscript:
#
#   Male sex:       2.872
#   IgG1:           1.263 per 1 g/L
#   Eczema:         6.234
#   Age at onset:   0.93 per year
#
# Small numerical differences can arise from rounding and the method used to
# calculate confidence intervals.


reported_or <- c(
  male = 2.872,
  igg1 = 1.263,
  eczema_binary = 6.234,
  age_onset = 0.93
)


observed_or <- primary_results %>%
  filter(term != "(Intercept)") %>%
  select(term, odds_ratio)


comparison <- tibble(
  term = names(reported_or),
  reported_or = unname(reported_or)
) %>%
  left_join(
    observed_or,
    by = "term"
  ) %>%
  mutate(
    absolute_difference =
      abs(odds_ratio - reported_or)
  )


cat("\nCOMPARISON WITH MANUSCRIPT TABLE 2\n")
cat("----------------------------------\n")

print(comparison)


# Do not require exact equality because published values are rounded.

stopifnot(
  all(comparison$absolute_difference < 0.02)
)


# ==============================================================================
# PART E. MODEL DIAGNOSTICS
# ==============================================================================


# ------------------------------------------------------------------------------
# 9. Check coefficient stability and possible separation
# ------------------------------------------------------------------------------

# Extremely large coefficients or standard errors can indicate sparse-data
# problems or separation. No automatic exclusion is performed here; this is
# an audit check.

diagnostic_table <- tibble(
  term = rownames(coef_table),
  coefficient = coef_table[, "Estimate"],
  standard_error = coef_table[, "Std. Error"]
) %>%
  mutate(
    abs_coefficient = abs(coefficient)
  )


cat("\nCOEFFICIENT DIAGNOSTIC\n")
cat("----------------------\n")

print(diagnostic_table)


if (
  any(
    !is.finite(diagnostic_table$coefficient) |
    !is.finite(diagnostic_table$standard_error)
  )
) {
  
  warning(
    "Non-finite coefficient or standard error detected."
  )
}


# ------------------------------------------------------------------------------
# 10. Examine correlations among continuous predictors
# ------------------------------------------------------------------------------

continuous_predictors <- primary_data %>%
  select(
    age_onset,
    igg1
  )


continuous_correlations <- cor(
  continuous_predictors,
  use = "complete.obs",
  method = "spearman"
)


cat("\nSPEARMAN CORRELATION BETWEEN CONTINUOUS PREDICTORS\n")
cat("--------------------------------------------------\n")

print(continuous_correlations)


# ------------------------------------------------------------------------------
# 11. Variance inflation factors
# ------------------------------------------------------------------------------

# Calculate VIFs directly to avoid adding another package dependency.
#
# For each predictor, VIF = 1 / (1 - R^2), where R^2 is obtained by regressing
# that predictor on all remaining predictors.
#
# This is straightforward here because all predictors in the final model are
# represented by a single numeric column.

x_vif <- model.matrix(
  ~ age_onset + male + igg1 + eczema_binary,
  data = primary_data
)[, -1, drop = FALSE]


calculate_vif <- function(x) {
  
  sapply(
    seq_len(ncol(x)),
    function(j) {
      
      response <- x[, j]
      others <- x[, -j, drop = FALSE]
      
      fit <- lm(
        response ~ others
      )
      
      r_squared <- summary(fit)$r.squared
      
      1 / (1 - r_squared)
    }
  )
}


vif_values <- calculate_vif(x_vif)

names(vif_values) <- colnames(x_vif)


vif_table <- tibble(
  variable = names(vif_values),
  vif = as.numeric(vif_values)
)


cat("\nVARIANCE INFLATION FACTORS\n")
cat("--------------------------\n")

print(vif_table)


# ==============================================================================
# PART F. SAVE RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 13. Create output directory
# ------------------------------------------------------------------------------

dir.create(
  "results/models",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 14. Save model object and result tables
# ------------------------------------------------------------------------------

saveRDS(
  primary_model,
  "data/derived/primary_model.rds"
)


write_csv(
  primary_results,
  "results/models/primary_model_results.csv"
)


write_csv(
  comparison,
  "results/models/primary_model_manuscript_comparison.csv"
)


write_csv(
  vif_table,
  "results/models/primary_model_vif.csv"
)


message("Primary multivariable analysis completed successfully.")
message("Complete-case sample: ", nrow(primary_data))