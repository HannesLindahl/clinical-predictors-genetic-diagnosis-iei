# ==============================================================================
# 06_sensitivity_analyses.R
#
# Study:
#   Clinical predictors of genetic diagnoses in adults evaluated for
#   inborn errors of immunity
#
# Purpose:
#   Reproduce the sensitivity analyses reported in Tables S6 and S7 by fitting
#   the primary multivariable logistic regression model separately among
#   participants with and without predominantly antibody deficiencies (PAD).
#
# Input:
#   data/derived/analysis_data.rds
#
# Outputs:
#   results/models/non_pad_model_results.csv
#   results/models/pad_model_results.csv
#
# Notes:
#   These are subgroup sensitivity analyses of the primary model:
#
#     molecular diagnosis ~ age at onset + male sex + IgG1 + eczema
#
#   Models are fitted using complete cases within each subgroup.
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
# PART A. CHECK PAD CLASSIFICATION
# ==============================================================================


# ------------------------------------------------------------------------------
# 3. Required variables
# ------------------------------------------------------------------------------

required_vars <- c(
  "molecular_diagnosis",
  "age_onset",
  "male",
  "igg1",
  "eczema_binary",
  "pad"
)


stopifnot(
  all(required_vars %in% names(dat))
)


# ------------------------------------------------------------------------------
# 4. Summarize PAD classification
# ------------------------------------------------------------------------------

pad_summary <- dat %>%
  count(
    pad,
    molecular_diagnosis,
    name = "n"
  )


cat("\nPAD CLASSIFICATION\n")
cat("------------------\n")

print(pad_summary)


cat("\nPAD GROUP SIZES\n")
cat("---------------\n")

print(
  dat %>%
    count(
      pad,
      name = "n"
    )
)


# ==============================================================================
# PART B. DEFINE SUBGROUP MODEL FUNCTION
# ==============================================================================


# ------------------------------------------------------------------------------
# 5. Function for subgroup logistic regression
# ------------------------------------------------------------------------------

fit_subgroup_model <- function(data, subgroup_name) {
  
  
  # --------------------------------------------------------------------------
  # Complete-case analysis
  # --------------------------------------------------------------------------
  
  model_data <- data %>%
    select(
      molecular_diagnosis,
      age_onset,
      male,
      igg1,
      eczema_binary
    ) %>%
    filter(
      complete.cases(.)
    )
  
  
  cat("\n", subgroup_name, "\n", sep = "")
  cat(
    paste0(
      rep("-", nchar(subgroup_name)),
      collapse = ""
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Complete cases:",
    nrow(model_data),
    "\n"
  )
  
  cat(
    "Molecular diagnoses:",
    sum(model_data$molecular_diagnosis == 1),
    "\n"
  )
  
  cat(
    "No molecular diagnosis:",
    sum(model_data$molecular_diagnosis == 0),
    "\n"
  )
  
  
  # Both outcome categories must be represented.
  
  stopifnot(
    length(
      unique(model_data$molecular_diagnosis)
    ) == 2
  )
  
  
  # --------------------------------------------------------------------------
  # Fit model
  # --------------------------------------------------------------------------
  
  fit <- glm(
    molecular_diagnosis ~
      age_onset +
      male +
      igg1 +
      eczema_binary,
    data = model_data,
    family = binomial(link = "logit")
  )
  
  
  if (!fit$converged) {
    warning(
      subgroup_name,
      " model did not converge."
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Extract results
  # --------------------------------------------------------------------------
  
  coef_table <- summary(fit)$coefficients
  
  
  ci_logit <- suppressMessages(
    confint(fit)
  )
  
  
  results <- tibble(
    term = rownames(coef_table),
    estimate_log_odds =
      coef_table[, "Estimate"],
    standard_error =
      coef_table[, "Std. Error"],
    p_value =
      coef_table[, "Pr(>|z|)"],
    odds_ratio =
      exp(coef_table[, "Estimate"]),
    conf_low =
      exp(ci_logit[, 1]),
    conf_high =
      exp(ci_logit[, 2])
  ) %>%
    mutate(
      variable = case_when(
        term == "(Intercept)" ~
          "Intercept",
        
        term == "age_onset" ~
          "Age at symptom onset (years)",
        
        term == "male" ~
          "Male sex",
        
        term == "igg1" ~
          "IgG1 (g/L)",
        
        term == "eczema_binary" ~
          "Eczema",
        
        TRUE ~ term
      ),
      
      subgroup = subgroup_name
    ) %>%
    select(
      subgroup,
      variable,
      term,
      estimate_log_odds,
      standard_error,
      odds_ratio,
      conf_low,
      conf_high,
      p_value
    )
  
  
  # --------------------------------------------------------------------------
  # Return model, data and results
  # --------------------------------------------------------------------------
  
  list(
    data = model_data,
    fit = fit,
    results = results
  )
}


# ==============================================================================
# PART C. NON-PAD SENSITIVITY ANALYSIS
# ==============================================================================


# ------------------------------------------------------------------------------
# 6. Non-PAD subgroup
# ------------------------------------------------------------------------------

non_pad_data <- dat %>%
  filter(
    pad == 0
  )


non_pad_model <- fit_subgroup_model(
  data = non_pad_data,
  subgroup_name = "Non-PAD"
)


cat("\nNON-PAD MODEL RESULTS\n")
cat("---------------------\n")

print(
  non_pad_model$results,
  n = Inf
)


# ==============================================================================
# PART D. PAD SENSITIVITY ANALYSIS
# ==============================================================================


# ------------------------------------------------------------------------------
# 7. PAD subgroup
# ------------------------------------------------------------------------------

pad_data <- dat %>%
  filter(
    pad == 1
  )


pad_model <- fit_subgroup_model(
  data = pad_data,
  subgroup_name = "PAD"
)


cat("\nPAD MODEL RESULTS\n")
cat("-----------------\n")

print(
  pad_model$results,
  n = Inf
)


# ==============================================================================
# PART E. COMBINE RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 8. Combined sensitivity-analysis table
# ------------------------------------------------------------------------------

sensitivity_results <- bind_rows(
  non_pad_model$results,
  pad_model$results
)


cat("\nCOMBINED SENSITIVITY ANALYSIS\n")
cat("-----------------------------\n")

print(
  sensitivity_results,
  n = Inf
)


# ==============================================================================
# PART F. COMPARE WITH REPORTED NON-PAD RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 9. Table S6 comparison
# ------------------------------------------------------------------------------

# Reported adjusted odds ratios from Table S6:
#
#   Male sex:       1.77
#   IgG1:           1.38
#   Eczema:         7.31
#   Age at onset:   0.92


reported_non_pad_or <- c(
  male = 1.77,
  igg1 = 1.38,
  eczema_binary = 7.31,
  age_onset = 0.92
)


observed_non_pad_or <- non_pad_model$results %>%
  filter(
    term != "(Intercept)"
  ) %>%
  select(
    term,
    odds_ratio
  )


non_pad_comparison <- tibble(
  term = names(reported_non_pad_or),
  reported_or =
    unname(reported_non_pad_or)
) %>%
  left_join(
    observed_non_pad_or,
    by = "term"
  ) %>%
  mutate(
    absolute_difference =
      abs(
        odds_ratio -
          reported_or
      )
  )


cat("\nCOMPARISON WITH TABLE S6\n")
cat("------------------------\n")

print(non_pad_comparison)


# ==============================================================================
# PART G. COMPARE WITH REPORTED PAD RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 10. Table S7 comparison
# ------------------------------------------------------------------------------

# Reported adjusted odds ratios from Table S7:
#
#   Male sex:       10.3
#   IgG1:           1.64
#   Eczema:         8.03
#   Age at onset:   0.92


reported_pad_or <- c(
  male = 10.3,
  igg1 = 1.64,
  eczema_binary = 8.03,
  age_onset = 0.92
)


observed_pad_or <- pad_model$results %>%
  filter(
    term != "(Intercept)"
  ) %>%
  select(
    term,
    odds_ratio
  )


pad_comparison <- tibble(
  term = names(reported_pad_or),
  reported_or =
    unname(reported_pad_or)
) %>%
  left_join(
    observed_pad_or,
    by = "term"
  ) %>%
  mutate(
    absolute_difference =
      abs(
        odds_ratio -
          reported_or
      )
  )


cat("\nCOMPARISON WITH TABLE S7\n")
cat("------------------------\n")

print(pad_comparison)


# ==============================================================================
# PART H. SAVE RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 11. Create output directory
# ------------------------------------------------------------------------------

dir.create(
  "results/models",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 12. Save model objects
# ------------------------------------------------------------------------------

saveRDS(
  non_pad_model$fit,
  "data/derived/non_pad_model.rds"
)


saveRDS(
  pad_model$fit,
  "data/derived/pad_model.rds"
)


# ------------------------------------------------------------------------------
# 13. Save numerical results
# ------------------------------------------------------------------------------

write_csv(
  non_pad_model$results,
  "results/models/non_pad_model_results.csv"
)


write_csv(
  pad_model$results,
  "results/models/pad_model_results.csv"
)


write_csv(
  sensitivity_results,
  "results/models/sensitivity_analysis_results.csv"
)


write_csv(
  non_pad_comparison,
  "results/models/non_pad_manuscript_comparison.csv"
)


write_csv(
  pad_comparison,
  "results/models/pad_manuscript_comparison.csv"
)


message("Sensitivity analyses completed successfully.")
message(
  "Non-PAD complete-case sample: ",
  nrow(non_pad_model$data)
)
message(
  "PAD complete-case sample: ",
  nrow(pad_model$data)
)