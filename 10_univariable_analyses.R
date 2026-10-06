# ==============================================================================
# 10_univariable_analyses.R
#
# Study:
#   Clinical predictors of genetic diagnoses in adults evaluated for
#   inborn errors of immunity
#
# Purpose:
#   Reproduce the univariable logistic regression analyses underlying
#   Supplementary Table S5.
#
#   The analysis is restricted to the 79 predictors that remained eligible
#   after the preprocessing, missingness, and no-variability exclusions
#   reconstructed in 03_variable_selection.R.
#
# Historical analysis:
#   - One predictor was evaluated at a time.
#   - Ordinary logistic regression was used for most predictors.
#   - Firth penalized logistic regression was used for celiac disease and
#     porto-sinusoidal vascular disease because of separation.
#   - Ventricular malignancy had only one positive observation and was
#     reported as non-estimable (NE) in the historical supplementary table.
#   - Each analysis used observations with non-missing outcome and predictor.
#
# Input:
#   data/derived/analysis_data.rds
#   data/derived/eligible_predictors.rds
#
# Outputs:
#   results/univariable/univariable_results.csv
#   results/univariable/univariable_model_diagnostics.csv
#
# Note:
#   Manuscript-facing labels and domain headings are handled separately from
#   the statistical analysis.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tibble)
library(readr)
library(logistf)


# ------------------------------------------------------------------------------
# 2. Load analysis data and eligible predictor list
# ------------------------------------------------------------------------------

dat <- readRDS(
  "data/derived/analysis_data.rds"
)

eligible_predictors <- readRDS(
  "data/derived/eligible_predictors.rds"
)


stopifnot(
  "molecular_diagnosis" %in% names(dat),
  length(eligible_predictors) == 79,
  length(unique(eligible_predictors)) == 79,
  all(eligible_predictors %in% names(dat))
)


# ==============================================================================
# PART A. HISTORICAL ANALYSIS RULES
# ==============================================================================


# ------------------------------------------------------------------------------
# 3. Historical handling of separation
#
# Review of the original analysis code showed that Firth penalized logistic
# regression was applied manually to celiac disease and porto-sinusoidal
# vascular disease.
#
# Ventricular malignancy had only one positive observation and was reported
# as non-estimable (NE) rather than being refitted using Firth regression.
# ------------------------------------------------------------------------------

firth_variables <- c(
  "Celiaki",
  "NRH/portosinusoidal vascular disease"
)

nonestimable_variables <- c(
  "Ventricular malignancy"
)


stopifnot(
  all(firth_variables %in% eligible_predictors),
  all(nonestimable_variables %in% eligible_predictors)
)


# ==============================================================================
# PART B. HELPER FUNCTIONS
# ==============================================================================


# ------------------------------------------------------------------------------
# 4. Identify binary predictors
# ------------------------------------------------------------------------------

is_binary_predictor <- function(x) {
  
  x_nonmissing <- x[!is.na(x)]
  
  length(unique(x_nonmissing)) == 2
}


# ------------------------------------------------------------------------------
# 5. Convert binary predictors to 0/1
#
# This function preserves numeric 0/1 variables and handles the encodings
# present in the historical dataset.
# ------------------------------------------------------------------------------

binary_to_01 <- function(x) {
  
  # Numeric 0/1
  
  if (is.numeric(x)) {
    
    vals <- sort(
      unique(
        x[!is.na(x)]
      )
    )
    
    if (
      length(vals) == 2 &&
      all(vals %in% c(0, 1))
    ) {
      return(x)
    }
  }
  
  
  x_chr <- trimws(
    as.character(x)
  )
  
  
  vals <- unique(
    x_chr[!is.na(x_chr)]
  )
  
  
  # Y / N
  
  if (
    all(vals %in% c("Y", "N"))
  ) {
    
    return(
      ifelse(
        is.na(x_chr),
        NA_real_,
        ifelse(
          x_chr == "Y",
          1,
          0
        )
      )
    )
  }
  
  
  # Yes / No
  
  if (
    all(vals %in% c("Yes", "No"))
  ) {
    
    return(
      ifelse(
        is.na(x_chr),
        NA_real_,
        ifelse(
          x_chr == "Yes",
          1,
          0
        )
      )
    )
  }
  
  
  # Sex: female = 0, male = 1
  
  if (
    all(vals %in% c("M", "F"))
  ) {
    
    return(
      ifelse(
        is.na(x_chr),
        NA_real_,
        ifelse(
          x_chr == "M",
          1,
          0
        )
      )
    )
  }
  
  
  # Character 0 / 1
  
  if (
    all(vals %in% c("0", "1"))
  ) {
    
    return(
      as.numeric(x_chr)
    )
  }
  
  
  stop(
    paste(
      "Unrecognized binary encoding:",
      paste(vals, collapse = ", ")
    )
  )
}


# ==============================================================================
# PART C. UNIVARIABLE ANALYSIS FUNCTION
# ==============================================================================


# ------------------------------------------------------------------------------
# 6. Analyze one predictor
# ------------------------------------------------------------------------------

run_univariable <- function(variable, data) {
  
  x <- data[[variable]]
  
  
  # --------------------------------------------------------------------------
  # Determine whether predictor is binary or continuous
  # --------------------------------------------------------------------------
  
  binary <- is_binary_predictor(x)
  
  
  if (binary) {
    
    x_analysis <- binary_to_01(x)
    predictor_type <- "binary"
    
  } else {
    
    x_analysis <- suppressWarnings(
      as.numeric(
        as.character(x)
      )
    )
    
    predictor_type <- "continuous"
    
    
    # Stop rather than silently introduce missing values through conversion.
    
    original_nonmissing <- sum(
      !is.na(x)
    )
    
    converted_nonmissing <- sum(
      !is.na(x_analysis)
    )
    
    
    if (
      converted_nonmissing != original_nonmissing
    ) {
      
      stop(
        paste(
          "Numeric conversion introduced missing values for:",
          variable
        )
      )
    }
  }
  
  
  # --------------------------------------------------------------------------
  # Complete-case data for this predictor
  # --------------------------------------------------------------------------
  
  analysis_data <- tibble(
    outcome = data$molecular_diagnosis,
    predictor = x_analysis
  ) %>%
    filter(
      !is.na(outcome),
      !is.na(predictor)
    )
  
  
  n_analyzed <- nrow(
    analysis_data
  )
  
  
  if (
    n_analyzed == 0
  ) {
    
    stop(
      paste(
        "No observations available for:",
        variable
      )
    )
  }
  
  
  # Number positive is reported only for binary characteristics.
  
  n_positive <- if (binary) {
    
    sum(
      analysis_data$predictor == 1
    )
    
  } else {
    
    NA_integer_
  }
  
  
  # --------------------------------------------------------------------------
  # Historical non-estimable variable
  #
  # Ventricular malignancy was reported as NE in the historical Table S5.
  # --------------------------------------------------------------------------
  
  if (
    variable %in% nonestimable_variables
  ) {
    
    return(
      tibble(
        variable = variable,
        predictor_type = predictor_type,
        n_analyzed = n_analyzed,
        n_positive = n_positive,
        method = "Not estimable in historical analysis",
        odds_ratio = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p_value = NA_real_,
        glm_coefficient = NA_real_,
        glm_standard_error = NA_real_
      )
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Historical Firth logistic regression
  #
  # Applied manually in the original analysis to:
  #   - Celiac disease
  #   - Porto-sinusoidal vascular disease
  # --------------------------------------------------------------------------
  
  if (
    variable %in% firth_variables
  ) {
    
    firth_fit <- logistf(
      outcome ~ predictor,
      data = analysis_data,
      pl = TRUE
    )
    
    
    beta <- firth_fit$coefficients[
      "predictor"
    ]
    
    
    return(
      tibble(
        variable = variable,
        predictor_type = predictor_type,
        n_analyzed = n_analyzed,
        n_positive = n_positive,
        method = "Firth logistic regression",
        odds_ratio = exp(beta),
        conf_low = exp(
          firth_fit$ci.lower[
            "predictor"
          ]
        ),
        conf_high = exp(
          firth_fit$ci.upper[
            "predictor"
          ]
        ),
        p_value = firth_fit$prob[
          "predictor"
        ],
        glm_coefficient = NA_real_,
        glm_standard_error = NA_real_
      )
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Standard logistic regression
  # --------------------------------------------------------------------------
  
  glm_fit <- glm(
    outcome ~ predictor,
    data = analysis_data,
    family = binomial(link = "logit")
  )
  
  
  coef_table <- summary(
    glm_fit
  )$coefficients
  
  
  beta <- coef(
    glm_fit
  )[["predictor"]]
  
  
  p_value <- coef_table[
    "predictor",
    "Pr(>|z|)"
  ]
  
  
  # Profile-likelihood confidence interval.
  
  ci_logit <- suppressMessages(
    confint(
      glm_fit,
      parm = "predictor",
      level = 0.95
    )
  )
  
  
  return(
    tibble(
      variable = variable,
      predictor_type = predictor_type,
      n_analyzed = n_analyzed,
      n_positive = n_positive,
      method = "Standard logistic regression",
      odds_ratio = exp(beta),
      conf_low = exp(ci_logit[1]),
      conf_high = exp(ci_logit[2]),
      p_value = p_value,
      glm_coefficient = beta,
      glm_standard_error = coef_table[
        "predictor",
        "Std. Error"
      ]
    )
  )
}


# ==============================================================================
# PART D. RUN ALL 79 UNIVARIABLE ANALYSES
# ==============================================================================


# ------------------------------------------------------------------------------
# 7. Run analyses
# ------------------------------------------------------------------------------

cat("\nRunning univariable analyses...\n")


univariable_results <- bind_rows(
  lapply(
    eligible_predictors,
    run_univariable,
    data = dat
  )
)


stopifnot(
  nrow(univariable_results) == 79
)


cat(
  "\nUnivariable analyses completed:",
  nrow(univariable_results),
  "predictors\n"
)


# ------------------------------------------------------------------------------
# 8. Verify that all 79 eligible predictors are represented
# ------------------------------------------------------------------------------

stopifnot(
  setequal(
    univariable_results$variable,
    eligible_predictors
  )
)


# ==============================================================================
# PART E. QUALITY-CONTROL CHECKS
# ==============================================================================


# ------------------------------------------------------------------------------
# 9. Regression methods used
# ------------------------------------------------------------------------------

cat("\nREGRESSION METHODS USED\n")
cat("-----------------------\n")

print(
  univariable_results %>%
    count(
      method
    )
)


# ------------------------------------------------------------------------------
# 10. Variables analyzed using Firth regression
# ------------------------------------------------------------------------------

firth_results <- univariable_results %>%
  filter(
    method ==
      "Firth logistic regression"
  )


cat("\nVARIABLES ANALYZED USING FIRTH REGRESSION\n")
cat("-----------------------------------------\n")

print(
  firth_results %>%
    select(
      variable,
      n_analyzed,
      n_positive,
      odds_ratio,
      conf_low,
      conf_high,
      p_value
    ),
  n = Inf
)


# ------------------------------------------------------------------------------
# 11. Variables reported as non-estimable
# ------------------------------------------------------------------------------

nonestimable_results <- univariable_results %>%
  filter(
    method ==
      "Not estimable in historical analysis"
  )


cat("\nVARIABLES REPORTED AS NON-ESTIMABLE\n")
cat("-----------------------------------\n")

print(
  nonestimable_results %>%
    select(
      variable,
      n_analyzed,
      n_positive
    ),
  n = Inf
)


# ------------------------------------------------------------------------------
# 12. Sample-size distribution
# ------------------------------------------------------------------------------

cat("\nSAMPLE SIZE DISTRIBUTION\n")
cat("------------------------\n")

print(
  univariable_results %>%
    count(
      n_analyzed
    ) %>%
    arrange(
      desc(n_analyzed)
    ),
  n = Inf
)


# ------------------------------------------------------------------------------
# 13. Statistically notable results
#
# Descriptive output only. Statistical significance is not used as an
# automatic variable-selection criterion.
# ------------------------------------------------------------------------------

cat("\nUNIVARIABLE RESULTS WITH P < 0.05\n")
cat("----------------------------------\n")

print(
  univariable_results %>%
    filter(
      !is.na(p_value),
      p_value < 0.05
    ) %>%
    arrange(
      p_value
    ) %>%
    select(
      variable,
      n_analyzed,
      n_positive,
      method,
      odds_ratio,
      conf_low,
      conf_high,
      p_value
    ),
  n = Inf
)


# ==============================================================================
# PART F. SPECIFIC REPRODUCIBILITY CHECKS
# ==============================================================================


# ------------------------------------------------------------------------------
# 14. Age at symptom onset
#
# Expected:
#   N = 124
#   OR approximately 0.92 per one-year increase
# ------------------------------------------------------------------------------

age_check <- univariable_results %>%
  filter(
    variable == "age_onset"
  )


cat("\nAGE AT SYMPTOM ONSET\n")
cat("--------------------\n")

print(
  age_check,
  n = Inf
)


# ------------------------------------------------------------------------------
# 15. Age at genetic testing
#
# This variable was included in the historical univariable analysis code and
# is therefore retained among the 79 eligible predictors. It was omitted from
# the manuscript version of Table S5.
# ------------------------------------------------------------------------------

age_genetic_test_check <- univariable_results %>%
  filter(
    variable == "age_genetic_test"
  )


cat("\nAGE AT GENETIC TESTING\n")
cat("----------------------\n")

print(
  age_genetic_test_check,
  n = Inf
)


# ------------------------------------------------------------------------------
# 16. Sex
#
# The historical Table S5 reports N positive = 67 for male sex, but the
# verified cohort contains 57 male participants. The historical analysis code
# separately modeled sex with female as the reference category.
#
# Expected:
#   N = 124
#   N positive = 57
#   OR approximately 2.14
# ------------------------------------------------------------------------------

sex_check <- univariable_results %>%
  filter(
    variable == "sex"
  )


cat("\nSEX\n")
cat("---\n")

print(
  sex_check,
  n = Inf
)


# ------------------------------------------------------------------------------
# 17. Eczema
# ------------------------------------------------------------------------------

eczema_check <- univariable_results %>%
  filter(
    variable == "eczema"
  )


cat("\nECZEMA\n")
cat("------\n")

print(
  eczema_check,
  n = Inf
)


# ------------------------------------------------------------------------------
# 18. IgG1
# ------------------------------------------------------------------------------

igg1_check <- univariable_results %>%
  filter(
    variable == "igg1"
  )


cat("\nIgG1\n")
cat("----\n")

print(
  igg1_check,
  n = Inf
)


# ------------------------------------------------------------------------------
# 19. Historical special cases
# ------------------------------------------------------------------------------

cat("\nHISTORICAL SPECIAL CASES\n")
cat("------------------------\n")

print(
  univariable_results %>%
    filter(
      variable %in% c(
        firth_variables,
        nonestimable_variables
      )
    ) %>%
    select(
      variable,
      n_analyzed,
      n_positive,
      method,
      odds_ratio,
      conf_low,
      conf_high,
      p_value
    ),
  n = Inf
)


# ------------------------------------------------------------------------------
# 20. Full 79-predictor audit output
# ------------------------------------------------------------------------------

cat("\nALL 79 ELIGIBLE PREDICTORS\n")
cat("--------------------------\n")

print(
  univariable_results %>%
    select(
      variable,
      n_analyzed,
      n_positive,
      method,
      odds_ratio,
      conf_low,
      conf_high,
      p_value
    ) %>%
    arrange(
      variable
    ),
  n = Inf
)


# ==============================================================================
# PART G. SAVE RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 21. Create output directory
# ------------------------------------------------------------------------------

dir.create(
  "results/univariable",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 22. Save statistical results
# ------------------------------------------------------------------------------

write_csv(
  univariable_results %>%
    select(
      variable,
      predictor_type,
      n_analyzed,
      n_positive,
      method,
      odds_ratio,
      conf_low,
      conf_high,
      p_value
    ),
  "results/univariable/univariable_results.csv"
)


# ------------------------------------------------------------------------------
# 23. Save additional model diagnostics
# ------------------------------------------------------------------------------

write_csv(
  univariable_results %>%
    select(
      variable,
      method,
      glm_coefficient,
      glm_standard_error
    ),
  "results/univariable/univariable_model_diagnostics.csv"
)


# ==============================================================================
# PART H. FINAL CHECKS
# ==============================================================================


stopifnot(
  nrow(univariable_results) == 79,
  nrow(firth_results) == 2,
  nrow(nonestimable_results) == 1
)


message("Univariable analyses completed successfully.")

message(
  "Predictors analyzed: ",
  nrow(univariable_results)
)

message(
  "Standard logistic regression models: ",
  sum(
    univariable_results$method ==
      "Standard logistic regression"
  )
)

message(
  "Firth logistic regression models: ",
  nrow(firth_results)
)

message(
  "Historically non-estimable variables: ",
  nrow(nonestimable_results)
)