# ==============================================================================
# 03_variable_selection.R
#
# Study:
#   Clinical predictors of genetic diagnoses in adults evaluated for
#   inborn errors of immunity
#
# Purpose:
#   Reproduce the clinical-variable prioritization workflow summarized in
#   Figure S1.
#
#   This script:
#     1. defines the historical candidate predictor set
#     2. excludes complement-related variables
#     3. excludes variables with >20% missingness
#     4. excludes variables with no variability
#     5. assigns eligible predictors to clinical/laboratory domains
#     6. performs domain-specific LASSO logistic regression
#     7. reports the complete LASSO results for audit
#     8. records the historically documented 22-variable candidate set
#
# IMPORTANT:
#   This was an exploratory variable-prioritization procedure. It was not
#   intended to develop or validate a clinical prediction model.
#
#   LASSO was used as one component of variable prioritization. The final
#   candidate-variable set was determined by reviewing the LASSO results
#   together with clinical relevance and statistical considerations.
#
# Input:
#   data/derived/analysis_data.rds
#
# Outputs:
#   data/derived/predictor_missingness.rds
#   data/derived/historical_predictor_screening.rds
#   data/derived/eligible_predictors.rds
#   data/derived/predictor_domains.rds
#   data/derived/lasso_summary.rds
#   data/derived/candidate_predictors_22.rds
#
#   results/variable_selection/lasso_domain_summary.csv
#   data/derived/historical_predictor_screening.csv
#   results/variable_selection/lasso_coefficients.csv
#   results/variable_selection/lasso_nonzero_lambda_min.csv
#   results/variable_selection/historical_candidate_comparison.csv
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tibble)
library(readr)
library(glmnet)


# ------------------------------------------------------------------------------
# 2. Load analysis dataset
# ------------------------------------------------------------------------------

dat <- readRDS(
  "data/derived/analysis_data.rds"
)

stopifnot(nrow(dat) == 124)


# ==============================================================================
# PART A. DEFINE REPOSITORY-AVAILABLE PREDICTOR SET
# ==============================================================================


# ------------------------------------------------------------------------------
# 3. Outcome
# ------------------------------------------------------------------------------

y <- dat$molecular_diagnosis

stopifnot(
  length(y) == 124,
  sum(y == 1) == 55,
  sum(y == 0) == 69,
  !anyNA(y)
)


# ------------------------------------------------------------------------------
# 4. Historical preprocessing audit
# ------------------------------------------------------------------------------

# The original exploratory variable-prioritization workflow began with
# 109 candidate predictors.
#
# Before domain-specific LASSO:
#
#   6 complement-related variables were excluded during preprocessing
#   22 variables were excluded because >20% of observations were missing
#   2 variables were excluded because they had no variability
#
# This left 79 eligible predictors.
#
# The deposited repository dataset intentionally contains only these 79
# eligible predictors, together with the molecular-diagnosis outcome and
# IUIS category required for other manuscript analyses.
#
# Consequently, the 109-to-79 screening step is documented here as historical
# preprocessing provenance; it cannot and should not be recomputed from the
# deposited dataset because the 30 excluded variables are not released.

historical_screening <- tibble(
  stage = c(
    "Initial candidate predictors",
    "Complement-related exclusions",
    "Excluded for >20% missingness",
    "Excluded for no variability",
    "Eligible predictors"
  ),
  n = c(
    109L,
    6L,
    22L,
    2L,
    79L
  )
)

print(historical_screening)

stopifnot(
  109L - 6L - 22L - 2L == 79L
)


# Variables known from the historical reconstruction to have had no
# variability.

historical_no_variability <- c(
  "HIV",
  "Autoimmune liver disease"
)


# ------------------------------------------------------------------------------
# 5. Define the 79 predictors available for reproducible analysis
# ------------------------------------------------------------------------------

# Source-variable names are used where possible. Variables renamed by
# 01_prepare_analysis_data.R use their analysis names.

eligible_vars <- c(
  "age_onset",
  "age_genetic_test",
  "sex",
  "consanguinity",
  "familial",
  "CVID",
  "Autoinflammation/periodic fever",
  "Splenomegaly",
  "Lymphproliferation",
  "GLILD",
  "NRH/portosinusoidal vascular disease",
  "Immunoglobulins",
  "bronchiectasis (X-ray)",
  "Nr. of Pneumonia (X-ray)",
  "1Resp inf",
  "2Resp inf",
  "3Resp virus",
  "4Resp PCJ/TBC",
  "Unusal microbes/localisation",
  "Gastrointestinal pathogens",
  "Blood stream infections",
  "CNS infection",
  "Skin/bone/joints",
  "Warts (HPV)",
  "HSV1",
  "Herpes zoster (VZV)",
  "Fungal-skin inf",
  "HSV genital",
  "Condyloma (HPV)",
  "Genital dysplasia (LSIL, HSIL)",
  "Malignancy",
  "Lymphoma",
  "Ventricular malignancy",
  "Other cancer",
  "Atopy",
  "asthma",
  "urticaria",
  "eczema",
  "AB allergy",
  "food/animal/pollen allergy",
  "anaphylaxis",
  "Autoimmunity",
  "Systemic disease: RA/SLE/Sjögren",
  "Other autoimmune disease",
  "CNS",
  "Skin",
  "IBD/microscopic colitis",
  "Thyroid dysregulation (thyroiditis, hyperthyroidism)",
  "Celiaki",
  "Pernicious anemia",
  "DM1",
  "Cytopenia",
  "CRP low",
  "CRPhigh",
  "SR low",
  "SR high",
  "Hb at admission",
  "WBC total",
  "neutrophils",
  "Lymfocyt",
  "Eosinophilia",
  "Trc",
  "AST",
  "ALT",
  "ALP",
  "GT",
  "Krea",
  "IgG before treatment/ or at admission",
  "IgA",
  "IgM",
  "igg1",
  "IgG2",
  "IgG3",
  "IgG4",
  "B CELL - before immunosuppression/trxpl",
  "CD16/56",
  "CD3",
  "CD4",
  "CD8"
)


stopifnot(
  length(eligible_vars) == 79,
  length(unique(eligible_vars)) == 79,
  all(eligible_vars %in% names(dat))
)


cat("\nREPOSITORY-AVAILABLE ELIGIBLE PREDICTORS\n")
cat("----------------------------------------\n")
cat("Number:", length(eligible_vars), "\n\n")


# ------------------------------------------------------------------------------
# 6. Missingness among the 79 released predictors
# ------------------------------------------------------------------------------

# This describes missingness in the deposited analysis dataset.
# It is not used to repeat the historical >20% screening.

missingness <- tibble(
  variable = eligible_vars,
  n_missing = colSums(
    is.na(dat[eligible_vars])
  ),
  pct_missing = 100 * colMeans(
    is.na(dat[eligible_vars])
  )
)

stopifnot(
  all(missingness$pct_missing <= 20)
)

print(
  missingness %>%
    arrange(desc(pct_missing)),
  n = Inf
)
# ==============================================================================
# PART C. DEFINE VARIABLE DOMAINS
# ==============================================================================


# ------------------------------------------------------------------------------
# 8. Domain definitions
# ------------------------------------------------------------------------------

# Candidate predictors were grouped into clinically related domains before
# LASSO prioritization. The definitions below reproduce the domains documented
# in the historical exploratory analysis.


# General clinical characteristics

general_vars <- c(
  "age_onset",
  "sex",
  "Autoinflammation/periodic fever",
  "familial",
  "CVID"
)


# CVID-related clinical manifestations

cvid_vars <- c(
  "CVID",
  "Splenomegaly",
  "Lymphproliferation",
  "GLILD",
  "NRH/portosinusoidal vascular disease",
  "Immunoglobulins",
  "bronchiectasis (X-ray)"
)


# Infectious manifestations

infection_vars <- c(
  "1Resp inf",
  "2Resp inf",
  "3Resp virus",
  "4Resp PCJ/TBC",
  "Unusal microbes/localisation",
  "Gastrointestinal pathogens",
  "Blood stream infections",
  "CNS infection",
  "Skin/bone/joints",
  "Warts (HPV)",
  "HSV1",
  "Herpes zoster (VZV)",
  "Fungal-skin inf",
  "HSV genital",
  "Condyloma (HPV)"
)


# Malignancy

malignancy_vars <- c(
  "Genital dysplasia (LSIL, HSIL)",
  "Malignancy",
  "Lymphoma",
  "Ventricular malignancy",
  "Other cancer"
)


# Atopy/allergy

allergy_vars <- c(
  "Atopy",
  "asthma",
  "urticaria",
  "eczema",
  "AB allergy",
  "food/animal/pollen allergy",
  "anaphylaxis"
)


# Autoimmunity

autoimmune_vars <- c(
  "Autoimmunity",
  "Systemic disease: RA/SLE/Sjögren",
  "Other autoimmune disease",
  "CNS",
  "Skin",
  "IBD/microscopic colitis",
  "Thyroid dysregulation (thyroiditis, hyperthyroidism)",
  "Celiaki",
  "Pernicious anemia",
  "DM1",
  "Cytopenia"
)


# Inflammatory laboratory markers

inflammation_lab_vars <- c(
  "CRPhigh",
  "SR high",
  "S-AA high",
  "Beta2m high",
  "IL-2 rec high"
)


# General laboratory variables

general_lab_vars <- c(
  "Hb at admission",
  "WBC total",
  "neutrophils",
  "Lymfocyt",
  "Eosinophilia",
  "Trc",
  "AST",
  "ALT",
  "ALP",
  "GT",
  "Krea"
)


# Immunoglobulins
#
# IgE was not included in the later historical domain model because its
# inclusion substantially reduced the complete-case sample.

immunoglobulin_vars <- c(
  "IgG before treatment/ or at admission",
  "IgA",
  "IgM",
  "igg1",
  "IgG2",
  "IgG3",
  "IgG4"
)


# Lymphocyte subsets

lymphocyte_vars <- c(
  "B CELL - before immunosuppression/trxpl",
  "Naive B-cells",
  "Switched memory B cells",
  "CD21 low",
  "CD16/56",
  "CD3",
  "CD4",
  "CD8",
  "Naive CD4 cells",
  "Naiva CD8",
  "T reg/CD4+CD5+CD127+",
  "Th1EM",
  "Th2EM",
  "Th17EM",
  "Th1CM",
  "Th2CM",
  "Th17CM"
)


# ------------------------------------------------------------------------------
# 9. Restrict domains to eligible variables
# ------------------------------------------------------------------------------

keep_eligible <- function(vars) {
  
  intersect(
    vars,
    eligible_vars
  )
}


domains <- list(
  General = keep_eligible(general_vars),
  CVID_related = keep_eligible(cvid_vars),
  Infection = keep_eligible(infection_vars),
  Malignancy = keep_eligible(malignancy_vars),
  Allergy = keep_eligible(allergy_vars),
  Autoimmunity = keep_eligible(autoimmune_vars),
  Inflammation_labs = keep_eligible(inflammation_lab_vars),
  General_labs = keep_eligible(general_lab_vars),
  Immunoglobulins = keep_eligible(immunoglobulin_vars),
  Lymphocytes = keep_eligible(lymphocyte_vars)
)


cat("\nDOMAIN SIZES\n")
cat("------------\n")

print(
  sapply(
    domains,
    length
  )
)


# Expected reconstruction:
#
# General              5
# CVID_related         7
# Infection           15
# Malignancy           5
# Allergy              7
# Autoimmunity        11
# Inflammation_labs    2
# General_labs        11
# Immunoglobulins      7
# Lymphocytes          5


expected_domain_sizes <- c(
  General = 5L,
  CVID_related = 7L,
  Infection = 15L,
  Malignancy = 5L,
  Allergy = 7L,
  Autoimmunity = 11L,
  Inflammation_labs = 2L,
  General_labs = 11L,
  Immunoglobulins = 7L,
  Lymphocytes = 5L
)


stopifnot(
  identical(
    sapply(domains, length),
    expected_domain_sizes
  )
)


# ==============================================================================
# PART D. CHECK DOMAIN COVERAGE
# ==============================================================================


# ------------------------------------------------------------------------------
# 10. Identify eligible variables not included in domain-specific LASSO
# ------------------------------------------------------------------------------

assigned_vars <- unique(
  unlist(domains)
)


unassigned_vars <- setdiff(
  eligible_vars,
  assigned_vars
)


# Five variables survived the initial eligibility filters but were not included
# in the domain-specific LASSO models in the historical analysis.

expected_unassigned <- c(
  "age_genetic_test",
  "consanguinity",
  "Nr. of Pneumonia (X-ray)",
  "CRP low",
  "SR low"
)


stopifnot(
  setequal(
    unassigned_vars,
    expected_unassigned
  )
)


cat(
  "\nELIGIBLE VARIABLES NOT INCLUDED IN DOMAIN-SPECIFIC LASSO\n"
)
cat(
  "--------------------------------------------------------\n"
)
cat(
  "Number:",
  length(unassigned_vars),
  "\n\n"
)

print(unassigned_vars)


# ==============================================================================
# PART E. PREPARE VARIABLES FOR LASSO
# ==============================================================================


# ------------------------------------------------------------------------------
# 11. Convert predictors to numeric representation
# ------------------------------------------------------------------------------

# glmnet requires a numeric design matrix.
#
# Binary Y/N variables are converted to 1/0.
# Numeric variables remain numeric.
# Other categorical variables are handled through model.matrix() below.

yn_to_binary <- function(x) {
  
  if (!is.character(x)) {
    return(x)
  }
  
  observed <- unique(
    na.omit(x)
  )
  
  if (
    length(observed) > 0 &&
    all(observed %in% c("Y", "N"))
  ) {
    
    return(
      ifelse(
        is.na(x),
        NA_real_,
        ifelse(
          x == "Y",
          1,
          0
        )
      )
    )
  }
  
  x
}


lasso_source <- dat %>%
  mutate(
    across(
      all_of(eligible_vars),
      yn_to_binary
    )
  )


# ==============================================================================
# PART F. DOMAIN-SPECIFIC LASSO
# ==============================================================================


# ------------------------------------------------------------------------------
# 13. LASSO function
# ------------------------------------------------------------------------------

# LASSO was used as an exploratory variable-prioritization tool.
#
# Separate cross-validated logistic LASSO models are fitted within each
# clinical/laboratory domain.
#
# The results are retained for all predictors entered into each model.
# Variables with non-zero coefficients at lambda.min are reported for
# transparency, but this was not used as a fully automated rule for defining
# the final candidate-variable set.
#
# The historical candidate-variable set was determined by reviewing the LASSO
# results together with clinical relevance, redundancy/collinearity, rarity,
# and suitability for subsequent regression modelling.


run_domain_lasso <- function(
    data,
    outcome,
    variables,
    seed = 123
) {
  
  if (length(variables) == 0) {
    return(NULL)
  }
  
  
  # --------------------------------------------------------------------------
  # Create domain-specific dataset
  # --------------------------------------------------------------------------
  
  domain_data <- data %>%
    select(all_of(variables)) %>%
    mutate(
      .outcome = outcome
    )
  
  
  # --------------------------------------------------------------------------
  # Complete-case analysis within this domain
  # --------------------------------------------------------------------------
  
  complete <- complete.cases(domain_data)
  
  domain_complete <- domain_data[
    complete,
    ,
    drop = FALSE
  ]
  
  
  y_complete <- domain_complete$.outcome
  
  predictor_data <- domain_complete %>%
    select(-.outcome)
  
  
  # glmnet requires a numeric design matrix.
  #
  # model.matrix() also handles categorical predictors such as sex.
  # The intercept is removed because glmnet supplies its own intercept.
  
  x_complete <- model.matrix(
    ~ .,
    data = predictor_data
  )[, -1, drop = FALSE]
  
  
  # --------------------------------------------------------------------------
  # Safety checks
  # --------------------------------------------------------------------------
  
  stopifnot(
    nrow(x_complete) == length(y_complete)
  )
  
  if (length(unique(y_complete)) < 2) {
    stop(
      "Outcome has no variability in complete cases for this domain."
    )
  }
  
  if (ncol(x_complete) < 2) {
    stop(
      "Domain contains fewer than two predictors after model-matrix creation."
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Cross-validated LASSO logistic regression
  # --------------------------------------------------------------------------
  
  set.seed(seed)
  
  cv_fit <- cv.glmnet(
    x = x_complete,
    y = y_complete,
    family = "binomial",
    alpha = 1,
    standardize = TRUE
  )
  
  
  # --------------------------------------------------------------------------
  # Extract coefficients at lambda.min and lambda.1se
  # --------------------------------------------------------------------------
  
  coef_min <- as.matrix(
    coef(
      cv_fit,
      s = "lambda.min"
    )
  )
  
  coef_1se <- as.matrix(
    coef(
      cv_fit,
      s = "lambda.1se"
    )
  )
  
  
  # Remove intercept from coefficient tables.
  
  coef_min <- coef_min[
    rownames(coef_min) != "(Intercept)",
    ,
    drop = FALSE
  ]
  
  coef_1se <- coef_1se[
    rownames(coef_1se) != "(Intercept)",
    ,
    drop = FALSE
  ]
  
  
  stopifnot(
    identical(
      rownames(coef_min),
      rownames(coef_1se)
    )
  )
  
  
  # --------------------------------------------------------------------------
  # Predictor-level results
  # --------------------------------------------------------------------------
  
  coefficient_table <- tibble(
    variable = rownames(coef_min),
    coefficient_lambda_min = as.numeric(coef_min[, 1]),
    coefficient_lambda_1se = as.numeric(coef_1se[, 1])
  ) %>%
    mutate(
      selected_lambda_min = coefficient_lambda_min != 0,
      selected_lambda_1se = coefficient_lambda_1se != 0
    )
  
  
  # --------------------------------------------------------------------------
  # Return results
  # --------------------------------------------------------------------------
  
  list(
    n_complete = nrow(domain_complete),
    n_predictors = ncol(x_complete),
    lambda_min = cv_fit$lambda.min,
    lambda_1se = cv_fit$lambda.1se,
    coefficient_table = coefficient_table,
    fit = cv_fit
  )
}


# ------------------------------------------------------------------------------
# 14. Run LASSO separately within each domain
# ------------------------------------------------------------------------------

lasso_results <- lapply(
  domains,
  function(vars) {
    
    run_domain_lasso(
      data = lasso_source,
      outcome = y,
      variables = vars,
      seed = 123
    )
  }
)


# ==============================================================================
# PART G. SUMMARIZE LASSO RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 15. Domain-level summary
# ------------------------------------------------------------------------------

lasso_summary <- bind_rows(
  lapply(
    names(lasso_results),
    function(domain_name) {
      
      result <- lasso_results[[domain_name]]
      
      if (is.null(result)) {
        
        return(
          tibble(
            domain = domain_name,
            n_complete = NA_integer_,
            n_predictors = 0L,
            lambda_min = NA_real_,
            lambda_1se = NA_real_,
            n_selected_lambda_min = 0L,
            n_selected_lambda_1se = 0L
          )
        )
      }
      
      
      tibble(
        domain = domain_name,
        n_complete = result$n_complete,
        n_predictors = result$n_predictors,
        lambda_min = result$lambda_min,
        lambda_1se = result$lambda_1se,
        n_selected_lambda_min =
          sum(result$coefficient_table$selected_lambda_min),
        n_selected_lambda_1se =
          sum(result$coefficient_table$selected_lambda_1se)
      )
    }
  )
)


cat("\nLASSO DOMAIN SUMMARY\n")
cat("--------------------\n")

print(
  lasso_summary,
  n = Inf
)


# ------------------------------------------------------------------------------
# 16. Predictor-level LASSO results
# ------------------------------------------------------------------------------

# Combine the complete coefficient tables from all domain-specific models.
#
# All predictors are retained in this table, including predictors whose
# coefficients were shrunk to zero.

lasso_coefficients <- bind_rows(
  lapply(
    names(lasso_results),
    function(domain_name) {
      
      result <- lasso_results[[domain_name]]
      
      if (is.null(result)) {
        return(NULL)
      }
      
      result$coefficient_table %>%
        mutate(
          domain = domain_name,
          n_complete = result$n_complete,
          lambda_min = result$lambda_min,
          lambda_1se = result$lambda_1se,
          .before = 1
        )
    }
  )
)


cat("\nPREDICTOR-LEVEL LASSO RESULTS\n")
cat("-----------------------------\n")

print(
  lasso_coefficients,
  n = Inf
)


# ------------------------------------------------------------------------------
# 17. Non-zero coefficients at lambda.min
# ------------------------------------------------------------------------------

# These results are provided as an audit of the reconstructed LASSO analyses.
#
# IMPORTANT:
# Non-zero coefficients at lambda.min were not treated as a deterministic
# selection rule for the subsequent analysis. LASSO was one component of the
# broader variable-prioritization procedure.

lasso_nonzero_lambda_min <- lasso_coefficients %>%
  filter(selected_lambda_min) %>%
  arrange(
    domain,
    desc(abs(coefficient_lambda_min))
  )


cat("\nNON-ZERO COEFFICIENTS AT LAMBDA.MIN\n")
cat("-----------------------------------\n")

print(
  lasso_nonzero_lambda_min,
  n = Inf
)


# ==============================================================================
# PART H. HISTORICALLY DOCUMENTED CANDIDATE VARIABLE SET
# ==============================================================================


# ------------------------------------------------------------------------------
# 18. Candidate variables retained after LASSO-informed review
# ------------------------------------------------------------------------------

# Domain-specific LASSO was used to prioritize potentially informative
# variables rather than as an automated variable-selection algorithm.
#
# The LASSO results were reviewed together with:
#
#   - clinical relevance and biological plausibility
#   - redundancy and correlation between predictors
#   - rarity of individual manifestations
#   - statistical issues such as separation
#   - suitability for subsequent multivariable modelling
#
# Contemporaneous documentation from the original analysis records that this
# review resulted in the following 22 candidate variables.
#
# These 22 variables are therefore explicitly recorded here rather than being
# reconstructed by imposing an arbitrary coefficient threshold on the current
# LASSO implementation.


candidate_vars_22 <- c(
  "age_onset",
  "CVID",
  "bronchiectasis (X-ray)",
  "NRH/portosinusoidal vascular disease",
  "igg1",
  "CD16/56",
  "B CELL - before immunosuppression/trxpl",
  "sex",
  "Immunoglobulins",
  "Splenomegaly",
  "GLILD",
  "eczema",
  "AB allergy",
  "food/animal/pollen allergy",
  "Celiaki",
  "ALT",
  "ALP",
  "Hb at admission",
  "IgG4",
  "IgG3",
  "CD8",
  "CD3"
)


# Verify that the documented candidate set contains exactly 22 unique
# variables and that all are present in the analysis dataset.

stopifnot(
  length(candidate_vars_22) == 22,
  length(unique(candidate_vars_22)) == 22,
  all(candidate_vars_22 %in% names(dat))
)


cat("\nHISTORICALLY DOCUMENTED CANDIDATE SET\n")
cat("-------------------------------------\n")
cat("Number of variables:", length(candidate_vars_22), "\n\n")

print(candidate_vars_22)


# ------------------------------------------------------------------------------
# 19. Show how the historical candidate set relates to current LASSO results
# ------------------------------------------------------------------------------

# This comparison is descriptive only.
#
# It is intentionally NOT used to redefine the historical candidate set.
# Its purpose is to make the reconstruction transparent and auditable.

historical_candidate_comparison <- lasso_coefficients %>%
  mutate(
    historical_candidate =
      variable %in% candidate_vars_22
  ) %>%
  arrange(
    domain,
    desc(historical_candidate),
    desc(abs(coefficient_lambda_min))
  )


cat("\nHISTORICAL CANDIDATES AND RECONSTRUCTED LASSO RESULTS\n")
cat("-----------------------------------------------------\n")

print(
  historical_candidate_comparison,
  n = Inf
)



# ==============================================================================
# CORRELATION ASSESSMENT AMONG THE 22 CANDIDATE VARIABLES
# ==============================================================================

# Correlation matrices were reviewed as part of the post-LASSO assessment of
# collinearity and redundancy among the 22 historically documented candidate
# variables.
#
# The reduction from 22 variables to the smaller candidate multivariable sets
# was NOT based on an algorithmic correlation threshold. Correlation structure
# was considered together with redundancy, clinical relevance, rarity of
# manifestations, and suitability for multivariable regression.


# ------------------------------------------------------------------------------
# Prepare candidate variables for correlation analysis
# ------------------------------------------------------------------------------

correlation_data <- dat %>%
  select(
    all_of(candidate_vars_22)
  )


# Convert variables to numeric form for correlation analysis.
#
# Numeric variables are retained as recorded.
# Binary Y/N variables are converted to 1/0.
# Sex is represented as male = 1, female = 0.

correlation_numeric <- correlation_data


for (v in names(correlation_numeric)) {
  
  if (is.numeric(correlation_numeric[[v]])) {
    next
  }
  
  values <- correlation_numeric[[v]]
  
  unique_values <- unique(
    na.omit(
      as.character(values)
    )
  )
  
  
  # Y/N variables
  
  if (
    length(unique_values) > 0 &&
    all(unique_values %in% c("Y", "N"))
  ) {
    
    correlation_numeric[[v]] <-
      ifelse(
        values == "Y",
        1,
        ifelse(
          values == "N",
          0,
          NA_real_
        )
      )
    
    next
  }
  
  
  # Sex
  
  if (
    all(unique_values %in% c("M", "F"))
  ) {
    
    correlation_numeric[[v]] <-
      ifelse(
        values == "M",
        1,
        ifelse(
          values == "F",
          0,
          NA_real_
        )
      )
    
    next
  }
  
  
  # Stop rather than silently assigning arbitrary numeric codes.
  
  stop(
    paste0(
      "Variable '",
      v,
      "' could not be converted safely for correlation analysis. ",
      "Observed values: ",
      paste(
        unique_values,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------------------------
# Calculate correlation matrix
# ------------------------------------------------------------------------------

# Spearman correlations are used because the candidate set contains a mixture
# of continuous, ordinal-like, and binary variables and several laboratory
# variables have non-normal distributions.
#
# Pairwise complete observations are used so that the correlation matrix does
# not unnecessarily exclude participants because of missingness in another
# candidate variable.

candidate_cor_matrix <- cor(
  correlation_numeric,
  use = "pairwise.complete.obs",
  method = "spearman"
)


cat("\nCORRELATION MATRIX: 22 CANDIDATE VARIABLES\n")
cat("-------------------------------------------\n")

print(
  round(
    candidate_cor_matrix,
    2
  )
)


# ------------------------------------------------------------------------------
# Create long-format correlation table
# ------------------------------------------------------------------------------

candidate_cor_long <- as.data.frame(
  as.table(candidate_cor_matrix)
) %>%
  as_tibble() %>%
  rename(
    variable_1 = Var1,
    variable_2 = Var2,
    correlation = Freq
  ) %>%
  mutate(
    variable_1 = as.character(variable_1),
    variable_2 = as.character(variable_2),
    absolute_correlation = abs(correlation)
  ) %>%
  filter(
    variable_1 != variable_2
  )


# Retain one copy of each variable pair.

candidate_order <- names(correlation_numeric)

candidate_cor_long <- candidate_cor_long %>%
  mutate(
    index_1 = match(
      variable_1,
      candidate_order
    ),
    index_2 = match(
      variable_2,
      candidate_order
    )
  ) %>%
  filter(
    index_1 < index_2
  ) %>%
  arrange(
    desc(absolute_correlation)
  ) %>%
  select(
    variable_1,
    variable_2,
    correlation,
    absolute_correlation
  )


cat("\nCANDIDATE VARIABLE CORRELATIONS, ORDERED BY MAGNITUDE\n")
cat("-----------------------------------------------------\n")

print(
  candidate_cor_long,
  n = Inf
)


# ==============================================================================
# PART I. SAVE VARIABLE-PRIORITIZATION OBJECTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 20. Create output directories
# ------------------------------------------------------------------------------

dir.create(
  "data/derived",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/variable_selection",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 21. Save internal R objects used by later repository scripts
# ------------------------------------------------------------------------------

saveRDS(
  missingness,
  "data/derived/predictor_missingness.rds"
)

saveRDS(
  eligible_vars,
  "data/derived/eligible_predictors.rds"
)

saveRDS(
  historical_screening,
  "data/derived/historical_predictor_screening.rds"
)

write_csv(
  historical_screening,
  "results/variable_selection/historical_predictor_screening.csv"
)

saveRDS(
  domains,
  "data/derived/predictor_domains.rds"
)

saveRDS(
  lasso_summary,
  "data/derived/lasso_summary.rds"
)

saveRDS(
  candidate_vars_22,
  "data/derived/candidate_predictors_22.rds"
)



write_csv(
  as.data.frame(candidate_cor_matrix) %>%
    rownames_to_column(
      "variable"
    ),
  "results/variable_selection/candidate_22_correlation_matrix.csv"
)


write_csv(
  candidate_cor_long,
  "results/variable_selection/candidate_22_correlations_long.csv"
)

# ------------------------------------------------------------------------------
# 22. Save human-readable audit tables
# ------------------------------------------------------------------------------

write_csv(
  lasso_summary,
  "results/variable_selection/lasso_domain_summary.csv"
)

write_csv(
  lasso_coefficients,
  "results/variable_selection/lasso_coefficients.csv"
)

write_csv(
  lasso_nonzero_lambda_min,
  "results/variable_selection/lasso_nonzero_lambda_min.csv"
)

write_csv(
  historical_candidate_comparison,
  "results/variable_selection/historical_candidate_comparison.csv"
)


message("Variable-prioritization analysis completed successfully.")
message("Eligible predictors after filtering: ", length(eligible_vars))
message("Historically documented candidate variables: ", length(candidate_vars_22))