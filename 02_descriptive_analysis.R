# ==============================================================================
# 02_descriptive_analysis.R
#
# Study:
#   Clinical predictors of genetic diagnoses in adults evaluated for
#   inborn errors of immunity
#
# Purpose:
#   Reproduce descriptive cohort results reported in the manuscript:
#
#   - cohort counts underlying Figure 1
#   - clinical manifestation domains
#   - Table 1
#   - clustering and heatmap underlying Figure 2
#
# Input:
#   data/derived/analysis_data.rds
#
# Output:
#   derived objects used by later figure/table scripts
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(tibble)


# ------------------------------------------------------------------------------
# 2. Load analysis dataset
# ------------------------------------------------------------------------------

dat <- readRDS("data/derived/analysis_data.rds")

stopifnot(nrow(dat) == 124)

# ------------------------------------------------------------------------------
# 3. Helper function - yes/no to binary
# ------------------------------------------------------------------------------


yn_to_binary <- function(x) {
  case_when(
    x == "Y" ~ 1L,
    x == "N" ~ 0L,
    is.na(x) ~ NA_integer_,
    TRUE ~ NA_integer_
  )
}

# ==============================================================================
# PART A. COHORT COUNTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 3. Cohort counts underlying Figure 1
# ------------------------------------------------------------------------------

# The first three counts originate from cohort ascertainment and are not
# derivable from the 124-participant repository dataset.

cohort_flow <- tibble(
  stage = c(
    "Patients attending adult IEI clinic",
    "Patients undergoing genetic testing",
    "Patients included in study",
    "Molecular diagnosis",
    "No molecular diagnosis"
  ),
  n = c(
    976L,
    153L,
    nrow(dat),
    sum(dat$molecular_diagnosis == 1),
    sum(dat$molecular_diagnosis == 0)
  )
)

print(cohort_flow)


# Verify counts represented in the manuscript.

stopifnot(
  cohort_flow$n[cohort_flow$stage == "Patients included in study"] == 124,
  cohort_flow$n[cohort_flow$stage == "Molecular diagnosis"] == 55,
  cohort_flow$n[cohort_flow$stage == "No molecular diagnosis"] == 69
)


# Diagnostic yield

diagnostic_yield <- mean(dat$molecular_diagnosis)

stopifnot(
  round(100 * diagnostic_yield, 1) == 44.4
)


# ==============================================================================
# PART B. CLINICAL MANIFESTATION DOMAINS
# ==============================================================================


# ------------------------------------------------------------------------------
# 4. Helper function for clinical manifestation variables
# ------------------------------------------------------------------------------

# Clinical manifestation fields in the source database were recorded as Y/N.
#
# IMPORTANT ASSUMPTION:
# For construction of the manifestation domains, a missing entry is treated
# as absence of that manifestation. This reproduces the coding used in the
# original analysis.
#
# This assumption applies ONLY to construction of the descriptive clinical
# domains and is not used as a general missing-data rule elsewhere.

manifestation_binary <- function(x) {
  
  case_when(
    x == "Y" ~ 1L,
    x == "N" ~ 0L,
    is.na(x) ~ 0L,
    TRUE ~ NA_integer_
  )
}


# ------------------------------------------------------------------------------
# 5. Define variables belonging to each clinical manifestation domain
# ------------------------------------------------------------------------------

# Definitions correspond to Table S2.

bacterial_fungal_vars <- c(
  "1Resp inf",
  "2Resp inf",
  "4Resp PCJ/TBC",
  "Unusal microbes/localisation",
  "Gastrointestinal pathogens",
  "Blood stream infections",
  "CNS infection",
  "Skin/bone/joints",
  "Fungal-skin inf",
  "bronchiectasis (X-ray)"
)


viral_vars <- c(
  "3Resp virus",
  "Warts (HPV)",
  "HSV1",
  "Herpes zoster (VZV)",
  "HSV genital",
  "Condyloma (HPV)",
  "Genital dysplasia (LSIL, HSIL)"
)


autoimmune_vars <- c(
  "Systemic disease: RA/SLE/Sjögren",
  "Other autoimmune disease",
  "CNS",
  "Skin",
  "DM1",
  "Thyroid dysregulation (thyroiditis, hyperthyroidism)",
  "IBD/microscopic colitis",
  "Celiaki",
  "Pernicious anemia"
)


lymphoproliferation_vars <- c(
  "Splenomegaly",
  "Lymphproliferation",
  "Cytopenia",
  "GLILD",
  "NRH/portosinusoidal vascular disease"
)


atopy_vars <- c(
  "asthma",
  "urticaria",
  "eczema",
  "AB allergy",
  "food/animal/pollen allergy",
  "anaphylaxis"
)


malignancy_vars <- c(
  "Lymphoma",
  "Ventricular malignancy",
  "Other cancer"
)


autoinflammation_vars <- c(
  "Autoinflammation/periodic fever"
)


# ------------------------------------------------------------------------------
# 6. Check that all required source variables exist
# ------------------------------------------------------------------------------

domain_source_vars <- c(
  bacterial_fungal_vars,
  viral_vars,
  autoimmune_vars,
  lymphoproliferation_vars,
  atopy_vars,
  malignancy_vars,
  autoinflammation_vars
)

missing_domain_vars <- setdiff(domain_source_vars, names(dat))

if (length(missing_domain_vars) > 0) {
  
  stop(
    paste(
      "Variables required for clinical domains are missing:",
      paste(missing_domain_vars, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------------------------
# 7. Function to calculate domain burden
# ------------------------------------------------------------------------------

domain_score <- function(data, variables) {
  
  data %>%
    select(all_of(variables)) %>%
    mutate(
      across(
        everything(),
        manifestation_binary
      )
    ) %>%
    rowSums()
}


# ------------------------------------------------------------------------------
# 8. Calculate domain burden scores
# ------------------------------------------------------------------------------

domains <- tibble(
  bacterial_fungal = domain_score(dat, bacterial_fungal_vars),
  viral = domain_score(dat, viral_vars),
  autoimmune = domain_score(dat, autoimmune_vars),
  lymphoproliferation = domain_score(dat, lymphoproliferation_vars),
  atopy = domain_score(dat, atopy_vars),
  malignancy = domain_score(dat, malignancy_vars),
  autoinflammation = domain_score(dat, autoinflammation_vars)
)


# Presence of a domain is defined as at least one manifestation.

domain_presence <- domains %>%
  mutate(
    across(
      everything(),
      ~ as.integer(.x > 0)
    )
  )


# ------------------------------------------------------------------------------
# 9. Verify domain prevalence against Table 1
# ------------------------------------------------------------------------------

observed_domain_counts <- colSums(domain_presence)

expected_domain_counts <- c(
  bacterial_fungal = 116L,
  viral = 75L,
  autoimmune = 40L,
  lymphoproliferation = 62L,
  atopy = 74L,
  malignancy = 20L,
  autoinflammation = 14L
)

print(observed_domain_counts)


stopifnot(
  identical(
    as.integer(observed_domain_counts[names(expected_domain_counts)]),
    as.integer(expected_domain_counts)
  )
)


# Add derived domain variables to descriptive dataset.

descriptive_data <- bind_cols(
  dat,
  domains %>%
    rename_with(~ paste0("domain_score_", .x)),
  domain_presence %>%
    rename_with(~ paste0("domain_present_", .x))
)


# ==============================================================================
# PART C. TABLE 1
# ==============================================================================


# ------------------------------------------------------------------------------
# 10. Helper functions for Table 1
# ------------------------------------------------------------------------------

median_iqr <- function(x) {
  
  q <- quantile(
    x,
    probs = c(0.25, 0.50, 0.75),
    na.rm = TRUE
  )
  
  sprintf(
    "%.2f [%.2f, %.2f]",
    q[2],
    q[1],
    q[3]
  )
}


n_percent <- function(x) {
  
  n <- sum(x == 1, na.rm = TRUE)
  
  sprintf(
    "%d (%.1f)",
    n,
    100 * n / length(x)
  )
}


# ------------------------------------------------------------------------------
# 11. Recode Table 1 binary variables
# ------------------------------------------------------------------------------

table1_data <- descriptive_data %>%
  mutate(
    familial_binary = yn_to_binary(familial),
    consanguinity_binary = yn_to_binary(consanguinity),
    
    # CVID remains a source clinical classification.
    cvid_binary = yn_to_binary(CVID)
  )


# ------------------------------------------------------------------------------
# 12. Construct Table 1
# ------------------------------------------------------------------------------

table1 <- tibble(
  variable = c(
    "Number of patients",
    "Age at onset, years, median [IQR]",
    "Age at genetic test, years, median [IQR]",
    "Male sex, n (%)",
    "Familial cases, n (%)",
    "Consanguinity, n (%)",
    "CVID, n (%)",
    
    "Bacterial/fungal, n (%)",
    "Viral, n (%)",
    "Atopy, n (%)",
    "Lymphoproliferation, n (%)",
    "Autoimmune disease, n (%)",
    "Malignancy, n (%)",
    "Autoinflammation, n (%)",
    
    "Genetic diagnosis achieved, n (%)"
  ),
  
  value = c(
    as.character(nrow(table1_data)),
    
    median_iqr(table1_data$age_onset),
    median_iqr(table1_data$age_genetic_test),
    
    n_percent(table1_data$male),
    n_percent(table1_data$familial_binary),
    n_percent(table1_data$consanguinity_binary),
    n_percent(table1_data$cvid_binary),
    
    n_percent(table1_data$domain_present_bacterial_fungal),
    n_percent(table1_data$domain_present_viral),
    n_percent(table1_data$domain_present_atopy),
    n_percent(table1_data$domain_present_lymphoproliferation),
    n_percent(table1_data$domain_present_autoimmune),
    n_percent(table1_data$domain_present_malignancy),
    n_percent(table1_data$domain_present_autoinflammation),
    
    n_percent(table1_data$molecular_diagnosis)
  )
)


print(table1, n = Inf)


# ------------------------------------------------------------------------------
# 13. Verify selected Table 1 values
# ------------------------------------------------------------------------------

stopifnot(
  median(table1_data$age_onset, na.rm = TRUE) == 16.5,
  median(table1_data$age_genetic_test, na.rm = TRUE) == 33,
  sum(table1_data$male == 1, na.rm = TRUE) == 57,
  sum(table1_data$familial_binary == 1, na.rm = TRUE) == 12,
  sum(table1_data$consanguinity_binary == 1, na.rm = TRUE) == 4,
  sum(table1_data$cvid_binary == 1, na.rm = TRUE) == 39
)


# ==============================================================================
# PART D. FIGURE 2: DOMAIN HEATMAP AND CLUSTERING
# ==============================================================================


# ------------------------------------------------------------------------------
# 14. Prepare heatmap matrix
# ------------------------------------------------------------------------------

# Figure 2 displays the burden (number of manifestations) within each domain,
# rather than simple domain presence.

heatmap_matrix <- domains %>%
  as.matrix()

rownames(heatmap_matrix) <- seq_len(nrow(heatmap_matrix))


# Patients are rows at this stage; transpose so that:
#
#   rows    = clinical domains
#   columns = participants

heatmap_matrix <- t(heatmap_matrix)


# More readable domain labels.

rownames(heatmap_matrix) <- c(
  "Bacterial/fungal",
  "Viral",
  "Autoimmune",
  "Lymphoproliferation",
  "Atopy",
  "Malignancy",
  "Autoinflammation"
)


# ------------------------------------------------------------------------------
# 15. Hierarchical clustering of participants
# ------------------------------------------------------------------------------

# Clustering reproduces the method used in Figure 2:
#
#   distance: Euclidean
#   linkage:  Ward's minimum-variance method (Ward.D2)

patient_distance <- dist(
  t(heatmap_matrix),
  method = "euclidean"
)

patient_clustering <- hclust(
  patient_distance,
  method = "ward.D2"
)


# ------------------------------------------------------------------------------
# 16. Molecular-diagnosis annotation
# ------------------------------------------------------------------------------

heatmap_annotation <- data.frame(
  `Molecular diagnosis` = factor(
    dat$molecular_diagnosis,
    levels = c(0, 1),
    labels = c("No", "Yes")
  )
)

rownames(heatmap_annotation) <- colnames(heatmap_matrix)


# ------------------------------------------------------------------------------
# 17. Save descriptive objects
# ------------------------------------------------------------------------------

dir.create(
  "data/derived",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  descriptive_data,
  "data/derived/descriptive_data.rds"
)

saveRDS(
  table1,
  "data/derived/table1.rds"
)

saveRDS(
  cohort_flow,
  "data/derived/cohort_flow.rds"
)

saveRDS(
  heatmap_matrix,
  "data/derived/heatmap_matrix.rds"
)

saveRDS(
  patient_clustering,
  "data/derived/patient_clustering.rds"
)

saveRDS(
  heatmap_annotation,
  "data/derived/heatmap_annotation.rds"
)


message("Descriptive analysis completed successfully.")
message("Clinical domain counts:")
print(observed_domain_counts)