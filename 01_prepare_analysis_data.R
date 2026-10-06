# ==============================================================================
# 01_prepare_analysis_data.R
#
# Study:
#   Clinical predictors of genetic diagnoses in adults evaluated for
#   inborn errors of immunity
#
# Purpose:
#   Import the repository dataset and create the common analysis dataset used
#   by all subsequent scripts.
#
# Input:
#   data/gimpid_clean.csv
#
# Output:
#   data/derived/analysis_data.rds
#
# Notes:
#   The repository dataset is assumed to have already undergone source-data
#   cleaning and de-identification before deposition.
#
#   This script performs analysis-level recoding only. It does not perform
#   variable selection, statistical modelling, or produce manuscript outputs.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(readr)


# ------------------------------------------------------------------------------
# 2. Import repository dataset
# ------------------------------------------------------------------------------

# TEMPORARY:
# During development, replace this path with the current full cleaned working
# dataset. The final repository will read the deposited dataset here.


dat <- read_csv("gimpid_clean.csv",
  show_col_types = FALSE
)


# ------------------------------------------------------------------------------
# 3. Basic integrity checks
# ------------------------------------------------------------------------------

# The published study cohort contains 124 participants.

stopifnot(nrow(dat) == 124)


# ------------------------------------------------------------------------------
# 4. Rename variables used throughout the analysis
# ------------------------------------------------------------------------------

# Source variable names are converted to concise analysis names here.
# Variables not listed remain unchanged.

dat <- dat %>%
  rename(
    genetic_diagnosis = `Relevant genetics`,
    age_onset = `Age at debut of symptoms`,
    age_genetic_test = `Age at genetic test`,
    sex = Sex,
    consanguinity = Consanguinity,
    familial = `Familial cases`,
    eczema = eczema,
    igg1 = `IgG1-subclass def`
  )


# ------------------------------------------------------------------------------
# 5. Define outcome
# ------------------------------------------------------------------------------

# Molecular diagnosis:
#   0 = no phenotype-concordant pathogenic/likely pathogenic variant
#   1 = molecular diagnosis achieved

stopifnot(
  all(
    na.omit(unique(dat$genetic_diagnosis)) %in% c("N", "Y")
  )
)

dat <- dat %>%
  mutate(
    molecular_diagnosis = case_when(
      genetic_diagnosis == "Y" ~ 1L,
      genetic_diagnosis == "N" ~ 0L,
      TRUE ~ NA_integer_
    )
  )


# Outcome should be available for every participant.

stopifnot(!anyNA(dat$molecular_diagnosis))


# Published cohort:
#   molecular diagnosis = 55
#   no molecular diagnosis = 69

stopifnot(
  sum(dat$molecular_diagnosis == 1) == 55,
  sum(dat$molecular_diagnosis == 0) == 69
)


# ------------------------------------------------------------------------------
# 6. Recode sex
# ------------------------------------------------------------------------------

# Inspect allowed source values before recoding.

stopifnot(
  all(
    na.omit(unique(dat$sex)) %in% c("M", "F")
  )
)

dat <- dat %>%
  mutate(
    male = case_when(
      sex == "M" ~ 1L,
      sex == "F" ~ 0L,
      TRUE ~ NA_integer_
    )
  )


# Published cohort contains 57 male participants.

stopifnot(sum(dat$male == 1, na.rm = TRUE) == 57)


# ------------------------------------------------------------------------------
# 7. Recode binary clinical variables
# ------------------------------------------------------------------------------

# Helper for variables recorded as Y/N.
#
# Missing values remain missing. Importantly, this function does NOT interpret
# missing values as absence of a manifestation.

yn_to_binary <- function(x) {
  
  case_when(
    x == "Y" ~ 1L,
    x == "N" ~ 0L,
    is.na(x) ~ NA_integer_,
    TRUE ~ NA_integer_
  )
}


# Variables required directly in the primary regression model.

dat <- dat %>%
  mutate(
    eczema_binary = yn_to_binary(eczema)
  )


# ------------------------------------------------------------------------------
# 8. Check primary continuous variables
# ------------------------------------------------------------------------------

stopifnot(is.numeric(dat$age_onset))
stopifnot(is.numeric(dat$igg1))


# Do not remove observations with missing values here.
# Complete-case populations are defined separately for each statistical
# analysis that requires them.


# ------------------------------------------------------------------------------
# 9. Define PAD status
# ------------------------------------------------------------------------------

# IUIS classification is already present in the repository dataset.
#
# PAD = predominantly antibody deficiencies.

if ("iuis_category" %in% names(dat)) {
  
  dat <- dat %>%
    mutate(
      pad = case_when(
        iuis_category == "Predominantly antibody deficiencies" ~ 1L,
        !is.na(iuis_category) ~ 0L,
        TRUE ~ NA_integer_
      )
    )
}


# ------------------------------------------------------------------------------
# 10. Create derived variables used for graphical presentation
# ------------------------------------------------------------------------------

# Age-at-onset categories used in Figure 3C.
#
# These categories are descriptive only. Age at onset remains continuous in
# the primary regression model.

dat <- dat %>%
  mutate(
    age_onset_group = cut(
      age_onset,
      breaks = c(-Inf, 5, 20, 30, Inf),
      labels = c("0–5", "6–20", "21–30", ">30"),
      right = TRUE
    )
  )


# ------------------------------------------------------------------------------
# 11. Final integrity checks
# ------------------------------------------------------------------------------

# Age-at-onset category totals reported in Figure 3C:
#   0–5 years: 38
#   6–20 years: 37
#   21–30 years: 21
#   >30 years: 28

expected_age_groups <- c(
  "0–5" = 38L,
  "6–20" = 37L,
  "21–30" = 21L,
  ">30" = 28L
)

observed_age_groups <- table(dat$age_onset_group)

stopifnot(
  identical(
    as.integer(observed_age_groups[names(expected_age_groups)]),
    unname(expected_age_groups)
  )
)


# Diagnostic yields reported in Figure 3C:
#   0–5:   32/38
#   6–20:  14/37
#   21–30:  4/21
#   >30:    5/28

expected_diagnoses <- c(
  "0–5" = 32L,
  "6–20" = 14L,
  "21–30" = 4L,
  ">30" = 5L
)

observed_diagnoses <- dat %>%
  group_by(age_onset_group) %>%
  summarise(
    diagnoses = sum(molecular_diagnosis),
    .groups = "drop"
  )

stopifnot(
  identical(
    observed_diagnoses$diagnoses,
    unname(expected_diagnoses)
  )
)


# ------------------------------------------------------------------------------
# 12. Save common analysis dataset
# ------------------------------------------------------------------------------

dir.create(
  "data/derived",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  dat,
  "data/derived/analysis_data.rds"
)


message("Analysis dataset created successfully.")
message("Participants: ", nrow(dat))
message("Variables: ", ncol(dat))