# ==============================================================================
# 08_supplementary_distributions.R
#
# Study:
#   Clinical predictors of genetic diagnoses in adults evaluated for
#   inborn errors of immunity
#
# Purpose:
#   Reproduce the descriptive distributions presented in:
#
#     Figure S3: Age at symptom onset according to molecular diagnostic status
#                in the full cohort, non-PAD subgroup, and PAD subgroup.
#
#     Figure S4: IgG1 concentration according to molecular diagnostic status
#                in the full cohort, non-PAD subgroup, and PAD subgroup.
#
# Input:
#   data/derived/analysis_data.rds
#
# Outputs:
#   results/descriptive/figureS3_age_onset_summary.csv
#   results/descriptive/figureS4_igg1_summary.csv
#   results/figures/figureS3_age_onset.pdf
#   results/figures/figureS4_igg1.pdf
#
# Notes:
#   These figures are descriptive. No additional hypothesis tests are introduced
#   here unless required to reproduce the published manuscript.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tibble)
library(readr)
library(ggplot2)


# ------------------------------------------------------------------------------
# 2. Load analysis dataset
# ------------------------------------------------------------------------------

dat <- readRDS(
  "data/derived/analysis_data.rds"
)


# ==============================================================================
# PART A. CHECK REQUIRED VARIABLES
# ==============================================================================


# ------------------------------------------------------------------------------
# 3. Required variables
# ------------------------------------------------------------------------------

required_vars <- c(
  "molecular_diagnosis",
  "age_onset",
  "igg1",
  "pad"
)


stopifnot(
  all(required_vars %in% names(dat))
)


# ------------------------------------------------------------------------------
# 4. Create common labels
# ------------------------------------------------------------------------------

dat_plot <- dat %>%
  mutate(
    diagnosis = factor(
      molecular_diagnosis,
      levels = c(0, 1),
      labels = c(
        "No molecular diagnosis",
        "Molecular diagnosis"
      )
    )
  )


# ==============================================================================
# PART B. FIGURE S3: AGE AT SYMPTOM ONSET
# ==============================================================================


# ------------------------------------------------------------------------------
# 5. Construct full-cohort and subgroup datasets
# ------------------------------------------------------------------------------

age_all <- dat_plot %>%
  filter(
    !is.na(age_onset),
    !is.na(molecular_diagnosis)
  ) %>%
  mutate(
    subgroup = "All"
  )


age_non_pad <- dat_plot %>%
  filter(
    pad == 0,
    !is.na(age_onset),
    !is.na(molecular_diagnosis)
  ) %>%
  mutate(
    subgroup = "Non-PAD"
  )


age_pad <- dat_plot %>%
  filter(
    pad == 1,
    !is.na(age_onset),
    !is.na(molecular_diagnosis)
  ) %>%
  mutate(
    subgroup = "PAD"
  )


age_plot_data <- bind_rows(
  age_all,
  age_non_pad,
  age_pad
) %>%
  mutate(
    subgroup = factor(
      subgroup,
      levels = c(
        "All",
        "Non-PAD",
        "PAD"
      )
    )
  )


# ------------------------------------------------------------------------------
# 6. Age-at-onset descriptive statistics
# ------------------------------------------------------------------------------

age_summary <- age_plot_data %>%
  group_by(
    subgroup,
    diagnosis
  ) %>%
  summarise(
    n = n(),
    median = median(age_onset),
    q1 = quantile(
      age_onset,
      0.25
    ),
    q3 = quantile(
      age_onset,
      0.75
    ),
    minimum = min(age_onset),
    maximum = max(age_onset),
    .groups = "drop"
  )


cat("\nFIGURE S3: AGE AT SYMPTOM ONSET\n")
cat("--------------------------------\n")

print(
  age_summary,
  n = Inf
)


# ------------------------------------------------------------------------------
# 7. Verify full-cohort age results
# ------------------------------------------------------------------------------

age_all_summary <- age_summary %>%
  filter(
    subgroup == "All"
  )


stopifnot(
  sum(age_all_summary$n) == 124
)


# Previously reproduced Figure 3 medians:
#
#   No molecular diagnosis: 24 years
#   Molecular diagnosis:     4 years

stopifnot(
  age_all_summary$median[
    age_all_summary$diagnosis ==
      "No molecular diagnosis"
  ] == 24,
  
  age_all_summary$median[
    age_all_summary$diagnosis ==
      "Molecular diagnosis"
  ] == 4
)


# ------------------------------------------------------------------------------
# 8. Figure S3
# ------------------------------------------------------------------------------

figureS3 <- ggplot(
  age_plot_data,
  aes(
    x = diagnosis,
    y = age_onset,
    fill = diagnosis
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.5,
    colour = NA
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    size = 1.4,
    alpha = 0.6
  ) +
  facet_wrap(
    ~ subgroup,
    nrow = 1
  ) +
  labs(
    x = NULL,
    y = "Age at symptom onset (years)"
  ) +
  guides(
    fill = "none"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )


# ==============================================================================
# PART C. FIGURE S4: IgG1
# ==============================================================================


# ------------------------------------------------------------------------------
# 9. Inspect IgG1 availability
# ------------------------------------------------------------------------------

igg1_availability <- dat_plot %>%
  summarise(
    total_n = n(),
    available_igg1 =
      sum(!is.na(igg1)),
    missing_igg1 =
      sum(is.na(igg1))
  )


cat("\nIgG1 AVAILABILITY\n")
cat("-----------------\n")

print(igg1_availability)


cat("\nIgG1 AVAILABILITY BY PAD STATUS\n")
cat("--------------------------------\n")

print(
  dat_plot %>%
    group_by(
      pad
    ) %>%
    summarise(
      n = n(),
      available_igg1 =
        sum(!is.na(igg1)),
      missing_igg1 =
        sum(is.na(igg1)),
      .groups = "drop"
    )
)


# ------------------------------------------------------------------------------
# 10. Construct full-cohort and subgroup datasets
# ------------------------------------------------------------------------------

igg1_all <- dat_plot %>%
  filter(
    !is.na(igg1),
    !is.na(molecular_diagnosis)
  ) %>%
  mutate(
    subgroup = "All"
  )


igg1_non_pad <- dat_plot %>%
  filter(
    pad == 0,
    !is.na(igg1),
    !is.na(molecular_diagnosis)
  ) %>%
  mutate(
    subgroup = "Non-PAD"
  )


igg1_pad <- dat_plot %>%
  filter(
    pad == 1,
    !is.na(igg1),
    !is.na(molecular_diagnosis)
  ) %>%
  mutate(
    subgroup = "PAD"
  )


igg1_plot_data <- bind_rows(
  igg1_all,
  igg1_non_pad,
  igg1_pad
) %>%
  mutate(
    subgroup = factor(
      subgroup,
      levels = c(
        "All",
        "Non-PAD",
        "PAD"
      )
    )
  )


# ------------------------------------------------------------------------------
# 11. IgG1 descriptive statistics
# ------------------------------------------------------------------------------

igg1_summary <- igg1_plot_data %>%
  group_by(
    subgroup,
    diagnosis
  ) %>%
  summarise(
    n = n(),
    median = median(igg1),
    q1 = quantile(
      igg1,
      0.25
    ),
    q3 = quantile(
      igg1,
      0.75
    ),
    minimum = min(igg1),
    maximum = max(igg1),
    .groups = "drop"
  )


cat("\nFIGURE S4: IgG1 CONCENTRATION\n")
cat("-----------------------------\n")

print(
  igg1_summary,
  n = Inf
)


# ------------------------------------------------------------------------------
# 12. Figure S4
# ------------------------------------------------------------------------------

figureS4 <- ggplot(
  igg1_plot_data,
  aes(
    x = diagnosis,
    y = igg1,
    fill = diagnosis
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.5,
    colour = NA
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    size = 1.4,
    alpha = 0.6
  ) +
  facet_wrap(
    ~ subgroup,
    nrow = 1
  ) +
  labs(
    x = NULL,
    y = "IgG1 (g/L)"
  ) +
  guides(
    fill = "none"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )


# ==============================================================================
# PART D. SAVE RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 13. Create output directories
# ------------------------------------------------------------------------------

dir.create(
  "results/descriptive",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 14. Save descriptive statistics
# ------------------------------------------------------------------------------

write_csv(
  age_summary,
  "results/descriptive/figureS3_age_onset_summary.csv"
)


write_csv(
  igg1_summary,
  "results/descriptive/figureS4_igg1_summary.csv"
)


write_csv(
  igg1_availability,
  "results/descriptive/figureS4_igg1_availability.csv"
)


# ------------------------------------------------------------------------------
# 15. Save figures
# ------------------------------------------------------------------------------

ggsave(
  filename =
    "results/figures/figureS3_age_onset.pdf",
  plot = figureS3,
  width = 10,
  height = 4.5,
  units = "in"
)


ggsave(
  filename =
    "results/figures/figureS4_igg1.pdf",
  plot = figureS4,
  width = 10,
  height = 4.5,
  units = "in"
)


message("Supplementary distribution analyses completed successfully.")
message(
  "Figure S3 full-cohort observations: ",
  nrow(age_all)
)
message(
  "Figure S4 full-cohort observations with IgG1: ",
  nrow(igg1_all)
)