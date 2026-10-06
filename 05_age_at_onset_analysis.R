# ==============================================================================
# 05_age_at_onset_analysis.R
#
# Study:
#   Clinical predictors of genetic diagnoses in adults evaluated for
#   inborn errors of immunity
#
# Purpose:
#   Reproduce the analyses and graphical results presented in Figure 3:
#
#     A. Distribution of age at symptom onset by molecular diagnostic status
#     B. Predicted probability of molecular diagnosis according to age at onset
#     C. Observed diagnostic yield across age-at-onset categories
#
# Input:
#   data/derived/analysis_data.rds
#
# Outputs:
#   data/derived/age_onset_model.rds
#   results/models/age_onset_model_results.csv
#   results/figures/figure3_age_at_onset.pdf
#
# Notes:
#   The logistic regression in this script evaluates age at symptom onset as
#   the sole predictor of molecular diagnosis. It is therefore distinct from
#   the adjusted multivariable model reported in Table 2.
#
#   Age at onset is modelled as a continuous variable in the regression.
#   Age categories are used only for descriptive presentation in Figure 3C.
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
# PART A. PREPARE AGE-AT-ONSET DATA
# ==============================================================================


# ------------------------------------------------------------------------------
# 3. Check required variables
# ------------------------------------------------------------------------------

required_vars <- c(
  "molecular_diagnosis",
  "age_onset",
  "age_onset_group"
)

stopifnot(
  all(required_vars %in% names(dat))
)


# ------------------------------------------------------------------------------
# 4. Analysis dataset
# ------------------------------------------------------------------------------

age_data <- dat %>%
  select(
    molecular_diagnosis,
    age_onset,
    age_onset_group
  ) %>%
  filter(
    !is.na(molecular_diagnosis),
    !is.na(age_onset)
  ) %>%
  mutate(
    diagnosis = factor(
      molecular_diagnosis,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )


cat("\nAGE-AT-ONSET ANALYSIS POPULATION\n")
cat("--------------------------------\n")
cat("Participants:", nrow(age_data), "\n")
cat(
  "Molecular diagnoses:",
  sum(age_data$molecular_diagnosis == 1),
  "\n"
)
cat(
  "No molecular diagnosis:",
  sum(age_data$molecular_diagnosis == 0),
  "\n"
)


# Age at onset should be available for the complete study cohort.

stopifnot(
  nrow(age_data) == 124,
  sum(age_data$molecular_diagnosis == 1) == 55,
  sum(age_data$molecular_diagnosis == 0) == 69
)


# ==============================================================================
# PART B. DESCRIPTIVE AGE DISTRIBUTION
# ==============================================================================


# ------------------------------------------------------------------------------
# 5. Summarize age at onset by molecular diagnostic status
# ------------------------------------------------------------------------------

age_summary <- age_data %>%
  group_by(diagnosis) %>%
  summarise(
    n = n(),
    median = median(age_onset),
    q1 = quantile(age_onset, 0.25),
    q3 = quantile(age_onset, 0.75),
    minimum = min(age_onset),
    maximum = max(age_onset),
    .groups = "drop"
  )


cat("\nAGE AT ONSET BY MOLECULAR DIAGNOSTIC STATUS\n")
cat("-------------------------------------------\n")

print(age_summary)


# ==============================================================================
# PART C. UNIVARIABLE LOGISTIC REGRESSION
# ==============================================================================


# ------------------------------------------------------------------------------
# 6. Fit age-only logistic regression model
# ------------------------------------------------------------------------------

age_model <- glm(
  molecular_diagnosis ~ age_onset,
  data = age_data,
  family = binomial(link = "logit")
)


stopifnot(age_model$converged)


cat("\nAGE-ONLY LOGISTIC REGRESSION\n")
cat("----------------------------\n")

print(
  summary(age_model)
)


# ------------------------------------------------------------------------------
# 7. Extract age coefficient
# ------------------------------------------------------------------------------

coef_table <- summary(age_model)$coefficients

age_beta <- coef(age_model)[["age_onset"]]

age_se <- coef_table[
  "age_onset",
  "Std. Error"
]

age_p <- coef_table[
  "age_onset",
  "Pr(>|z|)"
]


# Profile-likelihood confidence interval for the regression coefficient.

age_ci_logit <- suppressMessages(
  confint(
    age_model,
    parm = "age_onset"
  )
)


# ------------------------------------------------------------------------------
# 8. Express association per one-year increase
# ------------------------------------------------------------------------------

age_effect_1y <- tibble(
  scale = "Per 1-year increase",
  odds_ratio = exp(age_beta),
  conf_low = exp(age_ci_logit[1]),
  conf_high = exp(age_ci_logit[2]),
  p_value = age_p
)


cat("\nAGE EFFECT PER 1 YEAR\n")
cat("---------------------\n")

print(age_effect_1y)


# ------------------------------------------------------------------------------
# 9. Express association per 10-year increase
# ------------------------------------------------------------------------------

# In logistic regression, multiplying a continuous predictor by a change of k
# units corresponds to:
#
#   OR(k units) = exp(k * beta)
#
# The 10-year estimate therefore represents the same fitted age-only model,
# expressed on a 10-year rather than one-year scale.


age_effect_10y <- tibble(
  scale = "Per 10-year increase",
  odds_ratio = exp(10 * age_beta),
  conf_low = exp(10 * age_ci_logit[1]),
  conf_high = exp(10 * age_ci_logit[2]),
  p_value = age_p
)


cat("\nAGE EFFECT PER 10 YEARS\n")
cat("-----------------------\n")

print(age_effect_10y)


# ------------------------------------------------------------------------------
# 10. Compare with Figure 3B annotation
# ------------------------------------------------------------------------------

# Figure 3B reports:
#
#   OR per 10-year increase = 0.45
#   95% CI = 0.31–0.61
#
# Compare the reproduced age-only model with these published values.


figure3b_reported <- tibble(
  odds_ratio = 0.45,
  conf_low = 0.31,
  conf_high = 0.61
)


figure3b_comparison <- tibble(
  statistic = c(
    "OR",
    "95% CI lower",
    "95% CI upper"
  ),
  reported = c(
    figure3b_reported$odds_ratio,
    figure3b_reported$conf_low,
    figure3b_reported$conf_high
  ),
  reproduced = c(
    age_effect_10y$odds_ratio,
    age_effect_10y$conf_low,
    age_effect_10y$conf_high
  )
) %>%
  mutate(
    absolute_difference =
      abs(reported - reproduced)
  )


cat("\nCOMPARISON WITH FIGURE 3B\n")
cat("-------------------------\n")

print(figure3b_comparison)


# ==============================================================================
# PART D. PREDICTED PROBABILITY CURVE
# ==============================================================================


# ------------------------------------------------------------------------------
# 11. Generate prediction grid
# ------------------------------------------------------------------------------

prediction_data <- tibble(
  age_onset = seq(
    min(age_data$age_onset),
    max(age_data$age_onset),
    length.out = 300
  )
)


# Obtain predictions on the logit scale so that confidence intervals are
# calculated on the model's linear-predictor scale before transformation to
# probabilities.

prediction_link <- predict(
  age_model,
  newdata = prediction_data,
  type = "link",
  se.fit = TRUE
)


prediction_data <- prediction_data %>%
  mutate(
    fit_link = as.numeric(prediction_link$fit),
    se_link = as.numeric(prediction_link$se.fit),
    
    lower_link = fit_link - 1.96 * se_link,
    upper_link = fit_link + 1.96 * se_link,
    
    predicted_probability = plogis(fit_link),
    conf_low = plogis(lower_link),
    conf_high = plogis(upper_link)
  )


# ==============================================================================
# PART E. OBSERVED DIAGNOSTIC YIELD BY AGE CATEGORY
# ==============================================================================


# ------------------------------------------------------------------------------
# 12. Calculate observed yield
# ------------------------------------------------------------------------------

age_group_yield <- age_data %>%
  group_by(age_onset_group) %>%
  summarise(
    n = n(),
    diagnoses = sum(molecular_diagnosis),
    diagnostic_yield = mean(molecular_diagnosis),
    .groups = "drop"
  )


cat("\nOBSERVED DIAGNOSTIC YIELD BY AGE GROUP\n")
cat("--------------------------------------\n")

print(age_group_yield)


# ------------------------------------------------------------------------------
# 13. Verify Figure 3C counts
# ------------------------------------------------------------------------------

expected_n <- c(
  38L,
  37L,
  21L,
  28L
)

expected_diagnoses <- c(
  32L,
  14L,
  4L,
  5L
)


stopifnot(
  identical(
    age_group_yield$n,
    expected_n
  ),
  identical(
    age_group_yield$diagnoses,
    expected_diagnoses
  )
)


# Label used in Figure 3C.

age_group_yield <- age_group_yield %>%
  mutate(
    count_label = paste0(
      diagnoses,
      "/",
      n
    )
  )


# ==============================================================================
# PART F. FIGURE 3
# ==============================================================================


# ------------------------------------------------------------------------------
# 14. Panel A: age-at-onset distributions
# ------------------------------------------------------------------------------

figure3a <- ggplot(
  age_data,
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
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 1.7,
    alpha = 0.7
  ) +
  labs(
    x = "Molecular diagnosis",
    y = "Age at symptom onset (years)"
  ) +
  guides(
    fill = "none"
  ) +
  theme_classic(
    base_size = 13
  )


# ------------------------------------------------------------------------------
# 15. Panel B: predicted probability from age-only logistic regression
# ------------------------------------------------------------------------------

figure3b <- ggplot(
  prediction_data,
  aes(
    x = age_onset,
    y = predicted_probability
  )
) +
  geom_ribbon(
    aes(
      ymin = conf_low,
      ymax = conf_high
    ),
    alpha = 0.2
  ) +
  geom_line(
    linewidth = 1
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  labs(
    x = "Age at symptom onset (years)",
    y = "Predicted probability of molecular diagnosis"
  ) +
  theme_classic(
    base_size = 13
  )


# ------------------------------------------------------------------------------
# 16. Panel C: observed diagnostic yield by age category
# ------------------------------------------------------------------------------

figure3c <- ggplot(
  age_group_yield,
  aes(
    x = age_onset_group,
    y = diagnostic_yield
  )
) +
  geom_col(
    width = 0.7
  ) +
  geom_text(
    aes(
      label = count_label
    ),
    vjust = -0.5,
    size = 4
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  labs(
    x = "Age at symptom onset (years)",
    y = "Molecular diagnostic yield"
  ) +
  theme_classic(
    base_size = 13
  )


# ==============================================================================
# PART G. COMBINE FIGURE
# ==============================================================================


# ------------------------------------------------------------------------------
# 17. Assemble Figure 3
# ------------------------------------------------------------------------------

# patchwork is used only for figure assembly.

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop(
    "Package 'patchwork' is required to assemble Figure 3."
  )
}


figure3 <- figure3a +
  figure3b +
  figure3c +
  patchwork::plot_layout(
    widths = c(1, 1.25, 1)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )


# ==============================================================================
# PART H. SAVE RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 18. Create output directories
# ------------------------------------------------------------------------------

dir.create(
  "results/models",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 19. Save analysis objects
# ------------------------------------------------------------------------------

saveRDS(
  age_model,
  "data/derived/age_onset_model.rds"
)


# ------------------------------------------------------------------------------
# 20. Save numerical results
# ------------------------------------------------------------------------------

write_csv(
  age_summary,
  "results/models/age_onset_descriptive_summary.csv"
)

write_csv(
  age_effect_1y,
  "results/models/age_onset_effect_1y.csv"
)

write_csv(
  age_effect_10y,
  "results/models/age_onset_effect_10y.csv"
)

write_csv(
  figure3b_comparison,
  "results/models/figure3b_comparison.csv"
)

write_csv(
  age_group_yield,
  "results/models/age_group_diagnostic_yield.csv"
)


# ------------------------------------------------------------------------------
# 21. Save Figure 3
# ------------------------------------------------------------------------------

ggsave(
  filename = "results/figures/figure3_age_at_onset.pdf",
  plot = figure3,
  width = 11,
  height = 4.5,
  units = "in"
)


message("Age-at-onset analysis completed successfully.")
message("Participants: ", nrow(age_data))