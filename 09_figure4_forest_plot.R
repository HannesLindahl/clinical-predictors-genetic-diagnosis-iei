# ==============================================================================
# 09_figure4_forest_plot.R
#
# Study:
#   Clinical predictors of genetic diagnoses in adults evaluated for
#   inborn errors of immunity
#
# Purpose:
#   Reproduce Figure 4, showing adjusted odds ratios and 95% confidence
#   intervals from the primary multivariable logistic regression model.
#
# Input:
#   data/derived/primary_model.rds
#
# Outputs:
#   results/figures/figure4_forest_plot.pdf
#   results/figures/figure4_forest_plot.png
#   results/models/figure4_forest_plot_data.csv
#
# Notes:
#   The primary model includes:
#
#     molecular diagnosis ~ age at onset + male sex + IgG1 + eczema
#
#   Figure 4 displays male sex, IgG1, and eczema. Age at symptom onset is
#   included in the fitted model but omitted from the forest plot.
#
#   Odds ratios are obtained by exponentiating model coefficients.
#   Confidence intervals are profile-likelihood 95% confidence intervals,
#   consistent with the primary-model analysis.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tibble)
library(readr)
library(ggplot2)


# ------------------------------------------------------------------------------
# 2. Load primary model
# ------------------------------------------------------------------------------

primary_model <- readRDS(
  "data/derived/primary_model.rds"
)


stopifnot(
  inherits(
    primary_model,
    "glm"
  )
)


# ==============================================================================
# PART A. EXTRACT MODEL RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 3. Model coefficients and confidence intervals
# ------------------------------------------------------------------------------

coef_table <- summary(
  primary_model
)$coefficients


ci_logit <- suppressMessages(
  confint(
    primary_model
  )
)


model_results <- tibble(
  term = rownames(coef_table),
  
  estimate_log_odds =
    coef_table[, "Estimate"],
  
  standard_error =
    coef_table[, "Std. Error"],
  
  p_value =
    coef_table[, "Pr(>|z|)"],
  
  odds_ratio =
    exp(
      coef_table[, "Estimate"]
    ),
  
  conf_low =
    exp(
      ci_logit[, 1]
    ),
  
  conf_high =
    exp(
      ci_logit[, 2]
    )
)


cat("\nPRIMARY MODEL RESULTS USED FOR FIGURE 4\n")
cat("---------------------------------------\n")

print(
  model_results,
  n = Inf
)


# ==============================================================================
# PART B. PREPARE FIGURE 4 DATA
# ==============================================================================


# ------------------------------------------------------------------------------
# 4. Select variables displayed in Figure 4
# ------------------------------------------------------------------------------

figure4_data <- model_results %>%
  filter(
    term %in% c(
      "male",
      "igg1",
      "eczema_binary"
    )
  ) %>%
  mutate(
    variable = case_when(
      
      term == "male" ~
        "Male sex",
      
      term == "igg1" ~
        "IgG1 (per g/L)",
      
      term == "eczema_binary" ~
        "Eczema",
      
      TRUE ~ term
    )
  )


# ------------------------------------------------------------------------------
# 5. Set plotting order
# ------------------------------------------------------------------------------

figure4_data <- figure4_data %>%
  mutate(
    variable = factor(
      variable,
      levels = c(
        "Eczema",
        "Male sex",
        "IgG1 (per g/L)"
      )
    )
  ) %>%
  arrange(
    variable
  )


cat("\nFIGURE 4 PLOTTED VALUES\n")
cat("-----------------------\n")

print(
  figure4_data %>%
    select(
      variable,
      odds_ratio,
      conf_low,
      conf_high,
      p_value
    ),
  n = Inf
)


# ==============================================================================
# PART C. CHECK AGAINST REPORTED FIGURE/TABLE VALUES
# ==============================================================================


# ------------------------------------------------------------------------------
# 6. Reported adjusted odds ratios
# ------------------------------------------------------------------------------

# Figure 4 is based on the same primary multivariable model as Table 2.
#
# Reported adjusted odds ratios:
#
#   Male sex: 2.872
#   IgG1:     1.263
#   Eczema:   6.234

reported_or <- c(
  male = 2.872,
  igg1 = 1.263,
  eczema_binary = 6.234
)


figure4_comparison <- figure4_data %>%
  mutate(
    term = as.character(term),
    reported_or =
      unname(
        reported_or[term]
      ),
    absolute_difference =
      abs(
        odds_ratio -
          reported_or
      )
  ) %>%
  select(
    term,
    variable,
    reported_or,
    odds_ratio,
    absolute_difference
  )


cat("\nCOMPARISON WITH REPORTED FIGURE 4 / TABLE 2 VALUES\n")
cat("-------------------------------------------------\n")

print(
  figure4_comparison,
  n = Inf
)


# ==============================================================================
# PART D. CREATE FOREST PLOT
# ==============================================================================


# ------------------------------------------------------------------------------
# 7. Figure 4
# ------------------------------------------------------------------------------

figure4 <- ggplot(
  figure4_data,
  aes(
    x = odds_ratio,
    y = variable
  )
) +
  
  # Null value for an odds ratio.
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  # 95% confidence intervals.
  geom_errorbarh(
    aes(
      xmin = conf_low,
      xmax = conf_high
    ),
    height = 0.15,
    linewidth = 0.6
  ) +
  
  # Point estimates.
  geom_point(
    size = 3
  ) +
  
  # Logarithmic scale is appropriate for odds ratios.
  scale_x_log10(
    breaks = c(
      0.5,
      1,
      2,
      5,
      10,
      20,
      50
    )
  ) +
  
  labs(
    x = "Adjusted odds ratio (95% CI)",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  )


# ==============================================================================
# PART E. SAVE OUTPUTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 8. Create output directories
# ------------------------------------------------------------------------------

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/models",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 9. Save plotted numerical values
# ------------------------------------------------------------------------------

write_csv(
  figure4_data %>%
    mutate(
      variable =
        as.character(variable)
    ) %>%
    select(
      variable,
      term,
      odds_ratio,
      conf_low,
      conf_high,
      p_value
    ),
  "results/models/figure4_forest_plot_data.csv"
)


# ------------------------------------------------------------------------------
# 10. Save manuscript comparison
# ------------------------------------------------------------------------------

write_csv(
  figure4_comparison,
  "results/models/figure4_manuscript_comparison.csv"
)


# ------------------------------------------------------------------------------
# 11. Save Figure 4
# ------------------------------------------------------------------------------

ggsave(
  filename =
    "results/figures/figure4_forest_plot.pdf",
  plot = figure4,
  width = 6.5,
  height = 3.5,
  units = "in"
)


ggsave(
  filename =
    "results/figures/figure4_forest_plot.png",
  plot = figure4,
  width = 6.5,
  height = 3.5,
  units = "in",
  dpi = 300
)


message("Figure 4 forest plot completed successfully.")
message(
  "Variables displayed: ",
  nrow(figure4_data)
)