# ==============================================================================
# 07_iuis_diagnostic_yield.R
#
# Study:
#   Clinical predictors of genetic diagnoses in adults evaluated for
#   inborn errors of immunity
#
# Purpose:
#   Reproduce the molecular diagnostic yield by IUIS disease category
#   presented in Figure S2.
#
# Input:
#   data/derived/analysis_data.rds
#
# Outputs:
#   results/descriptive/iuis_diagnostic_yield.csv
#   results/figures/figureS2_iuis_diagnostic_yield.pdf
#
# Notes:
#   Diagnostic yield is calculated descriptively as the proportion of
#   participants with a molecular diagnosis within each IUIS category.
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
# PART A. CHECK IUIS VARIABLE
# ==============================================================================


# ------------------------------------------------------------------------------
# 3. Required variables
# ------------------------------------------------------------------------------

required_vars <- c(
  "molecular_diagnosis",
  "iuis_category"
)


stopifnot(
  all(required_vars %in% names(dat))
)


# ------------------------------------------------------------------------------
# 4. Inspect IUIS categories
# ------------------------------------------------------------------------------

cat("\nIUIS CATEGORIES IN ANALYSIS DATASET\n")
cat("-----------------------------------\n")

print(
  dat %>%
    count(
      iuis_category,
      sort = TRUE
    ),
  n = Inf
)


# Check for missing category assignments.

n_missing_iuis <- sum(
  is.na(dat$iuis_category) |
    dat$iuis_category == ""
)


cat(
  "\nParticipants without an IUIS category:",
  n_missing_iuis,
  "\n"
)


# ==============================================================================
# PART B. DIAGNOSTIC YIELD BY IUIS CATEGORY
# ==============================================================================


# ------------------------------------------------------------------------------
# 5. Calculate diagnostic yield
# ------------------------------------------------------------------------------

iuis_yield <- dat %>%
  filter(
    !is.na(iuis_category),
    iuis_category != ""
  ) %>%
  group_by(
    iuis_category
  ) %>%
  summarise(
    n = n(),
    diagnoses = sum(
      molecular_diagnosis == 1
    ),
    no_diagnosis = sum(
      molecular_diagnosis == 0
    ),
    diagnostic_yield =
      diagnoses / n,
    .groups = "drop"
  ) %>%
  arrange(
    desc(diagnostic_yield),
    desc(n)
  )


cat("\nDIAGNOSTIC YIELD BY IUIS CATEGORY\n")
cat("---------------------------------\n")

print(
  iuis_yield,
  n = Inf
)


# ------------------------------------------------------------------------------
# 6. Overall checks
# ------------------------------------------------------------------------------

stopifnot(
  sum(iuis_yield$n) ==
    nrow(dat) - n_missing_iuis,
  
  sum(iuis_yield$diagnoses) ==
    sum(
      dat$molecular_diagnosis == 1 &
        !is.na(dat$iuis_category) &
        dat$iuis_category != ""
    )
)


# ==============================================================================
# PART C. LABELS FOR FIGURE S2
# ==============================================================================


# ------------------------------------------------------------------------------
# 7. Create figure labels
# ------------------------------------------------------------------------------

iuis_yield <- iuis_yield %>%
  mutate(
    diagnostic_yield_percent =
      100 * diagnostic_yield,
    
    yield_label = paste0(
      round(
        diagnostic_yield_percent
      ),
      "%"
    ),
    
    n_label = paste0(
      "n = ",
      n
    )
  )


# Preserve ordering by diagnostic yield, but display Uncategorized last.

plot_order <- c(
  iuis_yield$iuis_category[
    iuis_yield$iuis_category != "Uncategorized"
  ],
  "Uncategorized"
)

iuis_yield <- iuis_yield %>%
  mutate(
    iuis_category_plot = factor(
      iuis_category,
      levels = rev(plot_order)
    )
  )


# ==============================================================================
# PART D. FIGURE S2
# ==============================================================================


# ------------------------------------------------------------------------------
# 8. Plot diagnostic yield
# ------------------------------------------------------------------------------

figureS2 <- ggplot(
  iuis_yield,
  aes(
    x = diagnostic_yield_percent,
    y = iuis_category_plot
  )
) +
  geom_col(
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(
        yield_label,
        " (",
        diagnoses,
        "/",
        n,
        ")"
      )
    ),
    hjust = -0.1,
    size = 3.5
  ) +
  scale_x_continuous(
    limits = c(0, 115),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    labels = function(x) {
      paste0(
        x,
        "%"
      )
    },
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    x = "Molecular diagnostic yield",
    y = NULL
  ) +
  theme_classic(
    base_size = 12
  )


# ==============================================================================
# PART E. SAVE RESULTS
# ==============================================================================


# ------------------------------------------------------------------------------
# 9. Create output directories
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
# 10. Save numerical results
# ------------------------------------------------------------------------------

write_csv(
  iuis_yield %>%
    select(
      iuis_category,
      n,
      diagnoses,
      no_diagnosis,
      diagnostic_yield,
      diagnostic_yield_percent
    ),
  "results/descriptive/iuis_diagnostic_yield.csv"
)


# ------------------------------------------------------------------------------
# 11. Save Figure S2
# ------------------------------------------------------------------------------

ggsave(
  filename =
    "results/figures/figureS2_iuis_diagnostic_yield.pdf",
  plot = figureS2,
  width = 8,
  height = 5.5,
  units = "in"
)


message("IUIS diagnostic-yield analysis completed successfully.")
message(
  "Participants represented: ",
  sum(iuis_yield$n)
)
message(
  "IUIS categories represented: ",
  nrow(iuis_yield)
)