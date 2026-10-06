# Clinical predictors of genetic diagnoses in adults evaluated for inborn errors of immunity

## Overview

This repository contains the analysis code and reproducibility documentation for the study:

**Clinical predictors of genetic diagnoses in adults evaluated for inborn errors of immunity**

The study investigated clinical and laboratory characteristics associated with obtaining a molecular diagnosis among adults evaluated for suspected inborn errors of immunity (IEI).

The study cohort comprised 124 individuals who underwent genetic testing and met the study inclusion criteria. Of these, 55 received a molecular diagnosis and 69 did not.

The analysis dataset contains the variables required to reproduce the principal descriptive analyses, the reproducible portion of the predictor-selection workflow, regression analyses, sensitivity analyses, and figures reported in the manuscript.

The analysis dataset is not included in this GitHub repository. The public analysis workflow begins with `data/gimpid_clean.csv`, which should be placed in the `data/` directory before running the analysis.

## Associated publication

**Article:** Manuscript under revision.

**Analysis code:**  
https://github.com/HannesLindahl/clinical-predictors-genetic-diagnosis-iei

## Analysis dataset

The analysis workflow uses:

`data/gimpid_clean.csv`

The analysis dataset contains:

- **124 observations** (study participants)
- **81 variables**
- **79 eligible clinical and laboratory predictor variables**
- **1 molecular diagnosis outcome variable**
- **1 IUIS classification variable**

Each row represents one study participant. No participant identifier is included in the analysis dataset.

The dataset contains detailed clinical information from a relatively small patient cohort and is therefore not distributed through this GitHub repository.

### Primary outcome

The primary outcome variable is:

`Relevant genetics`

This indicates whether genetic testing resulted in a molecular diagnosis considered relevant to the participant's clinical phenotype.

The analysis scripts create standardized analysis variables from the source-variable names where required.

### Predictor variables

The analysis dataset includes demographic, clinical, infectious, immunological, autoimmune, inflammatory, malignancy-related, and laboratory variables evaluated in the study.

A separate data dictionary documents the variables, coding, data types, laboratory units where applicable, and representation of missing values.

## Data provenance and processing

The analysis dataset was derived from secure study source data.

Source-data cleaning and preparation were performed in a private preprocessing workflow before creation of the analysis dataset. This included removal of direct and internal identifiers, data-type harmonization, handling of source-specific laboratory notation, and preparation of variables required for the analyses.

The private preprocessing code is not included in this repository because it operates on source data that are not part of the public research package.

The public reproducibility workflow begins with `data/gimpid_clean.csv`.

A detailed description of the processing history, predictor screening, variable-selection workflow, statistical analyses, and reproducibility boundaries is provided in:

`PROCESSING_AND_REPRODUCIBILITY.md`

## Analysis workflow

The analysis scripts are intended to be run sequentially:

1. `01_prepare_analysis_data.R`  
   Imports the analysis dataset, performs analysis-level recoding, defines the molecular diagnosis outcome and key analysis variables, and creates the derived analysis dataset.

2. `02_descriptive_analysis.R`  
   Produces cohort summaries and descriptive analyses, including clinical manifestation domains and the clustered clinical-manifestation heatmap.

3. `03_variable_selection.R`  
   Reconstructs the domain-based predictor prioritization from the 79 eligible predictors available in the analysis dataset and documents the historical predictor-screening stages.

4. `04_multivariable_analysis.R`  
   Fits the primary multivariable logistic regression model.

5. `05_age_at_onset_analysis.R`  
   Performs analyses of age at symptom onset and molecular diagnostic outcome.

6. `06_sensitivity_analyses.R`  
   Performs sensitivity analyses stratified by predominantly antibody deficiency (PAD) status.

7. `07_iuis_diagnostic_yield.R`  
   Calculates molecular diagnostic yield according to IUIS classification.

8. `08_supplementary_distributions.R`  
   Produces supplementary distributions of age at symptom onset and IgG1.

9. `09_figure4_forest_plot.R`  
   Produces the forest plot for the primary adjusted logistic regression model.

10. `10_univariable_analyses.R`  
    Performs the univariable predictor analyses, including Firth logistic regression for celiac disease and NRH/PSVD, consistent with the historical analysis.

The complete workflow (`01` through `10`) has been tested sequentially from a clean R session.

## Predictor screening and variable selection

The original predictor-screening workflow began with **109 candidate predictor variables**.

Before the domain-based variable-selection procedure:

- 6 complement-related variables were excluded;
- 22 variables were excluded because more than 20% of observations were missing; and
- 2 variables were excluded because they had no variability in the study cohort.

This resulted in **79 eligible predictors**.

The analysis dataset contains these 79 eligible predictors rather than all 109 variables considered during the historical preprocessing stage. Therefore, the exclusions that reduced the original 109 predictors to 79 are documented for provenance but cannot be independently recomputed from the analysis dataset.

The subsequent predictor-prioritization and statistical analyses based on the 79 eligible predictors are implemented in the public analysis code.

Further details are provided in `PROCESSING_AND_REPRODUCIBILITY.md`.

## Missing data

Missing values are stored as blank fields in the CSV file and are read as `NA` by the analysis workflow.

Missingness varies between variables. Predictor eligibility in the original analysis required no more than 20% missing observations.

Regression analyses use complete observations for the variables included in the respective model. Consequently, the analysis sample size may differ between analyses.

For descriptive construction of clinical manifestation domains, missing values in the relevant manifestation variables are treated as absence for domain construction, as documented in the analysis code and `PROCESSING_AND_REPRODUCIBILITY.md`.

## Reproducing the analyses

The analysis requires R and the packages loaded by the individual analysis scripts.

To reproduce the analyses:

1. Place `gimpid_clean.csv` in the `data/` directory.
2. Start a clean R session in the repository root.
3. Run scripts `01` through `10` sequentially.

Intermediate analysis objects are written to:

`data/derived/`

Figures and other generated outputs are written to the output locations specified by the individual scripts.

## Software environment

All analyses were performed in R.

The complete public analysis workflow (`01`–`10`) was tested sequentially from a clean R session using R version 4.5.0 (2025-04-11) on Windows 11 x64.

The principal R packages attached during the analysis workflow were `dplyr` (1.1.4), `readr` (2.1.5), `tidyr` (1.3.1), `tibble` (3.3.0), `Matrix` (1.7-3), `glmnet` (5.0), `ggplot2` (4.0.3), and `logistf` (1.26.1). Additional package dependencies were loaded through these packages as required.

A complete record of the R session and package versions after successful execution of the workflow is provided in:

`sessionInfo.txt`

## Reproducibility scope

This repository is designed to reproduce the analyses that can be performed from the analysis dataset.

Some aspects of the complete research workflow necessarily precede the analysis dataset and are therefore documented rather than computationally reproduced. These include source-data extraction and cleaning, removal of identifiers, exclusions based on variables not retained in the release dataset, and other processing requiring access to secure source data.

The public analysis code does not require access to individual gene or disease labels that are not necessary for reproducing the reported predictor analyses.

## Data confidentiality

The analysis dataset contains clinical research data from individuals evaluated for inborn errors of immunity. Although direct and internal identifiers have been removed, combinations of detailed clinical characteristics may remain sensitive.

For this reason, the clinical dataset is not distributed through this GitHub repository.

## Ethics

The study was approved by the Regional Ethical Review Board in Stockholm, Sweden, as part of the FUNGEN study (approval number 2011/116-31). An updated approval (2020-00125) was granted by the Swedish Ethical Review Authority.

Further details are provided in the associated publication.

## Citation

The associated manuscript is currently under revision. Publication and dataset citation information will be added when available.

When referring specifically to the analysis code, please use:

https://github.com/HannesLindahl/clinical-predictors-genetic-diagnosis-iei

## License and reuse

No separate license for the analysis code has yet been specified.

The absence of a code license should not be interpreted as permission to access or reuse the underlying clinical research data. The clinical dataset is not included in this repository.

## Contact

For questions regarding the study or analysis code, please contact the corresponding authors through the contact information provided in the associated publication.
