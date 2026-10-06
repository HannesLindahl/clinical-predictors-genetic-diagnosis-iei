# Clinical predictors of genetic diagnoses in adults evaluated for inborn errors of immunity

## Overview

This repository contains the analysis-ready dataset and supporting analysis code for the study:

**Clinical predictors of genetic diagnoses in adults evaluated for inborn errors of immunity**

The study investigated clinical and laboratory characteristics associated with obtaining a molecular diagnosis among adults evaluated for suspected inborn errors of immunity (IEI).

The study cohort comprised 124 individuals who underwent genetic testing and met the study inclusion criteria. Of these, 55 received a molecular diagnosis and 69 did not.

The deposited dataset contains the variables required to reproduce the principal descriptive analyses, the reproducible portion of the predictor-selection workflow, regression analyses, sensitivity analyses, and figures reported in the manuscript.

Direct and internal identifiers have been removed from the deposited dataset. Because the dataset contains detailed clinical information from a relatively small patient cohort, data access and reuse are subject to the conditions specified by the KI Data Repository.

## Associated publication

**Article:**  
[Full citation to be added after publication]

**DOI:**  
[Article DOI]

**Dataset DOI:**  
[KI Data Repository DOI]

**Analysis code:**  
[GitHub repository URL]

## Dataset

The primary analysis dataset is:

`data/gimpid_clean.csv`

It contains:

- **124 observations** (study participants)
- **81 variables**
- **79 eligible clinical and laboratory predictor variables**
- **1 molecular diagnosis outcome variable**
- **1 IUIS classification variable**

Each row represents one study participant.

No participant identifier is included in the deposited dataset.

### Primary outcome

The primary outcome variable is:

`Relevant genetics`

This indicates whether genetic testing resulted in a molecular diagnosis considered relevant to the participant's clinical phenotype.

The analysis scripts create standardized analysis variables from the deposited source-variable names where required.

### Predictor variables

The dataset includes demographic, clinical, infectious, immunological, autoimmune, inflammatory, malignancy-related, and laboratory variables evaluated in the study.

A detailed description of each variable, including coding, data type, units where applicable, and representation of missing values, is provided in the accompanying data dictionary:

`data_dictionary.csv`

[Update filename if the final data dictionary uses another name or format.]

## Data provenance and processing

The deposited dataset is derived from the study's secure source data.

Source-data cleaning and preparation were performed in a private preprocessing workflow before creation of the repository dataset. This included removal of direct/internal identifiers, data-type harmonization, handling of source-specific laboratory notation, application of study eligibility criteria, and preparation of variables required for the analyses.

The private preprocessing code is not included in the public repository because it operates on source data that are not part of the public research package.

The public analysis workflow begins with `data/gimpid_clean.csv`.

A detailed description of the processing history, predictor screening, variable-selection workflow, statistical analyses, and reproducibility boundaries is provided in:

`PROCESSING_AND_REPRODUCIBILITY.md`

## Analysis workflow

The analysis scripts are intended to be run sequentially:

1. `01_prepare_analysis_data.R`  
   Imports the deposited dataset, performs analysis-level recoding, defines the molecular diagnosis outcome and key analysis variables, and creates the derived analysis dataset.

2. `02_descriptive_analysis.R`  
   Produces cohort summaries and descriptive analyses, including clinical manifestation domains and the clustered clinical-manifestation heatmap.

3. `03_variable_selection.R`  
   Reconstructs the domain-based predictor prioritization from the 79 eligible predictors available in the deposited dataset and documents the historical predictor-screening stages.

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
    Performs the univariable predictor analyses, including prespecified use of Firth logistic regression for predictors requiring penalized estimation.

The scripts have been tested by running `01` through `10` sequentially from the deposited analysis dataset.

## Predictor screening and variable selection

The original predictor-screening workflow began with **109 candidate predictor variables**.

Before the domain-based variable-selection procedure:

- 6 complement-related variables were excluded;
- 22 variables were excluded because more than 20% of observations were missing; and
- 2 variables were excluded because they had no variability in the study cohort.

This resulted in **79 eligible predictors**.

The deposited dataset contains these 79 eligible predictors rather than all 109 variables considered during the historical preprocessing stage. Therefore, the exclusions that reduced the original 109 predictors to 79 are documented for provenance but cannot be independently recomputed from the deposited dataset.

The subsequent predictor-selection and statistical analyses based on the 79 eligible predictors are implemented in the public analysis code.

Further details are provided in `PROCESSING_AND_REPRODUCIBILITY.md`.

## Missing data

Missing values are represented as `NA`.

Missingness varies between variables. Predictor eligibility in the original analysis required no more than 20% missing observations.

Regression analyses use complete observations for the variables included in the respective model. Consequently, the analysis sample size may differ between analyses.

For descriptive construction of clinical manifestation domains, missing values in the relevant manifestation variables are handled according to the rules documented in the analysis code and `PROCESSING_AND_REPRODUCIBILITY.md`.

## Reproducing the analyses

The analysis requires R and the packages loaded by the individual analysis scripts.

To reproduce the analyses:

1. Obtain access to the deposited dataset according to the access conditions specified by the KI Data Repository.
2. Place `gimpid_clean.csv` in the `data/` directory.
3. Start a clean R session in the repository root.
4. Run scripts `01` through `10` sequentially.

Intermediate analysis objects are written to:

`data/derived/`

Figures and other generated outputs are written to the output directories specified by the scripts.

[Add exact R version/package environment instructions here if a lockfile, `sessionInfo()`, or package manifest is included in the final repository.]

## Reproducibility scope

The repository is designed to reproduce the analyses that can be performed from the deposited analysis dataset.

Some aspects of the complete research workflow necessarily precede the deposited dataset and are therefore documented rather than computationally reproduced. These include source-data extraction and cleaning, removal of identifiers, exclusions based on variables not retained in the release dataset, and other processing requiring access to the secure source data.

The public analysis code does not require access to individual gene or disease labels that are not necessary for reproducing the reported predictor analyses.

## Data access and confidentiality

The dataset contains clinical research data from individuals evaluated for inborn errors of immunity. Although direct and internal identifiers have been removed, combinations of detailed clinical characteristics may remain sensitive.

Data access, permitted reuse, and any applicable disclosure controls are determined by the KI Data Repository and the conditions associated with the dataset record.

**Access conditions:**  
[To be completed following KI Data Repository assessment]

Users of the dataset are responsible for complying with the applicable access conditions, ethical approvals, data-protection requirements, and terms of reuse.

## Ethics

The study was conducted under the applicable ethical approval(s).

**Ethical approval:**  
[Insert approving authority and approval/reference number(s)]

Further details are provided in the associated publication.

## Citation

If you use this dataset or analysis code, please cite:

**Article:**  
[Full article citation]

**Dataset:**  
[Dataset citation supplied by the KI Data Repository]

**Code:**  
[GitHub citation/DOI if applicable]

## License and reuse

**Dataset reuse terms:**  
[To be specified by the KI Data Repository]

**Code license:**  
[Insert selected code license]

The licensing of the analysis code does not override restrictions or conditions applying to the clinical dataset.

## Contact

For questions regarding the study or dataset:

[Name]  
[Department / research group]  
Karolinska Institutet  
[Institutional email]