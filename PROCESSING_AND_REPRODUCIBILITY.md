# Processing and reproducibility

## Purpose

This document describes the data-processing and statistical-analysis workflow for the study:

**Clinical predictors of genetic diagnoses in adults evaluated for inborn errors of immunity**

It supplements the repository `README.md` and data dictionary by documenting:

- the relationship between the secure source data and the deposited analysis dataset;
- preprocessing performed before creation of the deposited dataset;
- historical predictor screening and variable selection;
- the analyses that can be reproduced from the deposited dataset;
- handling of missing data and special statistical cases; and
- aspects of the original research workflow that cannot be independently reconstructed from the deposited dataset.

The public reproducibility workflow begins with:

`data/gimpid_clean.csv`

The deposited dataset contains 124 observations and 81 variables: 79 eligible clinical and laboratory predictors, the molecular diagnosis outcome, and the IUIS classification variable.

No participant identifier is included in the deposited dataset.

---

## 1. Study cohort and source data

This was a retrospective single-center study of patients evaluated at the Adult IEI Outpatient Clinic, Department of Infectious Diseases, Karolinska University Hospital.

Patients were eligible if genetic testing had been performed as part of the evaluation of a clinically diagnosed inborn error of immunity (IEI) before 2025.

A total of 153 patients were considered during the study period. After exclusion of patients who declined participation or could not be contacted, 124 patients were included in the final study cohort.

Clinical symptoms, laboratory results, microbiological findings, and radiological data were extracted retrospectively from electronic medical records. Investigations were performed as part of routine clinical care.

The deposited dataset represents the final analytic cohort. Source-data extraction, cohort assembly, and handling of identifiable clinical information occurred within the secure research environment and are not reproduced in the public repository.

---

## 2. Primary outcome

The primary study outcome was achievement of a molecular diagnosis.

In the deposited dataset, the source variable is:

`Relevant genetics`

A molecular diagnosis was defined as identification of a pathogenic or likely pathogenic variant considered sufficient to explain the patient's clinical phenotype.

Variants of uncertain significance and incidental findings were not considered diagnostic.

The public analysis scripts recode this variable into standardized analysis variables where required.

The final cohort contains:

- 55 participants with a molecular diagnosis; and
- 69 participants without a molecular diagnosis.

---

## 3. Secure preprocessing and creation of the deposited dataset

The deposited analysis dataset was generated from secure study source data using a preprocessing workflow that is not included in the public GitHub repository.

The private preprocessing workflow performed operations required before release of the analysis dataset, including:

- removal of direct and internal identifiers;
- harmonization of data types;
- conversion of laboratory variables stored as character values into numeric values where appropriate;
- harmonization of decimal notation and source-specific laboratory notation;
- handling of selected non-numeric laboratory entries;
- removal of variables not required for the reproducible analyses;
- incorporation of IUIS classification; and
- creation of the final analysis-ready release dataset.

The private preprocessing code is not required to reproduce the statistical analyses from the deposited dataset and is not distributed because it operates on secure source data that are not part of the research data package.

The resulting deposited dataset is:

`data/gimpid_clean.csv`

It contains:

- 124 rows;
- 81 columns;
- 79 eligible predictor variables;
- 1 molecular diagnosis outcome variable; and
- 1 IUIS classification variable.

Historical source-variable names have generally been retained in the deposited dataset to preserve correspondence with the original analyses. More descriptive definitions are provided in the accompanying data dictionary.

---

## 4. Historical predictor screening

### 4.1 Initial candidate variables

The original predictor-screening workflow began with 109 candidate predictor variables.

Before domain-based predictor prioritization, variables were screened according to predefined data-availability and analytic criteria.

The historical screening proceeded as follows:

| Screening stage | Number of variables |
| --- | ---: |
| Initial candidate predictors | 109 |
| Complement-related variables excluded | 6 |
| Variables excluded because of >20% missingness | 22 |
| Variables excluded because of no variability | 2 |
| Eligible predictors | 79 |

The two variables excluded because of no variability were:

- `HIV`
- `Autoimmune liver disease`

This resulted in 79 eligible predictors.

### 4.2 Reproducibility boundary for the 109-to-79 screening step

The deposited dataset contains the **79 eligible predictors**, not all 109 variables considered during the historical screening stage.

Consequently, the public analysis workflow documents the historical reduction from 109 to 79 predictors but cannot independently recompute all exclusions from the deposited dataset.

In particular, variables excluded because of complement status, excessive missingness, or absence of variability are generally not retained in the release dataset.

The public variable-selection script therefore treats the 109-to-79 screening as a documented historical preprocessing step.

As an integrity check, the public workflow verifies that all 79 expected eligible predictors are present in the deposited dataset and evaluates missingness among those released predictors.

The subsequent predictor-prioritization workflow beginning with the 79 eligible predictors can be reproduced from the deposited data.

---

## 5. Predictor domains

The 79 eligible predictors were organized into clinically defined domains used in the predictor-selection workflow.

These domains were:

- general characteristics;
- CVID-related characteristics;
- infectious manifestations;
- malignancy;
- allergy and atopy;
- autoimmunity;
- inflammation;
- general laboratory measurements;
- immunoglobulins; and
- lymphocyte subsets.

The domain definitions implemented in `03_variable_selection.R` reflect the historical analysis structure.

Five eligible predictors were not assigned to one of these domain-specific LASSO groups:

- age at genetic testing;
- consanguinity;
- number of radiologically documented pneumonias;
- lowest CRP; and
- lowest ESR.

These variables remain part of the set of 79 eligible predictors and are retained in the univariable analyses.

---

## 6. Domain-based predictor prioritization

Predictor prioritization was performed within clinical domains using LASSO logistic regression.

The purpose of this stage was variable prioritization rather than construction of the final inferential model solely through an automated selection procedure.

The public implementation reconstructs this domain-based prioritization using the 79 eligible predictors in the deposited dataset.

The historical prioritization identified 22 candidate predictors for further consideration:

- age at symptom onset;
- CVID;
- bronchiectasis;
- nodular regenerative hyperplasia / porto-sinusoidal vascular disease;
- IgG1;
- CD16/56;
- B-cell count;
- sex;
- immunoglobulin replacement therapy;
- splenomegaly;
- GLILD;
- eczema;
- antibiotic allergy;
- food/animal/pollen allergy;
- celiac disease;
- ALT;
- ALP;
- hemoglobin;
- IgG4;
- IgG3;
- CD8; and
- CD3.

The reconstructed public LASSO procedure and the historical selection record are both retained in the analysis script so that differences between historical and reconstructed selection can be identified rather than silently overwritten.

After domain-based prioritization, candidate predictors were further evaluated for collinearity, redundancy, interpretability, and clinical relevance before construction of the multivariable model.

The variable-selection workflow should therefore be interpreted as a combination of statistical prioritization and clinical review rather than a fully automated model-selection algorithm.

---

## 7. Primary multivariable model

The primary analysis used multivariable logistic regression with molecular diagnosis as the binary outcome.

The final model included:

- age at symptom onset;
- sex;
- serum IgG1 concentration; and
- eczema.

Age at symptom onset and IgG1 were modeled as continuous variables.

Sex was represented as male versus female, and eczema as present versus absent.

The model was fitted using complete observations for all variables included in the model.

The complete-case sample consisted of:

- 112 participants in total;
- 48 with a molecular diagnosis; and
- 64 without a molecular diagnosis.

Odds ratios and 95% confidence intervals were calculated from the fitted logistic regression model.

Profile-likelihood confidence intervals were used in the public analysis workflow.

Collinearity diagnostics are also calculated in the analysis script.

The corresponding analysis is implemented in:

`04_multivariable_analysis.R`

The forest plot displaying the principal adjusted associations is generated by:

`09_figure4_forest_plot.R`

Age at symptom onset is included in the adjusted model but is not displayed in the manuscript forest plot.

---

## 8. Age-at-onset analysis

Age at symptom onset was evaluated separately because of its strong association with molecular diagnostic outcome.

The continuous association was assessed using univariable logistic regression.

For presentation, age at symptom onset was additionally grouped into the following categories:

- 0–5 years;
- 6–20 years;
- 21–30 years; and
- >30 years.

These categories are used for descriptive visualization and do not replace the continuous age-at-onset variable in the primary multivariable model.

The analysis is implemented in:

`05_age_at_onset_analysis.R`

---

## 9. Sensitivity analyses by PAD status

Sensitivity analyses were performed according to predominantly antibody deficiency (PAD) status.

The same four predictors used in the primary multivariable model were evaluated:

- age at symptom onset;
- sex;
- serum IgG1 concentration; and
- eczema.

Separate complete-case logistic regression models were fitted for participants classified as PAD and non-PAD.

The sensitivity analyses are implemented in:

`06_sensitivity_analyses.R`

IUIS classification in the deposited dataset is used to derive the PAD/non-PAD grouping.

---

## 10. Univariable predictor analyses

Univariable associations between each of the 79 eligible predictors and molecular diagnosis are evaluated in:

`10_univariable_analyses.R`

Standard logistic regression is used where ordinary maximum-likelihood estimation is appropriate.

The historical analysis used Firth penalized logistic regression for two sparse predictors:

- celiac disease; and
- nodular regenerative hyperplasia / porto-sinusoidal vascular disease.

The public analysis reproduces this handling explicitly rather than selecting Firth regression automatically based on the observed results.

For the gastric cancer variable, the outcome association is not estimable because of the extremely sparse data configuration. It is therefore reported as not estimable rather than forcing a finite logistic-regression estimate.

Continuous predictors are analyzed on their original measurement scales. Accordingly, odds ratios for continuous predictors represent the change in odds associated with a one-unit increase in the corresponding variable unless otherwise stated.

Complete-case data are used separately for each univariable model, so the analysis sample size can vary between predictors.

---

## 11. Missing data

Missing data are represented as missing values when imported into R.

The historical predictor-screening criterion excluded candidate predictors with more than 20% missing observations.

All 79 predictors retained in the deposited dataset satisfy this eligibility criterion.

No statistical imputation is performed in the public analysis workflow.

Regression models use complete observations for the variables included in each respective model.

Consequently, sample size differs between analyses depending on the predictors included.

For construction of descriptive clinical manifestation domains used in the heatmap and related descriptive summaries, missing values in the relevant manifestation indicators are treated as absence of that manifestation. This convention is limited to the descriptive domain construction and should not be interpreted as a general missing-data rule for the regression analyses.

---

## 12. Laboratory variables

Laboratory measurements were obtained as part of routine clinical care.

Serum IgG, IgA, IgM, and IgG subclasses were measured using the Optilite turbidimetric platform. Immunoglobulin measurements were obtained before initiation of immunoglobulin replacement therapy whenever applicable.

The deposited dataset retains the original analysis measurement scales.

The data dictionary documents the units and interpretation of the released laboratory variables.

The historical source variable:

`IgG1-subclass def`

contains continuous serum IgG1 concentration values and is treated as a continuous variable in the analyses. The historical column name has been retained for reproducibility and should not be interpreted as a binary indicator of IgG1 subclass deficiency.

Peripheral blood lymphocyte subsets were assessed using multiparameter flow cytometry according to established clinical laboratory protocols.

---

## 13. Clinical manifestation domains and heatmap

Clinical manifestations were grouped into seven descriptive domains for visualization:

- bacterial/fungal infections;
- viral infections;
- atopy;
- lymphoproliferation;
- autoimmunity;
- malignancy; and
- autoinflammation.

These descriptive domains are distinct from the predictor domains used during variable selection.

The clinical-manifestation heatmap is generated in:

`02_descriptive_analysis.R`

Participants are clustered according to their clinical manifestation profiles using hierarchical clustering with Ward's method.

Molecular diagnosis status is included as an annotation.

The heatmap is descriptive and does not determine inclusion of variables in the final regression model.

---

## 14. IUIS classification and diagnostic yield

IUIS classification is included in the deposited dataset as:

`iuis_category`

Diagnostic yield according to IUIS classification is calculated in:

`07_iuis_diagnostic_yield.R`

This classification is also used to define PAD status for the sensitivity analyses.

Individual gene and disease labels are not required for these analyses and are not included in the deposited analysis dataset.

---

## 15. Supplementary distributions

Supplementary descriptive figures showing the distributions of age at symptom onset and serum IgG1 are generated by:

`08_supplementary_distributions.R`

Distributions are presented for the overall cohort and according to PAD/non-PAD status.

Missing IgG1 measurements are retained as missing and are not imputed.

---

## 16. Public analysis workflow

The public scripts should be run sequentially from a clean R session:

1. `01_prepare_analysis_data.R`
2. `02_descriptive_analysis.R`
3. `03_variable_selection.R`
4. `04_multivariable_analysis.R`
5. `05_age_at_onset_analysis.R`
6. `06_sensitivity_analyses.R`
7. `07_iuis_diagnostic_yield.R`
8. `08_supplementary_distributions.R`
9. `09_figure4_forest_plot.R`
10. `10_univariable_analyses.R`

The scripts have been tested sequentially in this order using the deposited analysis dataset.

The first script reads:

`data/gimpid_clean.csv`

and creates the derived analysis object used by subsequent scripts.

Intermediate analysis objects are stored under:

`data/derived/`

Generated figures and analysis outputs are written to the output locations defined in the corresponding scripts.

---

## 17. Reproducibility scope

The repository is intended to reproduce the analyses that can be performed from the deposited analysis dataset.

The following components are reproducible from the deposited data:

- descriptive cohort analyses;
- construction of the reported clinical manifestation domains;
- the clinical-manifestation heatmap;
- predictor prioritization beginning with the 79 eligible predictors;
- the primary multivariable logistic regression;
- age-at-onset analyses;
- PAD/non-PAD sensitivity analyses;
- diagnostic yield by IUIS category;
- supplementary age-at-onset and IgG1 distributions;
- the adjusted forest plot; and
- univariable analyses of the 79 eligible predictors.

The following components precede the deposited dataset and are therefore documented rather than independently reproducible from the public research package:

- extraction of information from electronic medical records;
- source-data cleaning and harmonization;
- cohort assembly within the secure research environment;
- removal of direct and internal identifiers;
- the complete historical screening from 109 candidate predictors to 79 eligible predictors, because excluded variables are not included in the deposited dataset;
- linkage of secure study identifiers to IUIS classification before removal of those identifiers; and
- analyses requiring individual gene, variant, or disease labels that are not necessary for reproduction of the predictor analyses.

This distinction is intentional. The deposited dataset contains the information required to audit and reproduce the reported statistical predictor analyses without distributing source variables or identifiers that are unnecessary for that purpose.

---

## 18. Confidentiality and data access

The deposited dataset contains detailed clinical research information from a relatively small cohort of individuals evaluated for IEI.

Direct and internal identifiers have been removed from the deposited dataset. However, combinations of clinical characteristics may remain sensitive.

The dataset should therefore not be described as anonymous solely on the basis of identifier removal.

Access conditions, disclosure controls, and permitted reuse are determined by the KI Data Repository and the conditions associated with the final dataset record.

The public analysis code does not require participant identifiers.

---

## 19. Data dictionary

The accompanying data dictionary provides one entry for each of the 81 variables in the deposited dataset.

For each variable, it documents where applicable:

- the original deposited variable name;
- a human-readable description;
- data type;
- measurement unit;
- coding;
- missing-value representation;
- observed missingness; and
- explanatory notes.

Historical source-variable names have been retained in the dataset where changing them could make correspondence with the original analysis less transparent.

---

## 20. Software environment

All analyses were performed in R.

The complete public analysis workflow (`01`–`10`) was tested sequentially from a clean R session using R version 4.5.0 (2025-04-11) on Windows 11 x64.

The principal R packages attached during the analysis workflow were `dplyr` (1.1.4), `readr` (2.1.5), `tidyr` (1.3.1), `tibble` (3.3.0), `Matrix` (1.7-3), `glmnet` (5.0), `ggplot2` (4.0.3), and `logistf` (1.26.1). Additional package dependencies were loaded through these packages as required.

A complete record of the R session and package versions after successful execution of the workflow is provided in `sessionInfo.txt`.

No manual editing of statistical results is required between scripts in the public reproducibility workflow.

---

## 21. Relationship between manuscript, dataset, and code

The research package consists of three related components:

**Manuscript**  
Describes the study design, clinical interpretation, statistical methods, and reported results.

**Deposited dataset**  
Contains the analysis-ready variables required for reproduction of the public statistical analyses, subject to the access conditions determined by the KI Data Repository.

**Public analysis repository**  
Contains the code required to transform the deposited dataset into analysis objects and reproduce the statistical analyses, tables, and figures supported by the released data.

The secure source dataset and private preprocessing workflow remain outside the public repository because they contain information and processing steps that are not necessary for reproduction of the released analyses.