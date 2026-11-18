# Early-Stage Diabetes Risk: Logistic Regression Analysis

Logistic regression analysis of the UCI *Early Stage Diabetes Risk Prediction* dataset
(520 patients, 16 predictors, binary response `class`: 1 = positive, 0 = negative; mean 0.615, fairly balanced).
Dataset: https://archive.ics.uci.edu/dataset/529/early+stage+diabetes+risk+prediction+dataset

## Contents

| File | Description |
|------|-------------|
| `Logistic_regression_report_EN.pdf` | 6-page report with results, tables and figures |
| `diabetes_analysis_en.R` | R script that produces all the analyses in the report (comments and labels in English) |

## How to run the script

1. Download `diabetes_data_upload.csv` from the UCI page above.
2. Edit the `setwd(...)` line at the top of the script so that it points to the folder containing the CSV.
3. Install the packages: `dplyr`, `car`, `ggplot2`, `gridExtra`, `ggrepel`, `sandwich`, `lmtest`, `caret`, **`corrplot`**
   (`corrplot` is used but not loaded in the script: add `library(corrplot)` before the correlation plot).
4. Run the script top to bottom. The bootstrap (R = 1999, `set.seed(123)`) is the slowest step.

## Report structure and matching script sections

| Report page | Script section | Content |
|-------------|----------------|---------|
| 1 | data loading, correlation | Variable description, class balance, correlation matrix, contingency tables for Polyuria / Polydipsia, Age boxplot |
| 2 | interactions, backward selection, Age dichotomization | Age x covariate plots, BIC backward selection, comparison of continuous / 3-quantile / binary Age, model summary |
| 3 | grouped diagnostics, overdispersion, rich vs poor | Grouped-data model, dispersion (1.406), HC3 robust standard errors, LRT and BIC, confusion matrices, accuracy and MCC on training data |
| 4 | leverage and separation | DIFCHISQ and Cook's distance, quasi-separation (Polyuria = 0 and Polydipsia = 0), coefficient variation after removing influential points |
| 5 | bootstrap, optimal threshold | Bootstrap (R = 1999) of the rich and poor models, MCC-based and cost-based threshold optimization |
| 6 | (report only) | Poor-model summary, coefficient interpretation, comment on robustness |

## Models

- **Selected by BIC backward selection** (`step(..., k = log(n))`; `k` is the per-parameter penalty, so `k = log(n)` gives the BIC):
  `class ~ Age * Itching + Gender + Polyuria + Polydipsia + Polyphagia + Genital.thrush + Irritability`
- **Age** is dichotomized as `Age_binary` (`Under_50` / `Over_50`), the specification with the best BIC and accuracy.
- **Rich model**: the formula above with `Age_binary`.
- **Poor model**: the rich model without `Polyphagia` and `Genital.thrush`.

## Main findings

- Rich model: accuracy 0.9461, MCC 0.886 on training data; poor model: 0.9346 and 0.861. BIC is practically identical (230.98 vs 230.90).
- Overdispersion (1.406) makes model-based standard errors optimistic. With HC3 robust standard errors, `Polyphagia` and `Genital.thrush` are no longer significant in the rich model.
- Bootstrap: the poor model has narrower intervals and low kurtosis; the rich model has heavy-tailed distributions (sign of instability).
- Influential observations all come from the subgroup with Polyuria = 0 and Polydipsia = 0. Removing the 9 with DIFCHISQ > 4 makes the model fail to converge (quasi-separation).
- The poor model is the more robust one, but all results are on training data: cross-validation or a test on unseen data is still needed.

## Known inconsistencies in the script (not changed)

- `library(corrplot)` is missing (see above).
- The comment on the robust covariance says "HC0", but the code uses `type = "HC3"`.
- The Age threshold is 50, while the comment mentions the median (47) and the report snippet shows `median(df$Age)`.
- Comment says "Model 3: DIFCHISQ > 10", but the code uses `DIFCHISQ > 20`. In `coef_comparison` the column order is Model 1, "Model 3 DIFCHISQ>20", "Model 2 DIFCHISQ>4", whereas in the report Model 2 is DIFCHISQ > 20 and Model 3 is DIFCHISQ > 4.
- The overdispersion check, HC3 and influence diagnostics were run on the rich model only; the threshold analysis was run on the poor model.
