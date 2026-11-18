# Dataset: Early Stage Diabetes Risk Prediction (UCI)
# https://archive.ics.uci.edu/dataset/529/early+stage+diabetes+risk+prediction+dataset

rm(list = ls())

library(dplyr)
library(car)        # Boot, Confint
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(sandwich)   # vcovHC
library(lmtest)     # coeftest
library(caret)      # confusionMatrixl

# ------------------------- data loading and preparation ----------------------------
setwd(" --- ")
df <- read.csv("diabetes_data_upload.csv")

# Categorical variables -> binary factors (1 = Yes / Positive / Male, 0 otherwise)
df[] <- lapply(df, function(x) {
  if (is.character(x)) {
    factor(ifelse(x %in% c("Yes", "Positive", "Male"), 1, 0))
  } else {
    x
  }
})
df$Age <- as.numeric(df$Age)

# The response variable is kept numeric (0/1)
df$class <- as.numeric(as.character(df$class))
str(df)

# ------------- correlation ----------------
# Descriptive statistics of the response class (mean = 0.615 -> sufficiently balanced class)
mean(df$class)
sd(df$class)

# Correlation matrix (binary variables are converted to 0/1)
df_corr <- df
binary_vars <- names(df_corr)[sapply(df_corr, is.factor)]
df_corr[binary_vars] <- lapply(df_corr[binary_vars], function(x) as.numeric(as.character(x)))

cor_matrix <- cor(df_corr, use = "complete.obs")

corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  addCoef.col = "black"
)

# Variables most correlated with the response: contingency tables
table(df$Polyuria, df$class)
table(df$Polydipsia, df$class)

# Boxplot of Age by class
boxplot(df$Age ~ df$class,
        main = "Age distribution by class",
        xlab = "Class",
        ylab = "Age",
        col = c("lightblue", "salmon"),
        border = "gray30",
        las = 1)

# ------------- possible interactions of Age with the other covariates ----------------
# Interesting situations: lines with opposite direction and not flat.
# The most interesting ones are Itching, Irritability and partial.paresis.
vars_int <- c("Gender", "Polyuria", "Polydipsia", "Irritability",
              "Itching", "Genital.thrush", "partial.paresis", "Obesity", "Alopecia")

plots <- lapply(vars_int, function(v) {
  ggplot(df, aes(x = Age, y = class, color = factor(.data[[v]]))) +
    geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
    labs(title = paste("Age:", v),
         x = "Age",
         y = "Prevalence of class = 1",
         color = v) +
    theme_minimal()
})

grid.arrange(grobs = plots, ncol = 3, nrow = 3)

# Boxplot of Age by class, conditioned on the three most interesting covariates
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (v in c("Itching", "Irritability", "partial.paresis")) {
  for (val in 0:1) {
    sel <- df[[v]] == val
    boxplot(df$Age[sel] ~ df$class[sel],
            main = paste(v, "=", val),
            xlab = "class", ylab = "Age",
            col = c("lightblue", "salmon"))
  }
}
par(mfrow = c(1, 1))

# ------------- model selection via backward elimination ----------------
# Backward selection with BIC criterion (k = log(n)) starting from the full model
scope_full <- class ~ Age * Itching + Gender + Polyuria + Polydipsia + sudden.weight.loss +
  weakness + Polyphagia + Genital.thrush + visual.blurring +
  Irritability + delayed.healing + partial.paresis +
  muscle.stiffness + Alopecia + Obesity

m_full <- glm(scope_full, data = df, family = binomial(link = "logit"))
m_step <- step(m_full, direction = "backward", k = log(nrow(df)))
summary(m_step)
# Resulting model:
# class ~ Age * Itching + Gender + Polyuria + Polydipsia +
#         Polyphagia + Genital.thrush + Irritability

# ------------- dichotomization of the Age variable ----------------
# Three age classes defined on the tertiles
df$Age_quantile <- cut(df$Age,
                       breaks = quantile(df$Age, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                       include.lowest = TRUE,
                       labels = c("Young", "Middle_aged", "Older"))
table(df$Age_quantile, df$class)

# Binary variable defined on the median = 47 (50 is used for easier interpretation)
df$Age_binary <- ifelse(df$Age < 50, "Under_50", "Over_50")
table(df$Age_binary, df$class)

# Comparison of the three Age specifications (continuous, 3 quantiles, binary)
m0 <- glm(class ~ Age * Itching + Gender + Polyuria + Polydipsia + Polyphagia +
            Genital.thrush + Irritability, data = df, family = binomial)
m1 <- glm(class ~ Age_quantile * Itching + Gender + Polyuria + Polydipsia + Polyphagia +
            Genital.thrush + Irritability, data = df, family = binomial)
m2 <- glm(class ~ Age_binary * Itching + Gender + Polyuria + Polydipsia + Polyphagia +
            Genital.thrush + Irritability, data = df, family = binomial)

models <- list(Age_continuous = m0, Age_3_quantiles = m1, Age_binary = m2)

accuracy_glm <- function(m, soglia = 0.5) {
  pred <- factor(as.integer(predict(m, type = "response") >= soglia), levels = c(0, 1))
  confusionMatrix(pred, factor(df$class, levels = c(0, 1)), positive = "1")$overall["Accuracy"]
}

results <- data.frame(
  Model    = names(models),
  Accuracy = sapply(models, accuracy_glm),
  BIC      = sapply(models, BIC)
)
results

# The model with binary Age (m2) is the best one (highest accuracy and lowest BIC)
summary(m2)
formula(m2)

# ------------- diagnostics on grouped data ----------------
# Since all variables are binary, the data are grouped by covariate combination
# (coefficient estimates and standard errors do not change).
vars <- setdiff(all.vars(formula(m2)), "class")

df_grouped <- df %>%
  group_by(across(all_of(vars))) %>%
  summarise(
    grouped_class = sum(class == 1),   # number of successes
    n = n(),                           # group size
    .groups = "drop"
  )
head(df_grouped[1:8])

m_grouped <- glm(
  cbind(grouped_class, n - grouped_class) ~ Age_binary * Itching + Gender + Polyuria +
    Polydipsia + Polyphagia + Genital.thrush + Irritability,
  data = df_grouped,
  family = binomial(link = "logit")
)
summary(m_grouped)

# Residuals vs fitted: residuals equal to 0 in the right half, the model predicts a
# group of observations very well -> possible heteroscedasticity problem.
par(mfrow = c(2, 2))
plot(m_grouped)
par(mfrow = c(1, 1))

# plot residuals weighted by number of grouped observations 
# Pearson residuals
res <- residuals(m_grouped, type = "pearson")

# Group sizes
n <- df_grouped$n

pred <- predict(m_grouped)

# Plot
plot(pred, res,
     pch = 19,
     cex = sqrt(n) / 2,
     xlab = "Predicted values",
     ylab = "Pearson residuals",
     main = "Pearson residuals weighted by group size")

abline(h = 0, lty = 2)
abline(h = c(-2, 2), col = "red", lty = 2)
plot(m_grouped, which = 1)

df_grouped[which(abs(res) > 2),]

# ------------- heteroscedasticity and overdispersion ----------------
# Theoretical variance of the model: Var(Yi) = mi * pi * (1 - pi), so it depends on the
# estimated probability: as pi approaches 0 or 1 the variance decreases.

# Estimate of the dispersion parameter (values > 1 indicate overdispersion,
# with potentially overly optimistic standard errors)
(dispersion <- sum(residuals(m_grouped, type = "pearson")^2) / df.residual(m_grouped))

# Check with quasibinomial model: unchanged estimates, standard errors corrected
# for the dispersion parameter
m_quasi <- glm(formula(m_grouped), data = df_grouped,
               family = quasibinomial(link = "logit"))
summary(m_quasi)

# Robust var-cov matrix (HC0): glm does not accept vcov as an argument, so coeftest is used.
# The p-values of Polyphagia and Genital.thrush become less significant.
coeftest(m_grouped, vcov. = vcovHC(m_grouped, type = "HC3"))

# ------------- comparison of rich model and poor model ----------------
# The models on non-grouped data are used (same estimates).
scope_rich <- class ~ Age_binary * Itching + Gender + Polyuria +
  Polydipsia + Polyphagia + Genital.thrush + Irritability
scope_poor <- class ~ Age_binary * Itching + Gender + Polyuria +
  Polydipsia + Irritability

m_rich <- glm(scope_rich, data = df, family = binomial)
m_poor <- glm(scope_poor, data = df, family = binomial)

summary(m_rich)
summary(m_poor)

# Nested models: likelihood ratio test (LRT) and BIC
anova(m_poor, m_rich, test = "LRT")
BIC(m_rich)
BIC(m_poor)

# ------------- performance metrics on training data ----------------
pred_prob_rich <- predict(m_rich, type = "response")
pred_prob_poor <- predict(m_poor, type = "response")

soglia <- 0.5
pred_rich <- ifelse(pred_prob_rich >= soglia, 1, 0)
pred_poor <- ifelse(pred_prob_poor >= soglia, 1, 0)

conf_matrix_rich <- confusionMatrix(data = as.factor(pred_rich), reference = as.factor(df$class))
conf_matrix_poor <- confusionMatrix(data = as.factor(pred_poor), reference = as.factor(df$class))

conf_matrix_rich$overall["Accuracy"]
conf_matrix_poor$overall["Accuracy"]

conf_matrix_rich$table
conf_matrix_poor$table

# Matthews Correlation Coefficient (MCC): takes TP, TN, FP and FN into account
mcc_from_counts <- function(TP, TN, FP, FN) {
  den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (den == 0) return(NA_real_)
  (TP * TN - FP * FN) / den
}

calc_MCC <- function(cm) {
  mcc_from_counts(TP = as.numeric(cm$table[2, 2]),
                  TN = as.numeric(cm$table[1, 1]),
                  FP = as.numeric(cm$table[2, 1]),
                  FN = as.numeric(cm$table[1, 2]))
}

# the rich model achieves a higher MCC
calc_MCC(conf_matrix_rich)
calc_MCC(conf_matrix_poor)   

# ------------- leverage points and separation problems ----------------
# Analysis on aggregated data: Cook > 0.5 would suggest removing an entire group
# (hundreds of observations); altering the dataset this much is not correct.
idx_cook_grouped <- which(cooks.distance(m_grouped) > 0.5)
df_grouped[idx_cook_grouped, 1:7]

# We therefore study the leverage points on the non-aggregated data (rich model)
plot(m_rich, which = 4)

df_model <- model.frame(m_rich)
id_model <- as.integer(rownames(df_model))

p_hat <- fitted(m_rich)
rP    <- residuals(m_rich, type = "pearson")
h     <- hatvalues(m_rich)

# Approximate change in Pearson's chi-square statistic
DIFCHISQ <- rP^2 / (1 - h)
cooks    <- cooks.distance(m_rich)

diag_df <- data.frame(
  id       = id_model,
  class    = df_model$class,
  fitted   = p_hat,
  pearson  = rP,
  leverage = h,
  DIFCHISQ = DIFCHISQ,
  cooks    = cooks
)

# Relationship between estimated probabilities and DIFCHISQ
ggplot(diag_df, aes(x = fitted, y = DIFCHISQ)) +
  geom_point(aes(size = cooks), shape = 1) +
  geom_text_repel(
    data = subset(diag_df, DIFCHISQ > 4 | cooks > quantile(cooks, 0.95)),
    aes(label = id),
    color = "red",
    size = 3
  ) +
  geom_hline(yintercept = 4, linetype = "dotted") +
  scale_size_area(max_size = 15) +
  labs(x = "Estimated probability",
       y = "Change in Pearson Chi-Square",
       size = "Cook's distance") +
  theme_classic()

# Observations with the highest DIFCHISQ (9 observations with DIFCHISQ > 4)
diag_df <- diag_df[order(-diag_df$DIFCHISQ), ]
head(diag_df, 15)
sum(diag_df$DIFCHISQ > 4)

# ------------- possible quasi-separation ----------------
# With Polydipsia = 1 and Polyuria = 1 the observations are all of class 1,
# but there is no complete separation: with both covariates at 0 there remain
# observations of class 1.
xtabs(~ class + Gender + Polydipsia + Polyuria, data = df)
# Polydipsia = 1, Polyuria = 1:
#         Gender
# class     0   1
#     0     0   0
#     1   113  80

# Mean values of fitted and DIFCHISQ by class, Polydipsia and Polyuria:
# the highest DIFCHISQ values are observed where Polyuria and Polydipsia are 0
diag <- data.frame(class = df_model$class,
                   Polydipsia = df_model$Polydipsia,
                   Polyuria = df_model$Polyuria,
                   fitted = p_hat,
                   DIFCHISQ = DIFCHISQ)

aggregate(cbind(fitted, DIFCHISQ) ~ class + Polydipsia + Polyuria,
          data = diag, FUN = mean)

diag <- diag[order(-diag$DIFCHISQ), ]
head(diag, 15)

# All the high DIFCHISQ values come from Polydipsia = 0 and Polyuria = 0
diag_PP00 <- diag[diag$Polydipsia == 0 & diag$Polyuria == 0, ]
head(diag_PP00, 15)

# More accurate analysis: 6 of the 9 observations with high DIFCHISQ have Gender = 1
diag_gender <- data.frame(class = df_model$class,
                          fitted = p_hat,
                          Polydipsia = df_model$Polydipsia,
                          Polyuria = df_model$Polyuria,
                          Gender = df_model$Gender,
                          DIFCHISQ = DIFCHISQ)

diag_PPG001 <- diag_gender[diag_gender$Polydipsia == 0 &
                             diag_gender$Polyuria == 0 &
                             diag_gender$Gender == 1, ]
head(diag_PPG001[order(-diag_PPG001$DIFCHISQ), ], 15)

# Some class 1 observations receive very low probabilities because
# Polydipsia = 1 and Polyuria = 1 push the estimated probability to 1.

# ------------- coefficient variation ----------------
# Model 1: all observations
# Model 2: without observations with DIFCHISQ > 4  (9 observations out of 520)
# Model 3: without observations with DIFCHISQ > 10 (5 observations out of 520)
# Warning: removing the observations poorly predicted by Polydipsia and Polyuria
# inflates the coefficients to fit the training data (risk of overfitting).
ids_4  <- diag_df$id[diag_df$DIFCHISQ > 4]
ids_20 <- diag_df$id[diag_df$DIFCHISQ > 20]

model_no_infl4  <- glm(scope_rich, data = df[-ids_4, ],  family = binomial)
model_no_infl20 <- glm(scope_rich, data = df[-ids_20, ], family = binomial)

coef_comparison <- cbind(`Model 1` = coef(m_rich),
                        `Model 3 DIFCHISQ>20` = coef(model_no_infl20),
                         `Model 2 DIFCHISQ>4`  = coef(model_no_infl4)
                         )
coef_comparison

# ------------------------ bootstrap ----------------------------------------
# Bootstrap (R = 1999) of the rich and poor models
set.seed(123)

# rich model
BOOT.MOD.rich <- Boot(m_rich, R = 1999)
summary(BOOT.MOD.rich, high.moments = TRUE)              # variability of the coefficients
Confint(BOOT.MOD.rich, level = 0.95, type = "perc")      # confidence intervals
hist(BOOT.MOD.rich, legend = "separate")                 # distribution of the coefficients

# poor model
BOOT.MOD.poor <- Boot(m_poor, R = 1999)
summary(BOOT.MOD.poor, high.moments = TRUE)
Confint(BOOT.MOD.poor, level = 0.95, type = "perc")
hist(BOOT.MOD.poor, legend = "separate")

# ---------------------------- optimal threshold -------------------------------------------
# The classification threshold is studied on the poor model
model_pred <- pred_prob_poor
y_true <- df$class
thresholds <- seq(0.01, 0.99, by = 0.01)

# Confusion matrix counts for a given threshold
conf_counts <- function(threshold, prob, y_true) {
  pred <- as.integer(prob >= threshold)
  cm <- table(factor(pred, levels = c(0, 1)),
              factor(y_true, levels = c(0, 1)))
  list(TN = as.numeric(cm[1, 1]), FN = as.numeric(cm[1, 2]),
       FP = as.numeric(cm[2, 1]), TP = as.numeric(cm[2, 2]))
}

# Optimization via MCC
MCC_threshold <- sapply(thresholds, function(t) {
  with(conf_counts(t, model_pred, y_true), mcc_from_counts(TP, TN, FP, FN))
})

results_mcc <- data.frame(threshold = thresholds, MCC = MCC_threshold)

best_threshold_mcc <- results_mcc$threshold[which.max(results_mcc$MCC)]
best_threshold_mcc
max(results_mcc$MCC, na.rm = TRUE)

plot(results_mcc$threshold, results_mcc$MCC,
     type = "l",
     xlab = "Classification threshold",
     ylab = "MCC",
     main = "MCC optimization")
abline(v = best_threshold_mcc, col = "red", lty = 2)

# Optimization via cost function
# In a medical setting a false negative (undetected diabetic patient) is worse
# than a false positive: cost_FN = 1.5, cost_FP = 1.
cost_function <- function(threshold, prob, y_true, cost_FN = 1.5, cost_FP = 1) {
  cc <- conf_counts(threshold, prob, y_true)
  cost_FN * cc$FN + cost_FP * cc$FP
}

costs <- sapply(thresholds, cost_function, prob = model_pred, y_true = y_true,
                cost_FN = 1.5, cost_FP = 1)

best_threshold_cost <- thresholds[which.min(costs)]
best_threshold_cost
min(costs)

plot(thresholds, costs,
     type = "l",
     xlab = "Threshold",
     ylab = "Classification cost",
     main = "Cost-sensitive threshold optimization")
abline(v = best_threshold_cost, col = "red", lty = 2)

