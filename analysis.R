library(ggplot2)
library(tidyverse)
library(GGally)
library(Hmisc)
library(car)
library(MuMIn)

library(VGAM)
require(nnet)
options(na.action = "na.fail")

# Data import
admissions <- read.csv("data/Admissions.csv")

# Data exploration

# Dependent variables
# age x admission: insignificant
anova_age <- aov(age ~ admission, data = admissions)
summary(anova_age)

# los x admission: significant
anova_los <- aov(los ~ admission, data = admissions)
summary(anova_los)

# age80 x admission: insignificant
chisq_age80 <- chisq.test(table(admissions$age80, admissions$admission))
chisq_age80

# white x admission: significant
chisq_white <- chisq.test(table(admissions$white, admissions$admission))
chisq_white

# died x admission: significant
chisq_died <- chisq.test(table(admissions$died, admissions$admission))
chisq_died

# provnum x admission: significant
chisq_hospital <- chisq.test(table(admissions$provnum, admissions$admission))
chisq_hospital

# Observe between-hospital trends
# Group by hospital
hospital_admissions <- admissions %>%
  group_by(provnum) %>%
  summarise(
    n = n(),
    prop_died = mean(as.numeric(as.character(died))), # Proportion of died
    prop_white = mean(as.numeric(as.character(white))), # Proportion of white
    prop_age80 = mean(as.numeric(as.character(age80))), # Proportion of age80
    mean_los = mean(los), # Mean length of stay
    # Proportions of each age level
    prop_age1 = mean(age == "1"),
    prop_age2 = mean(age == "2"),
    prop_age3 = mean(age == "3"),
    prop_age4 = mean(age == "4"),
    prop_age5 = mean(age == "5"),
    prop_age6 = mean(age == "6"),
    prop_age7 = mean(age == "7"),
    prop_age8 = mean(age == "8"),
    prop_age9 = mean(age == "9"),
    # Proportions of each admission level
    prop_ele = mean(admission == "Elective"),
    prop_urg = mean(admission == "Urgent"),
    prop_eme = mean(admission == "Emergency")
  )

summary(hospital_admissions)

# Issues:
# sample numbers vary a lot between hospitals
# Levels of mortality vary between hospitals
# Age varies a lot between hospitals
# Median hospital has no age1/age2 (i.e no pediatrics)
# Median hospital has no age9
# Median hospital has no emergency department

# If hospital has emergency ward
hospital_admissions <- hospital_admissions %>%
  mutate(has_emergency = ifelse(prop_eme > 0, 1, 0))

# If hospital has urgent ward
hospital_admissions <- hospital_admissions %>%
  mutate(has_urgent = ifelse(prop_urg > 0, 1, 0))

# If hospital has pediatric ward
hospital_admissions <- hospital_admissions %>%
  mutate(has_pediatric = ifelse(prop_age1 > 0, 1, 0))

# If hospital has youth ward
hospital_admissions <- hospital_admissions %>%
  mutate(has_youth = ifelse(prop_age2 > 0, 1, 0))

# If hospital has elderly ward
hospital_admissions <- hospital_admissions %>%
  mutate(has_elderly = ifelse(prop_age9 > 0, 1, 0))

# If hospital is specialist
hospital_admissions <- hospital_admissions %>%
  mutate(is_specialist = ifelse(prop_ele == 1 | prop_urg == 1 | prop_eme == 1, 1, 0))

admissions <- admissions %>%
  left_join(hospital_admissions %>% select(provnum, has_emergency, has_urgent, has_pediatric, has_youth, has_elderly, is_specialist), by = "provnum")
# Ignore hospital in main data
admissions$provnum <- NULL

# Remove specialist hospitals from analysis
admissions <- admissions %>%
  filter(is_specialist != 1)
admissions$is_specialist <- NULL

cor_matrix <- cor(admissions %>% select_if(is.numeric))
print(cor_matrix)

# Convert the correlation matrix into a data frame of pairwise correlations
cor_list <- as.data.frame(as.table(cor_matrix))
# Remove self-correlations (where Variable 1 = Variable 2)
cor_list <- cor_list[cor_list$Var1 != cor_list$Var2, ]
cor_list <- cor_list %>%
  mutate(abs_correlation = abs(Freq)) %>%
  arrange(desc(abs_correlation))

# age x age80: choose which is better predictor of age
# provnum x los: investigate variance between hospitals
# age x white: more older people are white
# age x died: more old people die
# age80 x died
# los x died: less time in hospital, less chance of death

# Data wrangling
admissions$provnum <- as.factor(admissions$provnum)
admissions$died <- as.factor(admissions$died)
admissions$white <- as.factor(admissions$white)
admissions$age <- as.factor(admissions$age)
admissions$age80 <- as.factor(admissions$age80)
admissions$admission <- as.factor(admissions$admission)
admissions$has_emergency <- as.factor(admissions$has_emergency)
admissions$has_urgent <- as.factor(admissions$has_urgent)
admissions$has_pediatric <- as.factor(admissions$has_pediatric)
admissions$has_youth <- as.factor(admissions$has_youth)
admissions$has_elderly <- as.factor(admissions$has_elderly)

# Model selection
# First model
init_model <- admission ~ age + age80 + white + los + died + has_emergency + has_urgent + has_pediatric + has_youth + has_elderly
model <- vglm(init_model, family=multinomial, data = admissions)

# Add interactions
model.mn <- multinom(init_model, dat = admissions)
model.nm.interac <- update(model.mn, .~. + age80 * age + has_youth * has_emergency + has_youth * has_urgent + has_youth * has_pediatric + has_youth * has_elderly + has_pediatric * has_emergency + has_urgent * has_emergency)
summary(model.nm.interac)

# Remove insignificant covariates

# Via ANOVA
Anova(model.nm.interac)
anova_model <- admission ~ los + died + has_emergency + has_urgent + has_elderly + has_youth + has_emergency:has_youth + has_youth:has_elderly
model.anova <- multinom(anova_model, dat = admissions)
AICc(model.anova)

# Via dredge
# Results of dredge in dredge.csv
dredge_model <- admission ~ los + died + white + has_elderly + has_emergency + has_urgent + has_pediatric + has_youth + has_emergency:has_youth + has_elderly:has_youth
model.dredge <- multinom(dredge_model, dat = admissions)
AICc(model.dredge)

# Although white and has_pediatric not selected in ANOVA, close to threshold and including them raises AICc
# Include these in final model
model.anova.dredge <- update(model.anova, .~. + white + has_pediatric)
AICc(model.anova.dredge)

# Final model combining ANOVA and dredge
model.final.mn <- model.anova.dredge
formula.final <- formula(model.final)
model.final.vglm <- vglm(formula.final, data = admissions, family = multinomial)


# Assumptions checking
# Assuming `formula.final` is the formula extracted from `model.final`
# Extract the right-hand side of the formula (covariates)
covariates <- as.character(formula.final)[3]

# Create formulas for the binary models
formula_elective_vs_rest <- as.formula(paste("binary_elective_vs_rest ~", covariates))
formula_emergency_vs_rest <- as.formula(paste("binary_emergency_vs_rest ~", covariates))
formula_urgent_vs_rest <- as.formula(paste("binary_urgent_vs_rest ~", covariates))

# Fit binomial logistic regression models using these formulas
model_elective_vs_rest <- glm(formula_elective_vs_rest, family = binomial, data = admissions)
model_emergency_vs_rest <- glm(formula_emergency_vs_rest, family = binomial, data = admissions)
model_urgent_vs_rest <- glm(formula_urgent_vs_rest, family = binomial, data = admissions)

# Summarize models
summary(model_elective_vs_rest)
summary(model_emergency_vs_rest)
summary(model_urgent_vs_rest)

AICc(model.final.mn, model_elective_vs_rest, model_emergency_vs_rest, model_urgent_vs_rest)

# Assess residuals
par(mfrow = c(3, 1))  # Set plotting area for three residual plots

plot_residuals <- function(model, title) {
  fitted_values <- fitted(model)
  residuals_deviance <- residuals(model, type = "deviance")

  plot(
    fitted_values, residuals_deviance,
    xlab = "Fitted Values",
    ylab = "Deviance Residuals",
    main = title,
    pch = 20, col = "blue"
  )
  abline(h = 0, lty = 2, col = "red")
}

# For each model
plot_residuals(model_elective_vs_rest, "Residuals vs Fitted (Elective vs Rest)")
plot_residuals(model_emergency_vs_rest, "Residuals vs Fitted (Emergency vs Rest)")
plot_residuals(model_urgent_vs_rest, "Residuals vs Fitted (Urgent vs Rest)")

coef_multinom <- coef(model.final)
print(coef_multinom)

coef_elective_vs_rest <- coef(model_elective_vs_rest)
coef_emergency_vs_rest <- coef(model_emergency_vs_rest)
coef_urgent_vs_rest <- coef(model_urgent_vs_rest)

binomial_coefficients <- data.frame(
  Elective_vs_Rest = coef_elective_vs_rest,
  Emergency_vs_Rest = coef_emergency_vs_rest,
  Urgent_vs_Rest = coef_urgent_vs_rest
)
print(binomial_coefficients)

# Create the comparison table
comparison_table <- data.frame(
  Multinomial_Elective = rep(0, ncol(coef_multinom)),  # Elective is baseline, so coefficients are all 0
  Multinomial_Emergency = coef_multinom["Emergency", ],
  Multinomial_Urgent = coef_multinom["Urgent", ],
  Binomial_Elective_vs_Rest = binomial_coefficients$Elective_vs_Rest,
  Binomial_Emergency_vs_Rest = binomial_coefficients$Emergency_vs_Rest,
  Binomial_Urgent_vs_Rest = binomial_coefficients$Urgent_vs_Rest
)

# Add row names for covariates
rownames(comparison_table) <- colnames(coef_multinom)


# Ensure covariates are aligned
covariates <- colnames(coef_multinom)

# Create a data frame with differences
differences <- data.frame(
  Emergency_Difference = coef_multinom["Emergency", covariates] - binomial_coefficients$Emergency_vs_Rest,
  Urgent_Difference = coef_multinom["Urgent", covariates] - binomial_coefficients$Urgent_vs_Rest,
  Elective_Difference = 0 - binomial_coefficients$Elective_vs_Rest  # Multinomial baseline (0) for Elective
)

# Print the results
print(differences)



# Observations come from multinomial distributions
# Residual plots are parallel lines in all three cases, but best for emergency vs rest

# AICc for all binomial models outperform multinomial model: by far, best is emergency vs rest
# model.final.mn          22 1584.4235
# model_elective_vs_rest  11 1396.0902
# model_emergency_vs_rest 11  365.8892
# model_urgent_vs_rest    11 1235.8663

# Coefficient differences in binomial vs multinomial are minimal with demographic and hospital stay info, but lots of coefficient differences in hospital structure
#                          Emergency_Difference Urgent_Difference Elective_Difference
# (Intercept)                        1.555994911      -0.525044435         -0.57200868
# los                                0.009673969       0.001039893          0.03975408
# died1                              0.078621827       0.019627569          0.42034458
# has_emergency1                    -1.475908153       0.043273349         -1.18270522
# has_urgent1                        0.129361595       0.506774631          0.43362102
# has_elderly1                      -0.098049243       0.009571130         -0.64731587
# has_youth1                         0.258583369       0.048766229         -1.20302056
# white1                            -0.121840658      -0.003091603         -0.47039122
# has_pediatric1                     0.123422552       0.040652210          0.52321164
# has_emergency1:has_youth1         -0.222762934       0.074960221          1.62914321
# has_elderly1:has_youth1            0.066094317      -0.082772564          0.26897299


# Independence of observations
# Do in report


# Linearity
# Extract the explanatory variables from the model's formula
explanatory_vars <- all.vars(formula(model.final.vglm)[-2])  # Get variables on the right-hand side of the formula

# Predicted probabilities from the multinomial model
predicted_probs <- fitted(model.final.vglm)

# Extract probabilities for each outcome
p_urgent <- predicted_probs[, "Urgent"]
p_emergency <- predicted_probs[, "Emergency"]
p_elective <- predicted_probs[, "Elective"]

# Calculate log-odds
log_odds_urgent_vs_elective <- log(p_urgent / p_elective)
log_odds_emergency_vs_elective <- log(p_emergency / p_elective)

# Set up a layout for multiple plots
num_vars <- length(explanatory_vars)
# par(mfrow = c(ceiling(num_vars / 2), 2))  # Adjust rows and columns dynamically
par(mfrow = c(1, 4))

# Loop through each explanatory variable and plot
for (var in explanatory_vars) {
  # Plot Urgent vs Elective
  plot(admissions[[var]], log_odds_urgent_vs_elective,
       main = paste("Log-Odds: Urgent vs Elective vs", var),
       xlab = var, ylab = "Log-Odds: Urgent vs Elective", pch = 19)

  # Plot Emergency vs Elective
  plot(admissions[[var]], log_odds_emergency_vs_elective,
       main = paste("Log-Odds: Emergency vs Elective vs", var),
       xlab = var, ylab = "Log-Odds: Emergency vs Elective", pch = 19)
}

# Length of stay: quarter of points clustered at -20 log odds than rest in urgent vs elective, 50:50 in emergency vs elective. Relationships themselves look linear.
# Clustering caused by systematic sepertion between two groups not explained by los, but linearity holds

# 


# IIA
# Do in report


# Graphical assessment
# Extract residuals and fitted values
mn_residuals <- residuals(model.final.vglm, type = "response")
mn_fitted <- fitted(model.final.vglm)

# Calculate observed proportions
mn_observed <- mn_fitted + mn_residuals

# Plot observed vs fitted proportions for each response category
response_categories <- colnames(mn_fitted)  # Get response levels: "Urgent", "Emergency", "Elective"
par(mfrow = c(1, 3))  # Set up plotting area for 3 plots

for (i in 1:3) {
  # Plot observed vs fitted for each category
  plot(mn_fitted[, i], mn_observed[, i],
       main = response_categories[i],
       xlab = "Fitted Probabilities",
       ylab = "Observed Proportions",
       pch = 19, col = "blue")

  # Add reference line (y = x)
  abline(a = 0, b = 1, col = "black", lty = 1)

  # Add a smoothed spline for trend
  lines(smooth.spline(mn_fitted[, i], mn_observed[, i], df = 4), col = "red", lty = 2)
}


# Numerical assessment
saturated_model <- vglm(admission ~ 1, family = multinomial(refLevel = 3), data = admissions)

# Chi^2
# Fit the saturated model (perfect fit)

# Deviance of the fitted model
fitted_deviance <- deviance(model.final.vglm)

# Deviance of the saturated model
saturated_deviance <- deviance(saturated_model)

# Calculate residual deviance
residual_deviance <- fitted_deviance - saturated_deviance

# Number of observations and parameters
n <- nrow(admissions)  # Number of observations
k <- length(unique(admissions$admission))  # Number of response levels
p <- length(coef(model.noage))  # Number of parameters in the fitted model

# Degrees of freedom
df <- (k - 1) * n - p

# Chi-squared test for goodness-of-fit
p_value <- 1 - pchisq(residual_deviance, df)

# Output results
cat("Residual Deviance:", residual_deviance, "\n")
cat("Degrees of Freedom:", df, "\n")
cat("P-value:", p_value, "\n")

residual_deviance <- deviance(model.noage)  # Extract residual deviance
df_residual <- df.residual(model.noage)     # Extract residual degrees of freedom
su
logLikNull <- logLik(saturated_model)
mcfadden_r2 <- function(model) {
  return(1-logLik(model)/logLikNull)
}
# Low R^2: 0.106
mcfadden_r2(model.final.vglm)
