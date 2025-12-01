# Install packages 

library(survival)   # survival objects, Cox model
library(survminer)  
library(dplyr)      
library(ggplot2)    
library(Hmisc)      
library(tidyr)      
library(knitr)      
library(broom)      
library(here)


#--- Load and Prepare Data ---
colon <- read.csv(here("colon.csv"))

set.seed(713)       # reproducibility

# Replace nodes = 0 with 1
colon %>%
  mutate(
    nodes = ifelse(nodes == 0, 1, nodes)
  )

# Remove missing data
colon %>%
  drop_na()

# Define survival object for recurrence-free survival
# status: 1 = recurrence, 0 = censored
surv_obj <- Surv(time = colon$time, event = colon$status)



#--- Variable Coding ---

# Recode treatment (rx) and sex to factors
colon %>%
  mutate(
    rx = factor(rx,
                levels = c(1, 2, 3),
                labels = c("Obs", "Lev", "Lev+5FU")),
    sex = factor(sex,
                 levels = c(0, 1),
                 labels = c("Female", "Male")),
    obstruct = factor(obstruct,
                      levels = c(0, 1),
                      labels = c("No obstruction", "Obstruction")),
    perfor  = factor(perfor,
                     levels = c(0, 1),
                     labels = c("No perforation", "Perforation")),
    adhere  = factor(adhere,
                     levels = c(0, 1),
                     labels = c("Adherent", "Not adherent"))
  )




#--- Descriptive Statistics ---

# Baseline characteristics by treatment group
baseline_table <- colon %>%
  group_by(rx) %>%
  summarise(
    n = n(),
    age_mean = mean(age, na.rm = TRUE),
    age_sd   = sd(age, na.rm = TRUE),
    nodes_median = median(nodes, na.rm = TRUE),
    nodes_IQR    = IQR(nodes, na.rm = TRUE),
    male_pct     = mean(sex == "Male") * 100,
    obstruct_pct = mean(obstruct == "Obstruction") * 100,
    perfor_pct   = mean(perfor == "Perforation") * 100,
    nonadherent_pct = mean(adhere == "Not adherent") * 100,
    .groups = "drop"
  )

print(baseline_table)



#--- Kaplan–Meier Curves by Treatment ---

fit_km_rx <- survfit(Surv(time, status) ~ rx, data = colon)

km_plot_rx <- ggsurvplot(
  fit_km_rx,
  data       = colon_rf,
  risk.table = TRUE,
  conf.int   = TRUE,
  legend.title = "Treatment",
  legend.labs  = levels(colon_rf$rx),
  xlab       = "Time since surgery (days)",
  ylab       = "Recurrence-free survival probability",
  ggtheme    = theme_minimal()
)

print(km_plot_rx)




########## 5. Cox Prediction Model ##########

# Full Cox model including treatment and potential predictors
cox_full <- coxph(
  Surv(surv_time, surv_event) ~ rx + sex + age +
    obstruct + perfor + adhere + nodes + extent + surg,
  data = colon_rf
)

# Summary (for your own reference)
summary(cox_full)

# Tidy hazard ratios (optional, you can refer to this but
# main focus is prediction rather than inference)
cox_full_tidy <- tidy(cox_full, exponentiate = TRUE, conf.int = TRUE)

kable(
  cox_full_tidy,
  digits = 3,
  caption = "Cox prediction model: hazard ratios (for reference)"
)


########## 6. Model Diagnostics ##########

## 6.1 Proportional Hazards Assumption

zph_full <- cox.zph(cox_full)

# Print PH tests
kable(
  tidy(zph_full),
  digits = 3,
  caption = "Test of proportional hazards assumption"
)

# Plot Schoenfeld residuals for visual check
ggcoxzph(zph_full)  # will open a multi-panel plot


## 6.2 Functional Form for Continuous Covariates (Martingale residuals)

martingale_resid <- residuals(cox_full, type = "martingale")

diagnostic_df <- colon_rf %>%
  mutate(martingale = martingale_resid)

# Age vs Martingale residuals
ggplot(diagnostic_df, aes(x = age, y = martingale)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  labs(
    title = "Martingale residuals vs. age",
    x     = "Age (years)",
    y     = "Martingale residuals"
  )

# Nodes vs Martingale residuals
ggplot(diagnostic_df, aes(x = nodes, y = martingale)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  labs(
    title = "Martingale residuals vs. nodes",
    x     = "Number of positive nodes",
    y     = "Martingale residuals"
  )

# Surgery duration vs Martingale residuals
ggplot(diagnostic_df, aes(x = surg, y = martingale)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  labs(
    title = "Martingale residuals vs. surgery duration",
    x     = "Surgery duration (minutes)",
    y     = "Martingale residuals"
  )


## 6.3 Influential Observations (Deviance residuals vs LP)

deviance_resid <- residuals(cox_full, type = "deviance")
lp <- predict(cox_full, type = "lp")

diagnostic_df <- diagnostic_df %>%
  mutate(
    deviance = deviance_resid,
    lp       = lp
  )

ggplot(diagnostic_df, aes(x = lp, y = deviance)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Deviance residuals vs. linear predictor",
    x     = "Linear predictor",
    y     = "Deviance residuals"
  )


########## 7. Predicted Survival Curves by Treatment ##########

# Create "typical patient" profile: mean/most common covariate values
newdata_template <- colon_rf %>%
  summarise(
    sex      = names(sort(table(sex), decreasing = TRUE))[1],
    age      = mean(age, na.rm = TRUE),
    obstruct = names(sort(table(obstruct), decreasing = TRUE))[1],
    perfor   = names(sort(table(perfor), decreasing = TRUE))[1],
    adhere   = names(sort(table(adhere), decreasing = TRUE))[1],
    nodes    = median(nodes, na.rm = TRUE),
    extent   = names(sort(table(extent), decreasing = TRUE))[1],
    surg     = mean(surg, na.rm = TRUE)
  )

# Create one row per treatment level
newdata_pred <- newdata_template %>%
  slice(rep(1, length(levels(colon_rf$rx)))) %>%
  mutate(rx = levels(colon_rf$rx))

# Fit predicted survival curves from Cox model
fit_pred <- survfit(cox_full, newdata = newdata_pred)

pred_plot_rx <- ggsurvplot(
  fit_pred,
  data         = newdata_pred,
  legend.title = "Treatment",
  legend.labs  = levels(colon_rf$rx),
  xlab         = "Time since surgery (days)",
  ylab         = "Predicted recurrence-free survival",
  ggtheme      = theme_minimal()
)

print(pred_plot_rx)

ggsave(
  filename = "figure_predicted_survival_by_treatment.png",
  plot     = pred_plot_rx$plot,
  width    = 7, height = 5, dpi = 300
)


########## 8. Internal Validation: K-fold CV ##########

# 10-fold cross-validation for:
#   - C-index (discrimination)
#   - Calibration slope

K <- 10
n <- nrow(colon_rf)

set.seed(713)
fold_id <- sample(rep(1:K, length.out = n))
colon_rf$fold_id <- fold_id

# Vector to store cross-validated linear predictor
lp_cv <- rep(NA_real_, n)

for (k in 1:K) {
  # Training and test split
  train_data <- colon_rf %>% filter(fold_id != k)
  test_data  <- colon_rf %>% filter(fold_id == k)
  
  # Fit Cox model in training data
  cox_k <- coxph(
    Surv(surv_time, surv_event) ~ rx + sex + age +
      obstruct + perfor + adhere + nodes + extent + surg,
    data = train_data
  )
  
  # Predict linear predictor (risk score) in test data
  lp_k <- predict(cox_k, newdata = test_data, type = "lp")
  
  # Store CV predictions
  lp_cv[colon_rf$fold_id == k] <- lp_k
}

colon_rf$lp_cv <- lp_cv

## 8.1 Cross-validated C-index

cindex_cv <- rcorr.cens(x = colon_rf$lp_cv,
                        S = Surv(colon_rf$surv_time, colon_rf$surv_event))

cindex_est <- as.numeric(cindex_cv["C Index"])

perf_table <- tibble(
  metric   = "C-index (10-fold CV)",
  estimate = cindex_est
)

kable(
  perf_table,
  digits = 3,
  caption = "Out-of-sample discrimination: C-index from 10-fold cross-validation"
)


## 8.2 Calibration slope from CV risk score

cox_calib <- coxph(
  Surv(surv_time, surv_event) ~ lp_cv,
  data = colon_rf
)

calib_tidy <- tidy(cox_calib, conf.int = TRUE, conf.level = 0.95) %>%
  filter(term == "lp_cv") %>%
  select(estimate, conf.low, conf.high, p.value)

kable(
  calib_tidy,
  digits = 3,
  caption = "Calibration slope for cross-validated risk score"
)


########## 9. Calibration Plot at Fixed Time (e.g., 3 years) ##########

# Choose time point (e.g., ~3 years = 1095 days)
t0 <- 1095

# Divide subjects into deciles of predicted risk (lp_cv)
colon_rf <- colon_rf %>%
  mutate(lp_decile = ntile(lp_cv, 10))

# Helper: KM survival at t0 in a subset
get_km_at_t0 <- function(data, t0) {
  fit <- survfit(Surv(surv_time, surv_event) ~ 1, data = data)
  s <- summary(fit, times = t0)$surv
  if (length(s) == 0) return(NA_real_)
  s
}

calib_df <- colon_rf %>%
  group_by(lp_decile) %>%
  summarise(
    mean_lp = mean(lp_cv, na.rm = TRUE),
    surv_t0 = get_km_at_t0(cur_data_all(), t0),
    .groups = "drop"
  ) %>%
  arrange(lp_decile)

# Calibration plot
ggplot(calib_df, aes(x = mean_lp, y = surv_t0)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(
    title = paste0("Calibration plot at t = ", t0, " days"),
    x     = "Mean cross-validated linear predictor (higher = higher risk)",
    y     = "Observed recurrence-free survival at t0"
  )

ggsave(
  filename = "figure_calibration_plot_t0.png",
  width    = 7, height = 5, dpi = 300
)

###############################################
# End of script
###############################################
