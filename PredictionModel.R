# Install / load packages -------------------------------------------------

library(survival)   # survival objects, Cox model
library(survminer)  # KM and Cox plots
library(dplyr)      # data wrangling
library(ggplot2)    # plotting
library(Hmisc)      # C-index (rcorr.cens)
library(tidyr)      # drop_na
library(knitr)      # kable
library(broom)      # tidy outputs
library(here)       # file paths

set.seed(713)       # reproducibility


#--- Load and Prepare Data -----------------------------------------------

colon <- read.csv(here("colon.csv"))

# Replace nodes = 0 with 1
colon <- colon %>%
  mutate(
    nodes = ifelse(nodes == 0, 1, nodes)
  )


# Define survival object for recurrence-free survival
surv_obj <- Surv(time = colon$time, event = colon$status)


#--- Variable Coding ------------------------------------------------------

# Recode treatment (rx) and sex to factors
colon <- colon %>%
  mutate(
    rx = factor(rx,
                levels = c("Obs", "Lev", "Lev+5FU")),
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



#--- Descriptive Statistics ----------------------------------------------

# Baseline characteristics by treatment group (on full colon, not colon_model)
baseline_table <- colon %>%
  group_by(rx) %>%
  summarise(
    n = n(),
    age_mean = mean(age, na.rm = TRUE),
    age_sd   = sd(age, na.rm = TRUE),
    nodes_median = median(nodes, na.rm = TRUE),
    nodes_IQR    = IQR(nodes, na.rm = TRUE),
    male_pct     = mean(sex == "Male", na.rm = TRUE) * 100,
    obstruct_pct = mean(obstruct == "Obstruction", na.rm = TRUE) * 100,
    perfor_pct   = mean(perfor == "Perforation", na.rm = TRUE) * 100,
    nonadherent_pct = mean(adhere == "Not adherent", na.rm = TRUE) * 100,
    .groups = "drop"
  )

print(baseline_table)


#--- Kaplan–Meier Curves by Treatment ------------------------------------

fit_km_rx <- survfit(Surv(time, status) ~ rx, data = colon)

km_plot_rx <- ggsurvplot(
  fit_km_rx,
  data       = colon,
  risk.table = TRUE,
  conf.int   = TRUE,
  pval       = FALSE,  # turned off to avoid "only 1 group" survdiff errors
  legend.title = "Treatment",
  legend.labs  = levels(colon$rx),
  xlab       = "Time since surgery (days)",
  ylab       = "Recurrence-free survival probability",
  ggtheme    = theme_minimal()
)

print(km_plot_rx)



#--- Cox Prediction Model -------------------------------------------------

cox <- coxph(
  Surv(time, status) ~ rx + sex + age +
    obstruct + perfor + adhere + nodes + extent + surg,
  data = colon
)

summary(cox)

cox_tidy <- tidy(cox, exponentiate = TRUE, conf.int = TRUE)
kable(cox_tidy, digits = 3,
      caption = "Cox prediction model: hazard ratios (for reference)")


#--- Model Diagnostics ----------------------------------------------------

## 1) Proportional Hazards Assumption

zph_full <- cox.zph(cox)
print(zph_full)

# Schoenfeld residual plots
ggcoxzph(zph_full)


## 2) Martingale residuals

martingale_resid <- residuals(cox, type = "martingale")

diagnostic_df <- colon %>%
  drop_na(time, status, rx, sex, age, obstruct, perfor, adhere, nodes, extent, surg) %>%
  mutate(martingale = martingale_resid)

# Age vs martingale
ggplot(diagnostic_df, aes(x = age, y = martingale)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  labs(
    title = "Martingale residuals vs. age",
    x     = "Age (years)",
    y     = "Martingale residuals"
  )

# Nodes vs martingale
ggplot(diagnostic_df, aes(x = nodes, y = martingale)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  labs(
    title = "Martingale residuals vs. nodes",
    x     = "Number of positive nodes",
    y     = "Martingale residuals"
  )

# Surgery duration vs martingale
ggplot(diagnostic_df, aes(x = surg, y = martingale)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  labs(
    title = "Martingale residuals vs. surgery duration",
    x     = "Surgery duration (minutes)",
    y     = "Martingale residuals"
  )


## 3) Influential observations: Deviance residuals

deviance_resid <- residuals(cox, type = "deviance")
lp <- predict(cox, type = "lp")

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


#--- Predicted Survival Curves by Treatment -------------------------------

# Build a "typical patient" profile using means / most common values
newdata_template <- colon_model %>%
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

newdata_pred <- newdata_template %>%
  slice(rep(1, length(levels(colon_model$rx)))) %>%
  mutate(rx = levels(colon_model$rx))

fit_pred <- survfit(cox, newdata = newdata_pred)

pred_plot_rx <- ggsurvplot(
  fit_pred,
  data         = newdata_pred,
  legend.title = "Treatment",
  legend.labs  = levels(colon_model$rx),
  xlab         = "Time since surgery (days)",
  ylab         = "Predicted recurrence-free survival",
  ggtheme      = theme_minimal()
)

print(pred_plot_rx)
ggsave("figure_predicted_survival_by_treatment.png",
       pred_plot_rx$plot, width = 7, height = 5, dpi = 300)


#--- Internal Validation: 10-fold Cross-Validation ------------------------

K <- 10
n <- nrow(colon_model)

set.seed(713)
fold_id <- sample(rep(1:K, length.out = n))
colon_model$fold_id <- fold_id

# store cross-validated linear predictor
lp_cv <- rep(NA_real_, n)

for (k in 1:K) {
  train_data <- colon_model %>% filter(fold_id != k)
  test_data  <- colon_model %>% filter(fold_id == k)
  
  cox_k <- coxph(
    Surv(time, status) ~ rx + sex + age +
      obstruct + perfor + adhere + nodes + extent + surg,
    data = train_data
  )
  
  lp_k <- predict(cox_k, newdata = test_data, type = "lp")
  lp_cv[colon_model$fold_id == k] <- lp_k
}

colon_model$lp_cv <- lp_cv

## C-index (discrimination, out-of-sample)
cindex_cv <- rcorr.cens(
  x = colon_model$lp_cv,
  S = Surv(colon_model$time, colon_model$status)
)

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


## Calibration slope
cox_calib <- coxph(
  Surv(time, status) ~ lp_cv,
  data = colon_model
)

calib_tidy <- tidy(cox_calib, conf.int = TRUE, conf.level = 0.95) %>%
  filter(term == "lp_cv") %>%
  select(estimate, conf.low, conf.high, p.value)

kable(
  calib_tidy,
  digits = 3,
  caption = "Calibration slope for cross-validated risk score"
)


#--- Calibration Plot at Fixed Time (e.g., 3 years) ----------------------

t0 <- 1095  # ~3 years

colon_model <- colon_model %>%
  mutate(lp_decile = ntile(lp_cv, 10))

get_km_at_t0 <- function(data, t0) {
  fit <- survfit(Surv(time, status) ~ 1, data = data)
  s <- summary(fit, times = t0)$surv
  if (length(s) == 0) return(NA_real_)
  s
}

calib_df <- colon_model %>%
  group_by(lp_decile) %>%
  summarise(
    mean_lp = mean(lp_cv, na.rm = TRUE),
    surv_t0 = get_km_at_t0(cur_data_all(), t0),
    .groups = "drop"
  ) %>%
  arrange(lp_decile)

ggplot(calib_df, aes(x = mean_lp, y = surv_t0)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(
    title = paste0("Calibration plot at t = ", t0, " days"),
    x     = "Mean cross-validated linear predictor (higher = higher risk) ",
    y     = "Observed recurrence-free survival at t0"
  )

ggsave("figure_calibration_plot_t0.png",
       width = 7, height = 5, dpi = 300)

#---------------------------- END OF SCRIPT ------------------------------
