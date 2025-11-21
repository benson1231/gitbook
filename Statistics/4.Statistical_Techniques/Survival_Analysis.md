# Survival Analysis

Survival analysis focuses on the time until an event of interest occurs, such as death, relapse, or device failure. It accounts for **censored data**, where the event has not occurred for some subjects during the study period.

---

## 1. Key Concepts
- **Survival Time**: Time from a defined starting point to the occurrence of the event
- **Censoring**: When the exact survival time is unknown (e.g., lost to follow-up)
- **Hazard Function**: The instantaneous rate of the event occurring at time t

---

## 2. Kaplan-Meier Estimator
Non-parametric method to estimate survival probabilities over time.

```r
# Kaplan-Meier curve in R
library(survival)
km_fit <- survfit(Surv(time, status) ~ group, data = data)
plot(km_fit, xlab = "Time", ylab = "Survival Probability")
```

- `time`: follow-up time
- `status`: 1 = event occurred, 0 = censored

---

## 3. Log-rank Test
Used to compare survival curves between groups.

```r
survdiff(Surv(time, status) ~ group, data = data)
```

---

## 4. Cox Proportional Hazards Model
Semi-parametric model that estimates hazard ratios while adjusting for covariates.

### Model
\[ h(t|X) = h_0(t) \cdot \exp(\beta_1 X_1 + \beta_2 X_2 + \cdots + \beta_k X_k) \]

```r
cox_model <- coxph(Surv(time, status) ~ age + sex + treatment, data = data)
summary(cox_model)
```

### Assumptions
- Proportional hazards: the ratio of hazards is constant over time
- Linearity in log-hazard

---

## 5. Visualization
- Kaplan-Meier plots for survival function
- Forest plots for hazard ratios

---

## 6. Summary
Survival analysis is essential for time-to-event data. Use Kaplan-Meier for visualization, log-rank test for group comparison, and Cox models for multivariable adjustment.
