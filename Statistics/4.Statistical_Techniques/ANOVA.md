# Analysis of Variance (ANOVA)

Analysis of Variance (ANOVA) is a statistical technique used to compare the means of three or more groups. It tests whether there are any statistically significant differences between the means of independent (unrelated) groups.

---

## 1. When to Use ANOVA
- Comparing mean blood pressure across three treatment groups
- Testing whether multiple teaching methods affect exam scores differently
- Analyzing differences in plant growth under different fertilizer types

---

## 2. Types of ANOVA

### a. One-way ANOVA
Used when comparing means across groups based on one independent variable (factor).

- **Example**: Comparing average test scores across four different schools

### b. Two-way ANOVA
Used when examining the effect of two independent variables and their interaction.

- **Example**: Testing the effects of diet (A/B) and exercise (yes/no) on weight loss

---

## 3. ANOVA Terminology
- **Between-group variance**: Variability due to differences between group means
- **Within-group variance**: Variability due to differences within each group
- **F-statistic**: Ratio of between-group variance to within-group variance

---

## 4. Assumptions of ANOVA
- Independence of observations
- Normally distributed residuals within each group
- Homogeneity of variances (equal variances across groups)

---

## 5. Interpreting ANOVA
- **Null hypothesis (H0)**: All group means are equal
- **Alternative hypothesis (Ha)**: At least one group mean is different
- If the p-value from the F-test < 0.05, reject H0

---

## 6. Post-hoc Tests
If ANOVA is significant, use post-hoc tests to identify which group means differ:
- Tukey’s HSD
- Bonferroni correction
- Scheffé’s method

---

## 7. Example (R Code)
```r
# One-way ANOVA
result <- aov(score ~ group, data = mydata)
summary(result)

# Post-hoc test
TukeyHSD(result)
```

---

## 8. Summary
ANOVA is a powerful method for testing differences among group means. It relies on specific assumptions and, when significant, should be followed by post-hoc comparisons to identify which groups differ.