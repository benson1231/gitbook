# Paired Comparisons

Paired comparisons are used when observations are collected in matched pairs or from the same subjects at two time points. These tests account for the dependency between measurements.

---

## 1. When to Use Paired Comparisons
- Before-and-after measurements (e.g., pre- and post-treatment)
- Matched subjects in a case-control study
- Repeated measures on the same subject

---

## 2. Parametric Test: Paired t-test
- Compares the mean difference between two related groups
- Assumes differences are normally distributed

```r
# Paired t-test in R
paired_result <- t.test(before, after, paired = TRUE)
summary(paired_result)
```

### Assumptions
- Data are continuous
- Differences follow a normal distribution
- Observations are dependent within pairs but independent between pairs

---

## 3. Non-parametric Alternative: Wilcoxon Signed-Rank Test
- Used when the assumption of normality is not met

```r
# Wilcoxon test in R
wilcox.test(before, after, paired = TRUE)
```

---

## 4. Visualization
- **Line plots** to show change per subject
- **Boxplots** of paired differences

```r
# Line plot example
matplot(t(cbind(before, after)), type = "l", lty = 1, col = 1:n)
```

---

## 5. Summary
Paired comparisons are essential for within-subject or matched designs. Use the paired t-test when assumptions are met and Wilcoxon test otherwise. Visualization helps reveal the direction and consistency of changes.