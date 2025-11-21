# T-tests and Z-tests

T-tests and Z-tests are inferential statistical tests used to compare means. They assess whether observed differences between groups are statistically significant.

---

## 1. When to Use
- Comparing mean blood pressure between treatment and control groups
- Evaluating changes in performance before and after intervention
- Determining if a sample mean differs from a known population mean

---

## 2. T-tests
Used when the population standard deviation is unknown and sample size is small.

### a. One-sample t-test
Compares the sample mean to a known or hypothesized population mean.

```r
# One-sample t-test
 t.test(sample_data, mu = 100)
```

### b. Independent two-sample t-test
Compares means between two independent groups.

```r
# Two-sample t-test
 t.test(group1, group2)
```

### c. Paired t-test
Used for paired or matched samples.

```r
# Paired t-test
 t.test(before, after, paired = TRUE)
```

### Assumptions
- Data are continuous and approximately normally distributed
- Homogeneity of variances (for independent samples)
- Observations are independent

---

## 3. Z-tests
Used when population variance is known or sample size is large (n > 30).

### a. One-sample Z-test
Compares a sample mean to a known population mean with known population standard deviation.

### b. Two-sample Z-test
Compares means from two independent large samples with known variances.

Z-tests are less commonly used in practice because population standard deviation is rarely known.

---

## 4. Interpreting Results
- **Null hypothesis (H₀)**: The means are equal
- **Alternative hypothesis (H₁)**: The means are different
- If p-value < 0.05, reject H₀ and conclude a significant difference exists

---

## 5. Summary
T-tests and Z-tests are fundamental tools in hypothesis testing for means. T-tests are preferred when standard deviation is unknown. Choose the test based on sample design and variance knowledge.