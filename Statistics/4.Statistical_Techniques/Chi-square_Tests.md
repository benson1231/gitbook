# Chi-square Tests

Chi-square (χ²) tests are used to determine whether there is a significant association between categorical variables. They compare the observed frequencies in each category to the frequencies expected under the null hypothesis.

---

## 1. When to Use Chi-square Tests
- Testing independence between gender and preference for a product
- Checking if a die is fair
- Analyzing survey results in contingency tables

---

## 2. Types of Chi-square Tests

### a. Test of Independence
Used to assess whether two categorical variables are associated.

- **Example**: Is there a relationship between smoking status (yes/no) and disease outcome (positive/negative)?

### b. Goodness-of-fit Test
Used to test whether the distribution of a single categorical variable fits a specified distribution.

- **Example**: Do observed coin toss results match the expected 50/50 distribution?

---

## 3. Assumptions of Chi-square Tests
- Observations are independent
- Categories are mutually exclusive
- Expected frequency in each cell should be ≥ 5 (for validity)

---

## 4. Chi-square Test Formula
\[ \chi^2 = \sum \frac{(O - E)^2}{E} \]
Where:
- \( O \): Observed frequency
- \( E \): Expected frequency

---

## 5. Interpreting Results
- **Null hypothesis (H0)**: No association between variables (or observed distribution matches expected)
- **Alternative hypothesis (Ha)**: There is an association (or distribution differs)
- If p-value < 0.05, reject H0

---

## 6. Example (R Code)
```r
# Test of independence
chisq.test(table(data$smoking, data$disease))

# Goodness-of-fit test
observed <- c(25, 30, 45)
expected <- c(33.3, 33.3, 33.3)
chisq.test(x = observed, p = expected/sum(expected))
```

---

## 7. Summary
Chi-square tests are non-parametric and widely used for analyzing categorical data. Ensure assumptions are met and interpret the results in the context of the data.
