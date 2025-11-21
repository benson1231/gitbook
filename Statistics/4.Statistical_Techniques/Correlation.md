# Correlation

Correlation measures the strength and direction of the linear relationship between two continuous variables.

---

## 1. When to Use Correlation
- Exploring the relationship between height and weight
- Assessing association between study hours and test scores
- Checking the link between temperature and electricity usage

---

## 2. Types of Correlation Coefficients

### a. Pearson Correlation (r)
Measures the strength and direction of a **linear** relationship between two continuous variables.
- Range: -1 to +1
- Assumes normal distribution and linearity

### b. Spearman Rank Correlation (ρ)
Non-parametric measure based on **ranked** values.
- Used when data is not normally distributed or has outliers

### c. Kendall’s Tau
Another rank-based correlation measure; more robust for small samples

---

## 3. Interpreting Correlation Coefficients
| Correlation (r) | Strength         |
|------------------|-------------------|
| 0.00–0.19        | Very weak         |
| 0.20–0.39        | Weak              |
| 0.40–0.59        | Moderate          |
| 0.60–0.79        | Strong            |
| 0.80–1.00        | Very strong       |

- Positive value: as one variable increases, the other tends to increase
- Negative value: as one increases, the other tends to decrease

---

## 4. Assumptions for Pearson Correlation
- Both variables are continuous
- Linearity
- Normal distribution
- No significant outliers

---

## 5. Example (R Code)
```r
# Pearson correlation
cor(data$height, data$weight, method = "pearson")

# Spearman correlation
cor(data$rank1, data$rank2, method = "spearman")

# Correlation test with p-value
cor.test(data$height, data$weight)
```

---

## 6. Visualizing Correlation
- **Scatter plots** to inspect linearity
- **Correlation matrices** for multiple variables
```r
# Base scatter plot
plot(data$height, data$weight)

# Correlation matrix (e.g., with corrplot package)
corrplot(cor(data[, c("var1", "var2", "var3")]))
```

---

## 7. Summary
Correlation quantifies the strength of association between variables. Choose the appropriate method based on data distribution and scale. Always visualize your data to check for patterns and anomalies.
