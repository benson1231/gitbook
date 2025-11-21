# Non-parametric Tests

Non-parametric tests are statistical methods that do not assume a specific distribution for the data. They are useful when data violate assumptions such as normality or equal variances.

---

## 1. When to Use Non-parametric Tests
- Data is ordinal or ranked
- Small sample size
- Non-normal distribution
- Presence of outliers

---

## 2. Common Non-parametric Tests

### a. Mann–Whitney U Test
- Alternative to the independent t-test
- Compares the ranks of two independent groups

```r
# Example in R
wilcox.test(group1, group2)
```

### b. Wilcoxon Signed-Rank Test
- Alternative to the paired t-test
- Compares matched or paired samples

```r
# Example in R
wilcox.test(before, after, paired = TRUE)
```

### c. Kruskal–Wallis Test
- Alternative to one-way ANOVA
- Compares ranks across three or more independent groups

```r
# Example in R
kruskal.test(score ~ group, data = data)
```

### d. Friedman Test
- Alternative to repeated-measures ANOVA
- Compares ranks in matched groups

```r
# Example in R
friedman.test(values ~ time | subject, data = data)
```

---

## 3. Advantages and Limitations
### Advantages
- Fewer assumptions
- Can be used with ordinal data
- More robust to outliers

### Limitations
- Less power than parametric tests
- Harder to interpret effect size

---

## 4. Summary
Non-parametric tests provide flexible alternatives when traditional parametric assumptions are violated. They are especially valuable in exploratory analysis, clinical trials, and small sample research.