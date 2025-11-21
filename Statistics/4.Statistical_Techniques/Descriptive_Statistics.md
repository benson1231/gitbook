# Descriptive Statistics

Descriptive statistics summarize and organize the features of a dataset, providing simple summaries about the sample and the measures.

---

## 1. Types of Descriptive Statistics

### a. Measures of Central Tendency
- **Mean**: Average of all values
- **Median**: Middle value when data is ordered
- **Mode**: Most frequently occurring value

### b. Measures of Dispersion
- **Range**: Difference between maximum and minimum values
- **Variance**: Average of squared deviations from the mean
- **Standard Deviation (SD)**: Square root of variance
- **Interquartile Range (IQR)**: Range between the 25th and 75th percentiles

### c. Shape of the Distribution
- **Skewness**: Measure of asymmetry
- **Kurtosis**: Measure of peakedness

---

## 2. Frequency Distribution
- Organizes data into classes or intervals
- Commonly visualized with histograms, bar charts, or frequency tables

---

## 3. Visualizing Descriptive Statistics
- **Histogram**: Displays distribution of continuous variables
- **Boxplot**: Shows median, quartiles, and outliers
- **Bar Chart**: Visualizes categorical variable counts

```r
# Example boxplot in R
boxplot(data$height)

# Histogram
hist(data$weight)
```

---

## 4. Summary Statistics in R
```r
# Basic summary
summary(data)

# Specific statistics
mean(data$score)
sd(data$score)
quantile(data$score, probs = c(0.25, 0.75))
```

---

## 5. Summary
Descriptive statistics provide foundational insights into the data before applying more complex analyses. They are essential for understanding data distribution, detecting outliers, and informing appropriate statistical methods.
