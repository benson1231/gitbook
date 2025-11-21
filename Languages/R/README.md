# R

R 是一種專門用於統計分析與資料視覺化的程式語言，廣泛應用於生物資訊、金融、社會科學等資料密集領域。它提供豐富的函數庫與繪圖工具，適合進行探索性資料分析與模型建構。

---

## 1. 基本語法

```r
# 變數指定
x <- 5
y = 10

# 向量與矩陣
vec <- c(1, 2, 3)
mat <- matrix(1:9, nrow=3)

# 函數定義與呼叫
square <- function(x) {
  return(x^2)
}
square(3)

# 條件與迴圈
if (x > 0) {
  print("正數")
}

for (i in 1:5) {
  print(i)
}
```

---

## 2. 資料處理套件（Tidyverse）

Tidyverse 是 R 中用於資料清理與分析的整合套件集合，包含：

* `dplyr`：資料篩選、排序、群組計算
* `tidyr`：資料轉置與整理
* `readr`：讀取 csv/tsv 檔案
* `ggplot2`：繪圖套件（語法風格優雅）

```r
library(dplyr)

data <- tibble(name=c("Amy", "Bob"), score=c(88, 95))

# 篩選與計算平均
filtered <- data %>% filter(score > 90)
avg <- data %>% summarize(mean_score = mean(score))
```

---

## 3. 資料視覺化（ggplot2）

```r
library(ggplot2)

ggplot(data, aes(x=name, y=score)) +
  geom_bar(stat="identity") +
  theme_minimal()
```

常見圖表類型：

* 長條圖（`geom_bar`）
* 折線圖（`geom_line`）
* 散點圖（`geom_point`）
* 箱型圖（`geom_boxplot`）

---

## 4. 機器學習與統計分析

R 提供大量套件支援回歸、分類、聚類與高維資料處理，如：

* `caret`：整合各種模型訓練與交叉驗證流程
* `glm()`：邏輯回歸與廣義線性模型
* `randomForest`, `xgboost`：高效分類器

```r
model <- glm(y ~ x1 + x2, data=df, family="binomial")
predict(model, newdata=df, type="response")
```

---

## 5. R 與生物資訊

R 在生物資訊學界有高度應用，以下是幾個常用套件：

* `Bioconductor`：基因表現、RNA-seq、微陣列分析等套件集合
* `DESeq2`、`edgeR`：差異基因分析
* `clusterProfiler`：功能性富集分析（GO、KEGG）
* `Seurat`：單細胞 RNA-seq 資料分析

```r
library(clusterProfiler)
result <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
```

---

R 是一種功能強大、社群活躍的分析語言。若以資料處理、統計建模或生物資料為核心工作，R 是非常值得學習與掌握的工具。
