# Linear Discriminant Analysis

Linear Discriminant Analysis（LDA）是一種監督式降維與分類方法，目標是在保留類別區分度的前提下，將高維資料投影至較低維空間。

---

### 1. 模型核心概念

LDA 假設每個類別的資料在同一個高斯分佈下，且具有相同的共變異矩陣。它的目標是最大化類間變異（between-class scatter）與最小化類內變異（within-class scatter）之比值，以尋找最佳投影方向。

---

### 2. 數學推導

給定 \$c\$ 個類別，每類有 \$n\_i\$ 筆樣本：

* 類內散佈矩陣：

$$
S_W = \sum_{i=1}^c \sum_{x \in D_i} (x - \mu_i)(x - \mu_i)^T
$$

* 類間散佈矩陣：

$$
S_B = \sum_{i=1}^c n_i(\mu_i - \mu)(\mu_i - \mu)^T
$$

最佳投影方向 \$w\$ 滿足：

$$
\hat{w} = \arg\max \frac{w^T S_B w}{w^T S_W w}
$$

---

### 3. 與 PCA 差異

| 比較項目 | LDA          | PCA      |
| ---- | ------------ | -------- |
| 類別資訊 | 有            | 無        |
| 目標   | 區分類別、最大化類間差異 | 最大化資料變異量 |
| 性質   | 監督式降維        | 非監督式降維   |

---

### 4. 優缺點

| 優點             | 缺點                 |
| -------------- | ------------------ |
| 考慮類別資訊，適合分類任務  | 假設共變異矩陣相同，現實中可能不滿足 |
| 可與其他模型結合使用     | 若特徵數大於樣本數，矩陣可能奇異   |
| 在高維資料下能有效降維與分群 | 僅限於線性邊界分類問題        |

---

### 5. Python 實作範例

```python
from sklearn.datasets import load_iris
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import matplotlib.pyplot as plt
import pandas as pd

# 載入資料
iris = load_iris()
X = iris.data
y = iris.target

# 建立 LDA 模型並轉換資料
lda = LinearDiscriminantAnalysis(n_components=2)
X_lda = lda.fit_transform(X, y)

# 畫圖
plt.figure()
colors = ['r', 'g', 'b']
for i, color in enumerate(colors):
    plt.scatter(X_lda[y == i, 0], X_lda[y == i, 1], label=iris.target_names[i], color=color)
plt.title('LDA: Iris Projection')
plt.xlabel('LD1')
plt.ylabel('LD2')
plt.legend()
plt.grid(True)
plt.show()
```

---

### 6. 應用場景

* 生物資訊：癌症分類、基因型鑑別
* 醫療影像分群
* 語音辨識、面部辨識
* 作為降維前處理方法以提升分類效能

---

LDA 是兼具降維與分類能力的工具，在類別分佈明確時尤其有效。與 PCA、SVM 等方法搭配，可組成強大的機器學習流程。
