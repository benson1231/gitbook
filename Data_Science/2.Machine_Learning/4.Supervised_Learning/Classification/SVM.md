# Support Vector Machines

Support Vector Machine（SVM）是一種監督式學習演算法，用於分類與回歸任務。它的核心目標是找到一條最佳超平面（optimal hyperplane）以最大化類別間的間隔（margin）。

---

### 1. 核心概念

SVM 尋找一條將不同類別資料分開的超平面，且兩側到最近資料點的距離（即 margin）最大。這些最接近邊界的點稱為**支持向量（support vectors）**。

* 若資料線性可分，則 SVM 可直接找到分隔超平面。
* 若資料非線性可分，則透過「核技巧（kernel trick）」將資料映射至高維空間，使其可分。

---

### 2. 數學形式（線性可分）

給定訓練資料 \$(x\_i, y\_i)\$，其中 \$y\_i \in {-1, 1}\$，SVM 解以下最佳化問題：

$$
\min_{\mathbf{w}, b} \frac{1}{2} \|\mathbf{w}\|^2 \\
\text{subject to } y_i(\mathbf{w}^T x_i + b) \geq 1 \quad \forall i
$$

當資料不可完全分開時，加入鬆弛變數與懲罰項形成 soft-margin SVM。

---

### 3. 常見核函數（Kernels）

| 核函數類型     | 表示                                        | 適用場景        |
| --------- | ----------------------------------------- | ----------- |
| 線性核       | \$K(x, x') = x \cdot x'\$                 | 特徵維度高、線性可分  |
| 多項式核      | \$K(x, x') = (x \cdot x' + c)^d\$         | 帶有交互特徵的中階資料 |
| 高斯 RBF 核  | \$K(x, x') = \exp(-\gamma \|x - x'\|^2)\$ | 資料非線性、形狀彎曲  |
| Sigmoid 核 | \$\tanh(\alpha x \cdot x' + c)\$          | 某些神經網路模擬應用  |

---

### 4. 優缺點

| 優點                 | 缺點                   |
| ------------------ | -------------------- |
| 適合高維稀疏資料           | 訓練時間與記憶體成本較高         |
| 可以有效處理非線性分類（使用核技巧） | 對超參數（如 C 與 kernel）敏感 |
| 有良好的理論基礎與泛化能力      | 難以直觀解釋模型             |

---

### 5. Python 實作範例（使用 RBF 核）

```python
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt

# 產生資料
X, y = datasets.make_classification(n_samples=100, n_features=2, n_informative=2,
                                    n_redundant=0, n_clusters_per_class=1, random_state=42)

# 分割資料
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 建立與訓練模型
model = SVC(kernel='rbf', C=1.0, gamma='scale')
model.fit(X_train, y_train)

# 預測與評估
y_pred = model.predict(X_test)
print(classification_report(y_test, y_pred))

# 視覺化資料點
plt.scatter(X[:, 0], X[:, 1], c=y, cmap='bwr', edgecolors='k')
plt.title('SVM Classification Result')
plt.show()
```

---

### 6. 應用範疇

* 文字與情感分類
* 圖像識別與臉部辨識
* 醫學資料分類與診斷
* 生物資訊：基因表現分類

---

SVM 為經典的強大分類工具，在高維與小樣本場景下尤其表現優異。其核技巧拓展了非線性資料處理的可能，是深度學習前的主力方法之一。
