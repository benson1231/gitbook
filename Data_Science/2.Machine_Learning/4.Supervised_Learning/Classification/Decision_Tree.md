# Decision Tree

決策樹是一種監督式學習演算法，適用於分類與回歸任務。它透過學習資料中的特徵條件來建立樹狀結構，做出決策或預測。

---

### 1. 演算法原理

決策樹模型以「遞迴分割（recursive partitioning）」方式進行：

* 由根節點開始，選擇能最佳區分目標變數的特徵做分割。
* 根據特徵值分割出子節點，重複進行直到終止條件滿足（如純度、深度限制）。
* 預測時根據輸入特徵從樹根一路走到葉節點。

---

### 2. 資訊衡量指標

分類任務中，常見的純度衡量指標如下：

* **Gini 指數（Gini Impurity）**：

$$
Gini = 1 - \sum p_i^2
$$

* **資訊增益（Information Gain）**（基於 Entropy）：

$$
Entropy = -\sum p_i \log_2(p_i)
$$

* **分類誤差（Classification Error）**：

$$
Error = 1 - \max(p_i)
$$

回歸任務中，則使用：

* 均方差（MSE）
* 絕對差（MAE）

---

### 3. 優缺點比較

| 優點         | 缺點            |
| ---------- | ------------- |
| 易於理解與視覺化   | 容易過擬合（尤其深樹）   |
| 可處理數值與類別資料 | 對資料微小變化敏感     |
| 不需特徵標準化    | 分裂偏好多值特徵，造成偏差 |
| 支援非線性決策邊界  | 單棵樹預測準確率不一定高  |

---

### 4. Python 實作

```python
from sklearn.tree import DecisionTreeClassifier
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

iris = load_iris()
X_train, X_test, y_train, y_test = train_test_split(
    iris.data, iris.target, test_size=0.3, random_state=42)

clf = DecisionTreeClassifier(criterion='gini', max_depth=3)
clf.fit(X_train, y_train)

y_pred = clf.predict(X_test)
print("準確率：", accuracy_score(y_test, y_pred))
```

---

### 5. 決策樹視覺化

```python
from sklearn.tree import plot_tree
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plot_tree(clf, feature_names=iris.feature_names, class_names=iris.target_names,
           filled=True, rounded=True)
plt.title("Decision Tree")
plt.show()
```

---

### 6. 常見參數說明

* `criterion`：分裂準則（'gini' 或 'entropy'）
* `max_depth`：樹的最大深度
* `min_samples_split`：內部節點再分裂所需的最小樣本數
* `min_samples_leaf`：葉節點的最小樣本數
* `max_features`：考慮的最大特徵數量

---

### 7. 應用場景

* 客戶分類、信用風險評估
* 疾病診斷、症狀推論
* 市場行為預測、客群細分

---

決策樹是構建集成模型（如隨機森林、梯度提升樹）的基礎模型，了解其架構對後續進階模型的理解非常重要。
