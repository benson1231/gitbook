# sklearn

`scikit-learn`（簡稱 `sklearn`）是 Python 中最常用的機器學習套件之一，提供分類、迴歸、分群、降維、模型選擇與資料前處理等功能。

---

## 一、安裝

```bash
pip install scikit-learn
```

---

## 二、主要功能模組

| 類別                        | 說明                                 |
| ------------------------- | ---------------------------------- |
| `sklearn.datasets`        | 提供內建資料集（如 iris、digits）             |
| `sklearn.model_selection` | 訓練/測試切分、交叉驗證、超參數搜尋                 |
| `sklearn.preprocessing`   | 特徵標準化、編碼、縮放                        |
| `sklearn.linear_model`    | 線性與羅吉斯迴歸等模型                        |
| `sklearn.tree`            | 決策樹與隨機森林                           |
| `sklearn.svm`             | 支援向量機                              |
| `sklearn.neighbors`       | KNN 與最近鄰搜尋                         |
| `sklearn.naive_bayes`     | 樸素貝氏分類器                            |
| `sklearn.cluster`         | 分群（KMeans、DBSCAN 等）                |
| `sklearn.decomposition`   | 降維（PCA、NMF 等）                      |
| `sklearn.metrics`         | 評估指標（accuracy、f1、confusion matrix） |

---

## 三、典型工作流程

```python
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

# 載入資料
X, y = load_iris(return_X_y=True)

# 分割訓練/測試集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 建立模型
model = RandomForestClassifier()
model.fit(X_train, y_train)

# 預測與評估
y_pred = model.predict(X_test)
print("Accuracy:", accuracy_score(y_test, y_pred))
```

---

## 四、常用前處理工具

```python
from sklearn.preprocessing import StandardScaler, OneHotEncoder

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
```

---

## 五、交叉驗證與模型選擇

```python
from sklearn.model_selection import cross_val_score, GridSearchCV

scores = cross_val_score(model, X, y, cv=5)
print("平均交叉驗證分數：", scores.mean())
```

---

## 六、應用場景

* 分類任務（疾病預測、影像辨識）
* 迴歸任務（價格預測、風險評估）
* 分群任務（顧客群分析）
* 降維與視覺化（PCA, t-SNE）

---

scikit-learn 提供一致的 API 與高效能實作，適合快速原型開發與教學使用，是進入機器學習實作的首選套件。
