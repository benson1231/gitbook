# Boosting

Boosting 是一種強大的集成學習方法，透過逐步訓練一系列弱學習器（如淺層決策樹），每個模型都試圖修正前一個模型犯下的錯誤，最終組合成強大的預測器。

---

## 1. 基本流程

Boosting 的基本流程如下：

1. 初始化模型並計算初始誤差。
2. 每一輪訓練一個新的弱模型，聚焦在前一輪預測錯誤的樣本上（給予較高權重）。
3. 將新的模型加入整體模型中，根據權重調整合併方式。
4. 重複直到達到預設迭代次數或誤差停止改善。

---

## 2. 代表演算法

### (1) AdaBoost（Adaptive Boosting）

* 每輪調整樣本權重，強化模型對錯誤分類樣本的學習。
* 適合用於二元分類問題。

### (2) Gradient Boosting

* 使用梯度下降法最小化任意損失函數。
* 可應用於分類與回歸問題。

### (3) XGBoost / LightGBM / CatBoost

* 高效能的 Gradient Boosting 實作，支援並行、剪枝、正則化與缺失值處理。
* 在資料科學競賽中廣泛使用。

---

## 3. 優點與缺點

| 優點            | 缺點             |
| ------------- | -------------- |
| 擬合能力強，能處理複雜關係 | 訓練時間長，需謹慎調參    |
| 可使用不同損失函數     | 對雜訊敏感，可能過擬合    |
| 支援特徵重要性分析     | 難以平行化（部分演算法例外） |

---

## 4. Python 實作：Gradient Boosting

```python
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

# 資料分割
X, y = load_iris(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 建立與訓練模型
clf = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3)
clf.fit(X_train, y_train)

# 預測與評估
y_pred = clf.predict(X_test)
print("Accuracy:", accuracy_score(y_test, y_pred))
```

---

## 5. 應用場景

* 信用風險建模
* 醫療預測與風險評估
* 顧客流失分析
* Kaggle 資料競賽

---

Boosting 是在建模精度與表現上最為強大的方法之一，適合對結果要求高的任務。透過對錯誤的逐步修正與模型疊加，能在保持模型可解釋性的同時，獲得優異的預測效能。
