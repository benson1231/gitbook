# Bagging（Bootstrap Aggregating）

Bagging（Bootstrap Aggregating）是一種集成學習技術，透過隨機取樣與多模型平均來降低模型的變異性與過擬合風險，特別適用於高方差模型如決策樹。

---

## 1. 基本概念

Bagging 的流程如下：

1. 從原始訓練資料中使用 bootstrap（有放回）方式抽樣，生成多個資料子集。
2. 在每個子集上分別訓練一個基模型（例如決策樹）。
3. 對於分類任務，最終預測使用多數投票；對於回歸任務，使用平均值。

這種方式能降低模型的變異性，進而提升泛化能力。

---

## 2. 優勢

* **降低過擬合風險**：即便基模型容易過擬合，Bagging 可透過集成減少變異。
* **可併行訓練模型**：每個基模型獨立訓練，可進行並行運算。
* **提升穩定性**：對於輸入資料的小變動不會造成結果大幅波動。

---

## 3. 典型實作：Random Forest

Random Forest 就是以 Bagging 為基礎的集成方法，並在訓練每棵樹時加入隨機特徵選擇以增強多樣性。

```python
from sklearn.ensemble import BaggingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

# 載入資料
X, y = load_iris(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Bagging 模型
base = DecisionTreeClassifier()
bag = BaggingClassifier(base_estimator=base, n_estimators=100, random_state=42)
bag.fit(X_train, y_train)

# 預測與評估
y_pred = bag.predict(X_test)
print("Accuracy:", accuracy_score(y_test, y_pred))
```

---

## 4. 常見應用

* 決策樹的穩定性提升
* 難以處理高變異資料的建模任務
* 提高小型樣本集模型的一致性

---

Bagging 是一種簡單卻強大的集成方法，尤其適合與不穩定學習器（如決策樹）搭配使用。透過多次抽樣與結果平均，它能有效減少模型對訓練資料的依賴，增強整體穩健性。
