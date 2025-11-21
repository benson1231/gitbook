# Stacking

Stacking 是一種集成學習方法，將多個基模型的預測結果當作輸入，餵給一個次級模型（meta-learner）來做最終預測。此方法能結合多種演算法的優點，有助於提升預測精度。

---

## 1. 基本架構

1. 建立多個基模型（第一層），可使用不同的機器學習演算法。
2. 使用交叉驗證，將每個基模型的預測結果作為新特徵。
3. 將這些預測結果餵給次級模型（第二層）進行學習與最終預測。

---

## 2. 優勢與劣勢

| 優勢          | 劣勢               |
| ----------- | ---------------- |
| 結合多種模型以提高表現 | 結構較複雜，訓練時間較長     |
| 可使用不同類型的基模型 | 難以除錯與解釋          |
| 有助於降低偏差與變異  | 若資料過少或基模型重疊，效益有限 |

---

## 3. Python 範例：Scikit-learn 實作

```python
from sklearn.ensemble import StackingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

# 資料分割
X, y = load_iris(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 定義基模型
base_models = [
    ('dt', DecisionTreeClassifier()),
    ('svc', SVC(probability=True))
]

# 定義堆疊模型
stack_model = StackingClassifier(
    estimators=base_models,
    final_estimator=LogisticRegression()
)

# 訓練與預測
stack_model.fit(X_train, y_train)
y_pred = stack_model.predict(X_test)
print("Accuracy:", accuracy_score(y_test, y_pred))
```

---

## 4. 應用場景

* 結合異質模型（例如 SVM + 樹模型 + 線性模型）
* 資料競賽（如 Kaggle）中提升最終模型表現
* 難以選擇單一最佳模型時

---

Stacking 能夠有效整合多種預測模型，尤其在資料結構複雜或基模型各有優勢時特別有用。不過其訓練成本與模型結構複雜度較高，需謹慎設計與驗證。
