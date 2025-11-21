# Random Forests

隨機森林是一種集成學習法，屬於 Bagging 方法的一種。它結合多棵決策樹進行預測，以提升模型的準確性與穩定性，並降低過擬合風險。

---

## 1. 基本概念

隨機森林是由多棵「隨機訓練的決策樹」所組成，每棵樹都是在不同的 bootstrap 資料子集中訓練得到，並在每次分裂時隨機選擇部分特徵。

* 分類任務：使用多數投票（majority voting）決定結果。
* 回歸任務：使用平均值作為預測結果。

---

## 2. 特點與優勢

* **降低過擬合風險**：透過隨機抽樣與特徵隨機選擇降低模型變異。
* **可處理高維資料**：適用於特徵數大於樣本數的問題。
* **提供特徵重要性指標**：可用來解釋哪些特徵對預測最具貢獻。
* **適用於分類與回歸**：在多類別分類與連續型預測都具優異表現。

---

## 3. Python 實作

```python
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

# 載入資料
X, y = load_iris(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 建立模型
clf = RandomForestClassifier(n_estimators=100, max_features='sqrt', random_state=42)
clf.fit(X_train, y_train)

# 預測與評估
y_pred = clf.predict(X_test)
print("Accuracy:", accuracy_score(y_test, y_pred))
```

---

## 4. 常見參數

| 參數名稱           | 說明             |
| -------------- | -------------- |
| `n_estimators` | 森林中樹的數量        |
| `max_features` | 每棵樹分裂時最多考慮的特徵數 |
| `max_depth`    | 樹的最大深度         |
| `bootstrap`    | 是否使用有放回抽樣      |
| `random_state` | 固定隨機性以利重現實驗結果  |

---

## 5. 特徵重要性視覺化

```python
import matplotlib.pyplot as plt
import pandas as pd

feat_importance = pd.Series(clf.feature_importances_, index=load_iris().feature_names)
feat_importance.sort_values().plot(kind='barh')
plt.title("Feature Importances")
plt.show()
```

---

## 6. 應用範疇

* 生物醫學特徵選擇（如基因篩選）
* 信用風險評估
* 圖像分類與文字分類
* 推薦系統中的用戶行為預測

---

隨機森林以其高效能與穩健性成為實務上常見的建模工具。尤其在資料較雜或特徵數眾多時，能提供良好的預測表現與模型解釋力。
