# Ensembling Methods

集成學習（Ensemble Learning）是將多個模型結合起來以提升預測表現的技術，透過多模型的整合能有效降低偏差（bias）、變異（variance），並提高泛化能力。

---

### 1. 類型分類

集成方法主要分為以下三大類：

| 方法類型     | 核心策略            | 特徵           |
| -------- | --------------- | ------------ |
| Bagging  | 並行訓練子模型，降低變異    | 高穩定性、可並行運算   |
| Boosting | 順序訓練子模型，矯正前模型錯誤 | 漸進式學習、能減偏差   |
| Stacking | 多層模型結合最終輸出      | 結合異質模型、提升靈活性 |

---

### 2. Bagging（Bootstrap Aggregating）

從訓練集中有放回地抽樣子資料集，分別訓練多個模型，最後使用平均（回歸）或投票（分類）方式整合預測結果。

#### 常見方法：

* **Random Forest**：基於 Bagging 的決策樹集成，透過隨機選擇特徵減少相關性。

```python
from sklearn.ensemble import RandomForestClassifier
clf = RandomForestClassifier(n_estimators=100)
clf.fit(X_train, y_train)
```

---

### 3. Boosting

每一輪模型學習前一輪的殘差，聚焦難分類樣本。

#### 常見方法：

* **AdaBoost**：依照錯誤率分配樣本權重。
* **Gradient Boosting**：以梯度下降優化損失函數。
* **XGBoost / LightGBM**：效能與準確率優化的梯度提升樹實作。

```python
from sklearn.ensemble import GradientBoostingClassifier
clf = GradientBoostingClassifier(n_estimators=100)
clf.fit(X_train, y_train)
```

---

### 4. Stacking（堆疊泛化）

利用多個基模型（level-0），並將其預測結果作為新特徵，輸入至最終的模型（level-1）進行預測。

```python
from sklearn.ensemble import StackingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC

estimators = [('svc', SVC(probability=True)), ('rf', RandomForestClassifier())]
stacking = StackingClassifier(estimators=estimators, final_estimator=LogisticRegression())
stacking.fit(X_train, y_train)
```

---

### 5. 優缺點比較

| 方法       | 優點          | 缺點          |
| -------- | ----------- | ----------- |
| Bagging  | 降低變異、提升穩定性  | 不易改善偏差      |
| Boosting | 減少偏差、擬合能力強  | 訓練耗時、對噪聲較敏感 |
| Stacking | 結合多種模型、擴展性佳 | 易過擬合、結構複雜   |

---

集成學習是提升模型效能的重要利器。根據資料特性與模型需求選擇適當策略，常可顯著提升預測表現，並廣泛應用於競賽、金融、生醫、推薦系統等領域。
