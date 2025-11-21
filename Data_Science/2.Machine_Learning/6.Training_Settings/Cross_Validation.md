# Cross-Validation

交叉驗證是一種模型驗證方法，將資料集劃分為多個子集，用於訓練與驗證模型效能，避免過擬合並提升泛化能力。

---

## 1. 為什麼要交叉驗證？

在建立機器學習模型時，如果僅依賴單一訓練／測試切割，可能因資料切割偏差導致結果不穩定。交叉驗證透過多次切割與評估平均效能，提供更可靠的評估依據。

---

## 2. 常見方法

### (1) K-Fold Cross-Validation

將資料隨機分為 K 個子集（folds）：

1. 每次使用 K-1 個子集作為訓練集，剩下的 1 個為驗證集。
2. 重複 K 次，每個子集都當過一次驗證集。
3. 將 K 次評估結果取平均。

```python
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.datasets import load_iris

X, y = load_iris(return_X_y=True)
model = LogisticRegression(max_iter=200)
scores = cross_val_score(model, X, y, cv=5)
print("Average Accuracy:", scores.mean())
```

### (2) Stratified K-Fold

在分類任務中，保留每個類別在每個 fold 中的比例一致，適合不平衡資料集。

### (3) Leave-One-Out (LOO)

每次只留一筆資料作為驗證，其餘作為訓練。適用於小型資料集，但計算成本高。

### (4) Repeated K-Fold

多次重複 K-Fold，可進一步降低樣本分配對模型評估的變異。

---

## 3. 模型選擇與調參結合

交叉驗證常與超參數搜尋結合，例如 GridSearchCV、RandomizedSearchCV，可同時進行模型評估與參數最佳化。

```python
from sklearn.model_selection import GridSearchCV

param_grid = {'C': [0.1, 1, 10]}
clf = GridSearchCV(LogisticRegression(), param_grid, cv=5)
clf.fit(X, y)
print("Best Parameters:", clf.best_params_)
```

---

## 4. 優點與限制

| 優點            | 限制               |
| ------------- | ---------------- |
| 可減少模型對切割偏誤敏感度 | 計算成本高（特別是在大型資料集） |
| 適用於模型比較與選擇    | 較難平行化（部分方法除外）    |
| 與調參可整合        | 小資料集仍可能高變異       |

---

交叉驗證是建立穩健模型不可或缺的一步，尤其在資料量有限或需評估多個模型／參數組合時，能大幅提升選擇正確模型的機會。
