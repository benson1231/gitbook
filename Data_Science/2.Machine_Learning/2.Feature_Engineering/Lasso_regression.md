# Lasso Regression（L1 正則化回歸）

Lasso（Least Absolute Shrinkage and Selection Operator）是一種線性回歸方法，加入了 L1 正則化項，使模型在減少過擬合的同時，還能執行特徵選擇。

---

## 1. 模型定義
對於線性回歸模型：
\[ \hat{y} = Xw + b \]

Lasso 的損失函數為：
\[ \text{Loss} = \sum (y_i - \hat{y}_i)^2 + \lambda \sum |w_j| \]

其中：
- 第一項是平方誤差（MSE）
- 第二項是權重係數的絕對值總和（L1 penalty）
- \( \lambda \) 控制正則化強度（越大代表壓縮越多）

---

## 2. 特點與優勢
- 可自動將不重要的特徵係數壓縮為 0 → 實現**特徵選擇**
- 與 Ridge Regression（L2）相比，更能產生稀疏模型
- 適合高維度、特徵冗餘的問題（如：文字分類、基因資料）

---

## 3. 實作（scikit-learn 範例）
```python
from sklearn.linear_model import Lasso
from sklearn.datasets import load_diabetes
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

# 載入資料
X, y = load_diabetes(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)

# 建立 Lasso 模型
model = Lasso(alpha=1.0)
model.fit(X_train, y_train)

# 模型評估
y_pred = model.predict(X_test)
print("MSE:", mean_squared_error(y_test, y_pred))
print("選中的特徵數量:", (model.coef_ != 0).sum())
```

---

## 4. 超參數 alpha 的調整
- `alpha` 越大 → 正則化越強，特徵被壓縮越多
- 常用 `LassoCV` 自動搜尋最佳 alpha：
```python
from sklearn.linear_model import LassoCV

model_cv = LassoCV(cv=5).fit(X_train, y_train)
print("最佳 alpha:", model_cv.alpha_)
```

---

## 5. Lasso vs Ridge vs ElasticNet
| 模型       | 正則化方式    | 特徵選擇能力 | 適用情境                         |
|------------|----------------|----------------|----------------------------------|
| Lasso      | L1             | 強（可為 0）   | 高維稀疏資料、希望特徵剃除        |
| Ridge      | L2             | 無（不為 0）   | 特徵都重要，但需抑制過擬合        |
| ElasticNet | L1 + L2        | 中             | 特徵多且具相關性（混合效果最佳）  |

---

## 6. 注意事項
- Lasso 偏好在特徵數量大於樣本數時使用
- 若特徵高度相關，Lasso 可能只保留其中一個 → 可考慮 ElasticNet
- 建議搭配標準化（StandardScaler）

---

## 7. 延伸閱讀
- `Ridge Regression.md`
- `ElasticNet.md`
- `FeatureSelection.md`
