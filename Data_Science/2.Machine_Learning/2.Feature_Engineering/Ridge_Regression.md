# Ridge Regression（L2 正則化回歸）

Ridge Regression 是一種在線性回歸模型中加入 L2 正則化的技巧，能有效解決過擬合問題，並對多重共線性具有良好處理能力。

---

## 1. 模型定義
對於基本線性模型：
\[ \hat{y} = Xw + b \]

Ridge Regression 的損失函數為：
\[ \text{Loss} = \sum (y_i - \hat{y}_i)^2 + \lambda \sum w_j^2 \]

其中：
- 第一項為 MSE（均方誤差）
- 第二項為權重平方和（L2 penalty）
- \( \lambda \) 控制正則化強度（越大越懲罰大權重）

---

## 2. 特點與優勢
- 抑制模型過度學習訓練資料（過擬合）
- 保留所有特徵（不會產生稀疏解）
- 對特徵間具有相關性（共線性）時能穩定估計

---

## 3. 實作（scikit-learn 範例）
```python
from sklearn.linear_model import Ridge
from sklearn.datasets import load_diabetes
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

# 載入資料
X, y = load_diabetes(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)

# 建立 Ridge 模型
model = Ridge(alpha=1.0)
model.fit(X_train, y_train)

# 模型評估
y_pred = model.predict(X_test)
print("MSE:", mean_squared_error(y_test, y_pred))
print("權重大小:", model.coef_)
```

---

## 4. 超參數 alpha 的調整
- `alpha` 越大 → 模型越平滑、但可能欠擬合
- 可使用 `RidgeCV` 交叉驗證選最佳 alpha：
```python
from sklearn.linear_model import RidgeCV

model_cv = RidgeCV(alphas=[0.1, 1.0, 10.0]).fit(X_train, y_train)
print("最佳 alpha:", model_cv.alpha_)
```

---

## 5. Ridge vs Lasso vs ElasticNet
| 模型       | 正則化方式    | 特徵選擇能力 | 適用情境                         |
|------------|----------------|----------------|----------------------------------|
| Lasso      | L1             | 強（可為 0）   | 特徵多且希望剃除部分變數         |
| Ridge      | L2             | 無（保留全部） | 特徵多但都重要，有共線性時穩定   |
| ElasticNet | L1 + L2        | 中             | 特徵多且高度相關，需平衡兩者     |

---

## 6. 注意事項
- 較不適合在特徵非常稀疏或不重要特徵很多的情境（改用 Lasso）
- 請搭配標準化（StandardScaler），以避免特徵尺度影響正則化效果

---

## 7. 延伸閱讀
- `Lasso Regression.md`
- `ElasticNet.md`
- `FeatureSelection.md`
