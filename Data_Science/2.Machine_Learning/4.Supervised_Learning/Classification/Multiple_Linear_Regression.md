# Multiple Linear Regression

多元線性回歸是線性回歸的延伸，當預測變數（自變數）不只一個時，可使用多元線性回歸模型來描述多個變數與應變數之間的線性關係。

---

### 1. 模型形式

$$
y = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \dots + \beta_p x_p + \epsilon
$$

其中：

* $y$：應變數
* $x_1, x_2, \dots, x_p$：多個自變數（共 $p$ 個）
* $\beta_0$：截距
* $\beta_1, \dots, \beta_p$：各變數的迴歸係數
* $\epsilon$：誤差項，代表隨機擾動

---

### 2. 最小平方法與向量表示

可將模型寫成向量矩陣形式：

$$
\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}
$$

* $\mathbf{y}$：$n \times 1$ 應變數向量
* $\mathbf{X}$：$n \times (p+1)$ 設計矩陣（包含常數項）
* $\boldsymbol{\beta}$：$(p+1) \times 1$ 參數向量
* $\boldsymbol{\epsilon}$：誤差向量

最小平方法估計式為：

$$
\hat{\boldsymbol{\beta}} = (\mathbf{X}^T\mathbf{X})^{-1}\mathbf{X}^T\mathbf{y}
$$

---

### 3. 假設與誤差項

與簡單線性回歸相同，誤差項 $\epsilon$ 須滿足：

* 常態分布：$\epsilon \sim N(0, \sigma^2)$
* 同方差性（homoscedasticity）：誤差變異數不隨 $x_i$ 而變化
* 獨立性：誤差彼此獨立
* 無多重共線性：自變數間不應高度相關

---

### 4. 評估指標

* **$R^2$**：模型對變異解釋能力
* **Adjusted $R^2$**：修正 $R^2$，考慮變數數量
* **AIC/BIC**：模型選擇指標
* **殘差分析**：檢查違反假設之情形（如殘差常態性）

---

### 5. 機器學習角度的損失函數（Loss Function）

多元線性回歸的損失函數仍為均方誤差（MSE）：

$$
\mathcal{L}(\boldsymbol{\beta}) = \frac{1}{n} \sum_{i=1}^n \left(y_i - (\beta_0 + \beta_1 x_{i1} + \dots + \beta_p x_{ip})\right)^2
$$

透過梯度下降法等演算法來最小化損失，學習出最佳參數。

---

### 6. Python 實作範例

```python
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

# 模擬資料
X = np.array([[1, 2], [2, 1], [3, 4], [4, 3], [5, 5]])
y = np.array([3, 3.5, 6, 6.5, 8])

model = LinearRegression()
model.fit(X, y)

print("截距：", model.intercept_)
print("係數：", model.coef_)
```

---

### 7. 延伸應用

* **特徵選擇（Feature Selection）**：例如逐步回歸、LASSO 等方法
* **正則化方法**：Ridge、Lasso 來處理多重共線性與過擬合
* **交叉驗證（Cross-validation）**：評估泛化能力

---

### 8. 結論

多元線性回歸能夠同時考慮多個變數的影響，是資料分析與機器學習中的重要基礎模型。了解其假設條件與限制，能協助我們建立更穩健且解釋力強的預測模型。
