## Simple Linear Regression

簡單線性回歸是一種統計方法，用來建構一條最佳擬合線，以解釋一個自變數（$x$）與一個應變數（$y$）之間的線性關係。

---

### 1. 模型形式

簡單線性回歸模型可表示為：

$$
y = \beta_0 + \beta_1 x + \epsilon
$$

其中：

* $y$：應變數（dependent variable）
* $x$：自變數（independent variable）
* $\beta_0$：截距（intercept）
* $\beta_1$：斜率（slope）
* $\epsilon$：誤差項（error term），假設為平均為 0 的常態分布

---

### 2. 最小平方法（Least Squares Method）

目標是最小化所有觀測值的殘差平方和（RSS）：

$$
RSS = \sum_{i=1}^{n} (y_i - \hat{y}_i)^2 = \sum_{i=1}^{n} (y_i - (\hat{\beta}_0 + \hat{\beta}_1 x_i))^2
$$

最佳估計量如下：

$$
\hat{\beta}_1 = \frac{\sum_{i=1}^{n}(x_i - \bar{x})(y_i - \bar{y})}{\sum_{i=1}^{n}(x_i - \bar{x})^2}
$$

$$
\hat{\beta}_0 = \bar{y} - \hat{\beta}_1 \bar{x}
$$

---

### 3. 解釋與應用

* **斜率 $\hat{\beta}_1$**：表示當 $x$ 增加一單位時，$y$ 平均變化多少。
* **截距 $\hat{\beta}_0$**：當 $x = 0$ 時 $y$ 的預測值。

---

### 4. 誤差項的意義與假設

誤差項 $\epsilon$ 代表模型無法解釋的變異，為真實觀測值與模型預測值之間的差。

常見假設如下：

* $\mathbb{E}[\epsilon] = 0$：誤差的期望值為 0。
* $\epsilon \sim N(0, \sigma^2)$：常態分布且變異數恆定。
* 誤差之間相互獨立，且與 $x$ 無關。

這些假設是線性回歸進行參數估計與統計推論的基礎。

---

### 5. 模型評估

* **$R^2$**（決定係數）：衡量模型對 $y$ 變異的解釋能力。

$$
R^2 = 1 - \frac{\sum (y_i - \hat{y}_i)^2}{\sum (y_i - \bar{y})^2}
$$

* **殘差圖（Residual Plot）**：可檢查線性假設與常態性。
* **假設檢定**：檢驗斜率是否顯著，例如：

$$
H_0: \beta_1 = 0 \quad vs. \quad H_1: \beta_1 \neq 0
$$

---

### 6. 損失函數與機器學習觀點（Loss Function in ML）

在線性回歸的機器學習實作中，會將最小平方法視為一種最小化「損失函數」的問題。

最常見的損失函數為**均方誤差（MSE）**：

$$
\mathcal{L}(\beta_0, \beta_1) = \frac{1}{n} \sum_{i=1}^{n} (y_i - (\beta_0 + \beta_1 x_i))^2
$$

在實作中，我們透過梯度下降法等最佳化演算法，最小化此損失函數來更新模型參數。

---

### 7. 實例（Python）

```python
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

x = np.array([1, 2, 3, 4, 5]).reshape(-1, 1)
y = np.array([2, 4, 5, 4, 5])

model = LinearRegression()
model.fit(x, y)

print("截距：", model.intercept_)
print("斜率：", model.coef_[0])

# 繪圖
plt.scatter(x, y, label="資料點")
plt.plot(x, model.predict(x), color='red', label="回歸線")
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.title("簡單線性回歸")
plt.show()
```

---

### 8. 限制與假設

1. **線性關係**：$x$ 與 $y$ 的關係為線性
2. **誤差常態分布**
3. **誤差變異數恆定（同質性）**
4. **觀測值獨立**

簡單線性回歸為許多進階統計與機器學習技術的基礎，理解其假設與應用至關重要。
