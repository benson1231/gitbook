# Optimizer

在機器學習與深度學習中，Optimizer（最佳化器）是用來調整模型參數以最小化損失函數（loss function）的演算法。透過反向傳播（backpropagation）計算梯度後，optimizer 根據這些梯度來更新模型的權重與偏差。

---

## 為什麼需要 Optimizer？

神經網路透過訓練資料來學習模式，其目標是讓模型的預測盡可能接近真實標籤。這個目標是透過損失函數來衡量的，optimizer 則用來最小化這個損失。

---

## 常見 Optimizer 演算法

### 1. Gradient Descent（梯度下降法）

* **原理**：沿著損失函數的梯度方向更新參數。
* **缺點**：對學習率（learning rate）非常敏感；可能陷入局部極小值。

### 2. Stochastic Gradient Descent（SGD，隨機梯度下降）

* **原理**：每次使用一筆樣本（或一小批 mini-batch）進行梯度更新。
* **優點**：速度快，可跳出局部極小值。
* **缺點**：更新不穩定、震盪大。

### 3. Momentum（動量法）

* **改進點**：引入過去梯度的指向，減少震盪。
* **公式**：

  ```
  v_t = γ * v_{t-1} + η * ∇L(θ)
  θ = θ - v_t
  ```

### 4. RMSProp（Root Mean Square Propagation）

* **改進點**：對每個參數維度自適應調整學習率。
* **適合**：非平穩目標問題，如 RNN。

### 5. Adam（Adaptive Moment Estimation）

* **綜合 Momentum + RMSProp**：同時考慮一階與二階矩估計。
* **優點**：收斂快、調參容易，是目前最常用的 optimizer。
* **公式概要**：

  ```
  m_t = β₁ * m_{t-1} + (1 - β₁) * ∇L(θ)
  v_t = β₂ * v_{t-1} + (1 - β₂) * (∇L(θ))²
  θ = θ - η * m_t / (√v_t + ε)
  ```

---

## 重要參數

* **Learning Rate（學習率）**：控制每次參數更新的幅度。
* **β₁, β₂（動量與 RMS 的衰減因子）**：影響更新的平滑程度（Adam/RMSProp）。
* **ε（防止除以 0 的小常數）**。

---

## 實作範例（TensorFlow）

```python
import tensorflow as tf

model = tf.keras.Sequential([
    tf.keras.layers.Dense(1, input_shape=(10,))
])

model.compile(
    optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
    loss='mean_squared_error'
)

# 假設 x_train, y_train 為訓練資料
model.fit(x_train, y_train, epochs=10, batch_size=32)
```

---

Optimizer 是訓練模型時不可或缺的元件，選擇合適的 optimizer 與調整其參數會大幅影響模型的訓練效率與最終表現。
