# Deep Learning

深度學習（Deep Learning）是機器學習的子領域，基於人工神經網路（Artificial Neural Networks）所發展，用於處理結構複雜、規模龐大的資料，尤其在圖像、語音與自然語言處理上展現出色表現。

---

## 核心概念

* **神經元（Neuron）**：模擬生物神經元的基本運算單元，接收輸入、加權計算後輸出。
* **層（Layer）**：由多個神經元組成，包含輸入層、隱藏層與輸出層。
* **激活函數（Activation Function）**：非線性轉換，使模型能學習複雜特徵（如 ReLU、Sigmoid、Tanh）。
* **前向傳播（Forward Propagation）**：計算每層輸出，直至模型輸出結果。
* **反向傳播（Backpropagation）**：計算損失函數的梯度並更新參數。

---

## 常見深度學習架構

| 架構          | 主要用途       | 說明                         |
| ----------- | ---------- | -------------------------- |
| DNN         | 表格、數值資料    | 基本多層全連接網路（fully connected） |
| CNN         | 圖像處理       | 卷積操作可抓取區域特徵                |
| RNN / LSTM  | 時序與語言模型    | 處理序列資料，具記憶性                |
| Transformer | NLP、語音、圖像等 | 採 attention 機制，擴展性強        |

---

## TensorFlow/Keras 模型範例

```python
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense

model = Sequential([
    Dense(64, activation='relu', input_shape=(10,)),
    Dense(64, activation='relu'),
    Dense(1, activation='sigmoid')
])

model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.fit(X_train, y_train, epochs=10, batch_size=32)
```

---

## 常見超參數（Hyperparameters）

* **學習率（learning rate）**：控制每次參數更新的幅度。
* **批次大小（batch size）**：每次訓練所取樣的資料數。
* **訓練次數（epochs）**：訓練資料被完整瀏覽的次數。
* **層數與神經元數量**：模型結構設計的關鍵。

---

## 模型評估指標

* **準確率（Accuracy）**
* **損失（Loss）**
* **混淆矩陣（Confusion Matrix）**
* **AUC / ROC 曲線**

---

## 應用範疇

* 醫療影像診斷（如腫瘤偵測）
* 自動駕駛（影像與雷達融合）
* 自然語言處理（如聊天機器人、翻譯）
* 金融風險預測（詐騙偵測、信用評分）

---

深度學習藉由大量資料與強大運算能力，能自主學習特徵並在多項任務中達成超越傳統方法的表現。理解基本架構與實作流程，是進入人工智慧世界的重要基礎。
