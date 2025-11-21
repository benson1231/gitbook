# Recommender Systems

推薦系統是一種機器學習應用，根據使用者的行為、偏好或特徵，預測並推薦其可能感興趣的項目（如商品、電影、音樂、文章等）。

---

## 1. 推薦系統的主要類型

### (1) 協同過濾（Collaborative Filtering）

依賴「使用者－物品」互動資料（如評分、點擊、購買紀錄），分為：

* **使用者式協同過濾（User-based CF）**：找出與目標使用者行為相似的其他使用者
* **物品式協同過濾（Item-based CF）**：推薦與使用者曾喜歡的物品相似的項目

```python
from sklearn.metrics.pairwise import cosine_similarity
from scipy.sparse import csr_matrix

# 假設 ratings 為使用者-物品矩陣
similarity = cosine_similarity(ratings)
```

### (2) 內容式過濾（Content-Based Filtering）

根據物品的特徵（如關鍵字、類別、描述）與使用者偏好進行配對。

* 對每個使用者建立其偏好向量
* 計算每個物品與偏好的相似度

### (3) 混合式方法（Hybrid Methods）

結合協同過濾與內容式方法，提升準確度與覆蓋率。

---

## 2. 評估指標

| 指標         | 說明                       |
| ---------- | ------------------------ |
| MAE / RMSE | 預測評分誤差                   |
| Precision  | 推薦項目中正確的比例               |
| Recall     | 所有應該推薦中正確推薦的比例           |
| F1-score   | Precision 與 Recall 的加權平均 |
| AUC / MAP  | 排名品質衡量，適用於排序任務           |

---

## 3. Python 套件與工具

* `Surprise`: 用於協同過濾與模型評估
* `LightFM`: 支援 hybrid、隱語義模型與冷啟動問題
* `implicit`: 適用於隱式回饋資料（如點擊）
* `TensorFlow Recommenders`: 支援深度學習推薦架構

---

## 4. 實作範例（使用 Surprise）

```python
from surprise import Dataset, SVD, accuracy
from surprise.model_selection import train_test_split

# 載入內建 movielens 資料
data = Dataset.load_builtin('ml-100k')
trainset, testset = train_test_split(data, test_size=0.25)

model = SVD()
model.fit(trainset)
predictions = model.test(testset)

print("RMSE:", accuracy.rmse(predictions))
```

---

## 5. 挑戰與議題

* **冷啟動問題（Cold Start）**：新使用者或新物品缺乏互動紀錄
* **稀疏性（Sparsity）**：使用者－物品矩陣常非常稀疏
* **多樣性與新穎性**：避免只推薦熱門內容，提升推薦多樣性
* **過擬合與隱私問題**：需平衡推薦準確性與使用者數據保護

---

推薦系統廣泛應用於電商、影音平台、社群網站與數位內容服務，是連結使用者與資訊的關鍵技術之一。
