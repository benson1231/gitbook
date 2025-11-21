# Theoretical Concepts

本章節彙整機器學習中關鍵的理論基礎，有助於理解模型訓練行為、預測誤差來源與正則化效果。這些概念雖不直接產出模型，但能提升分析、選模與除錯的深度。

---

## 1. Bias-Variance Tradeoff
- 說明模型在「學得太少 vs 學得太多」間的取捨
- 高偏差：模型過簡 → 欠擬合
- 高變異：模型太複雜 → 過擬合

---

## 2. Overfitting vs Underfitting
- 過擬合：模型對訓練集記憶過深，泛化能力差
- 欠擬合：模型學不到有效規則，訓練誤差也高
- 可透過交叉驗證、正則化、Early Stopping 改善

---

## 3. Parametric vs Nonparametric Models
| 類型             | 定義                                      | 範例                   |
|------------------|-------------------------------------------|------------------------|
| Parametric       | 模型結構固定，學習有限參數               | Linear Regression, SVM |
| Nonparametric    | 模型彈性高，無固定參數量                  | KNN, Decision Tree     |

---

## 4. Curse of Dimensionality
- 資料維度越高，資料越稀疏、距離度量失效、模型訓練變困難
- 常見對策：PCA/UMAP 降維、特徵選擇、正則化

---

## 5. Distance Metrics
- 用於衡量資料點間的相似度（常見於 KNN、Clustering）
- Euclidean（歐式距離）
- Manhattan（曼哈頓距離）
- Cosine Similarity（餘弦相似度）
- Mahalanobis（適合有相關性的特徵）

---

## 6. Loss Functions（損失函數）
- 衡量預測結果與實際標籤的差距，是模型訓練的目標
- 常見損失函數：
  - MSE, MAE：回歸任務
  - Cross-Entropy：分類任務
  - Hinge Loss：SVM 使用

---

這些理論章節可以搭配實驗或視覺化觀察加深理解，在模型評估與調參時特別有用。
