# Supervised Learning

Supervised Learning 是機器學習中最常見的一種學習類型，其核心特點是「有標籤資料」。模型透過已知的輸入（特徵）與對應的輸出（標籤）來學習預測規則，進而應用在未知資料上。

---

## 1. 主要任務類型

### A. 分類（Classification）
- 輸出為離散類別（例如：是/否、A/B/C）
- 常見應用：垃圾郵件判斷、疾病診斷、圖像辨識

### B. 回歸（Regression）
- 輸出為連續數值（例如：價格、溫度、機率）
- 常見應用：房價預測、銷售量預測、風險評估

---

## 2. 常見演算法

### 分類模型：
- Logistic Regression
- Decision Tree
- K-Nearest Neighbors (KNN)
- Support Vector Machine (SVM)
- Naive Bayes
- Random Forest
- Gradient Boosting (XGBoost, LightGBM)

### 回歸模型：
- Linear Regression（Simple/Multiple）
- Ridge / Lasso / ElasticNet
- Decision Tree Regression
- SVR（Support Vector Regression）
- Ensemble 回歸模型

---

## 3. 建模流程
1. 數據前處理與特徵工程
2. 選擇適合任務的模型類型
3. 拆分資料（train/validation/test）
4. 模型訓練與超參數調整（GridSearch / CV）
5. 模型評估（準確率、MSE、ROC-AUC 等）
6. 預測與部署

---

## 4. 評估指標（依任務不同）
| 任務類型   | 常見指標                        |
|------------|----------------------------------|
| 分類       | Accuracy, Precision, Recall, F1, ROC-AUC |
| 回歸       | MAE, MSE, RMSE, R²               |

---

## 5. 工具與套件
- `scikit-learn`: 主力工具庫，提供完整模型與 API
- `xgboost`, `lightgbm`: 高效能集成演算法
- `pandas`, `numpy`: 特徵與資料處理
- `matplotlib`, `seaborn`: 可視化結果

---

## 6. 延伸子章節
- `Classification/Logistic Regression.md`
- `Classification/KNN.md`
- `Regression/Linear Regression.md`
- `Regression/Lasso.md`
- `Ensembling Methods/Random Forest.md`

---

Supervised Learning 是所有 AI 模型應用的基礎，從金融到醫療、行銷到推薦系統都離不開這類方法。建議從簡單模型著手，逐步深入理解其假設與限制。
