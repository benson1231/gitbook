# Machine Learning

機器學習是一種讓電腦從資料中自動學習模式，並進行預測或決策的人工智慧技術。它在生物資訊、金融、醫療、製造與自然語言處理等領域皆有廣泛應用。

---

## 一、基本概念

* **資料集（Dataset）**：包含特徵（features）與標籤（labels）的資料樣本，用於訓練與測試模型。
* **模型（Model）**：從資料中學習規則的數學函數。
* **訓練（Training）**：讓模型從資料中學習。
* **預測（Prediction）**：使用訓練好的模型對未知資料進行推論。
* **評估（Evaluation）**：測量模型在新資料上的準確性與泛化能力。

---

## 二、學習類型

### 1. 監督式學習（Supervised Learning）

* 輸入與對應輸出已知。
* 常見任務：分類（Classification）、回歸（Regression）
* 範例：乳癌良惡性預測、股票價格預測

### 2. 非監督式學習（Unsupervised Learning）

* 輸出未知，探索資料內在結構。
* 常見任務：分群（Clustering）、降維（Dimensionality Reduction）
* 範例：族群基因變異分類、細胞亞群辨識

### 3. 強化學習（Reinforcement Learning）

* 通過與環境互動獲得回饋，學習最優策略。
* 範例：機器手臂控制、圍棋 AI、資源配置最佳化

---

## 三、常見演算法

| 類型 | 演算法範例                                  |
| -- | -------------------------------------- |
| 分類 | 決策樹、SVM、KNN、隨機森林、Logistic Regression   |
| 回歸 | 線性回歸、Lasso、Ridge Regression            |
| 分群 | K-means、階層式分群（Hierarchical Clustering） |
| 降維 | PCA、t-SNE、UMAP                         |

---

## 四、常見流程

1. 資料前處理（清洗、標準化、特徵工程）
2. 分割資料集（訓練集 / 驗證集 / 測試集）
3. 選擇模型與訓練
4. 模型調參（Hyperparameter Tuning）
5. 評估模型效能（Accuracy, AUC, F1-score）
6. 上線應用或進一步分析（Feature Importance, SHAP 等）

---

## 五、應用案例（生醫領域）

* 利用 Logistic Regression 預測疾病風險（如糖尿病、乳癌）
* 以 Random Forest 判別基因型與表現型之間的關係
* 使用深度學習模型分析病理影像或單細胞 RNA-seq
* 建立 Polygenic Risk Score（PRS）模型整合 GWAS 結果

---

機器學習是資料驅動科學的核心，結合統計學、演算法與實際資料，能提供強大的預測能力與資料洞察。在生物資訊與精準醫療時代，學習掌握機器學習方法是不可或缺的技能。
