# Unsupervised Learning

Unsupervised Learning 是在**沒有標籤資料**的情況下，探索資料內在結構的一種學習方法。其目的是發現樣本間的潛在模式、群組或表示方式，常應用於探索性資料分析與特徵壓縮。

---

## 1. 主要任務類型

### A. 聚類（Clustering）
- 將資料自動劃分為數個群集
- 常見應用：顧客分群、影像分類前處理、生物分型

### B. 降維（Dimensionality Reduction）
- 將高維資料壓縮為低維度空間，保留主要資訊
- 常見應用：視覺化、特徵預處理、雜訊過濾

### C. 關聯規則學習（Association Rule Learning）
- 發掘資料中的共現關係（如購物籃分析）
- 常見應用：推薦系統、關聯性行銷

---

## 2. 常見演算法

### 聚類
- K-Means
- Hierarchical Clustering（階層式）
- DBSCAN（密度為基礎）
- Gaussian Mixture Model (GMM)

### 降維
- PCA（主成分分析）
- t-SNE
- UMAP
- Autoencoder（神經網路型態）

### 關聯規則
- Apriori
- FP-Growth

---

## 3. 使用場景
- 沒有標籤資料可供訓練時（如新資料探索）
- 特徵太多，需先降維再餵入監督模型
- 需進行資料壓縮或視覺化呈現

---

## 4. 聚類模型評估指標
由於沒有「真實標籤」，常用以下內部評估指標：
- Silhouette Score
- Davies-Bouldin Index
- Calinski-Harabasz Index
- 可視化：降維後觀察分群分布

---

## 5. 工具與套件
- `scikit-learn`: 提供大部分非監督模型與評估工具
- `umap-learn`, `tsne`, `hdbscan`: 進階降維與分群工具
- `mlxtend`: 關聯分析（Apriori, FP-Growth）

---

## 6. 延伸子章節
- `Clustering/KMeans.md`
- `Clustering/DBSCAN.md`
- `Dimensionality Reduction/PCA.md`
- `Dimensionality Reduction/tSNE.md`

---

Unsupervised Learning 是資料分析的重要起點，適合在對資料未知或無標註的情況下挖掘潛在規律，也常作為監督式學習的前置處理步驟。
