# Feature Engineering

Feature Engineering 是機器學習中極為關鍵的一步，它指的是從原始資料中建構出有意義的特徵，來提升模型表現與泛化能力。優秀的特徵工程往往比更複雜的模型更能帶來明顯改善。

---

## 1. 主要目標
- 改善模型輸入品質
- 發掘隱藏的資料規律與模式
- 提高模型準確度與穩定性
- 降低維度或冗餘度，減少過擬合風險

---

## 2. 標準流程
1. **資料理解**：瞭解欄位意義、型別、資料分布與缺失情況
2. **特徵清理**：處理缺值、極端值、重複欄位與無效資訊
3. **特徵轉換**：包含標準化、正規化、類別編碼等
4. **特徵建構**：根據領域知識進行數學運算、交互項建立、新欄位衍生
5. **特徵選擇**：移除冗餘或無關欄位，提升模型效率
6. **特徵評估**：透過交叉驗證、重要性分析確認特徵價值

---

## 3. 常見方法
### 數值特徵處理
- Min-Max Scaling
- Standardization (Z-score)
- Binning（分箱）
- Log/Root 轉換（對偏態數據）

### 類別特徵處理
- Label Encoding
- One-Hot Encoding
- Frequency / Target Encoding
- Embedding（進階：NLP、推薦系統）

### 特徵建構
- 組合欄位（如價格/數量 = 單價）
- 時間特徵分解（年/月/星期幾）
- 群組統計（groupby mean/std）
- 多項式特徵（PolynomialFeatures）

### 特徵選擇（與降維相關）
- Variance Threshold
- Correlation Filtering
- L1 Regularization（Lasso）
- Recursive Feature Elimination（RFE）

---

## 4. 工具套件
| 類別       | 工具                      |
|------------|---------------------------|
| 資料處理   | `pandas`, `numpy`         |
| 編碼轉換   | `sklearn.preprocessing`   |
| 特徵建構   | `PolynomialFeatures`, 自定義函數 |
| 特徵選擇   | `sklearn.feature_selection` |
| 自動化     | `Feature-engine`, `tsfresh`（時間序列） |

---

## 5. 注意事項
- 特徵處理只能用「訓練集資訊」產生轉換器（防止資料洩漏）
- 記得搭配 `Pipeline` 封裝流程，避免錯誤
- 盡可能保留業務意義與可解釋性

---

## 6. 延伸主題（可擴充）
- Feature Importance 解釋（如 SHAP）
- 自動特徵工程（AutoML）
- 時間序列特徵工程

---

本資料夾適用於各種資料型態（表格、時間序列、類別資料），每個主題可進一步拆分成範例實作檔案。