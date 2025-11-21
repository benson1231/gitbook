# Evaluation

模型評估是機器學習流程中不可或缺的一環，用於判斷模型的表現是否符合預期。選擇正確的評估指標與方法，有助於比較不同模型、避免過擬合，並做出更合理的模型選擇。

---

## 1. 根據任務選擇評估指標

### 分類（Classification）
| 指標             | 說明                                 |
|------------------|--------------------------------------|
| Accuracy         | 整體正確率                          |
| Precision        | 真陽性占預測為正的比例（避免假陽性） |
| Recall (Sensitivity) | 真陽性占實際為正的比例（避免假陰性） |
| F1-score         | Precision 與 Recall 的調和平均      |
| ROC-AUC          | 衡量模型在各閾值下的分類能力         |
| Confusion Matrix | 真陽性/真陰性/假陽性/假陰性的統計    |

### 回歸（Regression）
| 指標     | 說明                                |
|----------|-------------------------------------|
| MAE      | 平均絕對誤差（Mean Absolute Error）   |
| MSE      | 平均平方誤差（Mean Squared Error）   |
| RMSE     | 均方根誤差（Root Mean Squared Error） |
| R² Score | 解釋變異數比例（1 為完美預測）         |

---

## 2. 交叉驗證評估（Cross-Validation）
- 避免因為單次切分資料導致評估不穩定
- 常見方法：K-Fold、Stratified K-Fold

```python
from sklearn.model_selection import cross_val_score
scores = cross_val_score(model, X, y, cv=5, scoring='accuracy')
print("平均正確率：", scores.mean())
```

---

## 3. 視覺化工具
### 分類：
- Confusion Matrix（熱圖）
- ROC 曲線 / AUC 區域
- Precision-Recall 曲線

### 回歸：
- 預測值 vs 實際值 散佈圖
- 殘差圖（residual plot）

---

## 4. 避免過擬合的評估方式
- **Train/Validation/Test 分離評估**：三方資料互不干擾
- **使用 Early Stopping 時以驗證集為準**
- **留意評估指標差距**：訓練 vs 驗證分數若差距大 → 過擬合可能

---

## 5. 模型比較技巧
- 使用統一交叉驗證流程
- 固定 random seed 確保公平比較
- 收集多個評估指標 → 做全方位比較

---

透過科學的評估方法，你才能對模型有合理的信心，也能找到最適合你資料與任務的演算法。