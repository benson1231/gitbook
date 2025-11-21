# Model Performance

機器學習模型在分類任務中，需要用多種指標來評估效能，尤其在資料不平衡（imbalanced dataset）時，單看準確率（Accuracy）可能會誤導。

---

## 🔢 混淆矩陣（Confusion Matrix）

|                  | 實際為 Positive | 實際為 Negative |
|------------------|----------------|----------------|
| 預測為 Positive  | TP（True Positive） | FP（False Positive） |
| 預測為 Negative  | FN（False Negative） | TN（True Negative） |

---

## ✅ 評估指標說明

| 指標         | 計算公式                             | 解釋 |
|--------------|--------------------------------------|------|
| Accuracy     | (TP + TN) / (TP + TN + FP + FN)      | 整體預測正確的比例 |
| Precision    | TP / (TP + FP)                       | 預測為 Positive 中，有多少是真的 |
| Recall（Sensitivity） | TP / (TP + FN)              | 實際為 Positive 中，被正確預測的比例 |
| F1-score     | 2 × (Precision × Recall) / (Precision + Recall) | Precision 與 Recall 的調和平均數，當資料不平衡時特別重要 |

---

## 🧠 使用情境

- **Accuracy 適合用於**：資料平衡的分類問題
- **Precision 適合用於**：錯誤預測 Positive 成本高（如垃圾郵件過濾）
- **Recall 適合用於**：漏掉 Positive 成本高（如癌症偵測）
- **F1-score 適合用於**：需要兼顧 Precision 和 Recall 的場景

---

## 📌 Python 計算 F1-score 範例

```python
from sklearn.metrics import precision_score, recall_score, f1_score

y_true = [1, 0, 1, 1, 0, 1, 0]
y_pred = [1, 0, 1, 0, 0, 1, 1]

precision = precision_score(y_true, y_pred)
recall = recall_score(y_true, y_pred)
f1 = f1_score(y_true, y_pred)

print(f"Precision: {precision:.2f}")
print(f"Recall: {recall:.2f}")
print(f"F1-score: {f1:.2f}")
