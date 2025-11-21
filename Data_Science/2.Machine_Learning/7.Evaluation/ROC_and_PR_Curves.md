# ROC and PR Curves

在分類模型的評估中，ROC 曲線與 PR 曲線是常用的可視化工具，用於衡量模型在不同閾值下的分類表現，特別適用於二元分類問題。

---

## 1. ROC 曲線（Receiver Operating Characteristic Curve）

ROC 曲線描繪 **真正率（True Positive Rate, TPR）** 與 **假陽性率（False Positive Rate, FPR）** 的變化：

* TPR = TP / (TP + FN)
* FPR = FP / (FP + TN)

### 特點：

* X 軸：FPR
* Y 軸：TPR
* 理想模型 ROC 曲線靠近左上角，表示高 TPR 與低 FPR
* **AUC（Area Under Curve）**：ROC 曲線下的面積，越接近 1 表示越好

```python
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt

fpr, tpr, _ = roc_curve(y_true, y_prob)
auc = roc_auc_score(y_true, y_prob)

plt.plot(fpr, tpr, label=f"ROC curve (AUC = {auc:.2f})")
plt.plot([0, 1], [0, 1], linestyle='--')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend()
plt.title("ROC Curve")
plt.show()
```

---

## 2. PR 曲線（Precision-Recall Curve）

PR 曲線展示 **精確率（Precision）** 與 **召回率（Recall）** 的權衡：

* Precision = TP / (TP + FP)
* Recall = TP / (TP + FN)

### 特點：

* X 軸：Recall
* Y 軸：Precision
* 對於正負樣本不均衡的資料，PR 曲線比 ROC 曲線更具區辨力

```python
from sklearn.metrics import precision_recall_curve, average_precision_score

precision, recall, _ = precision_recall_curve(y_true, y_prob)
ap = average_precision_score(y_true, y_prob)

plt.plot(recall, precision, label=f"PR curve (AP = {ap:.2f})")
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision-Recall Curve")
plt.legend()
plt.show()
```

---

## 3. 選擇建議

| 情境              | 建議使用      |
| --------------- | --------- |
| 正負樣本均衡          | ROC Curve |
| 正負樣本極度不均（如罕見事件） | PR Curve  |

---

## 4. 總結

* ROC 與 PR 曲線皆能視覺化分類模型在不同閾值下的表現
* AUC（ROC）與 AP（PR）作為整體評估指標
* 二者搭配使用能全面了解模型能力

這些指標能協助我們針對不同資料特性選擇最佳模型與適當閾值。
