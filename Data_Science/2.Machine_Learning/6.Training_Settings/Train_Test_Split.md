# Train-Test Split

在進行機器學習建模前，必須將資料集切分為訓練集與驗證集，目的是評估模型的泛化能力，防止過度擬合（overfitting）。常見的切分方式包括訓練集（training set）、驗證集（validation set）與測試集（test set）。

---

## 1. 為什麼要分訓練與驗證集？

* **訓練集**：用於模型的學習與參數更新。
* **驗證集**：用於調整超參數與評估訓練過程中的效能，幫助選擇最佳模型。
* **測試集（可選）**：在模型訓練完成後，用於最終評估其效能，模擬真實部署情境。

---

## 2. 分割方式與比例

常見分割比例為：

* 80% 訓練集 / 20% 驗證集
* 70% 訓練集 / 15% 驗證集 / 15% 測試集

### Python 實作（使用 Scikit-Learn）

```python
from sklearn.model_selection import train_test_split

X = data.drop("target", axis=1)
y = data["target"]

# 訓練集與驗證集分割
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)
```

---

## 3. 注意事項

* **隨機性**：使用 `random_state` 可確保實驗可重現。
* **資料平衡**：分類任務中，應考慮使用 `stratify=y` 以保持各類別比例一致。
* **時序資料**：應使用時間順序分割，避免資訊洩漏（data leakage）。

---

## 4. 延伸：交叉驗證（Cross-Validation）

為了更穩健地評估模型，可以採用 k-fold 交叉驗證，例如將資料分為 5 等份，輪流作為驗證集，其餘為訓練集，最終取平均效能。

```python
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier

model = RandomForestClassifier()
scores = cross_val_score(model, X, y, cv=5)
print("平均準確率:", scores.mean())
```

交叉驗證能有效減少單一切分結果的偶然性，是常見的模型驗證方法。

---

良好的資料切分策略是建立穩健模型的基石，能提升模型在未見資料上的表現。
