# SHAP

SHAP 是一種模型解釋方法，基於合作博弈理論中的 Shapley Value 概念，能精確地量化每個特徵對機器學習模型預測結果的貢獻。

---

## 一、SHAP 原理簡述

* **理論基礎**：來自博弈論的 Shapley Value，每個特徵視為貢獻預測的玩家。
* **主要目標**：給定一筆資料，計算每個特徵對預測的邊際貢獻。
* **優點**：公平、有理論保證、模型不可知（model-agnostic）、具備一致性與局部準確性。

---

## 二、SHAP 的應用場景

* 模型審核與可解釋性
* 建構可信任的 AI 系統（如醫療、金融）
* 特徵選擇與重要性排序
* Debug 黑箱模型預測錯誤

---

## 三、使用範例（以 XGBoost 為例）

```python
import shap
import xgboost
from sklearn.datasets import load_boston
from sklearn.model_selection import train_test_split

# 載入資料與模型訓練
X, y = load_boston(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y)
model = xgboost.XGBRegressor().fit(X_train, y_train)

# 建立解釋器與取得 SHAP 值
explainer = shap.Explainer(model, X_train)
shap_values = explainer(X_test)

# 視覺化
shap.summary_plot(shap_values, X_test)
```

---

## 四、SHAP 視覺化工具

| 圖形類型              | 說明                      |
| ----------------- | ----------------------- |
| `summary_plot`    | 顯示每個特徵的整體重要性（顏色代表特徵值高低） |
| `force_plot`      | 單筆預測的拆解貢獻圖              |
| `dependence_plot` | 特徵間交互關係視覺化              |
| `bar_plot`        | 特徵重要性長條圖                |
| `heatmap`         | 多筆樣本的特徵貢獻熱圖             |

---

## 五、SHAP 支援的模型類型

| 類型            | 支援方式                 |
| ------------- | -------------------- |
| Tree Models   | 使用 `TreeExplainer`   |
| Linear Models | 使用 `LinearExplainer` |
| Deep Learning | 使用 `DeepExplainer`   |
| 任意模型          | 使用 `KernelExplainer` |

---

## 六、進階提示

* 可以搭配 `shap.plots.waterfall()` 分析個別預測
* 可與 pandas 結合選出對特定預測影響最大的特徵
* `SHAP` 支援 `pandas.DataFrame` 與 `numpy.array` 資料結構

---

SHAP 是當今解釋機器學習模型預測最具理論與實務價值的工具之一，無論是用於研究還是產業應用，都能提供強大且清晰的特徵解釋力。
