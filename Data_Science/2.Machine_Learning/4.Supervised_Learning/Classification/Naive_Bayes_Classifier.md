# Naive Bayes Classifier

Naive Bayes 是一種基於機率的監督式分類演算法，假設特徵間條件獨立。它特別適用於文本分類（如垃圾郵件偵測）等高維度資料情境。

---

### 1. 理論基礎

根據貝葉斯定理：

$$
P(y\mid x_1, ..., x_n) = \frac{P(y) P(x_1, ..., x_n \mid y)}{P(x_1, ..., x_n)}
$$

在 Naive Bayes 中，假設特徵 \$x\_i\$ 在給定 \$y\$ 下彼此條件獨立，因此：

$$
P(x_1, ..., x_n \mid y) = \prod_{i=1}^{n} P(x_i \mid y)
$$

模型實際預測為最大後驗機率：

$$
\hat{y} = \arg\max_y P(y) \prod_{i=1}^{n} P(x_i \mid y)
$$

---

### 2. 常見 Naive Bayes 類型

| 類型            | 特徵類型       | 備註             |
| ------------- | ---------- | -------------- |
| GaussianNB    | 連續特徵       | 假設特徵服從高斯分布     |
| MultinomialNB | 計數型特徵（如詞頻） | 常用於文本分類        |
| BernoulliNB   | 二元特徵（0/1）  | 特別適合處理是否出現類型特徵 |

---

### 3. 優缺點

| 優點            | 缺點               |
| ------------- | ---------------- |
| 實作簡單、計算效率高    | 特徵間獨立假設在現實中不一定成立 |
| 對高維資料仍表現良好    | 對資料中出現 0 機率敏感    |
| 在小樣本下仍具良好泛化能力 | 無法學習特徵間的交互作用     |

---

### 4. Python 實作（以 MultinomialNB 為例）

```python
from sklearn.naive_bayes import MultinomialNB
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

texts = ["buy cheap now", "limited offer", "hello friend", "win big prize"]
labels = [1, 1, 0, 1]  # 1: spam, 0: ham

vec = CountVectorizer()
X = vec.fit_transform(texts)
X_train, X_test, y_train, y_test = train_test_split(X, labels, test_size=0.5, random_state=42)

model = MultinomialNB()
model.fit(X_train, y_train)

pred = model.predict(X_test)
print("準確率：", accuracy_score(y_test, pred))
```

---

### 5. 應用範疇

* 垃圾郵件偵測、情感分析
* 文件分類、文字探勘
* 疾病預測、風險預測

---

Naive Bayes 雖模型簡單，但在許多實務應用中效果驚人。特別適合做為基準模型（baseline）或用於高維稀疏特徵的情境。
