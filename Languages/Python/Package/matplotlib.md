# Matplotlib

`matplotlib` 是 Python 中最常用的資料視覺化套件之一，支援靜態、互動式與動畫圖表的產生。其核心模組為 `pyplot`，用法類似於 MATLAB。

---

## 一、安裝

```bash
pip install matplotlib
```

---

## 二、基本用法（折線圖）

```python
import matplotlib.pyplot as plt

x = [1, 2, 3, 4, 5]
y = [1, 4, 9, 16, 25]

plt.plot(x, y)
plt.title("平方數")
plt.xlabel("x 軸")
plt.ylabel("y 軸")
plt.show()
```

---

## 三、常見圖表類型

| 圖表類型    | 函數          | 說明           |
| ------- | ----------- | ------------ |
| 折線圖     | `plot()`    | 數據趨勢顯示       |
| 長條圖     | `bar()`     | 分類資料比較       |
| 散點圖     | `scatter()` | 數值分布與相關性     |
| 直方圖     | `hist()`    | 數值分布密度       |
| 圓餅圖     | `pie()`     | 占比視覺化        |
| 影像圖（熱圖） | `imshow()`  | 二維影像或矩陣資料視覺化 |

---

## 四、常用設定與功能

```python
plt.figure(figsize=(8, 5))        # 設定畫布大小
plt.grid(True)                    # 顯示格線
plt.legend(['實驗組'])            # 加入圖例
plt.savefig('output.png')         # 儲存圖片
```

---

## 五、子圖繪製（subplot）

```python
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
axs[0].plot([1, 2, 3], [1, 2, 3])
axs[1].bar(['A', 'B'], [3, 5])
plt.tight_layout()
plt.show()
```

---

## 六、中文顯示（Linux 範例）

```python
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
my_font = fm.FontProperties(fname=font_path)
plt.title("中文標題", fontproperties=my_font)
```

---

## 七、進階功能

* 動畫（`matplotlib.animation`）
* 三維繪圖（`mpl_toolkits.mplot3d`）
* 與 pandas 或 seaborn 搭配使用

---

Matplotlib 是資料科學與機器學習工作中不可或缺的圖表工具，從快速繪圖到進階自定義視覺化都能滿足需求。熟練它能提升資料探索與報告表達力。
