# skimage

`scikit-image`（簡稱 `skimage`）是 Python 中基於 NumPy、SciPy 的影像處理函式庫，提供了豐富的工具用於圖像讀取、過濾、邊緣偵測、區域分割、形態學處理等。

---

## 一、安裝

```bash
pip install scikit-image
```

---

## 二、常用模組與功能

| 模組             | 功能概述                   |
| -------------- | ---------------------- |
| `io`           | 圖片讀取與儲存                |
| `color`        | 顏色空間轉換（RGB, HSV, Gray） |
| `filters`      | 平滑化、邊緣偵測等濾波器功能         |
| `transform`    | 縮放、旋轉、仿射變換等            |
| `feature`      | 邊緣、角點、霍夫轉換等特徵偵測        |
| `morphology`   | 區域形態學處理（膨脹、腐蝕等）        |
| `segmentation` | 區塊切割與區域標記              |
| `measure`      | 物件性質測量（面積、周長等）         |

---

## 三、影像讀取與顯示

```python
from skimage import io
import matplotlib.pyplot as plt

image = io.imread('example.jpg')
plt.imshow(image)
plt.axis('off')
plt.show()
```

---

## 四、灰階轉換與邊緣偵測

```python
from skimage.color import rgb2gray
from skimage import filters

gray = rgb2gray(image)
edge = filters.sobel(gray)

plt.imshow(edge, cmap='gray')
plt.title('Sobel 邊緣')
plt.show()
```

---

## 五、二值化與形態學操作

```python
from skimage.filters import threshold_otsu
from skimage.morphology import dilation, disk

thresh = threshold_otsu(gray)
binary = gray > thresh
morph = dilation(binary, disk(3))

plt.imshow(morph, cmap='gray')
plt.title('二值化後膨脹')
plt.show()
```

---

## 六、區域分割與標記

```python
from skimage import segmentation, measure, color
from skimage.segmentation import clear_border

labels = measure.label(morph)
colored_labels = color.label2rgb(labels, bg_label=0)

plt.imshow(colored_labels)
plt.title('區域標記')
plt.show()
```

---

## 七、應用場景

* 醫學影像處理（細胞偵測、腫瘤邊界分析）
* 電腦視覺預處理（去雜訊、特徵擷取）
* 類神經網路前處理（圖像標準化、二值化）

---

`skimage` 提供高階抽象與與 NumPy 相容的影像操作接口，是進行影像分析與前處理不可或缺的工具之一，適合初學者與研究人員快速開發與實驗。
