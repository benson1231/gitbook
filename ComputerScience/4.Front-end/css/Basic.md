# Basic CSS

CSS（Cascading Style Sheets）是用來控制 HTML 元素的樣式表語言。透過 CSS，可以設定顏色、字型、排版、位置與動畫等視覺效果。

## 1. CSS 基本語法結構

```css
選擇器 {
  屬性: 值;
}
```

例如：

```css
h1 {
  color: blue;
  font-size: 24px;
}
```

## 2. 常見選擇器

### 元素選擇器

針對 HTML 標籤名稱指定樣式。

```css
p {
  line-height: 1.6;
}
```

### 類別選擇器（class）

以 `.` 開頭，選擇具有該 class 的元素。

```css
.card {
  border: 1px solid #ccc;
}
```

### ID 選擇器

以 `#` 開頭，選擇特定 id 的元素。

```css
#header {
  background-color: lightgray;
}
```

### 群組選擇器

一次設定多個元素。

```css
h1, h2, h3 {
  font-family: sans-serif;
}
```

---

## 3. 選擇器的專一性（Specificity）

當多個樣式同時套用到同一元素時，CSS 會依照「專一性」決定哪一個樣式優先。

### 專一性計算規則：

每種選擇器對應一個數值：

* 行內樣式（如 `style=""`）：1000 分
* ID 選擇器：100 分
* 類別、屬性與偽類選擇器：10 分
* 元素與偽元素選擇器：1 分

### 範例說明：

```css
h1 { color: black; }           /* specificity: 1 */
.title { color: green; }       /* specificity: 10 */
#mainTitle { color: red; }     /* specificity: 100 */
```

若同時套用，最終文字顏色為紅色，因為 ID 選擇器的專一性最高。

### 補充：專一性的實際運作原理

上述「1000 / 100 / 10 / 1」的分數是一種教學用比喻。實際上，CSS 會計算一個四位數向量 `(a, b, c, d)` 來比較不同選擇器的優先順序：

* `a`: 是否有行內樣式（有則為 1）
* `b`: ID 選擇器的數量
* `c`: 類別、屬性選擇器、偽類的數量
* `d`: 元素選擇器與偽元素的數量

這個向量的比較是逐位進行的。例如：

* `#main` → (0,1,0,0)
* `.title` → (0,0,1,0)
* `h1` → (0,0,0,1)

因此 ID 選擇器會優先於類別選擇器，而類別又優先於元素。

### 注意事項：

* `!important` 會強制覆蓋所有其他規則，但不建議常用。
* 維持選擇器簡潔清晰，可提升維護性。

---

## 4. 常見屬性與範例

### 文字樣式

```css
p {
  color: #333;
  font-size: 16px;
  font-weight: bold;
  text-align: center;
}
```

### 區塊樣式與邊框

```css
div {
  width: 300px;
  height: 200px;
  padding: 10px;
  margin: 20px auto;
  border: 1px solid black;
}
```

### 背景與顏色

```css
body {
  background-color: #f2f2f2;
}
```

### 連結樣式

```css
a {
  color: blue;
  text-decoration: none;
}

a:hover {
  text-decoration: underline;
}
```

---

## 5. 引入 CSS 的方式

### 行內樣式（Inline Style）

直接在 HTML 標籤中使用 `style` 屬性。

```html
<p style="color: red;">紅色文字</p>
```

### 內部樣式表（Internal CSS）

寫在 HTML 的 `<style>` 區塊中。

```html
<style>
  h1 { color: green; }
</style>
```

### 外部樣式表（External CSS）

將樣式寫在 `.css` 檔案中，並用 `<link>` 標籤引入。

```html
<link rel="stylesheet" href="styles.css">
```

---

CSS 是網頁設計與前端開發的基礎工具，學會它可以讓網頁更具吸引力與可用性。
