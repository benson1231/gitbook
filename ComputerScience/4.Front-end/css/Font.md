# Font CSS

字體（font）是網頁設計中呈現品牌風格與提高可讀性的重要元素。以下為 CSS 中與文字外觀、排版、樣式控制相關的屬性與觀念整理：

---

### 1. `font-family`：字體家族

定義文字的字型。可依序指定多個字體，當第一個無法顯示時，會使用備選字體。

```css
font-family: "Helvetica", Arial, sans-serif;
```

### 2. `font-size`：字體大小

設定文字的大小。常見單位有：

* `px`：像素，固定大小
* `em`：相對於父元素的字體大小
* `rem`：相對於根元素的字體大小
* `%`：相對於父元素的百分比

```css
font-size: 1.2rem;
```

### 3. `font-weight`：字體粗細

設定文字是否為粗體。
常見值：`normal`、`bold`、`lighter`、`bolder` 或數值（`100` 到 `900`，以 `400` 為常態）。

```css
font-weight: 700;
```

### 4. `font-style`：字體樣式

設定文字是否為斜體。

* `normal`（正常）
* `italic`（斜體）
* `oblique`（傾斜）

```css
font-style: italic;
```

### 5. `line-height`：行距

控制每行文字的垂直間距，有助於可讀性。

```css
line-height: 1.5;
```

### 6. `letter-spacing`：字母間距

設定字母之間的水平間距。

```css
letter-spacing: 0.05em;
```

### 7. `word-spacing`：單字間距

設定單字之間的距離。

```css
word-spacing: 4px;
```

### 8. `text-align`：文字對齊

設定文字的水平對齊方式。

* `left`：靠左
* `center`：置中
* `right`：靠右
* `justify`：左右對齊

```css
text-align: justify;
```

### 9. `color` 與 `background-color`

* `color`：文字顏色
* `background-color`：背景顏色

```css
color: #222;
background-color: #f5f5f5;
```

### 10. `opacity`：透明度

設定整體元素的透明度。

```css
opacity: 0.8;
```

---

## 字體分類與備援

* **Serif（襯線字體）**：字母末端有裝飾線條，如 Times New Roman。
* **Sans-serif（非襯線字體）**：字母簡潔無襯線，如 Arial、Helvetica。
* 建議設定多個備援字體，確保不同裝置都能正常顯示。

---

## 外部字體載入

### Google Fonts 使用方法：

```html
<link href="https://fonts.googleapis.com/css2?family=Roboto&display=swap" rel="stylesheet">
```

```css
font-family: 'Roboto', sans-serif;
```

### `@font-face` 自定義字型

可使用本機或遠端字體來源。

```css
@font-face {
  font-family: "MyFont";
  src: url("myfont.woff2") format("woff2");
}
```

---

這些屬性與技術組合起來可打造具備專業外觀、品牌風格與高可讀性的文字設計。
