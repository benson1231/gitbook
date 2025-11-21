# HTML屬性

HTML 屬性（attributes）用來補充標籤的資訊，控制其外觀、行為與識別方式。以下是常見的 HTML 屬性整理：

## 通用屬性

### `id`

為元素指定唯一識別碼，可用於 JavaScript 操作或 CSS 樣式設定。

```html
<div id="main">主要區塊</div>
```

### `class`

定義元素所屬類別，便於群組化控制樣式。

```html
<p class="note">提示文字</p>
```

### `style`

直接在元素上定義內嵌樣式（不建議大量使用）。

```html
<h1 style="color: green;">綠色標題</h1>
```

### `title`

提供額外資訊，滑鼠懸停時會顯示提示文字。

```html
<img src="icon.png" title="網站圖示">
```

---

## 連結與資源相關屬性

### `href`

用於 `<a>` 標籤，指定超連結的目標網址。

```html
<a href="https://example.com">點我</a>
```
錨點連結範例（跳至頁面內特定元素）：

```html
<a href="#section1">跳到 Section 1</a>
...
<h2 id="section1">Section 1</h2>
```

### `target`

設定連結開啟方式：`_self`（同頁）、`_blank`（新頁）。

```html
<a href="page.html" target="_blank">新分頁開啟</a>
```

### `src`

指定外部資源路徑，常見於 `<img>`、`<script>`、`<iframe>`。

```html
<img src="image.jpg">
```

### `alt`

圖片無法載入時顯示的替代文字，對於無障礙與 SEO 很重要。

```html
<img src="logo.png" alt="公司標誌">
```

---

## 表單相關屬性

### `type`

定義 `<input>` 的資料類型，如文字、密碼、電子郵件等。

```html
<input type="email">
```

### `value`

設定預設值。

```html
<input type="text" value="預設文字">
```

### `placeholder`

提示輸入內容的淡灰文字。

```html
<input type="text" placeholder="請輸入姓名">
```

### `disabled`

停用輸入元件。

```html
<input type="submit" value="送出" disabled>
```

### `checked`

指定預設被勾選的選項（適用於 radio/checkbox）。

```html
<input type="checkbox" checked>
```

### `readonly`

內容僅可讀不可修改。

```html
<input type="text" value="唯讀" readonly>
```

---

這些屬性是學習與撰寫 HTML 頁面不可或缺的基礎知識，熟悉它們能讓你更有效地控制標籤行為與樣式。
