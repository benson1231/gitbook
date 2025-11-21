# HTML標籤

HTML（超文件標記語言，HyperText Markup Language）是用來建立網頁的標記語言。以下是一些常見的基本 HTML 標籤介紹：

## 1. `<!DOCTYPE html>`

宣告這是一個 HTML5 文件。

```html
<!DOCTYPE html>
```

## 2. `<html>`

HTML 文件的根元素，包裹整個 HTML 頁面。

```html
<html>
  ...
</html>
```

## 3. `<head>` 與 `<body>`

`<head>` 包含網站的中繼資料，例如標題與 CSS；`<body>` 包含頁面主要的內容。

```html
<head>
  <title>我的第一個網頁</title>
</head>
<body>
  <h1>歡迎來到我的網站</h1>
</body>
```

## 4. `<title>`

定義瀏覽器標籤上顯示的標題。

```html
<title>首頁</title>
```

## 5. 標題標籤 `<h1>` \~ `<h6>`

用來定義標題層級，`<h1>` 為最大標題，`<h6>` 為最小標題。

```html
<h1>主標題</h1>
<h2>副標題</h2>
```

## 6. 段落 `<p>`

用來定義一段文字內容。

```html
<p>這是一段文字。</p>
```

## 7. 超連結 `<a>`

建立超連結，可連至外部網站或頁面內指定 `id` 位置。

### 外部連結範例：

```html
<a href="https://example.com">前往範例網站</a>
```

### 錨點連結範例（跳至頁面內特定元素）：

```html
<a href="#section1">跳到 Section 1</a>
...
<h2 id="section1">Section 1</h2>
```

## 8. 圖片 `<img>`

插入圖片，需指定 `src`（來源）與 `alt`（替代文字）。

```html
<img src="image.jpg" alt="示意圖">
```

## 9. 清單 `<ul>`、`<ol>`、`<li>`

* `<ul>`：無序清單
* `<ol>`：有序清單
* `<li>`：清單項目

```html
<ul>
  <li>蘋果</li>
  <li>香蕉</li>
</ul>
<ol>
  <li>第一步</li>
  <li>第二步</li>
</ol>
```

## 10. `<div>` 與 `<span>`

* `<div>`：區塊元素，用於布局或分組內容。
* `<span>`：行內元素，用於小範圍文字的樣式控制。

```html
<div class="container">
  <span style="color: red;">紅色文字</span>
</div>
```

## 11. 強調文字 `<em>` 與 `<strong>`

* `<em>`：表示語意上的強調，通常會以斜體顯示。
* `<strong>`：表示重要性，通常會以粗體顯示。

```html
<p>這是<em>非常</em>重要的提示。</p>
<p>請<strong>立刻</strong>回覆！</p>
```

## 12. 表格 `<table>`、`<tr>`、`<th>`、`<td>`

用來呈現結構化資料。

```html
<table border="1">
  <thead>
    <tr>
      <th>姓名</th>
      <th>年齡</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>小明</td>
      <td>18</td>
    </tr>
    <tr>
      <td>小華</td>
      <td>20</td>
    </tr>
  </tbody>
</table>
```

---

## 13. 語意化標籤（Semantic Elements）

語意化標籤使 HTML 結構更清晰，便於維護與搜尋引擎理解。

### 常見語意標籤：

* `<header>`：網站或區塊的標頭，通常包含 logo、標題或導覽列。
* `<nav>`：導覽列，內含網站主要連結或選單。
* `<main>`：頁面主體內容，整體唯一。
* `<footer>`：頁尾資訊，例如聯絡方式、版權聲明。
* `<section>`：具主題性的區塊，如章節、功能群組。
* `<article>`：獨立內容單位，例如部落格文章、留言、報導。
* `<aside>`：與主內容間接相關的補充資訊，如側欄、提示區塊。
* `<figure>`：包含圖片、圖表、影片等媒體。
* `<figcaption>`：用於描述 `<figure>` 的文字說明。
* `<video>` / `<audio>` / `<embed>`：用於嵌入影音或其他媒體檔案。

### 結構範例：

```html
<header>
  <h1>我的網站</h1>
  <nav>
    <a href="#home">首頁</a>
    <a href="#about">關於我們</a>
  </nav>
</header>

<main>
  <section>
    <article>
      <h2>最新文章</h2>
      <p>這是一篇文章的簡介內容。</p>
    </article>
  </section>
  <aside>
    <p>相關連結與推薦內容。</p>
  </aside>
</main>

<footer>
  <p>&copy; 2025 MySite 保留所有權利。</p>
</footer>
```

這些語意標籤有助於建立結構明確、可讀性佳的 HTML 文件。

## 14. 導覽錨點範例

以下是一個實際應用這些標籤與屬性的範例清單，可作為頁面內導覽：

```html
<ol>
  <li><a href="#top">Top</a></li>
  <li><a href="#bottom">Bottom</a></li>
</ol>
```

---

掌握這些 HTML 標籤是建立現代網頁的第一步，搭配 CSS 與 JavaScript 可以創造互動性與視覺豐富的網站。
