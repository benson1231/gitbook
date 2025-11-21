# CSS 切版與 Display 屬性教學

在網頁設計中，切版是指使用 HTML 與 CSS 將設計稿實作為具有正確結構與視覺的網頁。而其中最核心的屬性之一就是 `display`，它決定了元素的排版行為。

---

## 🔹 1. `display` 屬性總覽

| 值類型            | 說明                               |
| -------------- | -------------------------------- |
| `block`        | 占滿整行，高度寬度可設定。常見於 `<div>`、`<p>` 等 |
| `inline`       | 不會換行，無法設定寬高。常見於 `<span>`、`<a>`   |
| `inline-block` | 行內顯示但可設定寬高                       |
| `none`         | 元素完全隱藏（不佔空間）                     |
| `flex`         | 啟用彈性盒模型（flexbox）                 |
| `grid`         | 啟用網格排版模型                         |
| `table`        | 模擬表格行為                           |

---

## 🔹 2. 基本範例

### block vs inline

```html
<style>
  .block-box {
    display: block;
    width: 200px;
    background: lightblue;
  }
  .inline-box {
    display: inline;
    background: lightgreen;
  }
</style>

<div class="block-box">區塊元素</div>
<span class="inline-box">行內元素</span>
<span class="inline-box">行內元素2</span>
```

---

## 🔹 3. `inline-block`：行內不換行又可調大小

```html
<style>
  .inline-block-box {
    display: inline-block;
    width: 100px;
    height: 100px;
    background: pink;
    margin: 5px;
  }
</style>

<div class="inline-block-box"></div>
<div class="inline-block-box"></div>
```

---

## 🔹 4. `display: none` vs `visibility: hidden`

```css
.hidden {
  display: none;      /* 元素消失，不佔空間 */
}
.invisible {
  visibility: hidden; /* 元素隱藏，但仍佔空間 */
}
```

---

## 🔹 5. Flex 排版（常用於水平垂直置中）

```css
.flex-container {
  display: flex;
  justify-content: center; /* 主軸置中 */
  align-items: center;      /* 副軸置中 */
  height: 200px;
  background: #eee;
}
```

```html
<div class="flex-container">
  <div>置中內容</div>
</div>
```

---

## 🔹 6. Grid 排版（區域切割）

```css
.grid-container {
  display: grid;
  grid-template-columns: 1fr 2fr;
  gap: 10px;
}
```

```html
<div class="grid-container">
  <div>左側欄</div>
  <div>右側主內容</div>
</div>
```

---

## 🔹 7. 對齊方式（Alignment）

### 文字對齊（text-align）

```css
.text-center {
  text-align: center;
}
.text-left {
  text-align: left;
}
.text-right {
  text-align: right;
}
```

### 垂直置中（line-height 或 flex）

```css
.vertical-text {
  height: 100px;
  line-height: 100px;
  text-align: center;
}
```

或使用 Flex：

```css
.center-box {
  display: flex;
  justify-content: center;
  align-items: center;
  height: 200px;
}
```

### Grid 對齊

```css
display: grid;
place-items: center; /* 水平 + 垂直置中 */
```

---

## 📌 補充：常見排版組合

* `display: flex` + `gap`：用於橫向按鈕排列
* `display: inline-block`：可用於圖片或卡片等區塊的並排呈現
* `display: none`：用於隱藏 modal、下拉選單等

切版的核心不僅是視覺，也包含結構與語意，合理運用 `display` 搭配 `position`、`margin`、`padding` 可創造穩定的響應式排版。
