# RWD 響應式網頁設計教學（Responsive Web Design）

響應式網頁設計（RWD）是指網站能根據使用者裝置的螢幕尺寸，自動調整佈局與樣式，確保在手機、平板、電腦等裝置上都有良好的閱讀與操作體驗。

---

## 🔹 1. Media Queries 媒體查詢

使用媒體查詢可根據裝置寬度指定不同的 CSS 樣式。

```css
/* 手機（最大寬度 767px） */
@media (max-width: 767px) {
  body {
    background: lightyellow;
  }
}

/* 平板（768px 到 1023px） */
@media (min-width: 768px) and (max-width: 1023px) {
  body {
    background: lightblue;
  }
}

/* 桌機（1024px 以上） */
@media (min-width: 1024px) {
  body {
    background: lightgreen;
  }
}
```

---

## 🔹 2. 流動式寬度（Fluid Layout）

使用百分比（%）代替固定像素，讓區塊自動依螢幕調整。

```css
.container {
  width: 100%;
  max-width: 1200px;
  margin: 0 auto;
  padding: 0 20px;
}
```

---

## 🔹 3. 響應式圖片

```css
img {
  max-width: 100%;
  height: auto;
}
```

---

## 🔹 4. Viewport 設定（HTML）

讓瀏覽器根據裝置大小縮放畫面。

```html
<meta name="viewport" content="width=device-width, initial-scale=1.0">
```

---

## 🔹 5. Flex 與 Grid 響應式技巧

### Flex 範例：卡片寬度隨裝置變化

```css
.cards {
  display: flex;
  flex-wrap: wrap;
  gap: 16px;
}
.card {
  flex: 1 1 calc(33.33% - 16px);
}

@media (max-width: 768px) {
  .card {
    flex: 1 1 100%;
  }
}
```

### Grid 範例：自動欄數調整

```css
.grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
  gap: 16px;
}
```

---

## 🔹 6. 響應式工具與框架

* [Bootstrap](https://getbootstrap.com/)：內建 RWD 格線系統
* [Tailwind CSS](https://tailwindcss.com/)：實用類別快速響應
* [Media Query Cheatsheet](https://gist.github.com/gokulkrishh/242e68d1ee94ad05f488)

---

## ✅ 建議做法

* 所有圖片都應加上 `max-width: 100%`
* 使用相對單位（% / em / rem）替代 px
* 設定 viewport 避免畫面跑版
* 利用媒體查詢分階段調整樣式（Mobile-first 最佳）

響應式設計是現代網頁必備技能，讓使用者無論在哪種裝置上都獲得一致且良好的體驗。
