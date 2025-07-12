# Vite 介紹與使用教學

[Vite](https://vitejs.dev/) 是由 Evan You（Vue 的創作者）開發的現代前端建構工具，主打「極速啟動、模組原生支援、開發體驗極佳」，已成為取代 Webpack 的熱門選擇。

## Vite 的特點

* 🚀 **超快啟動速度**：利用原生 ES Module 載入，省去繁重打包步驟
* 🔥 **熱更新超快（HMR）**：變更檔案時即時回應
* 📦 **支援多框架**：Vue、React、Svelte、Lit 等皆可整合
* 🔧 **內建支援 TypeScript、PostCSS、CSS Modules 等**
* 🧩 **簡潔的設定檔**：使用 `vite.config.js`

## 安裝與初始化

### 使用官方指令建立專案

```bash
npm create vite@latest
```

接著依照提示選擇：

* 專案名稱
* 所使用的框架（如 Vue、React、Svelte）
* 是否使用 TypeScript

### 安裝依賴與啟動開發伺服器

```bash
cd my-vite-app
npm install
npm run dev
```

## 專案結構簡介

```
my-vite-app/
├── index.html
├── vite.config.js
├── src/
│   └── main.js
```

* `index.html` 是入口檔案，直接使用原生 HTML 管理
* `src/` 為主要開發目錄

## 建立正式環境打包（Build）

```bash
npm run build
```

輸出結果會在 `dist/` 資料夾中，可直接部署。

## 自訂 Vite 設定檔

```js
// vite.config.js
import { defineConfig } from 'vite';
import vue from '@vitejs/plugin-vue';

export default defineConfig({
  plugins: [vue()],
  server: {
    port: 3000,
    open: true
  }
});
```

## Vite 與 Webpack 比較

| 項目     | Vite                       | Webpack      |
| ------ | -------------------------- | ------------ |
| 啟動速度   | 極快（原生模組）                   | 慢（需整體打包）     |
| 熱更新    | 快速                         | 偏慢           |
| 預設支援框架 | Vue/React/TS/CSS Modules 等 | 需自行安裝 loader |
| 建置階段   | 使用 Rollup 打包               | 使用自己的打包流程    |
| 設定檔    | 簡單                         | 較繁雜          |

---

Vite 適合用於各種現代前端開發，從原型設計、個人專案到中大型應用皆可勝任。若你想要快速啟動、享受極佳開發體驗，Vite 是非常推薦的選擇。
