# Webpack 介紹與學習必要性分析

`Webpack` 是一個現代 JavaScript 應用程式的靜態模組打包工具（module bundler），能將專案中的 JavaScript、CSS、圖片等資源打包為可部署的靜態資產，廣泛用於前端開發流程中。

## Webpack 的核心功能

* 將多個模組（module）與資源打包成單一或多個 bundle
* 支援 ES6 模組語法（import/export）
* 使用 Loader 處理 CSS、SCSS、圖片等非 JavaScript 資源
* 使用 Plugin 擴充功能，如壓縮、環境注入、Html 產生等
* 支援開發伺服器（webpack-dev-server）與熱更新（HMR）

## 常見指令與設定

### 安裝 Webpack

```bash
npm install --save-dev webpack webpack-cli
```

### 建立基本設定檔 `webpack.config.js`

```js
const path = require('path');

module.exports = {
  entry: './src/index.js',
  output: {
    filename: 'bundle.js',
    path: path.resolve(__dirname, 'dist')
  },
  module: {
    rules: [
      {
        test: /\.css$/,
        use: ['style-loader', 'css-loader']
      }
    ]
  }
};
```

### 打包專案

```bash
npx webpack
```

## 是否仍需學 Webpack？

### ✅ 適合學習的情境：

* **需要高度自訂的前端建構流程**（如老專案、複雜依賴管理）
* **維護既有使用 Webpack 的大型專案**
* **進階 React/Vue 工程化開發需求**

### ❌ 可選擇跳過的情境：

* 使用 **Vite、Next.js、Nuxt、Create React App 等現代工具** 已隱藏 Webpack 細節
* 小型或原型開發可不必深入了解 Webpack 配置

## 替代工具簡介

| 工具      | 優勢                     |
| ------- | ---------------------- |
| Vite    | 更快的熱更新與模組載入、預設支援 ES 模組 |
| Parcel  | 零設定打包工具，適合快速啟動專案       |
| ESBuild | 超高速打包器，適合 CLI 或工具鏈整合   |

---

總結來說，**Webpack 不再是所有專案必學的工具**，但對需要自訂構建流程或維護老專案的開發者而言仍具備價值。若你熟悉 React/Vue 並希望深入工程化開發，Webpack 仍值得理解。
