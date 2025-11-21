# 🚀 使用 GitHub Pages 部署 Vite + React 網站

這是一份完整教學，介紹如何將使用 Vite 建立的 React 專案部署到 GitHub Pages，並搭配 `gh-pages` 套件。

---

## ✅ 前置條件

* 已完成 Vite + React 專案初始化
* 已建立 GitHub repository（例如：`benson1231.github.io`）
* 專案已使用 Git 管理

---

## 📦 安裝 `gh-pages`

```bash
npm install gh-pages --save-dev
```

---

## 🛠️ 修改 `package.json`

### 1️⃣ 新增 `homepage` 欄位（注意：此 repo 為根網域部署）

```json
"homepage": "https://benson1231.github.io"
```

### 2️⃣ 增加 `scripts`

```json
"scripts": {
  "dev": "vite",
  "build": "vite build",
  "preview": "vite preview",
  "lint": "eslint .",
  "predeploy": "npm run build",
  "deploy": "gh-pages -d dist"
}
```

---

## ⚙️ 設定 `vite.config.js`

若你使用的是 `benson1231.github.io` 這種根目錄部署，base 設定需為空字串 `''`：

```js
// vite.config.js
import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  base: '', // 根網域部署請設為空字串
  plugins: [react()],
})
```

---

## 🔧 初始化 Git 並推送（如果還沒做）

```bash
git init
git remote add origin https://github.com/benson1231/benson1231.github.io
git add .
git commit -m "Initial commit"
git push -u origin main
```

---

## 🚀 部署到 GitHub Pages

```bash
npm run deploy
```

這會將 `dist/` 目錄推送到 `gh-pages` 分支。

---

## 🌐 啟用 GitHub Pages

1. 前往 GitHub repo 頁面：`https://github.com/benson1231/benson1231.github.io`
2. 點選 `Settings` → `Pages`
3. 在 **Build and deployment** 區塊設定：

   * Source: `Deploy from a branch`
   * Branch: `gh-pages`
   * Folder: `/ (root)`
4. 點擊 `Save`

等待 30 秒～1 分鐘，網站會部署成功！

---

## ✅ 完成！

你可以透過以下網址瀏覽你的網站：

🔍 [https://benson1231.github.io](https://benson1231.github.io)

---

## 🩀 小提醒

* 若使用子目錄（非 `benson1231.github.io`），則須將 `vite.config.js` 中的 `base` 設為 `/repo-name/`
* 若使用 React Router，請記得處理 404 或 Hash 模式
