# React

React 是由 Meta（前 Facebook）開發的 JavaScript 函式庫，用於建立使用者介面（UI），特別是單頁應用程式（SPA）。它採用元件化設計，使得 UI 可重用、易於維護。

## 1. React 的核心概念

### 1.1 元件（Component）

React 的 UI 是由一個個元件組成。每個元件可以是函式或類別，並包含自己的狀態與生命週期。

```jsx
function Welcome(props) {
  return <h1>Hello, {props.name}</h1>;
}
```

### 1.2 JSX

JSX 是 JavaScript 的語法擴充，讓你可以在 JavaScript 中寫 HTML 標籤。

```jsx
const element = <h1>歡迎回來</h1>;
```

### 1.3 狀態與事件處理（State & Events）

```jsx
import { useState } from 'react';

function Counter() {
  const [count, setCount] = useState(0);
  return (
    <button onClick={() => setCount(count + 1)}>
      點我：{count}
    </button>
  );
}
```

### 1.4 父子元件傳遞資料（Props）

Props 是父元件傳遞給子元件的只讀資料。

```jsx
<Welcome name="Alice" />
```

---

## 2. React 優點

* 元件化結構：模組化開發，便於維護與重用
* 虛擬 DOM：提升渲染效能
* 單向資料流：數據管理可預測
* 廣大社群與生態系：如 React Router、Redux、Next.js 等

---

## 3. 建立 React 專案（使用 Vite）

Vite 是一個現代化的前端建構工具，啟動快、建置快，適合搭配 React 開發。

### 安裝步驟：

```bash
npm create vite@latest my-app -- --template react
cd my-app
npm install
npm run dev
```

### 專案結構簡介：

```
my-app/
├── index.html           ← HTML 入口
├── src/
│   ├── main.jsx         ← 入口 JS
│   └── App.jsx          ← 根元件
├── public/              ← 靜態資源資料夾
└── vite.config.js       ← Vite 設定檔
```

### 優點：

* 啟動速度快（原生 ESM + on-demand bundling）
* 適合現代框架如 React、Vue、Svelte 等
* 更簡化的設定與專案結構

---

React 是現代前端開發的重要基礎之一，無論是單頁應用、行動應用（React Native）或是 SSR（Next.js）開發都廣泛應用。
