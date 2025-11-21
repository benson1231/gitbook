# Vue

Vue 是一個漸進式的 JavaScript 前端框架，由尤雨溪（Evan You）開發，特別適合建構互動式網頁應用。它主打簡潔語法與雙向資料綁定，易於上手且具高度彈性。

## 1. Vue 的核心概念

### 1.1 宣告式渲染

Vue 使用模板語法（template syntax）來描述 UI。資料變動會自動反映在畫面上。

```html
<div id="app">
  <p>{{ message }}</p>
</div>
```

```js
const app = Vue.createApp({
  data() {
    return {
      message: 'Hello Vue!'
    };
  }
}).mount('#app');
```

### 1.2 雙向綁定（v-model）

```html
<input v-model="name">
<p>你的名字是：{{ name }}</p>
```

### 1.3 條件與迴圈渲染

```html
<p v-if="isLoggedIn">歡迎回來！</p>
<ul>
  <li v-for="item in items" :key="item.id">{{ item.text }}</li>
</ul>
```

### 1.4 元件化（Component）

Vue 鼓勵以元件為單位構建整個應用程式。每個元件都是一個自包含的可重用單位，包含 HTML 模板、JavaScript 邏輯與 CSS 樣式，並透過 `props` 接收資料，`emit` 傳遞事件。

#### 單檔元件（Single File Component, SFC）標準格式：

```vue
<!-- HelloWorld.vue -->
<template>
  <h1>Hello, {{ name }}</h1>
</template>

<script>
export default {
  name: 'HelloWorld',
  props: {
    name: {
      type: String,
      required: true
    }
  }
};
</script>

<style scoped>
h1 {
  color: teal;
}
</style>
```

#### 在父元件中引入與使用：

```vue
<template>
  <HelloWorld name="Vue 使用者" />
</template>

<script>
import HelloWorld from './HelloWorld.vue';

export default {
  components: {
    HelloWorld
  }
};
</script>
```

每個 `.vue` 檔案需包含：

* `<template>`：元件的 UI 結構
* `<script>`：元件邏輯與資料管理
* `<style>`：元件樣式，建議使用 `scoped` 限定範圍

---

## 2. Vue 優點

* 模板語法直觀，學習曲線平緩
* 雙向綁定強化表單互動
* 單檔元件（SFC）整合 HTML、JS、CSS
* 支援小型到大型應用的彈性架構
* 生態系完整（Vue Router、Pinia、Vite）

---

## 3. 建立 Vue 專案（使用 Vite）

### 安裝步驟：

```bash
npm create vite@latest my-vue-app -- --template vue
cd my-vue-app
npm install
npm run dev
```

### 專案結構簡介：

```
my-vue-app/
├── index.html
├── src/
│   ├── main.js         ← 入口檔案
│   └── App.vue         ← 根元件（單檔元件）
└── vite.config.js
```

---

Vue 是建立互動式網頁與單頁應用的強力工具，若你偏好 HTML+JS 模板式開發風格，Vue 是非常值得學習的選擇。
