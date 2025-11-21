# JavaScript 模組（Modules）教學

模組化可以讓 JavaScript 程式碼更具結構性、可重用與可維護。ES6 引入原生模組語法，使用 `export` 與 `import` 來分享與引用程式碼。

---

## 🔹 1. 模組檔案基本語法

### a. 導出（export）

```js
// utils.js
export const add = (a, b) => a + b;
export const sub = (a, b) => a - b;
```

### b. 匯入（import）

```js
// main.js
import { add, sub } from "./utils.js";

console.log(add(2, 3));  // 5
```

---

## 🔹 2. 匯出方式

### a. 命名匯出（Named Exports）

可匯出多個變數或函式，使用大括號引入：

```js
export const name = "Tom";
export function greet() {}
```

### b. 預設匯出（Default Export）

每個模組只能有一個 `default` 匯出：

```js
// user.js
export default function sayHi() {
  console.log("Hi!");
}
```

```js
// main.js
import sayHi from "./user.js";
sayHi();
```

---

## 🔹 3. 匯入別名（as）

```js
import { add as plus } from "./utils.js";
console.log(plus(5, 7));
```

---

## 🔹 4. 結合 default 與 named export

```js
// math.js
export default "MathUtils";
export const PI = 3.14;
```

```js
import name, { PI } from "./math.js";
```

---

## 🔹 5. 模組限制

* 模組檔案必須使用 `.js` 且執行於支援 ES Modules 的環境（如瀏覽器需 `<script type="module">`）。
* 相對路徑必須加上副檔名，如 `./utils.js`

---

## 🔹 6. 瀏覽器支援

```html
<script type="module" src="main.js"></script>
```

---

模組化是大型應用中不可或缺的一環，推薦搭配打包工具如 Webpack、Vite、ESBuild 等，實現更高效的模組管理與編譯。
