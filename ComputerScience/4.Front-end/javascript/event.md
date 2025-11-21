# JavaScript 事件處理教學

事件處理（Event Handling）是 JavaScript 與使用者互動的核心機制，可用來監聽點擊、輸入、鍵盤等動作。

---

## 🔹 1. 什麼是事件？

事件是一種發生在元素上的互動，例如：

* 使用者點擊（click）
* 滑鼠移動（mousemove）
* 鍵盤輸入（keydown）
* 表單提交（submit）
* 輸入變更（input/change）

---

## 🔹 2. 基本事件綁定方式

### 方法一：HTML 中直接綁定（不建議）

```html
<button onclick="alert('點到了')">點我</button>
```

### 方法二：JavaScript 使用 `addEventListener()`（推薦）

```js
const btn = document.querySelector("#myBtn");
btn.addEventListener("click", () => {
  alert("你點了按鈕！");
});
```

---

## 🔹 3. 常見事件類型

| 事件類型                     | 說明       |
| ------------------------ | -------- |
| `click`                  | 點擊元素     |
| `dblclick`               | 雙擊       |
| `mouseover` / `mouseout` | 滑鼠移入/移出  |
| `keydown` / `keyup`      | 鍵盤按下/放開  |
| `input` / `change`       | 表單輸入/變更值 |
| `submit`                 | 表單送出     |

---

## 🔹 4. 事件物件（Event Object）

事件處理函式會接收到一個 `event` 物件，提供事件相關資訊。

```js
btn.addEventListener("click", function(e) {
  console.log(e.target);  // 觸發事件的元素
});
```

---

## 🔹 5. 移除事件監聽器

```js
function handleClick() {
  console.log("點擊了");
}

btn.addEventListener("click", handleClick);
btn.removeEventListener("click", handleClick);
```

---

## 🔹 6. 阻止預設行為與冒泡

```js
link.addEventListener("click", function(e) {
  e.preventDefault();  // 阻止跳轉連結
});

child.addEventListener("click", function(e) {
  e.stopPropagation(); // 阻止事件冒泡
});
```

---

事件處理是建立互動網頁的關鍵基礎，配合 DOM 操作與條件判斷，可製作出豐富的使用者體驗。
