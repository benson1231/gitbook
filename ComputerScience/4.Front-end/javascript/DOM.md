# HTML DOM 與 JavaScript 操作教學

DOM（Document Object Model）是瀏覽器提供的 API，讓 JavaScript 可以操作 HTML 結構與內容。

---

## 🔹 1. DOM 是什麼？

DOM 是瀏覽器內部建立的文件樹狀結構，每個 HTML 元素都是一個節點（Node），JavaScript 可透過 DOM 操作這些節點。

```html
<body>
  <h1 id="title">Hello</h1>
  <button id="btn">Click me</button>
</body>
```

---

## 🔹 2. 取得 DOM 元素

```js
// 常用方式
const title = document.getElementById("title");
const button = document.querySelector("#btn");
const allDivs = document.querySelectorAll("div");
```

---

## 🔹 3. 修改內容與屬性

```js
title.textContent = "新標題";
title.style.color = "red";
button.setAttribute("disabled", true);
```

---

## 🔹 4. 建立與插入節點

```js
const newP = document.createElement("p");
newP.textContent = "我是新段落";
document.body.appendChild(newP);
```

---

## 🔹 5. 刪除節點

```js
newP.remove();
```

---

## 🔹 6. 事件監聽器（Event Listener）

```js
button.addEventListener("click", () => {
  alert("你點了按鈕！");
});
```

---

## 🔹 7. 表單與輸入框

```html
<form id="myForm">
  <input type="text" id="name">
  <button type="submit">送出</button>
</form>
```

```js
const form = document.getElementById("myForm");
form.addEventListener("submit", e => {
  e.preventDefault();
  const name = document.getElementById("name").value;
  console.log("輸入的名字：", name);
});
```

---

## ✅ 常用方法速查

| 方法名稱                    | 功能                  |
| ----------------------- | ------------------- |
| `getElementById(id)`    | 根據 id 選取元素          |
| `querySelector(sel)`    | 使用 CSS 選擇器選取一個元素    |
| `querySelectorAll(sel)` | 選取所有符合的元素（NodeList） |
| `createElement(tag)`    | 建立新元素節點             |
| `appendChild(node)`     | 將子節點加入到父節點後面        |
| `remove()`              | 將元素從 DOM 移除         |

DOM 操作是 JavaScript 互動網頁的基礎，熟練這些操作能建立出各種動態效果。
