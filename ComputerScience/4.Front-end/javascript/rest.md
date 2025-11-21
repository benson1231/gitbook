# JavaScript Rest 語法教學

`rest` 語法（...）用來收集剩餘的元素，常用於函式參數或解構賦值，是與 `spread` 語法相對應的用法。

---

## 🔹 1. 函式參數中使用 rest

```js
function sum(...numbers) {
  return numbers.reduce((total, num) => total + num, 0);
}

console.log(sum(1, 2, 3)); // 6
```

> `...numbers` 會收集所有傳入參數成為一個陣列。

---

## 🔹 2. 陣列解構中的 rest

```js
const [first, ...others] = [10, 20, 30, 40];
console.log(first);  // 10
console.log(others); // [20, 30, 40]
```

---

## 🔹 3. 物件解構中的 rest

```js
const user = {
  name: "Alice",
  age: 25,
  city: "Taipei"
};

const { name, ...restInfo } = user;
console.log(name);     // Alice
console.log(restInfo); // { age: 25, city: "Taipei" }
```

---

## ✅ 小結

| 用法位置 | 功能          |
| ---- | ----------- |
| 函式參數 | 收集多個傳入參數    |
| 陣列解構 | 收集剩下的陣列元素   |
| 物件解構 | 收集剩下的屬性為新物件 |

`rest` 是讓 JavaScript 變數處理更具彈性的重要語法，常與 `spread` 搭配一起使用。
