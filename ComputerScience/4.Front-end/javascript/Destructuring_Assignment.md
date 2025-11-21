# JavaScript 解構賦值（Destructuring Assignment）教學

解構賦值是一種簡潔的語法，能從陣列或物件中提取值並賦予變數。

---

## 🔹 1. 陣列解構

```js
const arr = ["Alice", 30, "Taipei"];
const [name, age, city] = arr;
console.log(name); // Alice
```

### 跳過元素

```js
const [first, , third] = [1, 2, 3];
console.log(third); // 3
```

### 剩餘元素

```js
const [head, ...rest] = [10, 20, 30, 40];
console.log(rest); // [20, 30, 40]
```

---

## 🔹 2. 物件解構

```js
const person = {
  name: "Ben",
  age: 25,
  city: "Kaohsiung"
};

const { name, age } = person;
console.log(name); // Ben
```

### 更名變數

```js
const { name: userName } = person;
console.log(userName); // Ben
```

### 預設值

```js
const { job = "學生" } = person;
console.log(job); // 學生
```

---

## 🔹 3. 解構用於函式參數

```js
function showUser({ name, age }) {
  console.log(`${name} - ${age}`);
}

const user = { name: "Cindy", age: 20 };
showUser(user); // Cindy - 20
```

---

## 🔹 4. 巢狀解構

```js
const data = {
  user: {
    name: "Derek",
    contact: {
      email: "derek@mail.com"
    }
  }
};

const {
  user: {
    contact: { email }
  }
} = data;
console.log(email); // derek@mail.com
```

---

解構賦值可以讓程式碼更簡潔直觀，是 ES6 中非常常用的語法，尤其在處理物件或 API 回傳資料時非常有用。
