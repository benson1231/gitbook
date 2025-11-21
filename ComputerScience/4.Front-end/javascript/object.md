# JavaScript 物件（Object）教學

JavaScript 中的物件是一種儲存鍵值對（key-value pairs）的資料結構，用來描述一個實體或集合，是核心語法之一。

---

## 🔹 1. 建立物件

### 方法一：物件字面量（Object Literal）

```js
const person = {
  name: "Alice",
  age: 30,
  isStudent: false
};
```

### 方法二：使用 `new Object()`

```js
const person = new Object();
person.name = "Alice";
person.age = 30;
```

---

## 🔹 2. 存取與修改屬性

### 點記號（dot notation）

```js
console.log(person.name);     // 讀取
person.age = 31;              // 修改
```

### 括號記號（bracket notation）

```js
console.log(person["name"]);
person["isStudent"] = true;
```

---

## 🔹 3. 新增與刪除屬性

```js
person.email = "alice@mail.com";   // 新增

delete person.isStudent;           // 刪除
```

---

## 🔹 4. 迭代物件屬性（for...in）

```js
for (let key in person) {
  console.log(`${key}: ${person[key]}`);
}
```

---

## 🔹 5. 巢狀物件（Nested Object）

```js
const student = {
  name: "Bob",
  scores: {
    math: 90,
    english: 85
  }
};

console.log(student.scores.math); // 90
```

---

## 🔹 6. 陣列中的物件

```js
const users = [
  { name: "Amy", age: 20 },
  { name: "Ben", age: 25 }
];

console.log(users[0].name); // Amy
```

---

## 🔹 7. 方法（物件中的函式）

```js
const car = {
  brand: "Toyota",
  start: function() {
    console.log("發動引擎");
  }
};

car.start();
```

### 簡寫方式

```js
const car = {
  brand: "Toyota",
  start() {
    console.log("啟動中...");
  }
};
```

---

## 🔹 8. `this` 關鍵字

`this` 代表目前的物件本身

```js
const user = {
  name: "Ken",
  greet() {
    console.log(`Hi, I'm ${this.name}`);
  }
};

user.greet(); // Hi, I'm Ken
```

---

物件是 JavaScript 的核心資料結構之一，能有效管理複雜資料與行為。
