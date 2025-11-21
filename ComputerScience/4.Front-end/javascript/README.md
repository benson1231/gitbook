# JavaScript 基礎語法教學

JavaScript 是一種在瀏覽器與伺服器上皆可運行的腳本語言，廣泛應用於網頁互動與應用開發。

---

## 🔹 1. 資料型態（Data Types）

### 原始型別（Primitive）

* `Number`：數字，例如 `10`、`3.14`
* `String`：字串，例如 `'Hello'`、`"World"`
* `Boolean`：布林值，`true` 或 `false`
* `undefined`：尚未賦值的變數
* `null`：空值
* `BigInt`：可表示超大整數
* `Symbol`：獨一無二的值，常用作物件屬性

### 物件型別（Object）

* `Object`、`Array`、`Function`、`Date`、`RegExp` 等皆為物件

---

## 🔹 2. 變數與常數（let / const / var）

```js
let age = 25;       // 可重新賦值
const PI = 3.14;    // 不可變動
var name = 'Tom';   // 舊語法（不建議）
```

---

## 🔹 3. 運算子（Operators）

### 算術運算子

`+` `-` `*` `/` `%` `**`（次方）

### 比較運算子

`==`、`===`、`!=`、`!==`、`>`、`<`、`>=`、`<=`

### 邏輯運算子

`&&`（且）、`||`（或）、`!`（非）

---

## 🔹 4. 流程控制（Control Flow）

### 條件判斷

```js
if (age >= 18) {
  console.log("成年人");
} else {
  console.log("未成年");
}
```

### switch 判斷

```js
switch (fruit) {
  case 'apple':
    console.log('蘋果');
    break;
  case 'banana':
    console.log('香蕉');
    break;
  default:
    console.log('未知');
}
```

### 迴圈（for / while / forEach）

```js
for (let i = 0; i < 5; i++) {
  console.log(i);
}

let j = 0;
while (j < 5) {
  console.log(j);
  j++;
}

["a", "b"].forEach(e => console.log(e));
```

---

## 🔹 5. 函式（Functions）

### 傳統函式

```js
function greet(name) {
  return `Hello, ${name}`;
}
```

### 匿名函式（表達式）

```js
const add = function (a, b) {
  return a + b;
};
```

### 箭頭函式（Arrow Function）

```js
const multiply = (x, y) => x * y;
```

---

JavaScript 是動態且彈性高的語言，熟悉以上語法是學習前端開發、互動網頁設計的第一步。
