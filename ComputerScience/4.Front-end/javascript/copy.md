# JavaScript 淺拷貝（Shallow Copy）與深拷貝（Deep Copy）教學

在 JavaScript 中，複製物件時需要注意是「淺拷貝」還是「深拷貝」，否則可能會導致資料同步異動的問題。

---

## 🔹 1. 淺拷貝（Shallow Copy）

淺拷貝只複製**第一層**屬性，若屬性值是物件或陣列，仍是**參考同一記憶體位址**。

### 常見方法：

```js
// 使用 Object.assign()
const obj1 = { a: 1, b: { c: 2 } };
const shallow = Object.assign({}, obj1);

// 使用展開運算子（spread）
const shallow2 = { ...obj1 };

shallow.b.c = 999;
console.log(obj1.b.c); // 999，原始物件也被改變
```

---

## 🔹 2. 深拷貝（Deep Copy）

深拷貝會**完整複製所有巢狀結構**，兩者資料互不影響。

### 方法一：JSON 方式（簡單但有限制）

```js
const obj2 = { a: 1, b: { c: 2 } };
const deep = JSON.parse(JSON.stringify(obj2));

deep.b.c = 888;
console.log(obj2.b.c); // 2，原始物件不變
```

> 限制：無法處理 `function`、`Date`、`Map`、`Set` 等特殊型別。

### 方法二：遞迴函式深拷貝

```js
function deepClone(obj) {
  if (obj === null || typeof obj !== 'object') return obj;

  if (Array.isArray(obj)) {
    return obj.map(deepClone);
  }

  const result = {};
  for (let key in obj) {
    result[key] = deepClone(obj[key]);
  }
  return result;
}
```

---

## 🔹 3. 使用第三方函式庫

如 [Lodash](https://lodash.com/)

```js
import cloneDeep from 'lodash/cloneDeep';
const newObj = cloneDeep(obj);
```

---

## ✅ 小結

| 類型  | 說明             | 是否複製巢狀物件 |
| --- | -------------- | -------- |
| 淺拷貝 | 複製第一層，巢狀仍共享記憶體 | ❌        |
| 深拷貝 | 所有層級皆複製新物件     | ✅        |

了解兩者差異有助於正確處理資料結構與狀態管理，避免意外同步修改的 Bug。
