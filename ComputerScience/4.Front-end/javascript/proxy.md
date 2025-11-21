# JavaScript Proxy（代理物件）教學

`Proxy` 是 ES6 引入的功能，允許你攔截對物件的存取行為，例如讀取屬性、設定值、刪除屬性等。

---

## 🔹 1. 基本語法

```js
const target = {
  name: "Alice",
};

const handler = {
  get: (obj, prop) => {
    console.log(`取得屬性 ${prop}`);
    return obj[prop];
  },
  set: (obj, prop, value) => {
    console.log(`設定 ${prop} 為 ${value}`);
    obj[prop] = value;
    return true;
  }
};

const proxy = new Proxy(target, handler);

console.log(proxy.name); // 觸發 get
proxy.age = 30;          // 觸發 set
```

---

## 🔹 2. 支援的陷阱（traps）方法

| 方法名稱                       | 說明                                          |
| -------------------------- | ------------------------------------------- |
| `get`                      | 讀取屬性時觸發                                     |
| `set`                      | 設定屬性時觸發                                     |
| `has`                      | 使用 `in` 運算子時觸發                              |
| `deleteProperty`           | 使用 `delete` 時觸發                             |
| `ownKeys`                  | 使用 `Object.keys()` 或 `for...in` 觸發          |
| `getOwnPropertyDescriptor` | 使用 `Object.getOwnPropertyDescriptor()` 時觸發  |
| ...                        | 其他如 `defineProperty`, `preventExtensions` 等 |

---

## 🔹 3. 實用案例

### a. 驗證輸入資料

```js
const user = new Proxy({}, {
  set(obj, prop, value) {
    if (prop === "age" && typeof value !== "number") {
      throw new TypeError("年齡必須為數字");
    }
    obj[prop] = value;
    return true;
  }
});

user.age = 25;   // ✅
user.age = "xx"; // ❌ TypeError
```

### b. 提供預設值

```js
const withDefault = new Proxy({}, {
  get(obj, prop) {
    return prop in obj ? obj[prop] : "預設值";
  }
});

console.log(withDefault.title); // 預設值
```

---

## 🔹 4. 搭配 Reflect 使用

```js
const proxy = new Proxy(obj, {
  get(target, prop, receiver) {
    return Reflect.get(target, prop, receiver);
  },
  set(target, prop, value, receiver) {
    return Reflect.set(target, prop, value, receiver);
  }
});
```

---

`Proxy` 是元程式設計（Meta-programming）的核心工具之一，能動態攔截並定義物件的基本操作，非常適合用於資料驗證、快取、偽裝物件等進階應用。
