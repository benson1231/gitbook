# JavaScript Fetch API 教學

`fetch()` 是一個用來向伺服器發送 HTTP 請求並取得回應的現代 JavaScript API。

---

## 🔹 1. 基本語法

```js
fetch(url)
  .then(response => response.json())
  .then(data => console.log(data))
  .catch(error => console.error("錯誤:", error));
```

---

## 🔹 2. GET 請求

```js
fetch("https://jsonplaceholder.typicode.com/posts/1")
  .then(res => res.json())
  .then(data => console.log(data));
```

---

## 🔹 3. POST 請求

```js
fetch("https://jsonplaceholder.typicode.com/posts", {
  method: "POST",
  headers: {
    "Content-Type": "application/json"
  },
  body: JSON.stringify({
    title: "Hello",
    body: "This is a test",
    userId: 1
  })
})
  .then(res => res.json())
  .then(data => console.log(data));
```

---

## 🔹 4. PUT 與 DELETE 請求

### PUT（更新整筆資料）

```js
fetch("https://jsonplaceholder.typicode.com/posts/1", {
  method: "PUT",
  headers: {
    "Content-Type": "application/json"
  },
  body: JSON.stringify({
    id: 1,
    title: "Updated Title",
    body: "New content",
    userId: 1
  })
});
```

### DELETE（刪除資料）

```js
fetch("https://jsonplaceholder.typicode.com/posts/1", {
  method: "DELETE"
});
```

---

## 🔹 5. async/await 寫法

```js
async function getPost() {
  try {
    const response = await fetch("https://jsonplaceholder.typicode.com/posts/1");
    const data = await response.json();
    console.log(data);
  } catch (error) {
    console.error("錯誤:", error);
  }
}

getPost();
```

---

## 🔹 6. 常見錯誤處理

* 伺服器未回應：`fetch()` 不會拋出錯誤，除非網路錯誤。
* 使用 `!response.ok` 判斷 HTTP 狀態：

```js
fetch("/api")
  .then(response => {
    if (!response.ok) throw new Error("伺服器錯誤");
    return response.json();
  })
  .then(data => console.log(data))
  .catch(err => console.error(err));
```

---

Fetch 是前端與後端互動的核心工具，搭配 `async/await` 可大幅簡化非同步流程與錯誤處理。
