# MongoDB 基礎語法教學（mongosh CLI）

以下為 MongoDB 使用 `mongosh` 指令列介面的基礎語法，包含資料庫操作、文件操作、查詢與更新等常用語句。

---

## 📁 資料庫操作

### 切換或建立資料庫

```js
use mydb  // 若 mydb 不存在，將在插入資料後建立
```

### 查看目前所有資料庫

```js
show dbs
```

---

## 📂 集合（Collection）操作

### 查看目前資料庫中的集合

```js
show collections
```

### 插入單筆文件

```js
db.users.insertOne({ name: "Alice", age: 30 })
```

### 插入多筆文件

```js
db.users.insertMany([
  { name: "Bob", age: 25 },
  { name: "Charlie", age: 35 }
])
```

---

## 🔍 查詢資料

### 查詢全部文件

```js
db.users.find()
```

### 查詢第一筆符合條件的文件

```js
db.users.findOne({ name: "Alice" })
```

### 計算符合條件的文件數量

```js
db.users.countDocuments({ age: { $gte: 30 } })
```

### 使用條件查詢

```js
db.users.find({ age: 30 })
```

### 條件運算符：

```js
// 大於等於、小於
{ age: { $gte: 30, $lt: 40 } }

// $in / $nin
{ name: { $in: ["Alice", "Bob"] } }
{ name: { $nin: ["Eve"] } }

// 欄位是否存在
{ email: { $exists: true } }

// 組合條件
{ $and: [ { age: { $gt: 20 } }, { name: "Alice" } ] }
{ $or:  [ { age: 25 }, { name: "Charlie" } ] }
{ age: { $not: { $gte: 30 } } }

// 正規表達式
{ name: /a/i }  // 模糊比對包含 a，不分大小寫
```

### 限制與排序

```js
// 只取前 2 筆
.find().limit(2)

// 跳過前 2 筆
.find().skip(2)

// 依年齡遞增排序
.find().sort({ age: 1 })

// 組合使用
.find().sort({ age: -1 }).limit(3)
```

---

## ✏️ 更新與刪除資料

### 更新單筆資料

```js
db.users.updateOne(
  { name: "Alice" },
  { $set: { age: 31 } }
)
```

### 更新多筆資料

```js
db.users.updateMany(
  { age: { $lt: 30 } },
  { $inc: { age: 1 } }
)
```

### 刪除單筆資料

```js
db.users.deleteOne({ name: "Alice" })
```

### 刪除多筆資料

```js
db.users.deleteMany({ age: { $gte: 35 } })
```

---

這些語法涵蓋 MongoDB 開發中最常使用的操作，適合初學者練習 CRUD 與複雜條件查詢。進階應用如聚合 (`aggregate`)、索引、驗證規則等可在後續學習補充。
