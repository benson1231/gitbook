# FastAPI for Database

FastAPI 可與多種資料庫整合，以下說明如何分別連接並操作 MySQL、PostgreSQL 與 MongoDB。

---

## MySQL / PostgreSQL（使用 SQLAlchemy + databases）

### 1. 安裝套件

```bash
pip install sqlalchemy databases asyncmy psycopg2-binary
```

* MySQL 使用 `asyncmy`
* PostgreSQL 使用 `psycopg2-binary`

### 2. 建立連線與模型

```python
from sqlalchemy import Column, Integer, String, create_engine, MetaData, Table
from databases import Database

DATABASE_URL = "mysql+asyncmy://user:password@localhost/testdb"
# DATABASE_URL = "postgresql://user:password@localhost/testdb"

database = Database(DATABASE_URL)
metadata = MetaData()

notes = Table(
    "notes",
    metadata,
    Column("id", Integer, primary_key=True),
    Column("text", String(100)),
)
```

### 3. 建立 FastAPI 應用

```python
from fastapi import FastAPI

app = FastAPI()

@app.on_event("startup")
async def startup():
    await database.connect()

@app.on_event("shutdown")
async def shutdown():
    await database.disconnect()

@app.get("/notes")
async def read_notes():
    query = notes.select()
    return await database.fetch_all(query)
```

---

## MongoDB（使用 Motor）

### 1. 安裝 Motor 套件

```bash
pip install motor
```

### 2. 建立連線與操作資料

```python
from fastapi import FastAPI
from motor.motor_asyncio import AsyncIOMotorClient

app = FastAPI()

client = AsyncIOMotorClient("mongodb://localhost:27017")
db = client.mydatabase

@app.post("/mongo")
async def insert_doc():
    result = await db.mycollection.insert_one({"name": "FastAPI", "type": "MongoDB"})
    return {"inserted_id": str(result.inserted_id)}

@app.get("/mongo")
async def get_docs():
    docs = await db.mycollection.find().to_list(100)
    return docs
```

---

## 小結

| 資料庫        | 套件                       | 特點與用途           |
| ---------- | ------------------------ | --------------- |
| MySQL      | `sqlalchemy`, `asyncmy`  | 傳統關聯式資料庫，支援異步連線 |
| PostgreSQL | `sqlalchemy`, `psycopg2` | 高穩定性，常用於後端系統    |
| MongoDB    | `motor`                  | NoSQL，適合彈性文件結構  |

搭配 ORM 工具與非同步連線，FastAPI 能輕鬆整合主流資料庫，打造完整 API 後端服務。
