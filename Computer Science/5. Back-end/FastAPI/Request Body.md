# Request Body as JSON

FastAPI 支援從前端或第三方服務以 JSON 格式傳送資料，並透過 `pydantic` 模型進行驗證與解析。

---

## 1. 定義 JSON 資料結構

使用 `pydantic.BaseModel` 定義接收格式：

```python
from fastapi import FastAPI
from pydantic import BaseModel

app = FastAPI()

class Item(BaseModel):
    name: str
    price: float
    in_stock: bool = True
```

---

## 2. 接收 JSON 請求內容

```python
@app.post("/items")
def create_item(item: Item):
    return {"received": item}
```

---

## 3. 呼叫端範例（JavaScript Fetch）

```js
fetch("/items", {
  method: "POST",
  headers: {"Content-Type": "application/json"},
  body: JSON.stringify({
    name: "Keyboard",
    price: 350.5,
    in_stock: false
  })
});
```

---

## 4. 使用別名與嵌套模型

```python
class Creator(BaseModel):
    username: str
    email: str

class Product(BaseModel):
    title: str
    description: str | None = None
    tags: list[str] = []
    owner: Creator

@app.post("/product")
def post_product(p: Product):
    return {"ok": True, "product": p.title, "owner": p.owner.username}
```

---

## 5. 驗證與錯誤回應

FastAPI 將自動產生 422 回應碼，當輸入不符合資料型別或缺少欄位時會給出結構化錯誤訊息。

```json
{
  "detail": [
    {
      "loc": ["body", "price"],
      "msg": "field required",
      "type": "value_error.missing"
    }
  ]
}
```

---

## 小結

| 概念      | 說明                                  |
| ------- | ----------------------------------- |
| JSON 格式 | 使用 `Content-Type: application/json` |
| 結構定義    | 使用 `BaseModel` 建立欄位型別               |
| 自動驗證與解析 | FastAPI 自動處理資料轉換與錯誤訊息               |
| 多層嵌套    | 模型中可包含其他模型（巢狀結構）                    |

JSON 是現代 API 傳輸資料的主流格式，搭配 FastAPI 與 Pydantic 可有效建立安全、嚴謹的資料驗證流程。
