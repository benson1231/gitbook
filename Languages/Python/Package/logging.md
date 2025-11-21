# logging

`logging` 是 Python 內建模組之一，用於紀錄程式執行過程中的重要事件與除錯資訊。相比 `print()`，它提供等級分類、格式化、輸出至檔案等更強大的功能。

---

## 一、基本使用

```python
import logging

logging.basicConfig(level=logging.INFO)
logging.info("這是一條資訊訊息")
logging.warning("這是一條警告")
logging.error("這是錯誤訊息")
```

---

## 二、訊息等級

| 等級       | 說明          |
| -------- | ----------- |
| DEBUG    | 偵錯用，最低層級    |
| INFO     | 一般運行資訊      |
| WARNING  | 警告，但程式仍可執行  |
| ERROR    | 錯誤，程式可能異常中止 |
| CRITICAL | 嚴重錯誤，導致系統中止 |

---

## 三、自訂輸出格式

```python
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
```

輸出範例：

```
2025-06-04 22:10:23,456 - INFO - 執行完成
```

---

## 四、輸出到檔案

```python
logging.basicConfig(
    filename='app.log',
    filemode='w',  # 'a' 為 append 模式
    level=logging.INFO,
    format='%(levelname)s:%(message)s'
)
```

---

## 五、與例外搭配使用

```python
try:
    1 / 0
except ZeroDivisionError as e:
    logging.exception("發生例外錯誤")
```

此時會自動包含 traceback 訊息。

---

## 六、進階用法：多 logger、handler 架構

```python
logger = logging.getLogger("myapp")
file_handler = logging.FileHandler("myapp.log")
stream_handler = logging.StreamHandler()

formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)
logger.addHandler(stream_handler)
logger.setLevel(logging.INFO)
```

---

`logging` 是撰寫大型應用程式、機器學習 pipeline 或 REST API 等程式的除錯與維護關鍵工具。建議在正式應用中以 log 取代 `print()`，便於除錯與追蹤。
