# Concurrent Programming

Python 提供多種方式實現並行（concurrent）與平行（parallel）運算，常見方法包括：

* `threading`：多執行緒，適合 I/O 密集型任務
* `multiprocessing`：多行程，適合 CPU 密集型任務
* `asyncio`：非同步協程，適合高延遲 I/O 操作
* `nest_asyncio`：允許在 Jupyter Notebook 或已啟動的 event loop 中重入 asyncio

---

## 一、threading：多執行緒

```python
import threading
import time

def worker(name):
    print(f"{name} 開始")
    time.sleep(2)
    print(f"{name} 結束")

thread1 = threading.Thread(target=worker, args=("執行緒 1",))
thread2 = threading.Thread(target=worker, args=("執行緒 2",))

thread1.start()
thread2.start()
thread1.join()
thread2.join()
```

適合用於：等待網路、磁碟 I/O，Web 爬蟲等場景。

---

## 二、multiprocessing：多行程

```python
from multiprocessing import Process

def worker():
    print("處理中...")

if __name__ == '__main__':
    p = Process(target=worker)
    p.start()
    p.join()
```

適合用於：圖像轉換、數值計算、深度學習訓練等高 CPU 任務。

---

## 三、asyncio：非同步協程

```python
import asyncio

async def say_hello():
    print("Hello")
    await asyncio.sleep(1)
    print("World")

asyncio.run(say_hello())
```

適合用於：API 串接、大量非同步網路請求、WebSocket 等。

---

## 四、nest\_asyncio：允許在 Notebook 中使用 asyncio

```python
import nest_asyncio
import asyncio

nest_asyncio.apply()  # 解決 event loop 已啟動錯誤

async def demo():
    await asyncio.sleep(1)
    print("協程執行完畢")

await demo()
```

適合用於：Jupyter Notebook、REPL 中測試 async 程式。

---

## 五、比較表

| 方法                | 適用情境          | 並行模式 | 特性               |
| ----------------- | ------------- | ---- | ---------------- |
| `threading`       | I/O 密集        | 並行   | 較簡易，受 GIL 限制     |
| `multiprocessing` | CPU 密集        | 平行   | 可多核運算，避免 GIL     |
| `asyncio`         | 非同步 I/O       | 協程   | 單執行緒高效佇列調度       |
| `nest_asyncio`    | Notebook 裡測試用 | 協程   | 解決 event loop 限制 |

---

善用這些並行/非同步工具可提升 Python 程式效能，選擇合適的模式與場景搭配能有效處理高效能或高併發任務。
