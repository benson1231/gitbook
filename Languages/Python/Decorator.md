# Decorator

裝飾器（Decorator）是 Python 中一種非常強大的語法結構，用來在**不修改原本函式或類別程式碼的情況下，動態增加其功能**。常見於日常開發中，例如登入驗證、日誌紀錄、權限管理等。

## 基本概念

裝飾器本質上是一個**接收函式並回傳新函式的函式**。

最簡單的裝飾器範例：

```python
def my_decorator(callback):
    def innerFunction():
        print("我會先出現")
        callback()
        print("我會在後面")
    return innerFunction

@my_decorator
def say_hello():
    print("Hello!")

say_hello()
```

執行結果：

```
我會先出現
Hello!
我會在後面
```

`@my_decorator` 等同於：

```python
say_hello = my_decorator(say_hello)
```

---

## 帶參數的裝飾器

如果原本的函式有參數，需要在 `wrapper` 中加上 `*args, **kwargs` 來接收不定數量的參數：

```python
def my_decorator(func):
    def wrapper(*args, **kwargs):
        print("開始執行...")
        result = func(*args, **kwargs)
        print("結束執行...")
        return result
    return wrapper

@my_decorator
def add(a, b):
    print(a + b)

add(3, 5)
```

---

## 裝飾器的常見應用

1. **權限驗證**
2. **日誌紀錄**（logging）
3. **計時器（performance profiling）**
4. **快取（caching）**

範例 - 計算函式執行時間：

```python
import time

def timer(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"執行時間：{end_time - start_time:.4f} 秒")
        return result
    return wrapper

@timer
def slow_function():
    time.sleep(2)
    print("慢速運算完成")

slow_function()
```

---

## 小提醒

- 使用 `functools.wraps` 可以保留原函式的名稱與說明文字（metadata）。
- 裝飾器也可以疊加（多層裝飾）。

```python
from functools import wraps

def decorator(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper
```

---

## 小結

裝飾器是 Python 中提高程式碼可重用性與清晰度的重要工具，理解其基本結構與使用場景後，可以大幅優化你的程式設計結構。

