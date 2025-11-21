# Exception

Python 提供多種內建例外類別，用於處理不同類型的錯誤狀況。理解各種例外的用途有助於撰寫更健壯的錯誤處理邏輯。

---

## 一、例外階層架構概觀

```
BaseException
├── SystemExit
├── KeyboardInterrupt
├── GeneratorExit
└── Exception
    ├── ArithmeticError
    │   ├── ZeroDivisionError
    │   └── OverflowError
    ├── AttributeError
    ├── BufferError
    ├── EOFError
    ├── ImportError
    │   └── ModuleNotFoundError
    ├── LookupError
    │   ├── IndexError
    │   └── KeyError
    ├── MemoryError
    ├── NameError
    │   └── UnboundLocalError
    ├── OSError
    │   ├── FileNotFoundError
    │   └── PermissionError
    ├── RuntimeError
    │   └── RecursionError
    ├── NotImplementedError
    ├── UnicodeError
    │   ├── UnicodeDecodeError
    │   └── UnicodeTranslateError
    └── SystemError
```

---

## 二、BaseException

* **BaseException**：所有異常的基礎類別，不建議直接使用（會捕捉到像是 KeyboardInterrupt）

---

## 三、Exception 衍生類別

| 類別名稱                  | 說明                    |
| --------------------- | --------------------- |
| ArithmeticError       | 數學錯誤的基類               |
| ZeroDivisionError     | 除以零                   |
| OverflowError         | 運算結果超出可表示範圍           |
| AttributeError        | 訪問不存在的屬性              |
| BufferError           | 緩衝區操作錯誤               |
| EOFError              | 到達檔案結尾（如 input() 無輸入） |
| ImportError           | 導入模組失敗                |
| ModuleNotFoundError   | 模組不存在（ImportError 子類） |
| IndexError            | 序列索引超出範圍              |
| KeyError              | 字典鍵不存在                |
| KeyboardInterrupt     | 使用者中斷（如 Ctrl+C）       |
| MemoryError           | 記憶體不足                 |
| NameError             | 名稱未定義                 |
| UnboundLocalError     | 局部變數使用前未賦值            |
| OSError               | 操作系統錯誤，如開檔失敗          |
| FileNotFoundError     | 檔案不存在（OSError 子類）     |
| PermissionError       | 權限不足的操作（OSError 子類）   |
| RuntimeError          | 其他無特定類型的執行期錯誤         |
| RecursionError        | 遞迴深度超限                |
| NotImplementedError   | 尚未實作的函式或類別            |
| UnicodeDecodeError    | 解碼錯誤                  |
| UnicodeTranslateError | 翻譯過程錯誤                |
| SystemError           | Python 解譯器內部錯誤        |

---

## 四、raise 語句：主動引發例外

在某些情境下，開發者可以透過 `raise` 明確引發例外，以中斷執行流程並提示錯誤原因。

### 語法範例：

```python
def open_register(employee_status):
    if employee_status == 'Authorized':
        print('Successfully opened cash register')
    else:
        raise TypeError("error: 員工未授權")
```

### 使用方式：

```python
raise ValueError("輸入值無效")
raise KeyError("找不到指定鍵")
```

### 自訂例外類別：

```python
class CustomError(Exception):
    pass

raise CustomError("自訂錯誤訊息")
```

---

## 五、實用建議

* 使用 `try...except` 處理特定錯誤，避免使用過度廣泛的 `except:` 或 `except BaseException:`。
* 可搭配 `finally` 或 `else` 區塊增進安全性與可讀性。

```python
try:
    result = 10 / 0
except ZeroDivisionError:
    print("不能除以零！")
finally:
    print("程式結束")
```

---

掌握這些常見例外類別與用途，能讓你撰寫出更穩健、安全且易於除錯的 Python 程式碼。
