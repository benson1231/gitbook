# 📘 Python 檔案操作（File I/O）教學

Python 提供了內建的檔案操作功能，可讓使用者輕鬆地讀取與寫入檔案。本文將介紹常用的檔案操作方法與使用範例。

---

## 🗂️ 一、開啟與關閉檔案

### `open()`

用來開啟檔案並返回檔案物件。

```python
f = open('example.txt', 'r')  # 以讀取模式開啟檔案
```

### 常見模式

| 模式    | 說明                        |
| ----- | ------------------------- |
| `'r'` | 讀取模式（檔案必須存在）              |
| `'w'` | 寫入模式（會覆蓋原檔案）              |
| `'a'` | 附加模式（寫入內容會加到檔尾）           |
| `'x'` | 建立新檔（若檔案存在會出錯）            |
| `'b'` | 二進位模式（與其他模式搭配使用，如 `'rb'`） |
| `'t'` | 文字模式（預設值）                 |

### `close()`

使用完檔案後應關閉，釋放系統資源。

```python
f.close()
```

📌 **建議用 `with` 語法自動關閉檔案：**

```python
with open('example.txt', 'r') as f:
    data = f.read()
```

---

## 📖 二、讀取檔案內容

### `read()`

一次讀取整個檔案內容（適合小檔案）。

```python
with open('example.txt', 'r') as f:
    content = f.read()
    print(content)
```

可指定讀取長度：

```python
f.read(10)  # 讀取前10個字元
```

### `readline()`

逐行讀取一行文字。

```python
with open('example.txt', 'r') as f:
    line = f.readline()
    print(line)
```

### `readlines()`

一次讀取所有行，回傳一個 **list**。

```python
with open('example.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        print(line.strip())
```

---

## ✏️ 三、寫入檔案內容

### `write()`

將字串寫入檔案。

```python
with open('output.txt', 'w') as f:
    f.write('Hello, Python!\n')
```

### `writelines()`

將多行字串（list）一次寫入檔案。

```python
lines = ['Line 1\n', 'Line 2\n', 'Line 3\n']
with open('output.txt', 'w') as f:
    f.writelines(lines)
```

### 附加模式 (`'a'`)

在檔尾新增內容而不覆蓋原資料。

```python
with open('output.txt', 'a') as f:
    f.write('Appended line.\n')
```

---

## 🔄 四、其他常用方法

| 方法                       | 說明                                 |
| ------------------------ | ---------------------------------- |
| `f.tell()`               | 回傳目前檔案指標位置（以位元組計）。                 |
| `f.seek(offset, whence)` | 移動檔案指標位置。<br>例如：`f.seek(0)` 會回到開頭。 |
| `f.flush()`              | 立即將緩衝區資料寫入檔案。                      |
| `f.truncate(size)`       | 截斷檔案到指定大小。                         |

```python
with open('output.txt', 'r+') as f:
    f.seek(0)
    print(f.tell())  # 顯示位置
```

---

## 📦 五、處理二進位檔案（binary mode）

讀取圖片、影片或其他非文字檔案時使用 `'rb'` 或 `'wb'`。

```python
with open('image.jpg', 'rb') as f:
    data = f.read()

with open('copy.jpg', 'wb') as f:
    f.write(data)
```

---

## 💡 六、檔案存在檢查

使用 `os` 模組檢查檔案是否存在：

```python
import os

if os.path.exists('example.txt'):
    print('檔案存在！')
else:
    print('檔案不存在！')
```

---

## 🧭 七、綜合範例

```python
# 建立、寫入、再讀取
with open('demo.txt', 'w') as f:
    f.write('Line 1\n')
    f.write('Line 2\n')

with open('demo.txt', 'r') as f:
    for line in f:
        print(line.strip())
```

---

✅ **重點整理：**

1. 使用 `with open()` 可自動管理檔案關閉。
2. 文字模式與二進位模式的區別取決於檔案類型。
3. 寫入模式 `'w'` 會覆蓋原檔案，附加模式 `'a'` 不會。
4. 使用 `seek()` 和 `tell()` 可精準控制讀寫位置。
