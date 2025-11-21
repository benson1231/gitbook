# Unit Testing

在軟體開發中，「單元測試」（Unit Testing）是用來檢驗個別函式或模組是否正確的重要機制。Python 提供內建的 `assert` 關鍵字以及 `unittest` 測試框架來輔助開發者自動化測試流程。

---

## 一、`assert`：簡單斷言語句

`assert` 是 Python 內建的關鍵字，用於快速檢查程式邏輯是否正確，條件為 `False` 時會拋出 `AssertionError`。

### 語法：

```python
assert 條件, "錯誤訊息（可選）"
```

### 範例：

```python
def divide(x, y):
    return x / y

assert divide(6, 2) == 3
assert divide(10, 5) == 2, "除法計算錯誤"
```

適合用於簡單測試或原型開發階段的自我檢查。

---

## 二、`unittest`：結構化的測試框架

Python 標準函式庫中的 `unittest` 模組提供完整的單元測試支援，可定義測試類別與多個測試方法。

### 基本結構：

```python
import unittest

def add(a, b):
    return a + b

class TestAddFunction(unittest.TestCase):

    def test_add_positive(self):
        self.assertEqual(add(2, 3), 5)

    def test_add_zero(self):
        self.assertEqual(add(0, 0), 0)

    def test_add_negative(self):
        self.assertEqual(add(-1, -1), -2)

if __name__ == '__main__':
    unittest.main()
```

### 常用方法一覽：

| 方法名稱                   | 說明                   |
| ---------------------- | -------------------- |
| `assertEqual(a, b)`    | 驗證 `a == b`          |
| `assertNotEqual(a, b)` | 驗證 `a != b`          |
| `assertTrue(x)`        | 驗證 `bool(x) 為 True`  |
| `assertFalse(x)`       | 驗證 `bool(x) 為 False` |
| `assertIs(a, b)`       | 驗證 `a is b`          |
| `assertIsNone(x)`      | 驗證 `x is None`       |
| `assertIn(a, b)`       | 驗證 `a in b`          |
| `assertRaises()`       | 測試是否拋出特定例外           |

---

## 三、進階應用技巧

### 測試前後動作：

```python
class TestExample(unittest.TestCase):

    def setUp(self):
        print("建立測試前的環境")

    def tearDown(self):
        print("清理測試後的資源")
```

### 模擬依賴（Mocking）：

```python
from unittest.mock import patch

@patch('module.function')
def test_with_mock(mock_func):
    mock_func.return_value = 10
    assert module.function() == 10
```

---

## 四、進階測試主題

### Parameterizing Tests（參數化測試）

使用迴圈或外部資料批次測試多組輸入：

```python
class TestMultiply(unittest.TestCase):

    def test_multiples(self):
        test_cases = [(2, 3, 6), (0, 10, 0), (-1, 5, -5)]
        for a, b, expected in test_cases:
            with self.subTest(a=a, b=b):
                self.assertEqual(a * b, expected)
```

### Test Fixtures（測試樣板）

`setUp()` 和 `tearDown()` 是常見樣板方法，在每個測試前後執行初始化與清理。

### Skipping Tests（跳過測試）

可使用 `@unittest.skip()` 或在條件下使用 `self.skipTest()`：

```python
class TestCondition(unittest.TestCase):

    @unittest.skip("跳過這個測試")
    def test_skip(self):
        self.assertEqual(1, 1)

    def test_maybe_skip(self):
        if True:  # 替換為實際條件
            self.skipTest("不符合條件，跳過")
        self.assertEqual(2, 2)
```

---

## 五、撰寫良好測試的原則

* 函數名稱以 `test_` 開頭
* 每個邏輯分支或邊界情況應有測試
* 單元測試應獨立，不依賴順序或外部狀態
* 測試資料固定、可重現（可配合 `fixtures`）

---

## 六、結語

單元測試與斷言是保障程式正確性與可維護性的基礎工具。小至 assert、大至完整的測試框架，都是提升開發品質不可或缺的工具。
