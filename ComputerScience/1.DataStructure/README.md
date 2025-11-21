# Data Structure

資料結構是程式設計與演算法的核心，決定資料如何在記憶體中儲存與操作。正確選擇資料結構能顯著提升程式效率與可維護性。

---

## 1. 基本資料結構

### 陣列（Array）

* 固定長度、連續記憶體配置
* 優點：快速索引（O(1)）
* 缺點：新增/刪除元素需搬移（O(n)）

```python
arr = [1, 2, 3]
print(arr[1])  # 輸出 2
```

### 鏈結串列（Linked List）

* 節點組成，每個節點包含值與下一個節點指標
* 優點：插入/刪除元素快速（O(1)）
* 缺點：索引需逐步尋訪（O(n)）

```python
class Node:
    def __init__(self, value):
        self.value = value
        self.next = None
```

### 堆疊（Stack）

* 後進先出（LIFO）
* 操作：push（壓入）、pop（彈出）

```python
stack = []
stack.append(1)
stack.pop()
```

### 佇列（Queue）

* 先進先出（FIFO）

```python
from collections import deque
queue = deque()
queue.append(1)
queue.popleft()
```

---

## 2. 樹與圖

### 二元樹（Binary Tree）

* 每個節點最多有兩個子節點
* 常見變種：二元搜尋樹（BST）、堆積樹（Heap）

```python
class TreeNode:
    def __init__(self, val):
        self.val = val
        self.left = None
        self.right = None
```

### 圖（Graph）

* 節點與邊組成，可為有向或無向圖
* 表示法：鄰接矩陣、鄰接串列
* 應用：社交網路、路徑搜尋

```python
graph = {
    'A': ['B', 'C'],
    'B': ['A', 'D'],
    'C': ['A'],
    'D': ['B']
}
```

---

## 3. 雜湊表與集合

### 雜湊表（Hash Table / Dictionary）

* 以鍵值對儲存資料，快速查找（平均 O(1)）

```python
d = {"apple": 3, "banana": 2}
print(d["apple"])
```

### 集合（Set）

* 儲存不重複元素

```python
s = set([1, 2, 2, 3])
print(s)  # {1, 2, 3}
```

---

## 4. 時間與空間複雜度（Big-O）

| 操作 | 陣列   | 鏈結串列 | 堆疊/佇列 | 雜湊表  | 搜尋樹      |
| -- | ---- | ---- | ----- | ---- | -------- |
| 插入 | O(n) | O(1) | O(1)  | O(1) | O(log n) |
| 刪除 | O(n) | O(1) | O(1)  | O(1) | O(log n) |
| 搜尋 | O(1) | O(n) | O(n)  | O(1) | O(log n) |

---

## 小結

掌握常見資料結構有助於撰寫高效能程式，並為學習演算法打下良好基礎。建議配合實作練習、LeetCode 題目與圖解資源深入理解各資料結構的適用場景與特性。
