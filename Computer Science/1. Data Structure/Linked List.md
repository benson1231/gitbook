# Linked List

`Linked List` 是一種基本的資料結構，由一系列節點（Node）組成，每個節點包含資料與指向下一個節點的指標（pointer）。

---

## 一、鍊結串列的分類

| 類型       | 描述            |
| -------- | ------------- |
| 單向鍊結串列   | 每個節點只指向下一個節點  |
| 雙向鍊結串列   | 每個節點有前後兩個指標   |
| 循環鍊結串列   | 最後一個節點指回第一個節點 |
| 雙向循環鍊結串列 | 同時具備雙向與循環特性   |

---

## 二、基本節點結構（Python 範例）

```python
class Node:
    def __init__(self, value):
        self.value = value
        self.next = None
```

---

## 三、單向鍊結串列操作

### 1. 建立鏈表與新增節點

```python
class LinkedList:
    def __init__(self):
        self.head = None

    def append(self, value):
        new_node = Node(value)
        if not self.head:
            self.head = new_node
        else:
            current = self.head
            while current.next:
                current = current.next
            current.next = new_node
```

### 2. 顯示所有節點

```python
    def display(self):
        current = self.head
        while current:
            print(current.value, end=" → ")
            current = current.next
        print("None")
```

### 3. 插入、刪除、搜尋

```python
    def insert(self, index, value):
        new_node = Node(value)
        if index == 0:
            new_node.next = self.head
            self.head = new_node
            return
        current = self.head
        for _ in range(index - 1):
            if current:
                current = current.next
        new_node.next = current.next
        current.next = new_node

    def delete(self, value):
        current = self.head
        if current and current.value == value:
            self.head = current.next
            return
        while current.next:
            if current.next.value == value:
                current.next = current.next.next
                return
            current = current.next

    def search(self, value):
        current = self.head
        while current:
            if current.value == value:
                return True
            current = current.next
        return False
```

---

## 四、Linked List 特性

* **插入與刪除快**：操作不需移動其他元素
* **隨機訪問慢**：需從頭開始搜尋
* **動態長度**：不需預先分配空間

---

## 五、Linked List vs Array 比較

| 特性    | Linked List | Array（List） |
| ----- | ----------- | ----------- |
| 記憶體配置 | 分散          | 連續          |
| 訪問速度  | O(n)        | O(1)        |
| 插入刪除  | O(1)（已知位置）  | O(n)        |
| 空間效率  | 多一點指標空間     | 純資料         |

---

Linked List 是所有資料結構的基礎之一，尤其適合插入與刪除操作頻繁的場景，例如實作 queue、stack、hash table 或圖結構等。
