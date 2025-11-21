# BFS（Breadth-First Search）

`BFS`，即廣度優先搜尋，是一種從起點開始，逐層探索資料結構（如樹或圖）的搜尋策略。它適合找出**最短路徑**或**最淺層符合條件的節點**，在實作上通常使用 **佇列（Queue）** 管理待搜尋節點。

---

## 一、BFS 的核心概念

- 從根節點開始，依層級順序訪問所有鄰接節點
- 優先探索「距離根最近」的節點
- 每訪問一個節點，就將它的所有子節點加入佇列中等待處理

---

## 二、節點結構設計（以樹為例）

```python
class TreeNode:
    def __init__(self, value):
        self.value = value
        self.children = []  # 可擴展為多叉樹結構
```

---

## 三、BFS 搜尋實作（傳回路徑）

```python
from collections import deque

def bfs(root_node, goal_value):
    path_queue = deque()
    path_queue.appendleft([root_node])  # 初始路徑包含根節點

    while path_queue:
        current_path = path_queue.pop()
        current_node = current_path[-1]
        print(f"Searching node with value: {current_node.value}")

        if current_node.value == goal_value:
            return current_path

        for child in current_node.children:
            new_path = current_path[:]
            new_path.append(child)
            path_queue.appendleft(new_path)

    return None  # 未找到則回傳 None
```

---

## 四、完整範例：樹上的 BFS

```python
# 建立樹節點
a = TreeNode("A")
b = TreeNode("B")
c = TreeNode("C")
d = TreeNode("D")
e = TreeNode("E")

# 建立節點之間的層級關係
a.children = [b, c]
b.children = [d, e]

# 呼叫 BFS 搜尋目標節點 "E"
path = bfs(a, "E")
if path:
    print("Path found:", " → ".join([node.value for node in path]))
else:
    print("Goal not found")
```

輸出結果：
```
Searching node with value: A
Searching node with value: B
Searching node with value: C
Searching node with value: D
Searching node with value: E
Path found: A → B → E
```

---

## 五、BFS 特性與適用情境

| 特性          | 說明                                     |
|---------------|------------------------------------------|
| 搜尋策略       | 逐層訪問，優先處理靠近根的節點                     |
| 適用資料結構    | 樹（Tree）、圖（Graph）                        |
| 典型應用場景    | 最短路徑、語意層級探索、AI 狀態空間搜尋、網路傳播等       |
| 實作資料結構    | 需用 `Queue` 儲存待探索節點路徑                     |
| 記憶體使用量    | 與層寬有關（最壞情況可能較 DFS 佔更多記憶體）           |

---

BFS 是一種結構明確、遍歷全面的搜尋方法，特別適合在資料階層明確、目標接近根節點的情境下使用。在處理圖形、決策樹、遊戲搜尋、自然語言語法分析等應用時非常常見。