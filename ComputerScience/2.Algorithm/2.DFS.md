# DFS（Depth-First Search）

`DFS`，即深度優先搜尋，是一種從起點出發，盡可能深入節點後才回退的搜尋策略。它適合用來探索所有可能的路徑、驗證某條路徑是否存在，或是進行拓樸排序與子結構搜尋等任務。

---

## 一、DFS 的核心概念

- 從根節點開始，優先走到底，再回頭探索其他分支
- 採用「遞迴」或「堆疊」記錄搜尋路徑
- 常見於圖與樹的遍歷或完全搜索任務

---

## 二、節點結構設計（以樹為例）

```python
class TreeNode:
    def __init__(self, value):
        self.value = value
        self.children = []  # 支援多叉樹
```

---

## 三、DFS 搜尋實作（傳回路徑）

```python
def dfs(node, goal_value, path=None):
    if path is None:
        path = []

    path.append(node)
    print(f"Visiting node: {node.value}")

    if node.value == goal_value:
        return path

    for child in node.children:
        result = dfs(child, goal_value, path[:])
        if result:
            return result

    return None
```

---

## 四、完整範例：樹上的 DFS

```python
# 建立樹節點
a = TreeNode("A")
b = TreeNode("B")
c = TreeNode("C")
d = TreeNode("D")
e = TreeNode("E")

# 建立層級關係
a.children = [b, c]
b.children = [d, e]

# 搜尋節點 "E"
path = dfs(a, "E")
if path:
    print("Path found:", " → ".join([node.value for node in path]))
else:
    print("Goal not found")
```

可能輸出結果：
```
Visiting node: A
Visiting node: B
Visiting node: D
Visiting node: E
Path found: A → B → E
```

---

## 五、DFS 特性與適用情境

| 特性          | 說明                                       |
|---------------|--------------------------------------------|
| 搜尋策略       | 儘可能深入節點，無解時才回退                           |
| 適用資料結構    | 樹（Tree）、圖（Graph）                            |
| 典型應用場景    | 遞迴處理、邊界追蹤、迷宮尋路、組合問題、結構遍歷               |
| 實作資料結構    | 可用遞迴（系統堆疊）或自行管理 `Stack` 結構                 |
| 記憶體使用量    | 與樹的深度成正比，適合層數深但每層節點不多的資料結構               |

---

DFS 提供一種靈活的探索方式，能夠快速深入結構底層並適用於路徑驗證、結構分析與邏輯推演。搭配剪枝或狀態紀錄可用於許多 AI、演算法競賽與遞迴型問題中。
