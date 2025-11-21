# Algorithm

演算法是解決問題的具體步驟與策略，結合資料結構能實現高效運算。常見類型如下：

---

## 1. 排序（Sorting）

### 常見演算法：

* **Bubble Sort**：兩兩比較交換，效率低 O(n²)
* **Selection Sort**：每輪選最小值 O(n²)
* **Insertion Sort**：插入排序 O(n²)，對幾乎排序資料有效
* **Merge Sort**：分治法，穩定排序 O(n log n)
* **Quick Sort**：分治法，平均 O(n log n)，最壞 O(n²)

---

## 2. 搜尋（Searching）

* **線性搜尋（Linear Search）**：從頭找起，O(n)
* **二分搜尋（Binary Search）**：排序陣列中快速搜尋，O(log n)

---

## 3. 遞迴與分治（Recursion & Divide and Conquer）

將問題拆解為子問題再解，適用於樹狀結構、排序等

```python
def factorial(n):
    return 1 if n == 0 else n * factorial(n-1)
```

---

## 4. 動態規劃（Dynamic Programming, DP）

儲存子問題解答以避免重複計算，常用於最短路徑、背包問題、最大子陣列和等

---

## 5. 貪婪演算法（Greedy）

每步選擇局部最優，期望導致全域最優，如最小生成樹、活動選擇問題

---

## 6. 回溯與 DFS / BFS

* **回溯法（Backtracking）**：用於排列組合、數獨、N 皇后等
* **深度優先搜尋（DFS）**：使用遞迴或堆疊遍歷圖或樹
* **廣度優先搜尋（BFS）**：使用佇列，自上而下或找最短路徑

---

## 7. 常見算法面試應用

* 雙指標（Two Pointers）
* 滑動視窗（Sliding Window）
* 前綴和（Prefix Sum）
* 堆與優先佇列（Heap / Priority Queue）
* 拓撲排序（Topological Sort）

---

掌握演算法設計與資料結構搭配，可提升解題效率與應對大型資料挑戰。建議透過實作、題庫訓練與圖解書籍建立直覺與策略感。
