# Heap Sort

`Heap Sort` 是一種基於 Heap 結構的排序演算法，其特色是總是操作在完全二元樹上，並使用 Max Heap (或 Min Heap)來達成排序目的。

---

## 一、演算法流程

1. 將數組轉成 Max Heap
2. 將根節點 (max value)和最後一個元素交換
3. 排除最後一個元素（已排完順）
4. 重新 heapify 剩下的樹
5. 重複步驟至全數排完

---

## 二、範例程式 (Python)

```python
def heapify(arr, n, i):
    largest = i
    left = 2 * i + 1
    right = 2 * i + 2

    if left < n and arr[left] > arr[largest]:
        largest = left
    if right < n and arr[right] > arr[largest]:
        largest = right

    if largest != i:
        arr[i], arr[largest] = arr[largest], arr[i]
        heapify(arr, n, largest)

def heap_sort(arr):
    n = len(arr)

    # 建立 Max Heap
    for i in range(n // 2 - 1, -1, -1):
        heapify(arr, n, i)

    # 排序
    for i in range(n - 1, 0, -1):
        arr[0], arr[i] = arr[i], arr[0]  # 最大值換至尾
        heapify(arr, i, 0)

# 範例
arr = [4, 10, 3, 5, 1]
heap_sort(arr)
print("Sorted array:", arr)
```

---

## 三、特性

| 項目    | 說明                    |
| ----- | --------------------- |
| 時間複雜度 | O(n log n) (最壞、最佳和平均) |
| 空間複雜度 | O(1) 【in-place排序】     |
| 穩定性   | 不穩定 (相等元素可能交換)        |
| 是否原地  | 是                     |

---

## 四、應用場景

* 對空間敏感的應用
* 需要 O(n log n) 性能但不依賴顏色性排序的地方

---

Heap Sort 雖然對空間有效，但總體性能經常不如 Quick Sort，但在不允許備用額外空間的場景下是非常好的選擇。
