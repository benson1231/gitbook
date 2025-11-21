## LaTeX

LaTeX 是一種強大的排版語言，常用於撰寫數學與科學文件。在 Markdown 中，許多平台（如 HackMD、Typora、Jupyter Notebook、Obsidian、GitHub Pages 等）支援 LaTeX 數學公式的渲染。

以下為常見 LaTeX 數學語法與範例，並搭配實際輸出結果：

---

### 1. 基本格式

* **行內公式**：使用 `$ ... $` 包裹

  ```md
  行內公式：$E = mc^2$
  ```

  輸出：$E = mc^2$

* **區塊公式**：使用 `$$ ... $$` 包裹

  ```md
  $$
  \int_0^\infty e^{-x} dx = 1
  $$
  ```

  輸出：

  $$
  \int_0^\infty e^{-x} dx = 1
  $$

---

### 2. 上標與下標

```md
$x^2$、$x_i$、$x_i^2$
```

輸出：$x^2$、$x_i$、$x_i^2$

---

### 3. 分數與根號

```md
$\frac{1}{2}$、$\sqrt{a^2 + b^2}$
```

輸出：$\frac{1}{2}$、$\sqrt{a^2 + b^2}$

---

### 4. 括號與自動大小

```md
$\left( \frac{a}{b} \right)^2$
```

輸出：$\left( \frac{a}{b} \right)^2$

---

### 5. 常見符號與函數

| 類別 | 指令        | 顯示                     |
| -- | --------- | ---------------------- |
| π  | `\pi`     | \$\pi\$                |
| ∞  | `\infty`  | \$\infty\$             |
| ≈  | `\approx` | \$\approx\$            |
| ≠  | `\neq`    | \$\neq\$               |
| ≤  | `\leq`    | \$\leq\$               |
| ≥  | `\geq`    | \$\geq\$               |
| ×  | `\times`  | \$\times\$             |
| ÷  | `\div`    | \$\div\$               |
| ±  | `\pm`     | \$\pm\$                |
| ∑  | `\sum`    | \$\sum\_{i=1}^n a\_i\$ |
| ∫  | `\int`    | \$\int\_0^1 f(x) dx\$  |

---

### 6. 邏輯與集合符號

```md
$A \cup B$, $A \cap B$, $A \subseteq B$
```

輸出：$A \cup B$, $A \cap B$, $A \subseteq B$

---

### 7. 矩陣表示（需平台支援 `amsmath`）

```md
$$
\begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix}
$$
```

輸出：

$$
\begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix}
$$

---

### 8. 多行對齊（如方程組）

```md
$$
\begin{align}
a + b &= c \\
d + e &= f
\end{align}
$$
```

輸出：

$$
\begin{align}
a + b &= c \\
d + e &= f
\end{align}
$$

---

### 9. 說明符號轉義

```md
\text{這是文字而非公式}、\_、\^、\\
```

輸出：\text{這是文字而非公式}、\_、\^、\\

---

### 10. LaTeX + Markdown 實例應用

```md
本文我們使用以下公式來表示能量與質量的關係：

$$
E = mc^2
$$
```

輸出：

本文我們使用以下公式來表示能量與質量的關係：

$$
E = mc^2
$$

---

掌握 LaTeX 可讓 Markdown 文件具備數學推導與專業格式，是學術與工程撰寫的重要工具。
