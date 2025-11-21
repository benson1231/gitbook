# argparse

`argparse` 是 Python 標準函式庫中用來解析命令列參數的模組，可用於構建 CLI（命令列介面）工具，讓使用者能夠在終端機中傳遞參數給 Python 腳本。

---

## 1. 基本用法

以下為基本的使用範例：

```python
import argparse

parser = argparse.ArgumentParser(description="加法計算器")
parser.add_argument("x", type=int, help="第一個數字")
parser.add_argument("y", type=int, help="第二個數字")

args = parser.parse_args()
print("結果：", args.x + args.y)
```

### 執行方式：

```bash
python script.py 3 5
# 輸出：結果： 8
```

---

## 2. 常用參數設定

| 方法                    | 功能說明                   |
| --------------------- | ---------------------- |
| `add_argument()`      | 定義命令列參數                |
| `type=int/float/str`  | 指定輸入資料型別               |
| `help="說明文字"`         | 顯示參數說明                 |
| `default=值`           | 指定參數預設值                |
| `required=True`       | 參數是否為必要                |
| `action='store_true'` | 遇到參數即設為 True（常用於 flag） |
| `choices=[...]`       | 限定參數可選擇的值範圍            |

---

## 3. 可選參數（Optional Argument）

```python
parser.add_argument("--verbose", action="store_true", help="顯示詳細輸出")
```

執行：

```bash
python script.py 3 5 --verbose
```

---

## 4. 多值輸入與列表解析

```python
parser.add_argument("--nums", nargs='+', type=int, help="整數清單")
```

執行：

```bash
python script.py --nums 1 2 3 4
```

---

## 5. 範例：進階使用

```python
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--mode", choices=["train", "test"], default="train")
parser.add_argument("--lr", type=float, default=0.001, help="學習率")
parser.add_argument("--epochs", type=int, default=10)

args = parser.parse_args()
print(f"模式: {args.mode}, 學習率: {args.lr}, 訓練輪數: {args.epochs}")
```

---

## 6. 顯示說明訊息（Help）

使用 `-h` 或 `--help` 參數會自動顯示說明文件：

```bash
python script.py -h
```

---

`argparse` 是建立靈活 Python CLI 的利器，熟練使用後可大幅提升腳本的可操作性與重用性。
