# MD5 驗證教學

## 什麼是 MD5？
MD5（Message-Digest Algorithm 5）是一種雜湊演算法，會將任意長度的資料轉換為一個 128-bit 的「指紋」。常用來驗證檔案完整性、確保下載過程未被破壞或篡改。

---

## 常見應用場景
- 驗證下載檔案是否正確（比對 MD5 值）
- 檔案傳輸後進行一致性校驗
- 軟體發佈提供 MD5 碼供使用者驗證

---

## Linux / macOS 驗證方法

### 1. 查看檔案 MD5 碼
```bash
md5sum filename
```

### 2. 比對 `.md5` 檔案內容
```bash
md5sum -c filename.md5
```

範例內容（filename.md5）：
```
d41d8cd98f00b204e9800998ecf8427e  filename
```
執行結果：
```
filename: OK
```

---

## Windows 驗證方法
使用 PowerShell：
```powershell
Get-FileHash filename -Algorithm MD5
```

---

## Python 驗證方式（進階）
```python
import hashlib

def compute_md5(filepath):
    md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5.update(chunk)
    return md5.hexdigest()

print(compute_md5("yourfile.fastq.gz"))
```

---

## 注意事項
- MD5 速度快但安全性不高，不適合密碼或數位簽章驗證
- 僅適用於一般檔案完整性驗證

---

## 常見錯誤排查
- 檔案名對不上（請確保 md5 檔案內檔名與實際檔案一致）
- 檔案換行符不同（Linux/Windows 換行符差異可能影響校驗）
- 傳輸過程中自動解壓或修改了內容（導致 MD5 不一致）
