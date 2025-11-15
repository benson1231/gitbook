# GitHub SSH 設定流程

以下為最標準、最安全、最佳實務的 GitHub SSH 設定流程。設定完成後：

* `git clone` / `git push` 不需輸入密碼或 Token
* 永久登入，不用每次驗證
* 適合 WSL、Linux、Mac、Windows 使用者

---

## 1. 產生 SSH Key

```bash
ssh-keygen -t ed25519 -C "your_email@example.com"
```

### `-C "your_email@example.com"` 是什麼？

這是一個「註解」（comment），通常放你的 GitHub Email，用來識別這把 key。

* **不影響登入機制**
* **純粹是識別用途**（讓你知道這把 key 是做什麼的）
* GitHub UI 上會顯示這段 comment

範例：

```
ed25519 AAAAC3Nz...  your_email@example.com
```

一路按 Enter 即可，會產生兩個檔：

```
~/.ssh/id_ed25519          # 私鑰（不能分享）
~/.ssh/id_ed25519.pub      # 公鑰（可以分享）
```

---

## 2. 輸出公鑰內容

```bash
cat ~/.ssh/id_ed25519.pub
```

複製整串內容。

---

## 3. 加入 GitHub SSH Keys

GitHub → **Settings → SSH and GPG keys → New SSH key**

貼上剛剛複製的公鑰。

---

## 4. 用 SSH 方式 Clone Repo（最重要）

### ✔ 正確做法（SSH）

```bash
git clone git@github.com:USERNAME/REPO.git
```

### ❌ 錯誤（HTTPS，會要求密碼或 token）

```bash
git clone https://github.com/USERNAME/REPO.git
```

---

## 5.（WSL 推薦）加入 ssh-agent 以免重開機後失效

### 啟動 ssh-agent

```bash
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
```

### 設為自動執行（建議）

```bash
echo 'eval "$(ssh-agent -s)"' >> ~/.bashrc
echo 'ssh-add ~/.ssh/id_ed25519' >> ~/.bashrc
```

---

## 6. 測試 SSH 與 GitHub 是否連線成功

```bash
ssh -T git@github.com
```

成功會看到：

```
Hi USERNAME! You've successfully authenticated, but GitHub does not provide shell access.
```

---

## 完成！

現在你所有 GitHub 動作都能：

* 不用密碼
* 不用 Token
* 更安全
* 更簡潔