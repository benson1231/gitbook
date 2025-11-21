# GitHub Flow 模型教學與實作流程

GitHub Flow 是一種簡潔、現代化的 Git 分支工作流程，特別適合持續整合（CI）與持續部署（CD）的敏捷開發團隊。

---

## 🔀 GitHub Flow 概念簡介

| 分支類型       | 說明                                                      |
| ---------- | ------------------------------------------------------- |
| `main`     | 永遠維持可部署狀態，為生產環境基準分支                                     |
| feature 分支 | 每個功能或修正由此分支開始，從 `main` 分出，完成後透過 Pull Request 合併回 `main` |

特點：

* 不使用 `develop` 分支
* 發佈與合併機制以 Pull Request 為主
* 每個新功能都從 `main` 分出，並盡快合併回去（保持精簡）

---

## 🛠️ 實作流程步驟

### 1. 從 `main` 建立新功能分支

```bash
git switch main
git pull origin main
git switch -c feature/signup-form
```

### 2. 開發並提交變更

```bash
git add .
git commit -m "Add signup form layout"
git push -u origin feature/signup-form
```

### 3. 發送 Pull Request（PR）

* 在 GitHub 上建立 Pull Request
* 指定 reviewer，進行程式碼審查
* 可自動觸發 CI/CD 流程（例如 GitHub Actions）

### 4. 通過審查並合併至 `main`

* 合併方式可選：`Merge commit`、`Squash and merge`、`Rebase and merge`

### 5. 自動或手動部署 `main`

* CI/CD 工具會在 PR 合併後自動部署到 staging/production 環境

---

## 📌 GitHub Flow 優點

* 適合快速發佈與 DevOps 團隊
* 每個 PR 都是單一功能，易於追蹤與回溯
* 簡化分支結構（僅需 `main` 與功能分支）

## ⚠️ 注意事項

* 不適合一次開發大量功能的長週期專案
* 所有變更都應透過 Pull Request 管理
* 搭配 CI/CD 工具使用最佳（如 GitHub Actions, Vercel, Netlify）

---

GitHub Flow 是一種輕量級但高度效率的工作流程，特別適合部署頻繁、快速疊代的現代化 Web 專案。
