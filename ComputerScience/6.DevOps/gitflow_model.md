# Git Flow 模型介紹與實作流程

Git Flow 是一種針對軟體專案版本控制的分支管理策略，適合具有發佈週期與多人協作的專案。

---

## 🧭 Git Flow 分支架構

Git Flow 模型將版本控制分為數個固定用途的分支：

| 分支          | 用途說明                                                 |
| ----------- | ---------------------------------------------------- |
| `main`      | 穩定的正式版本，部署到生產環境                                      |
| `develop`   | 整合開發分支，所有新功能會先合併到此                                   |
| `feature/*` | 功能開發分支，從 `develop` 分出，完成後合併回 `develop`               |
| `release/*` | 發佈準備分支，從 `develop` 分出，修復錯誤與優化，合併至 `main` 與 `develop` |
| `hotfix/*`  | 緊急修補分支，從 `main` 分出，修復後合併回 `main` 與 `develop`         |

---

## 🛠️ 分支操作範例

### 建立功能分支（feature）

```bash
git checkout develop
git switch -c feature/login-page
```

### 結束功能並合併回 develop

```bash
git checkout develop
git merge feature/login-page
git branch -d feature/login-page
```

### 建立發佈分支（release）

```bash
git checkout develop
git switch -c release/v1.0.0
```

### 發佈完成後：合併 release 到 main 與 develop

```bash
git checkout main
git merge release/v1.0.0

git checkout develop
git merge release/v1.0.0
git branch -d release/v1.0.0
```

### 緊急修復（hotfix）

```bash
git checkout main
git switch -c hotfix/fix-login-bug
```

修正後：

```bash
git checkout main
git merge hotfix/fix-login-bug

git checkout develop
git merge hotfix/fix-login-bug
git branch -d hotfix/fix-login-bug
```

---

## 🔖 版本號規則（Semantic Versioning）

在 Git Flow 中，發佈分支（`release/*`）與標籤（tag）通常遵循 **語意化版本** `MAJOR.MINOR.PATCH` 格式。

```
v1.2.3
│ │ └── PATCH  修訂版：僅修補錯誤，向下相容
│ └──── MINOR  次版本：新增向下相容的新功能
└────── MAJOR  主版本：重大變更，可能破壞相容性
```

### 建立標籤（Tag）

```bash
# 在完成 release 合併到 main 之後：
git checkout main
git tag v1.2.3 -m "Release v1.2.3"
git push origin v1.2.3
```

---

## 📌 Git Flow 優點

* 分支用途明確，結構穩定
* 適合多人協作與定期版本釋出
* 支援緊急修補流程而不影響正常開發

## ⚠️ 使用建議

* 不建議用於小型或原型專案（可能太複雜）
* 若使用 CI/CD，建議搭配 `main`/`release` 設定部署規則
* 可使用 `git-flow` 工具輔助操作（如 `brew install git-flow`）

---

Git Flow 是一套清晰的分支流程標準，適合中大型開發團隊進行版本控制與協作。理解各類分支的角色與建立/合併方式，有助於提升開發流程的可控性與穩定性。
