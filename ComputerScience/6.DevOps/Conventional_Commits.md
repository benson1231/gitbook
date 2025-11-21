# Conventional Commits 簡介

## 什麼是 Conventional Commits？

Conventional Commits 是一種 **提交訊息規範**，它透過簡單一致的格式，讓 commit 訊息更有結構性。這能幫助專案更容易產生變更日誌（CHANGELOG）、進行自動化版本管理、以及提升協作效率。

官方定義：[https://www.conventionalcommits.org/](https://www.conventionalcommits.org/)

---

## 基本格式

```
<type>(optional scope)!: <description>

[optional body]
[optional footer(s)]
```

### 例子

* `feat: add user login API`
* `fix(auth): correct password hash validation`
* `docs: update README with setup instructions`
* `chore(deps): bump snakemake to v9.1.10`
* `feat(api)!: change authentication method`

---

## 常見的 type 類型

* **feat**: 新功能 (feature)
* **fix**: 修正錯誤 (bug fix)
* **docs**: 只有文件變更
* **style**: 程式碼格式修正（不影響邏輯，例如縮排、空白）
* **refactor**: 重構（既不是修 bug，也不是加功能）
* **perf**: 改善效能的程式碼變更
* **test**: 新增或修正測試
* **chore**: 其他雜項（建置流程、工具設定、依賴套件更新）

---

## scope（可選）

用括號標示影響範圍，讓 commit 更精準：

* `fix(api): handle null response`
* `feat(ui): add dark mode toggle`
* `feat(core)!: refactor database schema`

---

## footer（選用）

用來放 **BREAKING CHANGE** 或 **issue reference**。

### BREAKING CHANGE

```
feat(auth)!: switch password hashing algorithm

BREAKING CHANGE: All existing user passwords must be reset.
```

### Issue reference

```
fix: correct typo in workflow config

Closes #123
```

---

## 如何影響版本號 (Semantic Versioning)

Conventional Commits 通常和 **Semantic Versioning (SemVer)** 搭配使用：

* **fix** → **Patch Release** (x.y.Z)

  * 例：`fix: correct typo in config` → 1.0.0 → 1.0.1

* **feat** → **Minor Release** (x.Y.z)

  * 例：`feat: add new API endpoint` → 1.0.1 → 1.1.0

* **BREAKING CHANGE 或 feat()!** → **Major Release** (X.y.z)

  * 例：`feat(api)!: update authentication method` → 1.1.0 → 2.0.0

### 總結

* `fix` = bug 修正 → patch 增加
* `feat` = 新功能 → minor 增加
* `BREAKING CHANGE` 或 `!` 標記 = 相容性破壞 → major 增加

---

## 為什麼要用 Conventional Commits？

* ✅ **一致性**：大家的 commit 訊息風格統一
* ✅ **自動化**：工具可以自動產生 CHANGELOG 或決定版本號 (Semantic Release)
* ✅ **溝通清晰**：快速理解每個 commit 的目的
* ✅ **支援自動化工作流**：像是 `release-please`、`semantic-release` 等工具

---

## 實務建議

1. **保持簡短**：標題通常不超過 72 字元
2. **英文動詞開頭**：統一風格（e.g., add, fix, update）
3. **描述要有意義**：避免 `update code` 這種模糊訊息

---

## 延伸閱讀

* 官方規範文件：[https://www.conventionalcommits.org/](https://www.conventionalcommits.org/)
* Semantic Release 工具：[https://semantic-release.gitbook.io/](https://semantic-release.gitbook.io/)
