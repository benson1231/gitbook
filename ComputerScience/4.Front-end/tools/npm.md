# NPM（Node Package Manager）介紹

`npm` 是 Node.js 的套件管理工具，能讓開發者快速安裝、分享與管理 JavaScript 函式庫與工具，是 Node.js 專案開發不可或缺的一環。

## 安裝與初始化

### 確認是否已安裝 npm

```bash
node -v
npm -v
```

### 初始化專案（產生 `package.json`）

```bash
npm init
# 或快速建立：
npm init -y
```

## 套件操作指令

### 安裝套件（local）

```bash
npm install <package-name>
# 範例：npm install express
```

### 全域安裝（global）

```bash
npm install -g <package-name>
# 範例：npm install -g nodemon
```

### 安裝開發環境使用套件（devDependencies）

```bash
npm install --save-dev <package-name>
# 範例：npm install --save-dev jest
```

### 升級套件

```bash
npm update <package-name>
```

### 移除套件

```bash
npm uninstall <package-name>
```

### 查看已安裝套件

```bash
npm list
npm list --depth=0    # 僅列出第一層套件
```

## 套件資訊查詢

```bash
npm info <package-name>
```

## 清除快取（避免安裝錯誤）

```bash
npm cache clean --force
```

## 套件版本鎖定與 package-lock.json

* `package.json`：記錄專案所需套件與其版本範圍
* `package-lock.json`：鎖定實際安裝的精確版本，確保團隊一致性

## 實用工具

* `npx <command>`：執行套件而不需全域安裝
* `.npmrc`：自訂 npm 設定檔（例如 registry、proxy）

---

NPM 提供了強大的套件管理功能，能協助 JavaScript 與 Node.js 專案高效協作與版本控制。建議熟悉常用指令與 `package.json` 結構，有助於日後維護與部署。
