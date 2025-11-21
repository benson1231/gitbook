# MySQL and MySQL Workbench

本章節說明如何安裝 MySQL 資料庫伺服器與 MySQL Workbench，讓使用者可以透過圖形化介面與指令操作資料庫。

---

## 1. 安裝 MySQL Server

### Windows/macOS：

前往官方網站下載安裝程式：

* [https://dev.mysql.com/downloads/mysql/](https://dev.mysql.com/downloads/mysql/)

選擇適合的作業系統後，依照指示進行安裝。

在安裝過程中請注意：

* **設定 root 密碼**：後續連線需要用到。
* **選擇預設埠號 3306**（若無特殊需求請勿更改）。

---

## 2. 安裝 MySQL Workbench

MySQL Workbench 是 MySQL 官方提供的圖形化操作工具。

### 下載網址：

* [https://dev.mysql.com/downloads/workbench/](https://dev.mysql.com/downloads/workbench/)

安裝完成後可用於：

* 管理資料庫與表格
* 撰寫與執行 SQL 指令
* 管理使用者與權限
* 查看資料庫內容與架構

---

## 3. 驗證連線

安裝完成後可透過 Workbench：

* 建立新連線，主機 `localhost`，帳號 `root`，密碼輸入安裝時設定的值
* 測試連線成功後，即可建立資料庫與表格

若使用終端機驗證：

```bash
mysql -u root -p
```

然後輸入密碼即可進入 MySQL Shell。

---

## 小結

| 工具              | 功能說明                        |
| --------------- | --------------------------- |
| MySQL Server    | 資料庫後端主程式，儲存與處理資料            |
| MySQL Workbench | 圖形化介面，進行 SQL 管理與測試          |
| 連線設定            | 主機 `localhost`，預設埠號為 `3306` |

完成這些步驟後，即可開始使用 FastAPI 或其他應用程式與 MySQL 整合。
