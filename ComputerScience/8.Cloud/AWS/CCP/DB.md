# AWS 資料庫核心服務概覽

AWS 提供多種**完全託管（Managed）、針對用途設計（Purpose-built）**的資料庫服務，幫助企業**現代化資料基礎架構**，以提高效能、彈性與安全性，詳見[官方文檔](https://aws.amazon.com/tw/products/databases/)。

---

## 1. AWS 資料庫的運行方式

* **在 Amazon EC2 上自行管理資料庫**

  * 完全自我維運，需自行安裝、配置、備份、擴容與高可用設計。
* **使用 AWS 託管資料庫服務**

  * 簡化設定與運維，由 AWS 處理備份、更新、擴容與高可用配置。
  * 降低維護成本，專注於應用程式與業務創新。

使用 **AWS Cloud 資料庫** 的好處：

* 針對用途選擇最佳資料庫引擎（Purpose-built）
* 簡化運維與自動擴展
* 高可用與高安全性
* 支援從傳統資料庫遷移、建構現代應用與擺脫舊式資料庫限制

---

## 2. AWS 資料庫類型與服務對應

| 資料庫類型                | 典型應用場景                   | 對應 AWS 服務                                      |
| -------------------- | ------------------------ | ---------------------------------------------- |
| **Relational（關聯式）**  | ERP、CRM、電商、傳統應用          | Amazon Aurora、Amazon RDS、Amazon Redshift（OLAP） |
| **Key-value（鍵值）**    | 高流量網站、電商、遊戲排行榜           | Amazon DynamoDB                                |
| **In-memory（內存型）**   | 快取、Session 管理、遊戲排行榜、地理應用 | Amazon ElastiCache、Amazon MemoryDB for Redis   |
| **Document（文件型）**    | 內容管理、目錄、用戶檔案             | Amazon DocumentDB（MongoDB 相容）                  |
| **Wide Column（寬列型）** | 高規模工業應用、車隊管理、路線最佳化       | Amazon Keyspaces（Cassandra 相容）                 |
| **Graph（圖形）**        | 欺詐偵測、社交網路、推薦系統           | Amazon Neptune                                 |
| **Time Series（時序）**  | IoT、DevOps、工業遙測          | Amazon Timestream                              |
| **Ledger（帳本）**       | 供應鏈、登記系統、金融交易            | Amazon QLDB                                    |

---

## 3. 主要資料庫服務簡介

### 3.1 Relational Database（關聯式資料庫）

* **Amazon Aurora**
  * [官方文檔](https://aws.amazon.com/tw/rds/aurora/)
  * 雲端原生、MySQL/PostgreSQL 相容
  * 支援自動擴展、Multi-AZ、跨區域複寫
  * 商業級效能與高可用，成本僅傳統 1/10

* **Amazon RDS**
  * [官方文檔](https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/Welcome.html)
  * 託管型關聯式資料庫，支援 7 種引擎（Aurora、MySQL、PostgreSQL、MariaDB、Oracle、SQL Server）
  * 提供 Multi-AZ 高可用、快照備份與自動擴容

* **Amazon Redshift**
  * [官方文檔](https://aws.amazon.com/tw/redshift/)
  * 雲端數據倉儲（OLAP），採用列式存儲，適合複雜查詢與分析
  * 適合大數據分析、商業智慧（BI）

### 3.2 NoSQL & 非關聯式資料庫

* **Amazon DynamoDB（Key-value）**
  * [官方文檔](https://aws.amazon.com/tw/dynamodb/)
  * 全託管、Serverless、單位毫秒級延遲
  * 自動擴展讀寫吞吐，適合高流量與彈性負載

* **Amazon ElastiCache / MemoryDB（In-memory）**

  * 超低延遲快取與 Session 管理
  * 相容 Redis 或 Memcached，MemoryDB 支援持久化與高可用

* **Amazon DocumentDB（Document）**

  * MongoDB 相容的雲端文件資料庫，適合 JSON 工作負載

* **Amazon Keyspaces（Wide Column）**

  * 託管 Cassandra，相容 API，適合高規模工業應用

* **Amazon Neptune（Graph）**
  * [官方文檔](https://aws.amazon.com/tw/neptune/)
  * 專為圖形資料設計，適用於社交網路、欺詐偵測、推薦系統

* **Amazon Timestream（Time Series）**

  * Serverless 時序資料庫，適合 IoT 與遙測數據分析

* **Amazon QLDB（Ledger）**

  * 提供透明、不可變且可加密驗證的交易日誌，適合金融與供應鏈系統

---

## 4. 選擇建議

* **傳統應用或交易型系統** → 關聯式（Aurora / RDS）
* **大數據分析、報表或 OLAP** → Redshift
* **高流量網站或遊戲排行榜** → DynamoDB / ElastiCache
* **文件或 JSON 內容管理** → DocumentDB
* **高連結關係資料** → Neptune
* **IoT 與時間序列數據** → Timestream
* **金融或帳本需求** → QLDB

---

這份筆記整理了 AWS **資料庫核心服務**、用途與對應場景，可作為設計雲端資料架構與選型的快速指南。
