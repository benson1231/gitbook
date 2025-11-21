# Docker Compose

本文件整理了 Docker Compose 的基本概念、常用語法與實務範例，適合用來管理多容器應用程式。

---

## 📖 基本概念

* **Docker Compose**：一個工具，用來定義和執行多容器應用程式。
* 使用 **YAML 檔** 之 `compose.yaml`，舊時命名 `docker-compose.yaml` 描述服務（service）、網路（network）、資料卷（volume）。
* 一個 `docker compose up` 就能啟動整個應用程式堆疊。

---

## 🔧 安裝與版本

```bash
# 查看版本
docker compose version
```

Docker Compose v2 已經整合到 Docker CLI 中，使用方式是：

```bash
docker compose <command>
```

而不是舊版的 `docker-compose`（建議升級）。

---

## 📑 compose.yaml 範例

### 1. Web + Redis 範例

```yaml
version: "3.9"
services:
  web:
    image: nginx:alpine
    ports:
      - "8080:80"
    volumes:
      - ./site:/usr/share/nginx/html:ro
    depends_on:
      - redis

  redis:
    image: redis:7
```

啟動：

```bash
docker compose up -d
```

---

### 2. Web + DB 範例（Flask + MySQL）

```yaml
version: "3.9"
services:
  db:
    image: mysql:8
    environment:
      MYSQL_ROOT_PASSWORD: example
      MYSQL_DATABASE: appdb
    volumes:
      - db_data:/var/lib/mysql

  web:
    build: ./web
    ports:
      - "5000:5000"
    depends_on:
      - db

volumes:
  db_data:
```

此範例示範了如何同時定義 **自建映像（web）** 與 **官方映像（mysql）**，並透過 volume 持久化資料。

---

## ⚙️ 常用指令

```bash
# 啟動服務（建議加 -d 背景執行）
docker compose up -d

# 停止並移除服務（保留 volumes）
docker compose down

# 查看服務狀態
docker compose ps

# 查看日誌
docker compose logs -f web

# 進入容器
docker compose exec web sh

# 執行一次性指令（結束後刪除容器）
docker compose run --rm web echo "hello"

# 重新建置服務
docker compose build --no-cache

# 驗證設定檔
docker compose config
```

---

## 🧩 常見設定鍵

| 鍵             | 功能              | 範例                             |
| ------------- | --------------- | ------------------------------ |
| `image`       | 指定映像            | `nginx:alpine`                 |
| `build`       | 從 Dockerfile 建置 | `build: ./app`                 |
| `ports`       | 埠口對映            | `"8080:80"`                    |
| `volumes`     | 掛載 volume 或目錄   | `./site:/usr/share/nginx/html` |
| `environment` | 環境變數            | `MYSQL_ROOT_PASSWORD: example` |
| `env_file`    | 從檔案載入環境變數       | `env_file: .env`               |
| `depends_on`  | 設定啟動順序          | `depends_on: [db]`             |
| `networks`    | 指定網路            | `networks: [appnet]`           |
| `restart`     | 重啟策略            | `restart: unless-stopped`      |

---

## 🔒 最佳實務

* 使用 `.env` 管理敏感資訊，避免寫死在 YAML。
* 資料庫一律掛載 volume，避免容器刪除後資料遺失。
* 不要過度依賴 `depends_on`，應在程式中實作「重試連線」。
* 建立專案專屬網路，讓服務能用名稱互通。
* 搭配 CI/CD，自動執行 `docker compose up -d` 部署。

---

## 📚 參考資源

* Docker 官方文件：[Docker Compose](https://docs.docker.com/compose/)
* Compose 範例集：[Awesome Compose](https://github.com/docker/awesome-compose)

> Awesome Compose 是 Docker 官方維護的範例集合，專門展示如何使用 Docker Compose 來快速建立常見應用程式架構。這些範例涵蓋前後端整合、資料庫、快取、訊息佇列以及完整應用堆疊，非常適合作為學習與實務參考。

---

### ✅ 結語

Docker Compose 是管理多容器應用的最佳工具，能讓你用一份 YAML 定義整個環境。掌握常用語法與最佳實務，可以讓開發、測試到部署都更加高效。
