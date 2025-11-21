# Docker

## 🐳 基本概念

* **Image（映像）**：唯讀模板，用於建立容器。
* **Container（容器）**：映像的執行個體，具有可寫層。
* **Volume**：持久化儲存或主機與容器間共享檔案。
* **Network**：容器彼此或對外的連線方式。
* **Dockerfile**：定義映像的建置步驟。

---

## 🔧 安裝與版本檢查

```bash
# 檢查版本
docker version

# 系統資訊
docker info

# Compose 版本
docker compose version
```

---

## 📦 映像管理

```bash
# 建立映像
docker build -t <image_name>:<tag> .

# 列出映像
docker images   # 或 docker image ls

# 標記映像
docker tag <image_name>:<tag> <docker_hub_user>/<repo>:<tag>

# 上傳/下載
docker push <docker_hub_user>/<repo>:<tag>
docker pull <docker_hub_user>/<repo>:<tag>

# 移除映像
docker rmi <image_name>:<tag>
```

---

## 🚀 容器管理

```bash
# 建立並啟動容器
docker run --name <container_name> -it --rm <image_name>:<tag> /bin/bash

# 背景執行
docker run -d --name <container_name> -p 8080:80 <image_name>:<tag>

# 查看容器
docker ps        # 執行中
docker ps -a     # 包含已停止

# 啟動/停止/刪除容器
docker start <container_name>
docker stop <container_name>
docker rm <container_name>

# 進入已創建容器
docker exec -it <container_name> bash
```

### `docker run` vs `docker exec`

| 指令 | 功能 | 使用時機 | 範例 |
|------|------|----------|------|
| `docker run` | 建立並啟動一個**新的容器** | 第一次從映像啟動服務，或想要一個新的環境 | `docker run -it ubuntu bash` |
| `docker exec` | 在**已存在且正在執行的容器**內執行命令 | 需要進入背景服務容器檢查或除錯 | `docker exec -it myweb bash` |

📌 小提醒：  
- 若容器背景模式啟動（`-d`），要進去查看時就用 `docker exec`。  
- 有些精簡映像（如 `alpine`）沒有 `bash`，需改用 `sh`：  
  ```bash
  docker exec -it <container_name> sh
  ```


### 常用旗標

| 旗標     | 說明   | 範例                   |
| ------ | ---- | -------------------- |
| `-it`  | 互動模式 | `-it bash`           |
| `--rm` | 自動刪除 | `--rm`               |
| `-d`   | 背景執行 | `-d`                 |
| `-p`   | 埠口對映 | `-p 8080:80`         |
| `-v`   | 掛載   | `-v $PWD/data:/data` |
| `-e`   | 環境變數 | `-e MODE=prod`       |

---

## 📂 Volume 與資料持久化

```bash
# 建立 volume
docker volume create <volume_name>

# 掛載 volume
docker run -d -v <volume_name>:/var/lib/mysql mysql

# 掛載主機目錄
docker run -d -v $PWD/site:/usr/share/nginx/html nginx
```
以上範例展示了 **兩種資料掛載方式**：

* 使用 `docker volume create` 建立的 **Volume** 由 Docker 自行管理，適合資料庫或需要持久化的資料。即使容器被刪除，Volume 仍會保留。
* 使用 `-v ./path:/container/path` 建立的 **Bind Mount** 會直接把主機目錄掛載到容器內，常用於開發環境同步程式碼或靜態檔案。

簡單來說：開發時方便修改 → 用 **Bind Mount**；生產環境需要穩定與可攜性 → 用 **Volume**。

---

## 🌐 網路

```bash
# 建立自訂網路
docker network create <network_name>

# 啟動並加入網路
docker run -d --name redis --network <network_name> redis
```
以上範例示範了如何建立自訂的 **bridge 網路** 並讓容器加入。好處是容器之間可以透過 **容器名稱** 互相通訊，而不用記 IP。這樣在專案或服務中會更方便管理。

---

## 📝 Dockerfile 建置

```dockerfile
FROM python:3.12-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .
CMD ["python", "app.py"]
```

```bash
# 建置
docker build -t <image_name>:<tag> .
# 執行
docker run --rm -p 8000:8000 <image_name>:<tag>
```

---

## 🧹 清理與診斷

```bash
# 清理未使用資源
docker system prune -f

# 清理未使用映像
docker image prune -a -f

# 查看資源使用
docker system df

# 查看容器日誌
docker logs -f <container_name>

# 即時監控
docker stats
```

---

## 🔒 最佳實務

* 使用最小化映像（`-slim` 或 `alpine`）。
* 盡量避免 root，Dockerfile 可建立專用使用者。
* 機密資訊用 `.env` 或 secrets，避免硬寫入映像。
* 合理使用 volume 持久化資料。
* 在 CI/CD 中加上建置快取與自動測試。

---

## 📚 參考資源

* Docker 官方文件：[https://docs.docker.com/manuals/](https://docs.docker.com/manuals/)
* Docker Hub 範例：[benson1231](https://hub.docker.com/u/benson1231)

---

### ✅ 結語

本指南整合了基礎指令與進階實務，能快速上手並支援專案開發。可依需求擴充至 CI/CD、自動化測試或 GPU/生物資訊等情境。


