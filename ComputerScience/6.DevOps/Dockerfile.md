# Dockerfile

本文件整理了 Dockerfile 的基本語法、常用指令與最佳實務，方便快速查閱與學習。

---

## 📝 基本語法

一個 Dockerfile 是一組指令的清單，Docker 在建置映像時會依序執行。常見語法如下：

```dockerfile
FROM <image>[:<tag>]          # 指定基底映像
WORKDIR <path>                # 設定工作目錄
COPY <src> <dest>             # 複製檔案或目錄到容器中
ADD <src> <dest>              # 類似 COPY，支援遠端 URL、壓縮檔解壓
RUN <command>                 # 在建置過程中執行指令
ENV <key>=<value>             # 設定環境變數
ARG <name>=<default>          # 建置參數（只能在建置時使用）
EXPOSE <port>                 # 文件性質，提示容器會使用的埠
CMD ["executable", "param"]   # 容器啟動時的預設指令，可被覆蓋
ENTRYPOINT ["executable"]     # 容器啟動時的主要指令，通常搭配 CMD
LABEL <key>=<value>           # 加上中繼資料
```

---

## 📦 最小範例

```dockerfile
FROM python:3.12-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
EXPOSE 8000
CMD ["python", "app.py"]
```

建置與執行：

```bash
docker build -t myapp:1.0 .
docker run -d -p 8000:8000 myapp:1.0
```

---

## 🏗️ 多階段建置

透過多階段建置，可以先用完整映像編譯，再將結果複製到輕量映像中，減少大小。

```dockerfile
# 建置階段
FROM node:20 AS build
WORKDIR /app
COPY package*.json ./
RUN npm ci
COPY . .
RUN npm run build

# 部署階段
FROM nginx:alpine
COPY --from=build /app/dist /usr/share/nginx/html
```

好處：最終映像更小、更安全。

---

## 🔧 常見技巧

### 1. 使用 `.dockerignore`

避免將不必要的檔案（如 `.git/`、`node_modules/`）打包進去：

```text
.git
node_modules
__pycache__/
*.log
```

### 2. 減少層數

把多個 RUN 指令合併：

```dockerfile
RUN apt-get update && apt-get install -y \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*
```

### 3. 指定非 root 使用者

```dockerfile
RUN useradd -m appuser
USER appuser
```

### 4. 健康檢查

```dockerfile
HEALTHCHECK --interval=30s --timeout=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1
```

---

## 📊 Dockerfile 範例

### Flask App

```dockerfile
FROM python:3.11-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .
EXPOSE 5000
CMD ["flask", "run", "--host=0.0.0.0"]
```

### R + Shiny App

```dockerfile
FROM rocker/shiny:latest
WORKDIR /srv/shiny-server/app
COPY . .
EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
```

---

## 🔒 最佳實務

* 使用輕量映像（`alpine`、`slim`）。
* 使用多階段建置降低最終映像大小。
* 避免硬編碼機密資訊，改用環境變數或 secrets。
* 使用 `--no-cache-dir` 或清理快取減少映像大小。
* 明確指定版本號（tag），避免使用 `latest`。

---

## 📚 參考資源

* Dockerfile 官方文件：[https://docs.docker.com/build/concepts/dockerfile/](https://docs.docker.com/build/concepts/dockerfile/)
* 最佳實務：[Docker Best Practices](https://docs.docker.com/develop/dev-best-practices/)

---

### ✅ 結語

Dockerfile 是建立映像的核心腳本。掌握基本語法與最佳實務，能幫助你建立更小、更快、更安全的映像。
