# DevOps

DevOps 是一種結合「開發（Development）」與「運維（Operations）」的文化與實踐方法，目標是提升軟體開發與部署的效率、穩定性與可持續性。

## 1. DevOps 的核心概念

### 1.1 自動化

涵蓋從程式碼提交到部署過程的自動化，包括建置、測試、部署與監控。

### 1.2 持續整合（CI）

持續地將開發者的程式碼整合到主分支，並自動進行建置與測試。

* 常用工具：GitHub Actions, Jenkins, Travis CI

### 1.3 持續部署（CD）

程式碼通過測試後，自動部署至 staging 或 production 環境。

* 常用工具：Docker, Kubernetes, Ansible

### 1.4 可觀察性（Observability）

透過日誌、指標與追蹤，了解系統在實際執行中的行為。

* 常見工具：Prometheus, Grafana, ELK stack

### 1.5 協作文化

強調開發與運維人員密切合作，透過自動化與版本控制縮短回饋週期，降低部署風險。

---

## 2. DevOps 工具鏈簡表

| 階段    | 目的        | 常用工具                            |
| ----- | --------- | ------------------------------- |
| 版本控制  | 程式碼管理     | Git, GitHub, GitLab             |
| 持續整合  | 自動建置與測試   | Jenkins, GitHub Actions, Travis |
| 容器化   | 打包與移植     | Docker, Podman                  |
| 編排    | 管理容器與資源配置 | Kubernetes, Docker Compose      |
| 部署    | 自動配置與交付   | Ansible, Helm                   |
| 監控與警報 | 系統健康檢查    | Prometheus, Grafana, ELK Stack  |

---

## 3. DevOps 的價值

* 減少人為錯誤與部署風險
* 加快開發至上線的時間（time-to-market）
* 提高系統穩定性與可回復性
* 建立透明且可追蹤的開發流程

DevOps 並非單一工具或流程，而是一套強調「文化、協作與自動化」的實踐方式，是現代軟體工程的重要支柱。
