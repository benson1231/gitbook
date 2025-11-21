# AWS ECS & EKS 教學筆記

AWS 提供兩種主要的容器管理服務：**Elastic Container Service (ECS)** 與 **Elastic Kubernetes Service (EKS)**，分別用於在雲端部署與管理容器化應用，詳見[官方文檔](https://aws.amazon.com/tw/containers/)。

## 1. ECS 與 EKS 的差異

| 服務      | 架構基礎                 | 適用場景                    | 優點                  |
| ------- | -------------------- | ----------------------- | ------------------- |
| **ECS** | AWS 自研的容器編排服務        | 想快速使用 AWS 原生容器服務        | 深度整合 AWS 服務，簡單易上手   |
| **EKS** | 托管的 Kubernetes (K8s) | 已熟悉 Kubernetes 或需多雲/混合雲 | 完整 Kubernetes 生態系支援 |

---

## 2. ECS 特色

* **無伺服器選項 (Fargate)**：可直接執行容器，無需管理 EC2 節點。
* **整合 AWS 服務**：原生支援 CloudWatch、ALB/NLB、IAM 權限控管。
* **Task & Service 架構**：以任務（Task）和服務（Service）為單位部署容器。

### ECS 運作流程

1. 建立 **Task Definition**（定義容器映像檔、資源、環境變數）
2. 建立 **Service**（指定運行副本數、負載平衡設定）
3. 選擇運行方式：

   * **EC2**：自行管理節點，彈性高
   * **Fargate**：完全無伺服器，由 AWS 管理基礎設施

---

## 3. EKS 特色

* **托管 Kubernetes**：AWS 負責管理控制平面，高可用且自動更新。
* **原生 K8s 體驗**：可使用 `kubectl`、Helm、Ingress Controller 等工具。
* **多雲/混合雲彈性**：可整合本地或其他雲的 Kubernetes 叢集。

### EKS 運作流程

1. 建立 **EKS Cluster**（控制平面由 AWS 管理）
2. 建立或註冊 **Worker Nodes**（EC2 或 Fargate）
3. 使用 `kubectl` 或 CD 工具部署應用（Deployment/Service/Ingress）

---

## 4. 選擇建議

* **使用 ECS**：

  * 對 Kubernetes 不熟悉
  * 想快速上雲，偏好 AWS 原生整合
  * 多數應用不需複雜多雲或混合架構
* **使用 EKS**：

  * 需要 Kubernetes 生態支援（Helm、Operator、CRD）
  * 有多雲或混合雲策略
  * 需要更高的可攜性與自訂彈性

---

## 5. 成本與管理建議

* **ECS Fargate** 適合中小型、事件驅動或不固定負載
* **EKS** 適合長期運行的大型應用，但需考慮節點成本
* 建議搭配 **Auto Scaling、CloudWatch、ALB/NLB** 以最佳化運維與成本

---

這份筆記總結了 **ECS 與 EKS** 的核心概念、運作流程與選擇建議，適合作為容器上雲的入門指南。
