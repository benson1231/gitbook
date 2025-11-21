# AWS 網路核心服務概覽

在雲端中運行資源時，需要確保這些資源能夠彼此通信，並且能安全地與使用者或其他網路互連。AWS 提供多種核心網路服務來幫助您建置與管理雲端網路。

---

## 1. Amazon VPC（Virtual Private Cloud）

* [官方文檔](https://aws.amazon.com/tw/vpc/)
* **定義**：在 AWS 雲端中建立的專屬虛擬網路，用來啟動您的資源。
* **特色**：

  * 提供 **邏輯隔離**，可將開發、測試與生產環境分開。
  * 控制進出 VPC 與資源的網路流量。
  * 支援與其他 VPC 或本地網路的連線（例如 VPN、VPC Peering、Direct Connect）。
* **使用情境**：

  * 在不同 VPC 啟動不同階段或類型的工作負載。
  * 自訂安全性與存取控制，決定封包如何在網路中流動。

---

## 2. Amazon Route 53

* [官方文檔](https://aws.amazon.com/tw/route53/)
* **定義**：高可用、可擴展的雲端 DNS 服務。
* **核心功能**：

  1. **網域名稱註冊**（Domain Registration）
  2. **DNS 路由**（DNS Routing）, [官方文檔](https://docs.aws.amazon.com/Route53/latest/DeveloperGuide/routing-policy.html)
  3. **健康檢查**（Health Check）
* **用途**：

  * 將網域名稱轉換為 IP 位址。
  * 將使用者請求導向至在 AWS 或本地運行的應用程式。
  * 提供多種路由策略（加權、延遲、地理位置等）。

---

## 3. Elastic Load Balancing（ELB）

*  [官網](https://aws.amazon.com/tw/elasticloadbalancing/)  [官方文檔](https://docs.aws.amazon.com/Route53/latest/DeveloperGuide/routing-to-elb-load-balancer.html)
* **定義**：負載平衡服務，可將進入的流量自動分配至多個目標（如 EC2 實例、容器或虛擬應用程式）。
* **特色**：

  * 提升應用程式的高可用性與容錯能力。
  * 支援跨多個可用區（AZ）的負載分配。
  * 提供單一存取點，隱藏後端實例數量與 IP 詳細資訊。
* **使用情境**：

  * 部署 Web 或 API 應用時分散使用者流量。
  * 與 Route 53 搭配實現全球流量導向與多區域高可用架構。

---

## 4. 網路安全與流量控制

* AWS 提供多層級的流量過濾與安全控制：

  * **安全群組（Security Groups）**：控制實例層級的入站與出站流量。
  * **網路 ACL（Network ACLs）**：控制子網層級的流量過濾。
  * **Route 53 + ELB**：結合 DNS 與負載平衡保護應用前端。

---

這份筆記整理了 AWS 的核心網路服務（VPC、Route 53、ELB），並說明其用途與安全控制，作為建置雲端網路架構的入門指南。
