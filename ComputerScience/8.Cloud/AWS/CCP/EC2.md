# AWS EC2 教學筆記

Amazon **Elastic Compute Cloud (EC2)** 是 AWS 提供的雲端虛擬伺服器服務，允許使用者彈性租用運算資源並根據需求自動擴展或縮減，詳見[官方文檔](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/concepts.html)。

## 1. EC2 特色

* **彈性運算**：可隨時啟動、停止或調整實例類型。
* **多種實例類型**：依照工作負載選擇，例如通用型、運算優化型、記憶體優化型或 GPU 型。
* **彈性網路**：可綁定彈性 IP、設定安全群組與子網路。
* **整合 AWS 服務**：可搭配 S3、RDS、EFS 等服務，方便資料存取。
* **付費模式多元**：提供按需、預留實例、節省方案與 Spot 定價。

## 2. 常見實例類型
[官方文檔](https://aws.amazon.com/tw/ec2/instance-types/)

| 類型               | 說明      | 使用情境            |
| ---------------- | ------- | --------------- |
| t 系列 (t3, t4g)   | 突發效能通用型 | 測試環境、小型網站       |
| m 系列 (m6i, m7i)  | 通用型     | 一般 Web 應用、後端服務  |
| c 系列 (c6g, c7g)  | 運算優化型   | 高運算需求、批次運算、科學模擬 |
| r 系列 (r6i, r7iz) | 記憶體優化型  | 資料庫、快取、記憶體運算    |
| g/p 系列 (g5, p4d) | GPU 型   | AI/ML 訓練、影像渲染   |

## 3. 啟動 EC2 實例流程

1. **登入 AWS Management Console** → 前往 **EC2** 服務
2. 點擊 **Launch Instance**
3. 選擇 **AMI (Amazon Machine Image)**，例如 Ubuntu、Amazon Linux 2
4. 選擇 **Instance Type**（如 t3.micro）
5. 設定 **Key Pair**（SSH 金鑰）
6. 設定 **Security Group**（開放必要的連接埠，如 22、80、443）
7. 選擇儲存空間與網路設定，然後啟動實例

## 4. 連線到 EC2

```bash
# 使用 SSH 連線範例
ssh -i my-key.pem ec2-user@<Public-IP>
```

> **注意**：確保 Security Group 已開放 TCP 22 連線，且 `.pem` 金鑰權限為 400。

## 5. 進階功能

* **Auto Scaling**：依流量自動增加或減少實例數量，詳見[官方文檔](https://aws.amazon.com/tw/ec2/autoscaling/)。
* **Load Balancer**：與 ALB/NLB 結合分散流量。
* **EBS 與 EFS**：彈性儲存與共享檔案系統。
* **CloudWatch**：監控實例 CPU、記憶體使用率。

## 6. 成本控制建議

* 測試用途可使用 **Free Tier t2.micro / t3.micro**
* 長期運行可考慮 **預留實例或 Savings Plans**
* 運算彈性工作可使用 **Spot Instances** 節省成本

---

這份筆記適合初學者快速了解 AWS EC2 的基本概念與使用流程，可依需求補充進階內容，例如 VPC 網路設計、AMI 自訂或使用 CloudFormation 自動化部署。
