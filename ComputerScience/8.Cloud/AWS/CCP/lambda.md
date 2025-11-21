# AWS Lambda 教學筆記

AWS **Lambda** 是一個無伺服器（Serverless）運算服務，允許開發者在不管理伺服器的情況下運行程式碼，並依執行次數與資源使用量計費，詳見[官方文檔](https://docs.aws.amazon.com/lambda/latest/dg/welcome.html)。

## 1. Lambda 特色

* **無伺服器管理**：無需啟動或維護 EC2 實例。
* **事件觸發運行**：支援來自 S3、DynamoDB、API Gateway、CloudWatch 等多種事件源。
* **自動擴展**：根據需求自動同時處理多個事件。
* **按使用付費**：依實際運行時間（毫秒級）與記憶體大小收費。
* **多語言支援**：支援 Python、Node.js、Java、Go、C# 等程式語言。

## 2. Lambda 架構概念

1. **Function（函數）**：部署的程式邏輯單位。
2. **Trigger（觸發器）**：啟動 Lambda 的事件來源，如 S3 上傳檔案或 API Gateway 請求。
3. **Execution Role（執行角色）**：Lambda 執行時使用的 IAM 角色，用於存取其他 AWS 資源。
4. **Layers（層）**：共用函式庫或依賴套件的封裝，供多個 Lambda 重複使用。
5. **Concurrency（併發）**：Lambda 同時處理多事件的能力，可設置並發上限。

## 3. 建立 Lambda 函數流程

1. 進入 **AWS Lambda Console** → 點擊 **Create function**
2. 選擇 **Author from scratch** 或使用範本
3. 選擇程式語言（如 Python 3.12）並設定 Execution Role
4. 上傳程式碼或在內建編輯器撰寫
5. 設定觸發器，例如 S3、API Gateway 或 EventBridge
6. 儲存並測試 Lambda 函數

## 4. 範例程式碼（Python）

```python
def lambda_handler(event, context):
    print("Received event:", event)
    return {
        'statusCode': 200,
        'body': 'Hello from Lambda!'
    }
```

## 5. 常見應用場景

* 後端 API：搭配 API Gateway 建立 Serverless REST API
* 事件驅動處理：S3 檔案上傳觸發影像壓縮或格式轉換
* 定時任務：使用 EventBridge 定期觸發 Lambda 處理資料
* 資料串接：整合 DynamoDB、S3、SNS/SQS 建立資料流管線

## 6. 成本控制建議

* 適用於不需長時間運算的輕量工作負載
* 依事件觸發與執行時間計費，避免閒置資源浪費
* 對頻繁觸發或長運行程式可考慮搭配 Step Functions 或 ECS

---

這份筆記適合初學者快速理解 AWS Lambda 的核心概念與使用流程，可依需求補充進階內容，如 Provisioned Concurrency、Lambda\@Edge 或與 Step Functions 的整合。
