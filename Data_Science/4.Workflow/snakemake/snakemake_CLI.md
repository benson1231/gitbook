# Snakemake CLI 常用指令整理

**Update: 2025.07.28**

Snakemake 提供強大的命令列介面 (CLI)，用來控制工作流程的執行、除錯與可視化。以下整理常見與實用的 CLI 指令。

---

## 🧭 基本執行指令

| 指令                        | 說明                |
| ------------------------- | ----------------- |
| `snakemake`               | 執行預設的 Snakefile   |
| `-s Snakefile`            | 指定 Snakefile 路徑   |
| `--cores N` 或 `-j N`      | 指定最多使用 N 個 CPU 核心 |
| `--dry-run` 或 `-n`        | 試跑，不實際執行命令        |
| `--printshellcmds` 或 `-p` | 顯示將執行的 shell 指令   |
| `--reason` 或 `-r`         | 顯示為何要執行該 rule     |
| `--directory DIR` 或 `-d`  | 切換執行目錄            |
| `--snakefile Snakefile`   | 明確指定 Snakefile 路徑 |

---

## 🗂️ 設定與參數覆寫

| 指令                                 | 說明                      |
| ---------------------------------- | ----------------------- |
| `--config key=value ...`           | 傳入 config dictionary 內容 |
| `--configfile config.yaml`         | 使用 YAML 格式的設定檔          |
| `--set-threads rule=4`             | 覆寫單一 rule 的 threads 數量  |
| `--set-resources rule:mem_mb=4000` | 覆寫 rule 使用的資源參數         |
| `--default-resources mem_mb=2000`  | 未定義資源的預設值               |

---

## 🏗️ 執行控制與錯誤容忍

| 指令                                  | 說明                         |
| ----------------------------------- | -------------------------- |
| `--rerun-incomplete`                | 重新執行不完整的 job               |
| `--keep-going`                      | 即使有 job 失敗也繼續執行其他 job      |
| `--force` 或 `--forceall`            | 強制重新執行指定/所有 rule           |
| `--forcerun rule1 rule2`            | 只重新執行特定 rule               |
| `--unlock`                          | 解鎖 workflow 狀態 (清除 lock 檔) |
| `--latency-wait N`                  | 等待檔案系統穩定時間 (秒)             |
| `--touch`                           | 觸碰輸出檔案，不執行實際命令             |
| `--rerun-triggers code input mtime` | 自訂觸發重新執行的條件                |

---

## 📦 環境與容器支援

| 指令                                      | 說明                    |
| --------------------------------------- | --------------------- |
| `--use-conda`                           | 啟用 conda 虛擬環境         |
| `--conda-prefix PATH`                   | 指定 conda 環境存放路徑       |
| `--conda-frontend conda/mamba`          | 使用 conda 或 mamba 安裝環境 |
| `--use-apptainer` / `--use-singularity` | 啟用容器執行環境              |
| `--apptainer-prefix`                    | 容器快取目錄                |
| `--apptainer-args`                      | 傳遞給容器的附加參數            |

---

## 📊 視覺化與報告

| 指令                     | 說明                           |
| ---------------------- | ---------------------------- |
| `--dag`                | 輸出 job DAG（dot/mermaid 格式）   |
| `--rulegraph`          | 輸出 rule DAG（不含檔案資訊）          |
| `--filegraph`          | 輸出 rule 與 input/output 的 DAG |
| `--summary`            | 顯示每個產出檔案的狀態摘要                |
| `--detailed-summary`   | 顯示輸出檔案、command 等細節           |
| `--report report.html` | 產生完整可讀的 HTML 報告              |
| `--lint`               | 檢查 Snakefile 可讀性與建議改進        |

---

## 🔁 大型工作流程支援

| 指令                   | 說明                        |
| -------------------- | ------------------------- |
| `--batch rule=1/3`   | 僅執行 rule 的第 1/3 批次輸入      |
| `--until target`     | 僅執行到指定 rule/檔案為止          |
| `--omit-from target` | 排除指定 rule/檔案及其 downstream |

---

## 📁 清理與狀態

| 指令                     | 說明                   |
| ---------------------- | -------------------- |
| `--cleanup-metadata`   | 移除檔案的 metadata 與版本資訊 |
| `--cleanup-shadow`     | 移除殘留的 shadow 資料夾     |
| `--delete-all-output`  | 刪除 workflow 產出的所有檔案  |
| `--delete-temp-output` | 僅刪除 temp() 指定的暫存檔    |
| `--keep-incomplete`    | 不刪除失敗 job 的不完整輸出     |
| `--list-untracked`     | 顯示 workflow 未使用的檔案   |

---

## 🧪 單元測試與除錯

| 指令                      | 說明                         |
| ----------------------- | -------------------------- |
| `--generate-unit-tests` | 自動產生每個 rule 的測試案例          |
| `--debug-dag`           | 顯示 DAG 中推論的 wildcard 與 job |
| `--verbose`             | 顯示除錯資訊                     |
| `--show-failed-logs`    | 顯示失敗 job 的 log             |

---

## 🧵 Profiles 支援

| 指令                    | 說明                       |
| --------------------- | ------------------------ |
| `--profile myprofile` | 指定全域配置檔，支援 cluster/cloud |
| `--workflow-profile`  | 指定工作流程專屬的設定檔             |

---

## 📚 參考連結

* CLI 文件：[Snakemake CLI Docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
* Profiles 範例：[snakemake-profiles GitHub](https://github.com/snakemake-profiles/doc)
* 所有 flags 一覽：執行 `snakemake -h`

