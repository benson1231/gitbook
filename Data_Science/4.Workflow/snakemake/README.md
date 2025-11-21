# 常見的 Snakemake 專案目錄結構

以下為一個模組化、可重用且適合版本控管與 CI/CD 的 Snakemake 專案架構建議：

```bash
my-snakemake-pipeline/
├── Snakefile                  # 主要流程入口，可 include 或 use 模組
├── config.yaml                # 使用者參數設定檔（建議放常數與樣本資訊）
├── envs/                      # Conda 環境定義檔
│   └── tool.yaml              # 每個工具一個環境
├── rules/                    # 各步驟規則拆分（使用 include）
│   ├── qc.smk
│   └── align.smk
├── modules/                  # Snakemake 模組（使用 use rule from）
│   └── fastp/
│       └── Snakefile
├── scripts/                  # 外部執行用 script（bash, R, Python）
│   └── helper.py
├── resources/                # 參考檔案（如 GTF, BED）
│   └── hg38.gtf
├── data/                     # 測試或示範資料
│   └── test_R1.fastq.gz
├── results/                  # 輸出資料目錄（可忽略於 git）
│   └── logs/
├── logs/                     # 執行 log 或錯誤輸出
├── .gitignore                # 忽略資料夾，如 results/, .snakemake/
├── README.md                 # 說明與執行方式文件
├── LICENSE                   # 授權條款（建議使用 MIT 或 GPL）
├── .github/                  # GitHub Actions CI（可選）
│   └── workflows/
│       └── test.yaml
└── docs/                     # 說明文件、圖表、報告（可搭 GitHub Pages）
    └── usage.md
```

---

## ✅ 建議事項

* 使用 `rules/` 管理各步驟邏輯，避免單一 Snakefile 過大。
* `envs/` 下每個 tool 一個 `.yaml`，方便 conda 管理與版本控管。
* `config.yaml` 儲存 sample ID、參數、路徑等變數，利於重現與部署。
* `results/` 應從 Git 中排除，保留資料夾結構供本地或遠端執行。
* 搭配 GitHub Actions 自動測試流程可加強 reproducibility。

---

如需進階設計（如 profile、cluster 執行、Snakemake 模組分發），建議搭配 [snakemake-wrapper](https://github.com/snakemake/snakemake-wrappers) 或 [cookiecutter-snakemake](https://github.com/audreyfeldroy/cookiecutter-snakemake) 架構範本使用。
