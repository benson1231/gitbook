# Snakemake 工作流程撰寫說明

Snakemake 使用 `Snakefile` 撰寫整個工作流程，其語法靈感來自 GNU Make，並透過規則（rule）來定義如何從輸入檔案生成輸出檔案。每個 rule 可包含輸入、輸出、資源限制、執行方式等元素，並可藉由萬用字元（wildcards）實現高度自動化與重複使用。

---

## 🧾 Snakefile 基本文法與結構

### 📌 一般語法結構：

```python
rule <名稱>:
    input:
        "輸入檔案1",
        "輸入檔案2"
    output:
        "輸出檔案"
    params:
        key1="value1",
        key2="value2"
    threads: 4                     # 指定使用 CPU threads 數量
    resources:
        mem_mb=4000               # 限定記憶體資源（單位 MB）
    conda:
        "envs/tool.yaml"          # 使用 conda 環境
    shell:
        "<shell 指令>"            # 使用 shell 指令進行處理
```

### 📌 萬用字元範例（wildcard）：

```python
rule fastqc:
    input:
        "raw/{sample}.fastq.gz"
    output:
        "qc/{sample}_fastqc.html"
    shell:
        "fastqc {input} -o qc"
```

此例中，`{sample}` 為萬用字元，Snakemake 將自動比對所有符合條件的輸入檔案並建立依賴。

---

## 🔗 特殊宣告說明

### `workdir:`

設定 workflow 的工作目錄：

```python
workdir: "/home/user/project"
```

### `include:`

引入其他 Snakefile：

```python
include: "rules/align.smk"
```

### `configfile:`

讀取 YAML 格式的設定檔：

```python
configfile: "config.yaml"
```

可透過 `config["sample"]` 方式在 rule 內部存取。

---

## 🔐 模組與重用 (模組化設計)

### `module:`

載入模組化的子流程：

```python
module foo:
    snakefile: "modules/foo/Snakefile"
    config:    "modules/foo/config.yaml"
```

### `use rule from` 語法：

從模組中導入某個 rule：

```python
use rule trim_reads from foo as trim_fastq
```

---

## ✅ 強制最低 Snakemake 版本

若需要限定最低 Snakemake 版本，可在開頭加上：

```python
from snakemake.utils import min_version
min_version("9.8")
```

---

## 📚 延伸閱讀

* 文法詳細說明：[https://snakemake.readthedocs.io/en/stable/snakefiles/writing\_snakefiles.html](https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html)
* Module 使用教學：[https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html)
* Wildcards 使用技巧：[https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards)
