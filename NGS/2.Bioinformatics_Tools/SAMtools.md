# samtools 速查表

**更新日期**: 2025.8.22
**參考來源**: [https://github.com/samtools/samtools](https://github.com/samtools/samtools)

---

samtools 是處理 SAM（Sequence Alignment/Map）與 BAM（Binary Alignment/Map）格式的核心工具套件。 SAM/BAM 格式是儲存 DNA/RNA 序列比對結果的標準格式，而 samtools 提供了高效能的轉換、排序、索引與查詢功能。其主要用途包括： 
- 格式轉換（SAM ↔ BAM）
- 比對檔案的篩選與擷取特定區域
- 統計比對結果（mapping 狀態）
- 對 BAM 進行排序與建立索引 合併、拆分與壓縮比對檔 

---

### 1. samtools view

`samtools view` 用於 **SAM/BAM 格式之間的轉換與篩選**。

```bash
# SAM 轉 BAM
samtools view -b -t example.fa.fai -o example.bam example.sam.gz

# BAM 轉 SAM
samtools view -h -o example.sam example.bam

# 查看所有比對結果（不含 header）
samtools view example.bam

# 查看所有比對結果（含 header）
samtools view -h example.bam

# 只輸出 BAM 的 header
samtools view -H example.bam

# 提取 chr1 1000000-2000000 的比對結果
samtools view example.bam 1:1000000-2000000 | head

# 根據 BED 檔指定區域提取比對
samtools view -L example.bed example.bam
```

---

### 2. samtools flagstat

`samtools flagstat` 用於快速統計 BAM 檔比對結果摘要（mapping 數據）。

```bash
samtools flagstat example.bam
```

---

### 3. samtools sort

`samtools sort` 用於對 BAM 檔排序，以利後續建立索引或快速查詢。

```bash
samtools sort -o example.sorted.bam example.bam
```

---

### 4. samtools index

`samtools index` 用於對排序後的 BAM 建立索引，方便快速隨機存取。

```bash
samtools index example.sorted.bam
```

---

### 5. samtools merge

`samtools merge` 用於合併多個 BAM 檔為單一檔案。

```bash
samtools merge -o example.bam example_bam1.bam example_bam2.bam
```
