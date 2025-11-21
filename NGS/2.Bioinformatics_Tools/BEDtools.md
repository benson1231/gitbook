# bedtools 使用速查表

**更新日期**: 2025.8.22
**參考來源**: [BEDtools GitHub](https://github.com/arq5x/bedtools2)

---

### 1. BEDtools Intersect

`bedtools intersect` 用於比較兩個基因體區間檔案 (A 與 B)，找出重疊、非重疊或具體的交集關係。

#### 常用參數

* **`-wa`** : 輸出有重疊的 A 區間。
* **`-wb`** : 輸出有重疊的 B 區間。
* **`-wo`** : 輸出 A 與 B 的原始欄位，並在最後加上重疊長度 (bp)。
* **`-wao`** : 輸出 A 與 B 的原始欄位，所有 A 都保留，若與 B 無重疊，B 欄位以"."表示，最後一欄為 0。
* **`-loj`** : Left outer join，所有 A 都保留，若沒有重疊則 B 欄位填 `.`。
* **`-u`** : 只輸出與 B 有重疊的 A。
* **`-v`** : 只輸出沒有與 B 重疊的 A。
* **`-c`** : 合併計算所有 B 的重疊數，一個 A 只輸出一行。
* **`-C`** : 對每個 B 檔案分別計算，一個 A 會輸出多行（每個 B 一行）。
* **`-s`** : 要求同股 (strand)，需有 strand 欄位。
* **`-S`** : 要求反股 (opposite strand)，需有 strand 欄位。
* **`-f 浮點數`** : A 至少有多少比例需被 B 覆蓋。
* **`-F 浮點數`** : B 至少有多少比例需被 A 覆蓋。
* **`-r`** : 搭配 `-f` 與 `-F`，要求互相覆蓋比例皆滿足。
* **`-split`** : 把 BED12 的 exon blocks 或 BAM 的 spliced exons 分開來看，每個 block 視為獨立區間來算交集。
* **`-filenames`** : 在輸出中加上 B 的來源檔名。
* **`-names`** : 自訂多個 B 檔的名稱。

#### 使用範例

```bash
# 找出 A 與 B 有重疊的區間 (輸出 A)
bedtools intersect -a A.bed -b B.bed -u

# 找出 A 中完全沒有重疊的區間
bedtools intersect -a A.bed -b B.bed -v

# 同時輸出 A 與 B 的原始欄位
bedtools intersect -a A.bed -b B.bed -wa -wb

# 輸出 A、B 原始欄位並加上重疊長度
bedtools intersect -a A.bed -b B.bed -wo

# 保留所有 A，若無重疊則 overlap=0
bedtools intersect -a A.bed -b B.bed -wao

# Left outer join，若無重疊則 B 欄位為 .
bedtools intersect -a A.bed -b B.bed -loj

# 計算每個 A 與 B 的重疊次數
bedtools intersect -a A.bed -b B.bed -c

# 與上相同，但若無重疊則輸出 0
bedtools intersect -a A.bed -b B.bed -C

# 要求至少 50% 的 A 被 B 覆蓋
bedtools intersect -a A.bed -b B.bed -f 0.5 -u

# 要求 A 與 B 互相至少有 50% 的重疊比例
bedtools intersect -a A.bed -b B.bed -f 0.5 -F 0.5 -r

# 只保留同股的重疊
bedtools intersect -a A.bed -b B.bed -s

# 對 BED12 格式，以 exon 為單位計算交集
bedtools intersect -a A.bed -b B.bed -wo -split

# 使用多個 B 檔並輸出來源檔名
bedtools intersect -a A.bed -b B1.bed B2.bed -wa -wb -filenames
```

---

### 2. BAM 轉 BED

`bedtools bamtobed` 用於將 BAM 格式轉換為 BED 格式。

#### 使用範例

```bash
# 轉換 BAM 為 BED，並將 CIGAR 字串放在第七欄
bedtools bamtobed -i sample.bam -cigar > sample.bed

# 將有剪接的比對拆分為 exon 區段
bedtools bamtobed -i sample.bam -split > sample_split.bed
```

---

### 3. BED 轉 BAM

`bedtools bedtobam` 用於將 BED 格式轉換為 BAM 格式。

#### 使用範例

```bash
# 使用 genome header 檔 (染色體長度) 轉換為 BAM
bedtools bedtobam -i sample.bed -g genome_header.txt > sample.bam

# 保留 BED12 的 exon 結構並轉為 BAM
bedtools bedtobam -i sample.bed -g genome_header.txt -bed12 > sample_exon.bam
```

---

### 4. 從 BED 區間提取 FASTA 序列

`bedtools getfasta` 用於從基因組序列中擷取 BED 區間對應的 FASTA 序列。

#### 使用範例

```bash
# 從 BED 區間提取 FASTA 序列
bedtools getfasta -fi genome.fa -bed regions.bed -fo output.fasta

# 將 exon 區段擷取並合併為一條序列
bedtools getfasta -fi genome.fa -bed exons.bed -fo exons_output.fasta -split
```

---

此速查表提供常用 BEDtools 指令的精簡範例，適合作為生物資訊流程的快速參考。更多詳細說明可參考 [官方文件](https://bedtools.readthedocs.io/)。
