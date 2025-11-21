# Linux Command

在生物資訊領域中，熟練使用 Linux 指令是進行高通量定序（NGS）分析不可或缺的技能。本篇筆記整理了常用的 Linux 基礎指令，涵蓋資料操作、文件處理、資料搜尋與壓縮等，提供日常分析流程的快速參考。

---

## 基礎操作指令

### echo
輸出文字到終端機或其他命令中。

```bash
# 輸出文字到螢幕上
echo hello
```

### pwd
顯示當前所在的工作目錄。

```bash
# 顯示目前所在的資料夾路徑
pwd
```

### cd
切換目錄。

```bash
cd ~           # 回到使用者家目錄
cd .           # 保持在目前目錄
cd ..          # 回到上層目錄
cd /Users/     # 以絕對路徑切換目錄
```

### ls
列出目錄下的檔案與子目錄。

```bash
ls         # 列出目錄下的所有檔案
ls -l      # 以詳細列表格式列出
ls -t      # 按修改時間排序
ls -a      # 顯示所有檔案（包括隱藏檔）
ls -lh     # 以人類易讀的格式列出（如 KB, MB）
ls -lt     # 結合時間排序與詳細列表
```

### man
查詢指令的使用手冊。

```bash
# 查看 ls 指令的說明文件
man ls
```

---

## 檔案與目錄管理

### mkdir
建立新目錄。

```bash
# 建立單一新目錄
mkdir new_dir

# 建立巢狀目錄（父目錄不存在也一併建立）
mkdir -p dir1/dir2
```

### cp
複製檔案或目錄。

```bash
# 複製 file1 為 file2
cp file1 file2

# 將 file1 複製到目錄 dir/
cp file1 dir/

# 遞迴複製整個目錄 dir1 到 dir2
cp -r dir1 dir2
```

### mv
移動或重新命名檔案/目錄。

```bash
# 將 file1 重新命名為 file2
mv file1 file2

# 將 file1 移動到目錄 dir/
mv file1 dir/
```

### rm
刪除檔案或目錄。

```bash
# 刪除單一檔案
rm file1

# 遞迴刪除整個目錄及其內容
rm -r dir/
```

### rmdir
刪除**空**目錄。

```bash
# 刪除空的目錄
rmdir empty_dir
```

---

## 文字檔案查看與處理

### more / less
分頁瀏覽大檔案。

```bash
# 使用 more 逐頁查看檔案內容
more file.txt

# 使用 less 查看檔案內容，可上下滾動
less file.txt
```

### head / tail
查看檔案的開頭或結尾內容。

```bash
# 顯示檔案的前 10 行
head file.txt

# 顯示檔案的前 5 行
head -n 5 file.txt

# 顯示檔案的最後 10 行
tail file.txt

# 顯示檔案的最後 5 行
tail -n 5 file.txt

# 持續監控檔案新增內容（如 log 檔）
tail -f file.txt
```

### wc
統計檔案的行數、字數、字節數。

```bash
wc file.txt     # 顯示行數、字數、字節數
wc -l file.txt  # 只顯示行數
wc -w file.txt  # 只顯示字數
wc -c file.txt  # 只顯示字節數
```

### cat
連接並顯示檔案內容。

```bash
# 顯示 file1.txt 的內容
cat file1.txt

# 顯示兩個檔案的內容
cat file1.txt file2.txt

# 將 file1.txt 的內容寫入 file3.txt
cat file1.txt > file3.txt
```

---

## 資料流與資料處理

### 標準輸入、輸出、錯誤 (stdin, stdout, stderr)

```bash
# 將檔案內容作為指令輸入 (command為任何指令，如cat file.txt)
command < input.txt

# 將指令輸出寫入檔案
command > output.txt

# 將錯誤訊息寫入檔案
command 2> error.txt

# 將輸出追加到檔案末尾
command >> output.txt
```

### 管線 ( | )

```bash
# 將一個指令的輸出，作為另一個指令的輸入
command1 | command2
```

### sort / uniq
排序並處理重複行。

```bash
sort file.txt     # 文字排序
sort -r file.txt  # 反向排序
sort -n file.txt  # 數字排序
sort -u file.txt  # 排序後去除重複行

uniq file.txt     # 去除重複行（需先 sort）
uniq -c file.txt  # 顯示每行出現次數
```

### cut
從檔案中抽取指定欄位。

```bash
cut -f 1,3 file.txt           # 以 tab 分隔，選取第 1 和第 3 欄
cut -d ',' -f 2 file.csv      # 指定逗號作為分隔符，選取第 2 欄
```

### comm
比較兩個排序過的檔案。

```bash
comm file1.txt file2.txt     # 顯示共同與差異行
comm -1 file1.txt file2.txt  # 只顯示 file2.txt 特有的行
```

---

## 壓縮與解壓縮

### gzip / gunzip

```bash
gzip file.txt       # 壓縮成 file.txt.gz
gunzip file.txt.gz  # 解壓縮 file.txt.gz
```

### bzip2 / bunzip2

```bash
bzip2 file.txt        # 壓縮成 file.txt.bz2
bunzip2 file.txt.bz2  # 解壓縮 file.txt.bz2
```

### tar

```bash
# 建立 tar 打包檔案
tar -cvf archive.tar file1 file2

# 建立並 gzip 壓縮 tar 檔案
tar -czvf archive.tar.gz file1 file2

# 解壓 tar 檔案
tar -xvf archive.tar

# 解壓 tar.gz 壓縮檔
tar -xzvf archive.tar.gz

# 查看 tar 檔案內容
tar -tvf archive.tar
```

### zcat

```bash
# 直接查看壓縮檔案內容
zcat file.gz
```

---

## 進階搜尋與管理

### 通配符 (Wildcard)

```bash
ls *.txt         # 列出所有 txt 檔
ls file?.txt     # 匹配 file1.txt, fileA.txt 等
ls file[123].txt # 匹配 file1.txt, file2.txt, file3.txt
ls file[!123].txt # 排除 1,2,3 開頭的檔案
```

### 正則表達式 (Regex)

```bash
grep "f.l" file.txt     # 任意單字母匹配
grep "^start" file.txt  # 匹配以 start 開頭的行
grep "end$" file.txt    # 匹配以 end 結尾的行
grep "ab*c" file.txt    # 匹配 ac、abc、abbc 等
grep "[aeiou]" file.txt # 匹配元音字母行
grep "\." file.txt     # 匹配包含 '.' 的行
```

### grep

```bash
grep "pattern" file.txt       # 搜尋指定文字

grep -i "pattern" file.txt    # 忽略大小寫

grep -n "pattern" file.txt    # 顯示行號

grep -v "pattern" file.txt    # 顯示不包含文字的行

grep -r "pattern" /path/to/dir # 遞迴搜尋目錄

grep -E "pattern1|pattern2" file.txt # 多模式搜尋
```

### find

```bash
find /path/to/dir -name "filename"    # 依名稱搜尋
find /path/to/dir -type f              # 搜尋檔案
find /path/to/dir -type d              # 搜尋目錄
find /path/to/dir -size +1M            # 找出大於 1MB 的檔案
find /path/to/dir -mtime -7            # 最近 7 天內修改的檔案
find /path/to/dir -name "*.log" -exec rm {} \; # 刪除所有 .log 檔
```

### chmod

```bash
chmod 755 file.txt   # rwxr-xr-x
chmod 644 file.txt   # rw-r--r--
chmod u+x file.txt   # 增加使用者執行權限
chmod g-w file.txt   # 移除群組寫入權限
chmod o+r file.txt   # 增加其他人讀取權限
```

---

## 資料傳輸工具

### curl

```bash
curl -O https://example.com/file.txt    # 下載文件
curl -X POST -d "key=value" https://example.com/api   # 發送 POST 請求
curl -I https://example.com              # 查看 HTTP 頭訊息
curl -H "Authorization: Bearer token" https://example.com/api # 加入授權標頭
```

### wget

```bash
wget https://example.com/file.txt   # 下載文件
wget -c https://example.com/file.txt # 斷點續傳
wget -r https://example.com         # 遞迴下載整個網站
wget --limit-rate=200k https://example.com/file.txt # 限速下載
```

---

這份指令筆記涵蓋了生物資訊日常分析中最常用的 Linux 工具與技巧，未來可根據需求持續擴充與深化。

