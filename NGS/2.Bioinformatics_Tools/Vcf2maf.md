awk '!/^#/ {gsub(/^chr/, "", $1)} 1' \
  /home/benson/project/WGS-report/WGS_data/test/NGS2401027.vcf \
  > /home/benson/project/WGS-report/WGS_data/test/NGS2401027.nochr.vcf



zcat /home/benson/project/WGS-report/WGS_data/test/NGS2401027.hard-filtered.vcf.gz \
  | awk '!/^#/ {gsub(/^chr/, "", $1)} 1' \
  > /home/benson/project/WGS-report/WGS_data/test/NGS2401027.nochr.vcf


# 取出 header
grep "^#" /home/benson/project/WGS-report/WGS_data/test/NGS2401027.nochr.vcf > header.vcf

# 取出前 5000 筆非註解（變異）行
grep -v "^#" /home/benson/project/WGS-report/WGS_data/test/NGS2401027.nochr.vcf | head -n 5000 > body.vcf

# 合併為測試用 VCF
cat header.vcf body.vcf > /home/benson/project/WGS-report/WGS_data/test/NGS2401027.head5000.vcf


awk 'BEGIN{FS=OFS="\t"} /^#/ || $1 ~ /^(MT|[1-9]|1[0-9]|2[0-2]|X|Y)$/ {print}' \
  /home/benson/project/WGS-report/WGS_data/test/NGS2401027.head5000.vcf \
  > /home/benson/project/WGS-report/WGS_data/test/NGS2401027.head5000.primary.vcf
