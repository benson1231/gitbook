## ✅ Read the Docs 設定檔與建置說明

### 📁 `.readthedocs.yaml`

```yaml
version: 2

build:
  os: ubuntu-22.04
  tools:
    python: "3.10"

mkdocs:
  configuration: mkdocs.yml

python:
  install:
    - requirements: requirements.txt
```

---

### 🛠️ `mkdocs.yml`

```yaml
site_name: My Documentation
theme:
  name: readthedocs

nav:
  - Home: README.md
  - Resource: Resource.md
  - Computer Science:
      - Introduction: Computer Science/intro.md
  - NGS:
      - 16S Analysis: NGS/16s.md
      - Shotgun: NGS/shotgun.md

```

---

### 📦 `requirements.txt`

```txt
mkdocs
mkdocs-rtd-theme
```

---

### 🧪 本地測試指令（選擇性）

```bash
pip install -r requirements.txt
mkdocs build
mkdocs serve  # 預覽網站 http://localhost:8000
```
