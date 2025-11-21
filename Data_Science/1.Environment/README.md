# Environment

良好的環境管理是進行資料科學與機器學習專案的基礎，能夠確保實驗的可重現性、版本一致性與開發效率。以下介紹常見的環境管理工具與實踐方式。

---

## 一、Python 環境管理工具

### 1. Conda

* 支援 Python 與 R
* 管理套件與虛擬環境
* 常見於機器學習與深度學習開發

```bash
conda create -n myenv python=3.10
conda activate myenv
conda install numpy pandas scikit-learn
```

### 2. venv（標準 Python 工具）

* 適合輕量專案，搭配 pip 安裝套件

```bash
python -m venv env
source env/bin/activate  # Windows: env\Scripts\activate
pip install numpy pandas
```

---

## 二、套件管理

### 1. pip + requirements.txt

* 以文字記錄所有安裝套件與版本

```bash
pip freeze > requirements.txt
pip install -r requirements.txt
```

### 2. conda environment.yml

* YAML 格式描述環境，便於團隊共用

```yaml
name: myenv
channels:
  - defaults
  - conda-forge
dependencies:
  - python=3.10
  - numpy
  - pandas
  - scikit-learn
```

---

## 三、專案管理建議

* 使用 `.gitignore` 排除環境資料夾（如 `env/` 或 `__pycache__/`）
* 為每個專案建立獨立環境，避免套件衝突
* 善用 Jupyter Notebook 或 JupyterLab 開發原型
* 整合 VSCode、PyCharm 等 IDE 以提升工作效率

---

## 四、Docker 容器化（進階）

若需在不同機器或雲端執行一致環境，推薦使用 Docker：

```Dockerfile
FROM python:3.10-slim
WORKDIR /app
COPY requirements.txt ./
RUN pip install -r requirements.txt
COPY . .
CMD ["python", "main.py"]
```

搭配 `docker-compose.yml` 管理多容器專案（如 Jupyter + PostgreSQL）也非常實用。

---

## 五、整合範例

```bash
conda create -n prj_env python=3.9
conda activate prj_env
pip install jupyterlab matplotlib seaborn
jupyter lab
```

---

良好的環境控制能幫助資料科學家減少部署與開發過程中的錯誤，並促進團隊協作與研究再現性，是數據導向專案不可或缺的一環。
