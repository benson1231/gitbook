# RDKit

本篇筆記整理 RDKit 的常用功能，並透過與 PubChem 結合，展示如何查詢結構式、產生分子圖像與取得 SMILES/InChI 等資訊。

## 1. 安裝所需套件

```bash
pip install rdkit-pypi pubchempy
```

> 注意：`rdkit` 建議於 conda 環境中安裝以避免相依問題。

```bash
conda install -c rdkit rdkit
```

## 2. 透過 PubChem 查詢化合物並取得 RDKit 分子物件

```python
from rdkit import Chem
from rdkit.Chem import Draw

# 幾個常見分子的 SMILES
mol = Chem.MolFromSmiles('CCO')          # ethanol

Draw.MolToImage(mol)
```

## 3. 範例：查詢「caffeine」

```python
mol, smiles, inchi = get_rdkit_mol_from_name("caffeine")
print("SMILES:", smiles)
print("InChI:", inchi)
Draw.MolToImage(mol).show()
```

## 4. 常見 RDKit 操作

- `Chem.MolFromSmiles(smiles)`：從 SMILES 建立分子
- `AllChem.Compute2DCoords(mol)`：計算 2D 座標
- `Draw.MolToImage(mol)`：輸出 PIL 圖像
- `Chem.MolToInchi(mol)`：轉換為 InChI

## 5. 延伸應用

- 計算分子指紋與相似度
- 產生分子 3D 構型（`AllChem.EmbedMolecule`）
- 與機器學習模型結合（分子屬性預測）

---

如需查詢更多化合物，可搭配 `pubchempy` 使用 `cid`, `name`, `formula` 等欄位靈活檢索。

```python
pcp.get_compounds("aspirin", namespace="name")
```
