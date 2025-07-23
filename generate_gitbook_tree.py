import os

def generate_tree_md(root_dir, base_path=".", level=0):
    output = ""
    indent = "  " * level

    # 忽略以 "." 開頭的檔案或資料夾
    entries = sorted([
        entry for entry in os.listdir(root_dir)
        if not entry.startswith(".")
    ])

    for entry in entries:
        full_path = os.path.join(root_dir, entry)
        rel_path = os.path.relpath(full_path, base_path).replace("\\", "/")

        if os.path.isdir(full_path):
            output += f"{indent}- 📁 [{entry}]({rel_path}/README.md)\n"
            output += generate_tree_md(full_path, base_path, level + 1)
        else:
            if entry.endswith(".md") and entry.lower() != "readme.md":
                output += f"{indent}- 📄 [{entry}]({rel_path})\n"
    return output

# 🔧 指定你的 GitBook 根目錄
root_directory = "./"  # 例如 "./content" 或 "./src"

tree_markdown = generate_tree_md(root_directory, base_path=root_directory)

# 輸出結果存成 markdown
with open("gitbook_tree.md", "w", encoding="utf-8") as f:
    f.write("# 📚 GitBook 導覽目錄\n\n")
    f.write(tree_markdown)

print("✅ 已產生 gitbook_tree.md")
