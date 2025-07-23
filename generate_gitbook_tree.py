import os

def generate_tree_md(root_dir, base_path=".", level=0):
    output = ""
    indent = "  " * level

    # Ignore files or directories starting with "."
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

# 🔧 Set the root directory of your GitBook
root_directory = "./"  # e.g., "./content" or "./src"

# Generate the directory tree in Markdown format
tree_markdown = generate_tree_md(root_directory, base_path=root_directory)

# Write output to a markdown file
with open("gitbook_tree.md", "w", encoding="utf-8") as f:
    f.write("# 📚 GitBook Navigation Tree\n\n")
    f.write(tree_markdown)

print("✅ gitbook_tree.md has been generated.")
