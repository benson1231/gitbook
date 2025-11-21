import os
import argparse

def generate_catalog(root_dir, base_path=".", level=0):
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
            output += generate_catalog(full_path, base_path, level + 1)
        else:
            if entry.endswith(".md") and entry.lower() != "readme.md":
                output += f"{indent}- 📄 [{entry}]({rel_path})\n"
    return output


def main():
    parser = argparse.ArgumentParser(
        description="📚 Generate a GitBook Catalog in Markdown."
    )
    parser.add_argument(
        "root_directory",
        help="Root directory of your GitBook (e.g., ./content or ./src)"
    )
    parser.add_argument(
        "-o", "--output",
        default="catalog.md",
        help="Output markdown file (default: catalog.md)"
    )
    args = parser.parse_args()

    # Generate tree
    tree_markdown = generate_catalog(args.root_directory, base_path=args.root_directory)

    # Write to file
    with open(args.output, "w", encoding="utf-8") as f:
        f.write("# 📚 目錄\n\n")
        f.write(tree_markdown)

    print(f"✅ {args.output} has been generated.")


if __name__ == "__main__":
    main()
