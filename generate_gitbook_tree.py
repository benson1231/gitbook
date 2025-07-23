import os

def generate_summary_md(root_dir, base_path=".", level=0):
    output = ""
    indent = "  " * level

    # List all entries excluding hidden ones (starting with ".")
    entries = sorted([
        entry for entry in os.listdir(root_dir)
        if not entry.startswith(".")
    ])

    for entry in entries:
        full_path = os.path.join(root_dir, entry)
        # Normalize path and make it relative (use forward slashes)
        rel_path = os.path.relpath(full_path, base_path).replace("\\", "/")

        if os.path.isdir(full_path):
            # Link to the folder's README.md
            output += f"{indent}- [{entry}]({rel_path}/README.md)\n"
            # Recursively process subdirectories
            output += generate_summary_md(full_path, base_path, level + 1)
        else:
            # Add Markdown file to summary (excluding README.md itself)
            if entry.endswith(".md") and entry.lower() != "readme.md":
                output += f"{indent}- [{entry}]({rel_path})\n"
    return output

# 🔧 Set the root directory of your GitBook content
root_directory = "./"  # Change this to your content folder if needed

# Generate the content for SUMMARY.md
summary_markdown = "# Summary\n\n" + generate_summary_md(root_directory, base_path=root_directory)

# Write the result to a file
with open("SUMMARY.md", "w", encoding="utf-8") as f:
    f.write(summary_markdown)

print("✅ SUMMARY.md has been generated successfully.")
