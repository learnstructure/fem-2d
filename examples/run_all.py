#run all examples in struct_draw/examples
import subprocess
from pathlib import Path

# Define the directory containing Python files (use "." for the current folder)
script_dir = Path("./examples/dynamic")
print(f"Running all Python scripts in: {script_dir.resolve()}")

# Loop through and execute every .py file
for file in script_dir.glob("*.py"):
    # Avoid running this controller script recursively if it's in the same folder
    if file.name == "run_all.py":
        continue

    print(f"Executing: {file.name}")
    try:
        subprocess.run(["python", str(file)], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing {file.name}: {e}")
