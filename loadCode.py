import subprocess
import os

REPO_NAME = "GrandPotentialPython"
REPO_URL = f"https://github.com/tagtog12000/GrandPotentialPython.git"
# --- 2. Clone the GitHub repo ---
if not os.path.exists(REPO_NAME):
    print(f"🔽 Cloning repository: {REPO_URL}")
    subprocess.run(["git", "clone", REPO_URL])
else:
    print("✅ Repository already cloned.")
# --- 3. Install requirements ---
if os.path.exists("requirements.txt"):
    print("📦 Installing Python dependencies ...")
    subprocess.run(["pip", "install", "-r", "requirements.txt"], check=True)
else:
    print("⚠️ No requirements.txt found — skipping dependency install.")

!git clone https://github.com/tagtog12000/GrandPotentialPython.git
%cd GrandPotentialPython

# STEP 2: Compile single C++ file
!g++ GP.cpp -o code

# STEP 3: Run it
!./code

# --- Run Section ---

nMax = int(input("Enter the maximum order nMax: "))


eps = int(input("Enter the particle statistics, eps (-1 for fermions and 1 for bosons): "))

try:
    # Run the compiled executable with both arguments
    result = subprocess.run(
        ["./code", str(nMax), str(eps)],
        capture_output=True, text=True, check=True
    )
    print("✅ Program output:")
    print(result.stdout)
except subprocess.CalledProcessError as e:
    print("❌ Error while running executable:")
    print(e.stderr)
