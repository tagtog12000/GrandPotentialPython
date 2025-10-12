# loadCode.py
import subprocess
import os

def setup_repo():
    REPO_NAME = "GrandPotentialPython"
    REPO_URL = f"https://github.com/tagtog12000/GrandPotentialPython.git"

    if not os.path.exists(REPO_NAME):
        print(f"🔽 Cloning repository: {REPO_URL}")
        subprocess.run(["git", "clone", REPO_URL])
    else:
        print("✅ Repository already cloned.")

    if os.path.exists("requirements.txt"):
        print("📦 Installing Python dependencies ...")
        subprocess.run(["pip", "install", "-r", "requirements.txt"], check=True)
    else:
        print("⚠️ No requirements.txt found — skipping dependency install.")

def compile_and_run(nMax, eps):
    subprocess.run(["g++", "GP.cpp", "-o", "code"], check=True)
    result = subprocess.run(["./code", str(nMax), str(eps)],
                            capture_output=True, text=True)
    print(result.stdout)
