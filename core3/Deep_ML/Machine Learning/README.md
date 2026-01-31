# Machine Learning Projects

This repository contains **Python solutions** for various **machine learning and deep learning coding problems**, covering topics such as **linear algebra, statistics, ML, and deep learning**.  
Each project is organized in its own folder (`PXX`) with:

- A **README.md** describing the problem, dataset, or task, including examples and the solution approach.  
- A **main.py** file demonstrating the implementation using **Python**.

---

## Example: How to Add a New Project

```bash
# Navigate to your ML projects folder
cd ~/01-R-and-D/core3/Deep_ML/Machine\ Learning

# Create a new project folder
mkdir P1

# Create README and main.py files
echo "# P1" > P1/README.md
touch P1/main.py

# Stage and commit changes
git add P1
git commit -m "Add P1 directory with README and main.py"

# Pull latest changes and rebase
git pull origin main --rebase

# Push your new project to GitHub
git push origin main
