# Machine Learning Projects


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
