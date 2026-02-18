# Machine Learning Projects


cd ~/01-R-and-D/core3/Deep_ML/Machine\ Learning

mkdir P1

echo "# P1" > P1/README.md

touch P1/main.py

git add P1

git commit -m "Add P1 directory with README and main.py"

git pull origin main --rebase

git push origin main


# Go to the Machine Learning folder
cd ~/01-R-and-D/core3/Deep_ML/Machine\ Learning

# Create folders P4 to P20 with README.md and main.py
for i in {4..20}; do
    mkdir P$i
    echo "# P$i" > P$i/README.md
    touch P$i/main.py
done

# Go back to repo root
cd ~/01-R-and-D

# Stage all new folders
git add core3/Deep_ML/Machine\ Learning/P{4..20}

# Commit all new folders
git commit -m "Add P4 to P20 directories with README and main.py"

# Push to GitHub
git push origin main

