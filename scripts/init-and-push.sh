#!/bin/bash

set -e

REPO_NAME="cdi-domain-template"
GITHUB_USER="datainsights"

if [ ! -f "scripts/cdi-install-packages.R" ]; then
  echo "❌ Please run this script from the root of the $REPO_NAME directory."
  exit 1
fi

echo "🚀 Initializing Git repo for $REPO_NAME..."
git init

# Create commit if not already done
if ! git rev-parse --quiet --verify HEAD; then
  git add .
  git commit -m "Initial commit"
  git branch -M main  # <-- Important: ensure branch is renamed
fi

# Create GitHub repo and push
if ! gh repo view "$GITHUB_USER/$REPO_NAME" &> /dev/null; then
  echo "📡 Creating GitHub repo $GITHUB_USER/$REPO_NAME..."
  gh repo create "$GITHUB_USER/$REPO_NAME" --private --source=. --remote=origin --push
else
  echo "🔗 GitHub repo already exists. Connecting to remote..."
  git remote add origin "https://github.com/$GITHUB_USER/$REPO_NAME.git"
  git push -u origin main
fi

echo "✅ Repo successfully pushed to GitHub!"