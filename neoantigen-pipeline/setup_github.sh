#!/bin/bash
# ============================================================================
# GitHub Repository Setup Script
# ============================================================================
# Run this script after extracting the pipeline to upload to GitHub
#
# Usage:
#   1. Extract neoantigen-pipeline.zip
#   2. cd neoantigen-pipeline
#   3. chmod +x setup_github.sh
#   4. ./setup_github.sh YOUR_GITHUB_USERNAME
#
# ============================================================================

set -e

# Check if username provided
if [ -z "$1" ]; then
    echo "Usage: ./setup_github.sh YOUR_GITHUB_USERNAME"
    echo "Example: ./setup_github.sh chukwuma-io"
    exit 1
fi

GITHUB_USER=$1
REPO_NAME="neoantigen-pipeline"

echo "============================================"
echo "  Neoantigen Pipeline - GitHub Setup"
echo "============================================"
echo ""

# Check if git is installed
if ! command -v git &> /dev/null; then
    echo "Error: git is not installed"
    exit 1
fi

# Check if gh CLI is installed (optional)
HAS_GH=false
if command -v gh &> /dev/null; then
    HAS_GH=true
fi

echo "Step 1: Initializing git repository..."
git init

echo ""
echo "Step 2: Adding files..."
git add .

echo ""
echo "Step 3: Creating initial commit..."
git commit -m "Initial commit: Neoantigen prediction pipeline

Features:
- HLA Class I typing (OptiType)
- HLA Class II typing (HLA-HD)
- Neoantigen prediction (pVACtools + NetMHCpan/IIpan)
- Priority scoring algorithm
- Clinical HTML reports
- Docker/Singularity support"

echo ""
echo "Step 4: Setting up remote..."

if [ "$HAS_GH" = true ]; then
    echo "GitHub CLI detected. Creating repository..."
    gh repo create $REPO_NAME --public --source=. --remote=origin --push
    echo ""
    echo "✓ Repository created and pushed!"
    echo "  URL: https://github.com/$GITHUB_USER/$REPO_NAME"
else
    echo ""
    echo "GitHub CLI not installed. Please complete these steps manually:"
    echo ""
    echo "  1. Go to https://github.com/new"
    echo "  2. Create a new repository named: $REPO_NAME"
    echo "  3. Run these commands:"
    echo ""
    echo "     git branch -M main"
    echo "     git remote add origin https://github.com/$GITHUB_USER/$REPO_NAME.git"
    echo "     git push -u origin main"
    echo ""
fi

echo ""
echo "============================================"
echo "  Setup Complete!"
echo "============================================"
echo ""
echo "Next steps:"
echo "  1. Update README.md with your GitHub username"
echo "  2. Add any additional documentation"
echo "  3. Configure GitHub Actions if desired"
echo ""
echo "To update the repository:"
echo "  git add ."
echo "  git commit -m 'Your commit message'"
echo "  git push"
echo ""
