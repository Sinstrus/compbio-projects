# Git Sync Workflow for Two Laptops

This guide explains how to keep your computational biology projects synchronized between your home laptop and work laptop using Git and GitHub.

## Repository Information
- **GitHub Repository**: https://github.com/Sinstrus/compbio-projects
- **Home Laptop Path**: `/home/cnguy/projects/`
- **Work Laptop Path**: `/home/cnguyen/projects/`

---

## Initial Setup on Work Laptop

You'll need to do this ONE TIME on your work laptop:

### 1. Install and Configure Git

```bash
# Check if Git is installed
git --version

# If not installed, install it:
sudo apt update && sudo apt install git

# Configure Git with your information
git config --global user.name "Cnguy"
git config --global user.email "cnguyen258@gmail.com"
```

### 2. Set Up SSH Keys

```bash
# Generate SSH key
ssh-keygen -t ed25519 -C "cnguyen258@gmail.com" -f ~/.ssh/id_ed25519 -N ""

# Start SSH agent and add key
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519

# Display your public key
cat ~/.ssh/id_ed25519.pub
```

**Add the SSH key to GitHub:**
1. Copy the output from the `cat` command above
2. Go to https://github.com/settings/keys
3. Click "New SSH key"
4. Title: "Work Laptop - Ubuntu WSL"
5. Paste the key and click "Add SSH key"

**Test the connection:**
```bash
ssh -T git@github.com
# Should see: "Hi Sinstrus! You've successfully authenticated..."
```

### 3. Clone the Repository

```bash
# Navigate to your projects directory (create if needed)
cd ~
mkdir -p projects
cd projects

# Clone the repository
git clone git@github.com:Sinstrus/compbio-projects.git

# Navigate into the repository
cd compbio-projects
```

---

## Daily Workflow

### Starting Work (Pull Latest Changes)

**Always start by pulling the latest changes from GitHub:**

```bash
cd /home/cnguyen/projects  # or /home/cnguy/projects on home laptop

# Check current status
git status

# Pull latest changes from GitHub
git pull origin main
```

This ensures you have the most recent version of your code before you start working.

### During Work (Make Changes)

Work on your projects as normal. You can check what's changed at any time:

```bash
# See what files have been modified
git status

# See detailed changes in files
git diff

# See changes for a specific file
git diff path/to/file.py
```

### Ending Work (Push Changes to GitHub)

**When you're done working and want to sync your changes:**

```bash
# 1. Check what's changed
git status

# 2. Add changed files
git add .                    # Add all changes
# OR
git add specific_file.py     # Add specific files only

# 3. Create a commit with a descriptive message
git commit -m "Brief description of what you changed"

# 4. Push to GitHub
git push origin main
```

**Example workflow:**
```bash
git status
git add DeepRisdiplam/train_deep_risdiplam_v16.py
git commit -m "Add new training script with improved architecture"
git push origin main
```

---

## Common Scenarios

### Scenario 1: You Made Changes on Both Laptops

If you forgot to pull before working and made changes on both laptops:

```bash
# Try to pull
git pull origin main

# If there are conflicts, Git will tell you
# Edit the conflicted files to resolve conflicts
# Then:
git add .
git commit -m "Resolved merge conflicts"
git push origin main
```

**Best practice**: Always pull before starting work to avoid conflicts!

### Scenario 2: Check if You Have Unpushed Changes

```bash
# See if you have local commits not yet pushed
git log origin/main..HEAD

# If empty, you're synced. If not, you have unpushed commits.
```

### Scenario 3: Undo Local Changes

```bash
# Discard changes to a specific file (careful! can't undo)
git checkout -- filename.py

# Discard all local changes (careful! can't undo)
git reset --hard HEAD
```

### Scenario 4: View Commit History

```bash
# See recent commits
git log --oneline -10

# See detailed history
git log
```

---

## Quick Reference Commands

```bash
# Status and information
git status                    # See current changes
git log --oneline -10        # View recent commits
git diff                     # See detailed changes

# Syncing
git pull origin main         # Get latest changes from GitHub
git push origin main         # Send your changes to GitHub

# Making commits
git add .                    # Stage all changes
git add filename.py          # Stage specific file
git commit -m "message"      # Create commit
git push origin main         # Push to GitHub

# Branches (advanced)
git branch                   # List branches
git branch feature-name      # Create new branch
git checkout feature-name    # Switch to branch
git merge feature-name       # Merge branch into current
```

---

## Important Notes

### What's Excluded from Git

The `.gitignore` file automatically excludes:
- Large data files: `.fa`, `.fastq`, `.bam`, `.vcf`, `.tsv`, `.csv`
- Model files: `.pt`, `.pth`, `.h5`, `.pkl`
- Data directories: `data/`, `results/`, `outputs/`
- Python artifacts: `__pycache__/`, `*.pyc`, `venv/`

**This means:**
- Large datasets and model checkpoints won't sync (you'll need to regenerate or transfer separately)
- Only your code, scripts, and small result files will sync
- This keeps the repository fast and within GitHub's limits

### Transferring Large Files

For large data files that aren't in Git:
```bash
# From home to work (or vice versa)
# Option 1: Use scp between laptops if on same network
scp -r ~/projects/DeepRisdiplam/data/ user@work-laptop:~/projects/DeepRisdiplam/

# Option 2: Use cloud storage (Google Drive, Dropbox, etc.)
# Option 3: Use external drive/USB
```

### Repository Size

Your repository is currently ~4.3MB, which is very manageable. GitHub has a 100MB file size limit and recommends repositories under 1GB.

---

## Troubleshooting

### "Permission denied (publickey)"

Your SSH key isn't set up correctly:
```bash
# Check if key exists
ls -la ~/.ssh/id_ed25519*

# If missing, regenerate (see "Set Up SSH Keys" above)
# Make sure to add to GitHub

# Test connection
ssh -T git@github.com
```

### "Updates were rejected because the remote contains work"

You have changes on GitHub that aren't on your laptop:
```bash
git pull origin main
# Resolve any conflicts
git push origin main
```

### "Your branch is ahead of 'origin/main'"

You have local commits not pushed to GitHub:
```bash
git push origin main
```

### Need to See Full List of Changed Files

```bash
git status
git diff --name-only
```

---

## Best Practices

1. **Always pull before starting work**: `git pull origin main`
2. **Commit often with descriptive messages**: Small, focused commits are better than large ones
3. **Push at the end of each work session**: Don't leave unpushed changes
4. **Check status frequently**: `git status` is your friend
5. **Don't commit large data files**: The `.gitignore` handles this, but be aware
6. **Use meaningful commit messages**: Future you will thank you

---

## Example Daily Workflow

**Morning (starting work):**
```bash
cd ~/projects
git pull origin main
# Start working...
```

**During the day (checkpoint):**
```bash
git add .
git commit -m "Implemented new feature for sequence analysis"
git push origin main
```

**Evening (end of work):**
```bash
git status
git add .
git commit -m "Final updates for the day - fixed bugs in training script"
git push origin main
```

**Next day on other laptop:**
```bash
cd ~/projects
git pull origin main
# You now have all yesterday's changes!
```

---

## Getting Help

- Git documentation: https://git-scm.com/doc
- GitHub guides: https://guides.github.com
- Check this file: `cat SYNC_WORKFLOW.md`
