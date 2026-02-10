# BluePebble HPC Setup Instructions

## Recommended Workflow

**Local → GitHub → BluePebble**

1. **Make changes locally** (edit code, update Docker image)
2. **Commit and push** to GitHub
3. **On BluePebble**: `git pull` to get latest code/scripts
4. **Rebuild Docker image** (if code changed) → push to GHCR
5. **Re-pull container** on BluePebble (if image updated)

## Storage Strategy on BluePebble

| What | Where | Why |
|------|-------|-----|
| **Scripts** (repo clone) | `~/TissueMorphology/` | Home directory (small files, backed up) |
| **Container cache** | `/user/work/$USER/.apptainer/` | Scratch (large temp files during pull) |
| **Container image** (`.sif`) | `/user/work/$USER/containers/` | Scratch (large, regeneratable) |
| **Output data** | `/user/work/$USER/chaste_output/` | Scratch (large, fast I/O) |
| **Job logs** | `/user/work/$USER/logs/` | Scratch (timestamped execution logs) |

**Why scratch storage for containers?**
- Home directory quota: **50 GB** (fills up fast with one container)
- Scratch space quota: **2-5 TB** per user
- Containers are ~2-10 GB compressed, ~20-40 GB cache during pull
- Scratch has faster I/O for reading/writing large simulation outputs

## Prerequisites (Local)

1. **Build and push Docker image**:
   ```bash
   cd Chaste/projects/TissueMorphology
   ./docker/build_and_push.sh
   ```

2. **Authentication** (if not already done):
   ```bash
   docker login ghcr.io -u lusCombeNCut
   # Paste GitHub Personal Access Token
   ```

## On BluePebble HPC

### 1. Configure scratch storage (one-time)

SSH to BluePebble:
```bash
ssh sv22482@bp1-login.acrc.bris.ac.uk
```

Add to your `~/.bashrc`:
```bash
echo 'export APPTAINER_CACHEDIR=/user/work/$USER/.apptainer' >> ~/.bashrc
source ~/.bashrc
```

Create container directory:
```bash
mkdir -p /user/work/$USER/containers
mkdir -p /user/work/$USER/chaste_output
mkdir -p /user/work/$USER/logs
```

### 2. Clone the TissueMorphology repository

```bash
cd ~
git clone https://github.com/lusCombeNCut/TissueMorphology.git
cd TissueMorphology
```

Now you have all the scripts! When you make changes locally and push to GitHub, just run:
```bash
git pull
```

### 3. Authenticate with GitHub Container Registry

**One-time authentication** (required to pull private container images):

```bash
apptainer remote login -u lusCombeNCut docker://ghcr.io
```

When prompted for password, enter your **GitHub Personal Access Token** (not GitHub password).

**Create a PAT** (if you don't have one):
1. Go to https://github.com/settings/tokens
2. Click "Generate new token" → "Generate new token (classic)"
3. Select scope: **`read:packages`** (required)
4. Copy the token and paste when prompted

**Verify authentication:**
```bash
apptainer remote status docker://ghcr.io
```

### 4. Pull the container image

Run the setup script:
```bash
cd ~/TissueMorphology
bash hpc/setup_hpc.sh
```

This will:
- Load Apptainer module
- Set cache directory to scratch space (if not in `.bashrc`)
- Check GHCR authentication
- Pull the Docker image as a `.sif` file to `/user/work/$USER/containers/`
- Verify the container
- Create output and logs directories

**First pull takes 10-30 minutes** (downloads ~5 GB, extracts to ~20 GB cache, creates ~10 GB .sif)

### 5. Submit jobs

```bash
cd ~/TissueMorphology
sbatch hpc/submit_vertex_organoid.sh
```

Monitor:
```bash
squeue -u $USER                                    # Check job status
tail -f /user/work/$USER/logs/vertex_organoid_*.log  # Watch latest log
ls -lht /user/work/$USER/logs/ | head              # List recent logs
ls -lh /user/work/$USER/chaste_output/             # Check output files
```

**Output organization:**
- Log files: `/user/work/$USER/logs/vertex_organoid_YYYY-MM-DD_HH-MM-SS_jobID.log`
- Simulation output: `/user/work/$USER/chaste_output/vertex_organoid_YYYY-MM-DD_HH-MM-SS_jobID/`
- All output is timestamped for easy tracking of multiple runs

## Updating workflow

When you change code locally:

```bash
# Local machine
cd Chaste/projects/TissueMorphology
git add .
git commit -m "Updated simulation parameters"
git push

# If you changed C++ code, rebuild & push Docker image:
./docker/build_and_push.sh
```

On BluePebble:

```bash
cd ~/TissueMorphology
git pull

# If Docker image was updated, re-pull container:
bash hpc/setup_hpc.sh

# Submit updated job
sbatch hpc/submit_vertex_organoid.sh
```

---

## Key BluePebble-specific details

- **Module loading required**: `module load apptainer` (loads default version)
- **Cache directory**: Set to `/user/work/$(whoami)/.apptainer` to avoid filling home directory
- **Account code**: `semt036404` (already set in submit script)
- **Apptainer command**: Use `apptainer` (not `singularity`)

## Troubleshooting

**"The following module(s) are unknown: apptainer/X.X.X"**

Check available versions and load the default:
```bash
module avail apptainer
module load apptainer
```

**"apptainer: command not found"**
```bash
module load apptainer
```

**"Cannot create directory /home/chaste/build/Testing/Temporary"**

This is fixed in the latest version of the submit script. Pull the latest changes:
```bash
cd ~/TissueMorphology
git pull
```

**"UID: readonly variable" errors**

The submit script has been updated to remove `--cleanenv` flag and properly bind writable directories for CTest. Make sure you're using the latest version.

**"Permission denied" on cache directory**
```bash
mkdir -p /user/work/$(whoami)/.apptainer
export APPTAINER_CACHEDIR=/user/work/$(whoami)/.apptainer
```

**Check available modules**
```bash
module avail apptainer
```

**Home directory full?**
Container files are large. Make sure:
```bash
# Check your home quota
du -sh ~/*
quota -s

# Make sure .sif is in scratch, not home:
ls -lh /user/work/$USER/containers/tissuemorphology.sif

# Cache should also be in scratch:
echo $APPTAINER_CACHEDIR
# Should show: /user/work/$USER/.apptainer
```

**Make APPTAINER_CACHEDIR permanent** (add to `~/.bashrc`):
```bash
echo 'export APPTAINER_CACHEDIR=/user/work/$USER/.apptainer' >> ~/.bashrc
source ~/.bashrc
```
