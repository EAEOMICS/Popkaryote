#!/bin/bash

# ==============================================================================
# SCRIPT SETUP V6: TRUST YML + R SMART INSTALL
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. & 2. CONDA VERIFICATION AND INSTALLATION (INTERACTIVE)
# ------------------------------------------------------------------------------
CONDA_BASE_PATH=""

if ! command -v conda &> /dev/null; then
    echo ">>> Conda not found. Starting installation..."
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    INSTALLER_FILE="miniconda.sh"
    
    if [ ! -f "$INSTALLER_FILE" ]; then
        wget $MINICONDA_URL -O $INSTALLER_FILE
    fi
    
    INSTALL_DIR="$HOME/miniconda3"
    INSTALLED_SUCCESSFULLY=false

    while [ "$INSTALLED_SUCCESSFULLY" = false ]; do
        echo ">>> Attempting to install Conda in: $INSTALL_DIR"
        bash $INSTALLER_FILE -b -p "$INSTALL_DIR" 
        
        if [ $? -eq 0 ]; then
            echo ">>> Installation successful."
            INSTALLED_SUCCESSFULLY=true
            CONDA_BASE_PATH="$INSTALL_DIR"
        else
            echo "!!! ERROR: Failed to install in '$INSTALL_DIR'"
            read -p ">>> Please enter a new full path (e.g., /home/user/conda_env): " NEW_PATH
            INSTALL_DIR="${NEW_PATH/#\~/$HOME}"
            if [ -z "$INSTALL_DIR" ]; then INSTALL_DIR="$HOME/miniconda3"; fi
        fi
    done
    rm $INSTALLER_FILE
    source "$CONDA_BASE_PATH/etc/profile.d/conda.sh"
else
    echo ">>> Conda is already installed."
    CONDA_BASE=$(conda info --base)
    source "$CONDA_BASE/etc/profile.d/conda.sh"
fi

# From this point on, stop the script if any errors occur
set -e 

# ------------------------------------------------------------------------------
# 3. CREATE ENVIRONMENT BASED ON YOUR YML
# ------------------------------------------------------------------------------
YML_PATH="./envs/popkaryote3.yml"
ENV_NAME="popkaryote"

echo ">>> Creating/Updating environment '$ENV_NAME' using $YML_PATH..."
if [ ! -f "$YML_PATH" ]; then
    echo "Error: File $YML_PATH does not exist."
    exit 1
fi

# Your YML already contains compilers and graphics libraries (systemfonts, etc.)
conda env create -f "$YML_PATH" -n "$ENV_NAME" || conda env update -f "$YML_PATH" -n "$ENV_NAME"

# ------------------------------------------------------------------------------
# 4. INSTALL MISSING R PACKAGES (CSV)
# ------------------------------------------------------------------------------
echo ">>> Activating environment '$ENV_NAME'..."
conda activate "$ENV_NAME"

echo ">>> Generating R installation script..."

cat <<EOF > install_packages_R.R
# --- START R SCRIPT ---

# 1. Initial Setup
if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools", repos="https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")

library(devtools)

# 2. Read CSV
if (!file.exists("r_packages.csv")) {
    stop("Error: r_packages.csv not found.")
}
ip <- read.csv("r_packages.csv")

# 3. EXCLUSION LIST
# These packages are already in base R or part of the core; attempting to reinstall them usually fails.
base_packages <- c(
    "base", "compiler", "datasets", "graphics", "grDevices", "grid",
    "methods", "parallel", "splines", "stats", "stats4", "tcltk", 
    "tools", "utils"
)

# 4. Smart Installation Function
install_pkg_smart <- function(pkg, ver) {
    
    # If it is a base package, skip
    if (pkg %in% base_packages) {
        return(TRUE)
    }

    # CRUCIAL CHECK:
    # Your YML already installed packages like 'systemfonts', 'textshaping', 'ragg'.
    # R will detect here that they already exist and will NOT attempt to reinstall them, avoiding errors.
    if (requireNamespace(pkg, quietly=TRUE)) {
        message(paste(">>> [OK]", pkg, "is already installed."))
        return(TRUE)
    }

    message(paste(">>> Attempting to install:", pkg, "Desired version:", ver))

    # --- A. GITHUB (Specific packages) ---
    if (pkg == "BactDating") {
        try(devtools::install_github("xavierdidelot/BactDating"))
        return()
    }
    if (pkg == "LDWeaver") {
        # Using the corrected repo
        try(devtools::install_github("sudaraka88/LDWeaver")) 
        return()
    }

    # --- B. BIOCONDUCTOR (Priority for ggtree, etc.) ---
    # If the package is in Bioconductor, use BiocManager.
    # This handles ggtree, Biostrings, etc.
    if (pkg %in% BiocManager::available()) {
         message(paste(">>> Installing via Bioconductor:", pkg))
         try_bioc <- try({
             BiocManager::install(pkg, update = FALSE, ask = FALSE)
         }, silent = TRUE)
         
         if (!inherits(try_bioc, "try-error")) return()
    }

    # --- C. CRAN (Specific version) ---
    message(paste(">>> Installing via CRAN (Exact version):", pkg))
    try_cran <- try({
        install_version(pkg, version = ver, repos = "https://cloud.r-project.org")
    }, silent = TRUE)

    if (inherits(try_cran, "try-error")) {
        warning(paste("!!! FAILED to install:", pkg))
    }
}

# 5. Execute Loop
for(i in 1:nrow(ip)) {
    install_pkg_smart(ip\$Package[i], ip\$Version[i])
}

print(">>> R package installation completed.")
# --- END R SCRIPT ---
EOF

echo ">>> Running R script..."
Rscript install_packages_R.R

echo "=============================================================================="
echo ">>> Setup Completed Successfully."
echo ">>> Environment ready: popkaryote"
echo "=============================================================================="