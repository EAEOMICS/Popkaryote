if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager") # Install BiocManager
library(devtools)

# Read the package list
ip <- read.csv("r_packages.csv")

for(i in 1:nrow(ip)) {
    pkg <- ip$Package[i]
    ver <- ip$Version[i]
    
    # Check if the package is already installed
    if (!requireNamespace(pkg, quietly=TRUE)) {
        message(paste("Attempting to install:", pkg, "version:", ver))
        
        # First attempt: standard method (CRAN) for specific versions
        try_install <- try({
            install_version(pkg, version = ver, repos = "https://cloud.r-project.org")
        }, silent = TRUE)
        
        # If CRAN fails (e.g., the package is from Bioconductor), use BiocManager
        if (inherits(try_install, "try-error")) {
            message(paste("Not found on CRAN. Attempting Bioconductor for:", pkg))
            
            # Note: BiocManager installs the version compatible with your current R version. 
            # Forcing an exact old version in Bioc is complex and error-prone.
            BiocManager::install(pkg, update = FALSE, ask = FALSE)
        }
    }
}