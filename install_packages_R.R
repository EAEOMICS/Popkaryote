if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager") # Instalar BiocManager
library(devtools)

ip <- read.csv("r_packages.csv")

for(i in 1:nrow(ip)) {
    pkg <- ip$Package[i]
    ver <- ip$Version[i]
    
    if (!requireNamespace(pkg, quietly=TRUE)) {
        message(paste("Intentando instalar:", pkg, "versión:", ver))
        
        # Intentamos primero con el método estándar (CRAN)
        try_install <- try({
            install_version(pkg, version = ver, repos = "https://cloud.r-project.org")
        }, silent = TRUE)
        
        # Si falla (ej. es de Bioconductor), usamos BiocManager
        if (inherits(try_install, "try-error")) {
            message(paste("No encontrado en CRAN. Intentando Bioconductor para:", pkg))
            # Nota: BiocManager instala la versión compatible con tu versión de R, 
            # forzar una versión antigua exacta en Bioc es complejo y propenso a errores.
            BiocManager::install(pkg, update = FALSE, ask = FALSE)
        }
    }
}