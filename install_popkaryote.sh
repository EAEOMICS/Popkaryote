#!/bin/bash

# ==============================================================================
# SCRIPT SETUP V6: TRUST YML + R SMART INSTALL
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. & 2. VERIFICACIÓN E INSTALACIÓN DE CONDA (INTERACTIVO)
# ------------------------------------------------------------------------------
CONDA_BASE_PATH=""

if ! command -v conda &> /dev/null; then
    echo ">>> Conda no encontrado. Iniciando instalación..."
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    INSTALLER_FILE="miniconda.sh"
    
    if [ ! -f "$INSTALLER_FILE" ]; then
        wget $MINICONDA_URL -O $INSTALLER_FILE
    fi
    
    INSTALL_DIR="$HOME/miniconda3"
    INSTALLED_SUCCESSFULLY=false

    while [ "$INSTALLED_SUCCESSFULLY" = false ]; do
        echo ">>> Intentando instalar Conda en: $INSTALL_DIR"
        bash $INSTALLER_FILE -b -p "$INSTALL_DIR" 
        
        if [ $? -eq 0 ]; then
            echo ">>> Instalación exitosa."
            INSTALLED_SUCCESSFULLY=true
            CONDA_BASE_PATH="$INSTALL_DIR"
        else
            echo "!!! ERROR: No se pudo instalar en '$INSTALL_DIR'"
            read -p ">>> Ingresa una nueva ruta completa (ej: /home/usuario/conda_env): " NEW_PATH
            INSTALL_DIR="${NEW_PATH/#\~/$HOME}"
            if [ -z "$INSTALL_DIR" ]; then INSTALL_DIR="$HOME/miniconda3"; fi
        fi
    done
    rm $INSTALLER_FILE
    source "$CONDA_BASE_PATH/etc/profile.d/conda.sh"
else
    echo ">>> Conda ya está instalado."
    CONDA_BASE=$(conda info --base)
    source "$CONDA_BASE/etc/profile.d/conda.sh"
fi

# A partir de aquí, si hay errores, paramos el script
set -e 

# ------------------------------------------------------------------------------
# 3. CREAR ENTORNO BASADO EN TU YML
# ------------------------------------------------------------------------------
YML_PATH="./envs/popkaryote3.yml"
ENV_NAME="popkaryote"

echo ">>> Creando/Actualizando entorno '$ENV_NAME' usando $YML_PATH..."
if [ ! -f "$YML_PATH" ]; then
    echo "Error: El archivo $YML_PATH no existe."
    exit 1
fi

# Tu YML ya contiene los compiladores y librerías gráficas (systemfonts, etc.)
conda env create -f "$YML_PATH" -n "$ENV_NAME" || conda env update -f "$YML_PATH" -n "$ENV_NAME"

# ------------------------------------------------------------------------------
# 4. INSTALACIÓN DE PAQUETES R FALTANTES (CSV)
# ------------------------------------------------------------------------------
echo ">>> Activando entorno '$ENV_NAME'..."
conda activate "$ENV_NAME"

echo ">>> Generando script de instalación de R..."

cat <<EOF > install_packages_R.R
# --- INICIO SCRIPT R ---

# 1. Configuración inicial
if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools", repos="https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")

library(devtools)

# 2. Leer CSV
if (!file.exists("r_packages.csv")) {
    stop("Error: r_packages.csv no encontrado.")
}
ip <- read.csv("r_packages.csv")

# 3. LISTA DE EXCLUSIÓN
# Estos paquetes ya vienen en R base o son parte del núcleo, intentar reinstalarlos suele dar error.
base_packages <- c(
    "base", "compiler", "datasets", "graphics", "grDevices", "grid",
    "methods", "parallel", "splines", "stats", "stats4", "tcltk", 
    "tools", "utils"
)

# 4. Función de instalación inteligente
install_pkg_smart <- function(pkg, ver) {
    
    # Si es paquete base, saltar
    if (pkg %in% base_packages) {
        return(TRUE)
    }

    # CHEQUEO CRUCIAL:
    # Tu YML ya instaló paquetes como 'systemfonts', 'textshaping', 'ragg'.
    # R detectará aquí que ya existen y NO intentará reinstalarlos, evitando errores.
    if (requireNamespace(pkg, quietly=TRUE)) {
        message(paste(">>> [OK]", pkg, "ya está instalado."))
        return(TRUE)
    }

    message(paste(">>> Intentando instalar:", pkg, "Versión deseada:", ver))

    # --- A. GITHUB (Paquetes específicos) ---
    if (pkg == "BactDating") {
        try(devtools::install_github("xavierdidelot/BactDating"))
        return()
    }
    if (pkg == "LDWeaver") {
        # Usando el repo corregido
        try(devtools::install_github("sudaraka88/LDWeaver")) 
        return()
    }

    # --- B. BIOCONDUCTOR (Prioritario para ggtree, etc.) ---
    # Si el paquete está en Bioconductor, usamos BiocManager.
    # Esto manejará ggtree, Biostrings, etc.
    if (pkg %in% BiocManager::available()) {
         message(paste(">>> Instalando vía Bioconductor:", pkg))
         try_bioc <- try({
             BiocManager::install(pkg, update = FALSE, ask = FALSE)
         }, silent = TRUE)
         
         if (!inherits(try_bioc, "try-error")) return()
    }

    # --- C. CRAN (Versión específica) ---
    message(paste(">>> Instalando vía CRAN (Versión exacta):", pkg))
    try_cran <- try({
        install_version(pkg, version = ver, repos = "https://cloud.r-project.org")
    }, silent = TRUE)

    if (inherits(try_cran, "try-error")) {
        warning(paste("!!! FALLÓ la instalación de:", pkg))
    }
}

# 5. Ejecutar Loop
for(i in 1:nrow(ip)) {
    install_pkg_smart(ip\$Package[i], ip\$Version[i])
}

print(">>> Instalación de paquetes R completada.")
# --- FIN SCRIPT R ---
EOF

echo ">>> Ejecutando script de R..."
Rscript install_packages_R.R

echo "=============================================================================="
echo ">>> Setup Completado Exitosamente."
echo ">>> Entorno listo: popkaryote"
echo "=============================================================================="