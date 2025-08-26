################################################################################
################################################################################
# ==============================================================================
# Script: PCA de variáveis bioclimáticas (WordClim v2.1)
# Autor: Dr. Leandro José da Silva
# Data: [25/08/2025]
#
# Descrição:
# Este script realiza a análise de Componentes Principais (PCA) 
# aplicada às variáveis bioclimáticas do WorldClim v2.1, 
# permitindo reduzir a dimensionalidade e eliminar colinearidade 
# entre variáveis ambientais.
#
# Etapas principais:
# 1. Leitura das camadas do clima atual (baseline, 1970–2000) e 
#    projeções futuras (CMIP6, SSP585, 2061–2080).
# 2. Recorte e máscara para os limites do Brasil.
# 3. Amostragem de pixels e execução do PCA no baseline.
# 4. Seleção automática do número de PCs suficientes para explicar 
#    ≥95% da variância climática.
# 5. Projeção dos cenários futuros nos mesmos eixos PCA do baseline.
# 6. Exportação dos rasters dos PCs (baseline e futuros).
# 7. Geração de produtos gráficos e analíticos:
#    - Scree plot (variância explicada e acumulada)
#    - Heatmap dos loadings (pesos das variáveis em cada PC)
#    - Biplot (PC1 x PC2)
# 8. Exportação de tabelas com os loadings em CSV.
#
# Utilidade:
# - Cria variáveis ambientais sintéticas (PCs) não correlacionadas,
#   que podem ser usadas em modelagem de nicho e análises climáticas,
#   garantindo comparabilidade entre clima atual e projeções futuras.
################################################################################


################################################################################
# ==============================================================================
# Pacotes necessários
# ==============================================================================

# Lista de pacotes
pkgs <- c(
  "terra",        
  "raster",         
  "sf",             
  "rnaturalearth",
  "dplyr",        
  "ggplot2",      
  "reshape2"        
)

# Instalar (se necessário) e carregar
invisible(lapply(pkgs, function(p){
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

# Ativar cálculos geodésicos precisos em sf
sf::sf_use_s2(TRUE)

# Definir semente para reprodutibilidade
set.seed(123)

################################################################################
# ==============================================================================
# 1) Escolher pastas/arquivos
# ==============================================================================
# Escolha a pasta do projeto
proj_dir <- choose.dir(caption = "Escolha a PASTA do projeto")
setwd(proj_dir)

# Escolha a pasta onde estão os 19 arquivos do baseline (wc2.1_30s_bio_1..19.tif)
base_dir <- choose.dir(caption = "Escolha a pasta do BASELINE (wc2.1_30s_bio_1..19.tif)")

# Escolha os arquivos FUTUROS (GeoTIFF multibanda: 1 arquivo = 19 bandas)
future_files <- c(
  file.choose(new = FALSE), # MIROC6
  file.choose(new = FALSE), # EC-Earth3-Veg
  file.choose(new = FALSE)  # IPSL-CM6A-LR
)

################################################################################
# ==============================================================================
# 2) Brasil: recorte + máscara
# ==============================================================================
bra_sf <- rnaturalearth::ne_countries(scale = 10, country = "Brazil", returnclass = "sf")
bra_v  <- terra::vect(bra_sf)

################################################################################
# ==============================================================================
# 3) Carregar BASELINE (19 variaveis bioclimática - clima atual) ex:(wc2.1_30s_bio_1)
# ==============================================================================
base_files <- list.files(base_dir, pattern = "bio_?\\d+\\.tif$|bio_?\\d+$", full.names = TRUE, ignore.case = TRUE)
# garantir que estão na ordem 1..19
ord <- order(as.numeric(gsub("\\D", "", basename(base_files))))
base_files <- base_files[ord]

if (length(base_files) < 19) stop("Não encontrei as 19 camadas do baseline.")

base_all <- terra::rast(base_files)
# recorte/máscara Brasil
base_br  <- terra::mask( terra::crop(base_all, bra_v), bra_v)

# renomear garantido: bio1..bio19
names(base_br) <- paste0("bio", 1:terra::nlyr(base_br))

################################################################################
# ==============================================================================
# 4) PCA no BASELINE
# ==============================================================================
# Amostrar pixels para PCA (para não carregar tudo na memória)
samp_pts <- terra::spatSample(base_br[[1]], size = 20000, method = "random", na.rm = TRUE, as.points = TRUE)
vals     <- terra::extract(base_br, samp_pts, ID = FALSE) |> as.data.frame()
vals     <- na.omit(vals)

# PCA (centrado e escalado)
pca <- prcomp(vals, center = TRUE, scale. = TRUE)

# Quantos PCs para >=95% variância?
eig  <- pca$sdev^2
cumv <- cumsum(eig) / sum(eig)
npcs <- which(cumv >= 0.95)[1]
cat(sprintf("Manter %d PCs (cobre %.1f%% da variância)\n", npcs, 100*cumv[npcs]))

################################################################################
# ==============================================================================
# 5) Projetar BASELINE nos PCs (rasters PC1..PCn)
# ==============================================================================
# Nota: predict(SpatRaster, prcomp) projeta cada pixel nos loadings do PCA
base_pcs <- terra::predict(base_br, pca, index = 1:npcs, na.rm = TRUE)
names(base_pcs) <- paste0("PC", 1:npcs)

# garantir pasta de saída
dir.create(file.path("out","pca_vars","baseline"), recursive = TRUE, showWarnings = FALSE)
terra::writeRaster(base_pcs, filename = file.path("out","pca_vars","baseline","baseline_PC.tif"),
                   overwrite = TRUE)
cat("✓ PCs do baseline salvos em out/pca_vars/baseline/baseline_PC.tif\n")

################################################################################
# ==============================================================================
# 6) Função para carregar FUTURO multibanda, alinhar, e projetar nos PCs
# ==============================================================================
proj_future_to_pcs <- function(fut_tif, pca_model, mask_poly, npcs, out_dir){
  r <- terra::rast(fut_tif)
  names(r) <- paste0("bio", 1:terra::nlyr(r))   # garantir bio1..bio19
  
  # recorte/máscara Brasil
  r_br <- terra::mask( terra::crop(r, mask_poly), mask_poly)
  
  # reamostrar para a grade do baseline, se necessário (WorldClim 30s costuma bater)
  if (!compareGeom(r_br, base_br, stopOnError = FALSE)){
    r_br <- terra::project(r_br, base_br)  # garante mesma grade/CRS
  }
  
  # projetar nos PCs do baseline
  pcs_fut <- terra::predict(r_br, pca_model, index = 1:npcs, na.rm = TRUE)
  names(pcs_fut) <- paste0("PC", 1:npcs)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, "future_PC.tif")
  terra::writeRaster(pcs_fut, out_file, overwrite = TRUE)
  return(out_file)
}

################################################################################
# ==============================================================================
# 7) Projetar cada FUTURO nas bases PCA do baseline
# ==============================================================================
# criar nomes de saída por arquivo
fut_names <- basename(future_files)
# tente extrair um rótulo amigável (modelo + ssp + período)
friendly <- gsub("^wc2\\.1_30s_bioc_|\\.tif$","", fut_names)

out_paths <- c()
for(i in seq_along(future_files)){
  out_dir_i <- file.path("out","pca_vars","futures", friendly[i])
  out_paths[i] <- proj_future_to_pcs(future_files[i], pca, bra_v, npcs, out_dir_i)
  cat("✓ PCs futuros salvos em:", out_paths[i], "\n")
}

################################################################################
# ==============================================================================
# 8) (Opcional) Visualizar PC1 (baseline x 1 futuro)
# ==============================================================================
plot(base_pcs[[1]], main = "Baseline PC1")
fut_example <- terra::rast(out_paths[1])
plot(fut_example[[1]], main = paste0("Futuro PC1 – ", friendly[1]))


################################################################################
# ==============================================================================
# 9) - (Opcional) Verificar resolução dos arquivos gerados
# ==============================================================================

# baseline
base_pc_file <- "out/pca_vars/baseline/baseline_PC.tif"
if (file.exists(base_pc_file)) {
  r_base <- terra::rast(base_pc_file)
  cat("Resolução do baseline_PC.tif:\n")
  print(res(r_base))
}

# futuros
fut_dirs <- list.dirs("out/pca_vars/futures", recursive = FALSE)
for (fd in fut_dirs) {
  fut_file <- file.path(fd, "future_PC.tif")
  if (file.exists(fut_file)) {
    r_fut <- terra::rast(fut_file)
    cat("\nResolução de", fut_file, ":\n")
    print(res(r_fut))
  }
}

################################################################################
# ==============================================================================
# 10 - (Opicional) Gráficos
# ==============================================================================
# ------------------------------------------------------------------------------
# 1) Variância explicada (scree plot)
# ------------------------------------------------------------------------------

eig  <- pca$sdev^2
var_exp <- eig/sum(eig)*100
cumvar  <- cumsum(var_exp)

plot(var_exp, type="b", pch=19,
     xlab="Componente Principal", ylab="% Variância Explicada",
     main="Scree Plot")


# Destacar até o PC6
abline(v=6, lty=3, col="gray50")
points(1:6, var_exp[1:6], col="blue", pch=19, cex=1.2)
text(6.3, cumvar[6], labels=paste0(round(cumvar[6],1),"%"), col="blue", pos=4)

# Linha da variância acumulada
lines(cumvar, type="b", col="red", pch=19)

# Linhas horizontais nos 6 primeiros PCs
abline(h=var_exp[1:6], lty=3, col="gray50")

legend("topright", legend=c("Individual","Acumulada","Cutoff PC6"),
       col=c("black","red","blue"), lty=c(1,1,2), pch=c(19,19,19))

# ------------------------------------------------------------------------------
# 2) Loadings (pesos das variáveis em cada eixo)
# ------------------------------------------------------------------------------
loadings <- pca$rotation[,1:6]  # só os 6 primeiros PCs que você manteve
round(loadings, 3)              # arredondar para ver melhor

# salvar como CSV se quiser
write.csv(round(loadings,3), "out/pca_vars/loadings_PCs.csv")

# ------------------------------------------------------------------------------
# 3) Heatmap dos loadings
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)

df_load <- as.data.frame(loadings)
df_load$var <- rownames(df_load)
df_m <- melt(df_load, id.vars="var")

ggplot(df_m, aes(x=variable, y=var, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
  theme_minimal() +
  labs(title="Contribuição das variáveis em cada PC",
       x="Componente Principal", y="Variável bioclimática")

# ------------------------------------------------------------------------------
# 4) Biplot clássico (PC1 x PC2)
# ------------------------------------------------------------------------------
biplot(pca, scale=0.5, cex=0.6, main="Biplot PCA - Baseline (PC1 x PC2)")


# scree plot
# Autovalores e variância explicada
eig  <- pca$sdev^2
var_exp <- eig/sum(eig)*100
cumvar  <- cumsum(var_exp)

# Plot básico
plot(var_exp, type="b", pch=19,
     xlab="Componente Principal", ylab="% Variância Explicada",
     main="Scree Plot")

# Linha da variância acumulada
lines(cumvar, type="b", col="red", pch=19)

# Destacar até o PC6
abline(v=6, lty=2, col="blue")
points(1:6, var_exp[1:6], col="blue", pch=19, cex=1.2)
text(6.3, cumvar[6], labels=paste0(round(cumvar[6],1),"%"), col="blue", pos=4)

# Linha da variância acumulada
lines(cumvar, type="b", col="red", pch=19)

legend("topright", legend=c("Individual","Acumulada","Cutoff PC6"),
       col=c("black","red","blue"), lty=c(1,1,2), pch=c(19,19,NA))