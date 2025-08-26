# PCA WordClim (Brasil)

Este repositório contém um script em **R** para realizar a **Análise de Componentes Principais (PCA)** das variáveis bioclimáticas do **WorldClim v2.1**, recortadas para o território brasileiro.  
O pipeline gera **variáveis sintéticas (PCs)** que eliminam colinearidade e retêm ≥95% da variância climática, permitindo comparabilidade entre clima atual e projeções futuras (CMIP6).

---

## Funcionalidades principais
1. Leitura das **19 variáveis bioclimáticas** do WorldClim v2.1 (baseline: 1970–2000).  
2. Leitura de cenários **futuros CMIP6 (2061–2080, SSP5-8.5)**:  
   - MIROC6  
   - EC-Earth3-Veg  
   - IPSL-CM6A-LR  
3. **Recorte e máscara** para os limites do Brasil (Natural Earth).  
4. Execução da **PCA** no baseline, com seleção de PCs até atingir ≥95% da variância explicada.  
5. **Projeção dos cenários futuros** no espaço PCA definido pelo baseline.  
6. Geração de **produtos analíticos**:  
   - Scree plot (variância individual e acumulada).  
   - Heatmap dos loadings (contribuição de cada BIO por eixo).  
   - Biplot (PC1 × PC2).  
7. **Exportação dos resultados**:  
   - `baseline_PC.tif` → PCs do baseline.  
   - `future_PC.tif` → PCs dos futuros, um por GCM.  
   - `loadings_PCs.csv` → pesos das variáveis em cada PC.  

---

