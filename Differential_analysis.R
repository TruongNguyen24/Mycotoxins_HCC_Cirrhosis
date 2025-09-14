# loading library

library(limma)
library(sva)
library(SummarizedExperiment)
library(tidyr)
library(ggplot2)

# Dataset of metabolomics: SummarizedExperiment object

se_50_log #this is SE object of metabolomics data including rowData of metabolite annotation, colData of sample annotation 

log_mat <- assay(se_50_log, "abundance")
assay(se_50, "logcounts") <- log_mat
se_50
scaled_mat <- t( scale(t(log_mat), center = TRUE, scale = TRUE) )

assay(se_50, "log_zscore") <- scaled_mat

hist(assay(se_50, "log_zscore") [,2])
hist(assay(se_50, "logcounts") [,2])
hist(assay(se_50, "abundance") [,2])

metadata <- colData(se_50)

# PCA to test the batch effects
mat_pca <- assay(se_50, "log_zscore")         # features × samples

pca <- prcomp(t(mat_pca), center = FALSE, scale. = FALSE)


scores <- as.data.frame(pca$x[, 1:3])                  # PC1–PC3
scores <- cbind(scores, as.data.frame(colData(se_50)))    # add batch, group, etc.

var_explained <- pca$sdev^2 / sum(pca$sdev^2)
barplot(100 * var_explained[1:10],
        names.arg = paste0("PC", 1:10),
        ylab = "% variance", las = 2)

library(ggplot2)
scores$run_id <- as.factor(scores$run_id)

ggplot(scores, aes(PC1, PC2, colour = run_id, shape = liver_control)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = run_id),
               type   = "norm",   # assumes roughly elliptical clusters
               level  = 0.95,     # 95 % “circle”
               linetype = 2,      # dashed outline
               size     = 0.6) +
  labs(
    title = "PCA of metabolomics – coloured by batch",
    x = sprintf("PC1 (%.1f%%)", 100 * var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * var_explained[2])
  ) +
  theme_classic()

# Using ComBat for batch correction

design <- model.matrix(~ hcc_cir_ctrl + sex + age + smoked_ever + hbsag_result + af_bin + run_id, metadata)
mod0 <- model.matrix(~ sex + age + smoked_ever + hbsag_result + af_bin + run_id, data = metadata)
common  <- intersect(rownames(metadata), rownames(design))
length(common)
metadata2 <- metadata[rownames(metadata) %in% common, ]
identical(rownames(metadata2), rownames(modcombat))

batch <- as.factor(metadata2$run_id)

combat_edata <- ComBat(dat= mat_filter, batch=batch, mod=design, par.prior=TRUE, prior.plots= FALSE, BPPARAM = bpparam("SerialParam"))

pca <- prcomp(t(combat_edata), center = FALSE, scale. = FALSE)

scores <- as.data.frame(pca$x[, 1:3])                  # PC1–PC3

metadatanew <- as.data.frame(colData(se_50))

intersect <- intersect(rownames(scores), rownames(metadatanew))

length(intersect)
metadatanew <- metadatanew[rownames(metadatanew) %in% intersect, ]
identical(rownames(scores), rownames(metadatanew))

scores <- cbind(scores, metadatanew)    # add batch, group, etc.

var_explained <- pca$sdev^2 / sum(pca$sdev^2)
barplot(100 * var_explained[1:10],
        names.arg = paste0("PC", 1:10),
        ylab = "% variance", las = 2)

library(ggplot2)
scores$run_id <- as.factor(scores$run_id)

scores$hcc_cirr_control <- ifelse(scores$pcat == 1, "Cirrhosis", ifelse(scores$pcat == 2, "HCC", "control"))

ggplot(scores, aes(PC1, PC2, colour = run_id, shape = hcc_cirr_control)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = run_id),
               type   = "norm",   # assumes roughly elliptical clusters
               level  = 0.95,     # 95 % “circle”
               linetype = 2,      # dashed outline
               size     = 0.6) +
  labs(
    title = "PCA of metabolomics – coloured by batch",
    x = sprintf("PC1 (%.1f%%)", 100 * var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * var_explained[2])
  ) +
  theme_classic()

# Differential analysis with limma

design2 <- model.matrix(~ hcc_cir_ctrl + sex + age + smoked_ever + hbsag_result + af_bin + run_id, metadata)
common  <- intersect(colnames(mat), rownames(design2))
missing_rows <- setdiff(colnames(mat), rownames(design2))
missing_rows

fit2 <- lmFit(mat, design2, na.action = na.exclude)
fit2 <- eBayes(fit2, robust = TRUE)
fit2 

dt2 <- decideTests(fit2, adjust.method = "BH", p.value = 0.05) 
tab2 <- summary(dt2)
tab2

res2_cirr  <- topTable(fit2, coef = "hcc_cir_ctrlCirrhosis", adjust = "BH", number = Inf)
res2_cirr

res2_hcc  <- topTable(fit2, coef = "hcc_cir_ctrlHCC", adjust = "BH", number = Inf)
res2_hcc

