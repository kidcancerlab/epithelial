library(ggplot2)
library(gridExtra)
library(rrrSingleCellUtils)
library(Seurat)
library(tidyverse)

# Load and label data
ctl <- tenx_load_qc("~/lab/Counts/S0066/filtered_feature_bc_matrix/")
ctl <- ctl <- subset(ctl, subset = nFeature_RNA > 300 & nCount_RNA < 45000 & percent.mt <16)
ctl$src <- "S0066"
ctl$cond <- "Unaffected"
ctl$

met1 <- tenx_load_qc("~/lab/Counts/S0067/filtered_feature_bc_matrix/")
met1 <- subset(met1, subset = nFeature_RNA > 300 & nCount_RNA < 45000 & percent.mt <16)
met1$src <- "S0067"
met1$cond <- "Met-associated"

met2 <- tenx_load_qc("~/lab/Counts/S0068/filtered_feature_bc_matrix/")
met2 <- subset(met2, subset = nFeature_RNA > 300 & nCount_RNA < 45000 & percent.mt <16)
met2$src <- "S0068"
met2$cond <- "Met-associated"

# Merge datasets and plot
lung <- merge(ctl, y = c(met1, met2), add.cell.ids = c("ctl", "met1", "met2"), project = "F420")

lung <- NormalizeData(lung) %>%
  ScaleData()  %>%
  FindVariableFeatures()
lung <- RunPCA(lung, features = VariableFeatures(lung), verbose = F) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.1) %>%
  RunUMAP(dims = 1:20)

DimPlot(lung, group.by = "src")
DimPlot(lung)

FeaturePlot(lung, features = "Krt8")
FeaturePlot(lung, features = "Pdpn")

# Subset out the epithelial cells
epi <- subset(lung, idents = 8)
p1 <- DimPlot(epi, group.by = "cond") +
  coord_fixed()
print(p1)

epi <- NormalizeData(epi) %>%
  ScaleData()  %>%
  FindVariableFeatures()
epi <- RunPCA(epi, features = VariableFeatures(lung), verbose = F) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.1) %>%
  RunUMAP(dims = 1:20)

DimPlot(epi)
epi$cond <- factor(epi$cond, levels = c("Unaffected", "Met-associated"))
p2 <- DimPlot(epi, split.by = "cond") +
  coord_fixed()
print(p2)

# Create the expression dot plot
gl <- c("Lcn2", "Lrg1", "Ptges", "Glrx", "Krt8", "Krt19", "Sprr1a",
  "Nupr1", "Sfn", "Fn1", "Cldn4", "Cdkn1a", "Itgb6", "Ager", "Hopx",
  "Igfbp2", "Cav1", "Cavin2", "Itgb1", "Fabp5", "Sftpb", "Sftpc",
  "Lamp3", "Etv5", "Abca3", "Fgfr2", "Mki67", "Prc1", "Cdk1", "Cenpe")
p3 <- DotPlot(epi, feature = gl) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p3)

# Quantify cell types in met-associated vs unaffected lung
pops <- table(epi$cond, epi$seurat_clusters)
pops <- pops / rowSums(pops)
pops <- as_tibble(pops, .name_repair = "minimal")
colnames(pops) <- c("Location", "Cluster", "Cells")

# Reorder the labels for the x-axis and plot it
pops$Location <- factor(pops$Location, 
  levels = c("Unaffected", "Met-associated"))

ggplot(pops, aes(x = Location, y = Cells, fill = Cluster)) +
  geom_bar(stat = "identity")

