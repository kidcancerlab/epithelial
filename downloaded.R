library(tidyverse)
library(Seurat)
library(colorspace)
if(!"rrrSingleCellUtils" %in% installed.packages()) {
  devtools::install_github("kidcancerlab/rrrSingleCellUtils")
}
if(!"GEOquery" %in% installed.packages()) {
  BiocManager::install("GEOquery")
}
library(rrrSingleCellUtils)
library(GEOquery)

# Download the single cell data sets from GSE145031
getGEOSuppFiles("GSE145031", fetch_files = T)
untar("GSE145031/GSE145031_RAW.tar", exdir = "GSE145031")
file.remove("GSE145031/GSE145031_RAW.tar")

days <- c("D0", "D14", "D28")
cells <- c("AT2", "Other")

# Create directories, organize the downloaded raw data, and create Seurat object
for(j in 1:2) {
  dir.create(str_c("GSE145031", "/", cells[j]))
  for(i in 1:3) {
    dir.create(str_c("GSE145031", "/", cells[j], "/", days[i]))
  }
} 

sobjs <- list()
ser <- 4304609
k <- 1
for(i in 1:3) {
  for(j in 1:2) {
    orig <- str_c('GSE145031\\GSM', ser, '_*')
    dest <- str_c('GSE145031\\', cells[j], '\\', days[i])
    shell(str_c('move ', orig, " ", dest))
    
    shell(str_c('ren ', dest, "\\*mtx.gz matrix.mtx.gz"))
    shell(str_c('ren ', dest, "\\*barcodes.tsv.gz barcodes.tsv.gz"))
    shell(str_c('ren ', dest, "\\*gene.tsv.gz features.tsv.gz"))
    
    sobjs[[k]] <- Read10X(dest) %>%
      CreateSeuratObject(min.cells = 5, min.features = 800)
    print(VlnPlot(sobjs[[k]], 
                  features = c("nFeature_RNA", "nCount_RNA"), 
                  ncol = 2))
    sobjs[[k]] <- subset(sobjs[[k]], nCount_RNA < 12000)
    sobjs[[k]]$cell <- cells[j]
    sobjs[[k]]$time <- days[i]
    
    k = k+1
    ser = ser+1
  }
}

# Basic processing and clustering, cell type assignment, graph
comb <- merge(sobjs[[1]],
              y = c(sobjs[[2]], sobjs[[3]], sobjs[[4]], sobjs[[5]], sobjs[[6]]),
              add.cell.ids = c("1", "2", "3", "4", "5", "6"),
              project = "Lung Injury")

comb <- comb %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures()

comb <- comb %>%
  RunPCA(features = VariableFeatures(comb)) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.3) %>%
  RunUMAP(dims = 1:20)

DimPlot(comb, pt.size = 1, label = T) +
  coord_fixed() +
  ggtitle("All Cells Combined") 

DimPlot(comb, pt.size = 1, group.by = "cell", split.by = "time") +
  coord_fixed() +
  ggtitle("All Cells by Sort")

# Identify the different cell types using a mouse lung single cell atlas
# Automatic prediction of cell types using SingleR
library(SingleR)
library(SingleCellExperiment)

# Load the mouse lung reference library
lung_ref <- readRDS("R:/RESRoberts/Bioinformatics/GenRef/scRNA-Altas/MurLung/mLungRef-simple.rds")

# Convert the Seurat object to a SCE object
comb_sce <- as.SingleCellExperiment(comb)

# Make cell type predictions using SingleR
comb_pred <- SingleR(test = comb_sce,
                     ref = lung_ref,
                     labels = lung_ref$label)

# Transfer labels back to the Seurat object and plot
comb$singleR <- comb_pred$labels
Idents(comb) <- comb$singleR

DimPlot(comb, pt.size = 1, label = T) +
  coord_fixed() +
  theme(legend.position = "none") +
  ggtitle("All Cells - Referenced")


# Generate feature plots with key genes
FeaturePlot(comb,
            features = "Cd14",
            pt.size = 1,
            split.by = "time",
            order = T, 
            cols = c("lightgoldenrodyellow", "red3"),
            min.cutoff = 1) +
  coord_fixed() +
  theme(legend.position = "none")

FeaturePlot(comb,
            features = "Krt8",
            pt.size = 1,
            split.by = "time",
            order = T, 
            cols = c("lightgoldenrodyellow", "red3"),
            min.cutoff = 1) +
  coord_fixed() +
  theme(legend.position = "none")

FeaturePlot(comb,
            features = "Il1b",
            pt.size = 1,
            split.by = "time",
            order = T, 
            cols = c("lightgoldenrodyellow", "red3"),
            min.cutoff = 1) +
  coord_fixed() +
  theme(legend.position = "none")

FeaturePlot(comb,
            features = "Pdpn",
            pt.size = 1,
            split.by = "time",
            order = T, 
            cols = c("lightgoldenrodyellow", "red3"),
            min.cutoff = 1) +
  coord_fixed() +
  theme(legend.position = "none")
