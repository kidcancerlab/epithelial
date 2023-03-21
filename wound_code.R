library(cowplot) 
library(cluster)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork) 
library(reticulate)
library(Seurat)
library(sctransform)
 
###load raw
setwd("/Applications/scRNA_seq raw files/sc k7m2 f420/S0066-B6-Ctl-lung/filtered_feature_bc_matrix")
raw.B6=Read10X("/Applications/scRNA_seq raw files/sc k7m2 f420/S0066-B6-Ctl-lung/filtered_feature_bc_matrix")
setwd("/Applications/scRNA_seq raw files/sc k7m2 f420/S0068-F420-Epcam-enriched-mets use this/filtered_feature_bc_matrix")
raw.F420=Read10X("/Applications/scRNA_seq raw files/sc k7m2 f420/S0068-F420-Epcam-enriched-mets use this/filtered_feature_bc_matrix")
setwd("/Applications/scRNA_seq raw files/sc k7m2 f420/S0074-Balb-C-Ctl-lung/filtered_feature_bc_matrix")
raw.balbc=Read10X("/Applications/scRNA_seq raw files/sc k7m2 f420/S0074-Balb-C-Ctl-lung/filtered_feature_bc_matrix")
setwd("/Applications/scRNA_seq raw files/sc k7m2 f420/S0075-Balb-C-K7M2-Epcam-enriched/filtered_feature_bc_matrix")
raw.k7m2=Read10X("/Applications/scRNA_seq raw files/sc k7m2 f420/S0075-Balb-C-K7M2-Epcam-enriched/filtered_feature_bc_matrix")
C57BL6<-CreateSeuratObject(counts = raw.B6, min.cells = 3, min.features =200 )
C57BL6$sample <- "C57BL6"
C57BL6[["percent.mt"]] <- PercentageFeatureSet(C57BL6, pattern = "^mt-")
F420<-CreateSeuratObject(counts = raw.F420, min.cells = 3, min.features =200 )
F420$sample <- "F420"
F420[["percent.mt"]] <- PercentageFeatureSet(F420, pattern = "^mt-")
###merge control and met datasets 
B6.f420.combined=merge(C57BL6, y=F420, add.cell.ids = c("C57BL6","F420"))
###scTransform merged data sets to include all genes in scaled data 
B6.f420.combined<- SCTransform(B6.f420.combined, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE)
###dim reduction
B6.f420.combined <- RunPCA(B6.f420.combined, verbose = FALSE)
B6.f420.combined <- RunUMAP(B6.f420.combined, dims = 1:30, verbose = FALSE)
###cluster low resolution "birds eye view" of lung microenvironment
B6.f420.combined<- FindNeighbors(B6.f420.combined, dims = 1:30, verbose = FALSE)
B6.f420.combined<- FindClusters(B6.f420.combined,resolution=0.2, verbose = FALSE)
###initial DimPlot
DimPlot(B6.f420.combined, label = TRUE, split.by = "sample", pt.size = 1,cols = scales::alpha(plot_cols, 0.6)) +coord_fixed()+theme(axis.text =element_text(size = 5), axis.title =element_text(size = 8), plot.title = element_text(size=10, face = "bold", hjust = 0.5), strip.text = element_text(size=8, face = "bold"), plot.subtitle = element_text(size=8, hjust = 0.5),legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"), strip.background =element_blank())
###Marker genes for each cluster 
initial.markers <- FindAllMarkers(B6.f420.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
###Remove cycling and dead cells clusters 3 and 13 noted by high expression of mitochondrial and cell cycling genes respectively
B6.f420.combined.subset<-subset(B6.f420.combined, idents = c("3", "13"), invert=T)
###Cluster assignments based on canonical markers and marker genes. merge clustes 
B6.f420.combined.subset<-RenameIdents(B6.f420.combined.subset, `0`="Macrophage", `1`="Alveolar macrophage",`2`="Endothelial cell",`4`="T cell", `5`="Dendritic cell", `6`="Proximal airway cell", `7`="Distal airway cell", `8`="Fibroblast", `9`="B cell", `10`="Endothelial cell", `11`="Dendritic cell", `12`="Granulocyte",  `14`=" F420", `15`="Erythrocyte", `16`="Smooth muscle cell", `17`="Adipocyte", `18`="Mesothelial cell")
###add cluster identity to metadata
B6.f420.combined.subset$celltype.sample <- paste(Idents(B6.f420.combined.subset), B6.f420.combined.subset$sample, sep = "_") 
B6.f420.combined.subset$celltype <- Idents(B6.f420.combined.subset)
###markers for plotting 
B6.F420.combined.plotmarkers<-c("Itgam", "Fcgr1", "Cd68", "Adgre1", "Ccr2", "Cx3cr1", "Ly6c2",  "Itgax", "Mertk", "Marco","Siglecf","Egr2", "Bhlhe41", "Fabp1", "Olr1", "Alox5" ,  "Vwf","Tek","Acer2", "Cavin3","Pecam1",  "Emcn", "Lyve1", "Ptprb","Plvap" ,"Cd3d", "Cd4", "Cd8b1", "Tcf7","Il7r", "Lef1","Dusp10", "Itk","Cd207", "Cd209a", "Ccr7", "Ccl17", "Ccl22", "Wdfy4", "Clec4b1", "H2-DMb1", "Etv3", "Itgae", "Xcr1","Epcam", "Scgb1a1", "Tuba1a","Foxj1", "Sftpd","Lamp3", "Etv5", "Ager", "Aqp5", "Col1a2", "Col3a1", "Col1a1", "Loxl1", "Pdgfra","Pdgfrb", "Fn1", "Ltbp4", "Bgn", "Cd19", "Cd79a", "Cd79b", "Iglc2", "Ighm", "Pax5" , "Cd33", "Csf3r","S100a9", "Cxcr2","Mxd1","Ptafr" , "Steap1","Crabp2","Wisp1","Dclk1", "Twist1","Il1rl1","Ctgf","Hbb-bs", "Hba-a1","Myh11","Cox4i2", "Pde5a", "Cfd", "Car3", "Adipoq","Msln")
###Addmodules GSEA with msigdbr
library(msigdbr)
all_gene_sets = msigdbr(species = "Mus musculus")
Inflammation_list<-filter(all_gene_sets, gs_name=="HALLMARK_INFLAMMATORY_RESPONSE")
Inflammation_list<-list(c(Inflammation_list$gene_symbol))
B6.f420.combined.subset<-AddModuleScore(B6.f420.combined.subset, features =Inflammation_list, name = 'Inflammation_features')
Lung_fibrosis_list<-filter(all_gene_sets, gs_name=="WP_LUNG_FIBROSIS")
Lung_fibrosis_list<-list(c(Lung_fibrosis_list$gene_symbol))
B6.f420.combined.subset<-AddModuleScore(B6.f420.combined.subset, features =Lung_fibrosis_list, name = 'Lung_fibrosis_features')
TNFA_list<-filter(all_gene_sets, gs_name=="HALLMARK_TNFA_SIGNALING_VIA_NFKB")
TNFA_list<-list(c(TNFA_list$gene_symbol))
B6.f420.combined.subset<-AddModuleScore(B6.f420.combined.subset, features =TNFA_list, name = 'TNFA_features')
ECM_degrad_list<-filter(all_gene_sets, gs_name=="REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX")
ECM_degrad_list<-list(c(ECM_degrad_list$gene_symbol))
B6.f420.combined.subset<-AddModuleScore(B6.f420.combined.subset, features =ECM_degrad_list, name = 'ECM_degrad_features')
ECM_org_list<-filter(all_gene_sets, gs_name=="REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION")
ECM_org_list<-list(c(ECM_org_list$gene_symbol))
B6.f420.combined.subset<-AddModuleScore(B6.f420.combined.subset, features =ECM_org_list, name = 'ECM_org_features')
Wound_response_list<-filter(all_gene_sets, gs_name=="GOBP_RESPONSE_TO_WOUNDING")
Wound_response_list<-list(c(Wound_response_list$gene_symbol))
B6.f420.combined.subset<-AddModuleScore(B6.f420.combined.subset, features =Wound_response_list, name = 'Wound_response_features')
###subset epithelial cells and re-normalize
b6.f420.combined.AEC=subset(B6.f420.combined.subset, idents = "Distal airway cell")
b6.f420.combined.AEC<- SCTransform(b6.f420.combined.AEC, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE)
###dim reduction
b6.f420.combined.AEC<- RunPCA(b6.f420.combined.AEC, verbose = FALSE)
b6.f420.combined.AEC<- RunUMAP(b6.f420.combined.AEC, dims = 1:30, verbose = FALSE)
###cluster
b6.f420.combined.AEC<- FindNeighbors(b6.f420.combined.AEC, dims = 1:30, verbose = FALSE)
b6.f420.combined.AEC<- FindClusters(b6.f420.combined.AEC,resolution=0.3, verbose = FALSE)
###remove contaminating non-epithelial cells-immune cells
b6.f420.combined.AEC=subset(b6.f420.combined.AEC, idents = c("1","2", "3", "6"))
###re-cluster with pure epithelial cell populations 
b6.f420.combined.AEC_1<- SCTransform(b6.f420.combined.AEC, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE)
b6.f420.combined.AEC_1<- RunPCA(b6.f420.combined.AEC_1, verbose = FALSE)
b6.f420.combined.AEC_1<- RunUMAP(b6.f420.combined.AEC_1, dims = 1:30, verbose = FALSE)
b6.f420.combined.AEC_1<- FindNeighbors(b6.f420.combined.AEC_1, dims = 1:30, verbose = FALSE)
b6.f420.combined.AEC_1<- FindClusters(b6.f420.combined.AEC_1,resolution=0.6, verbose = FALSE)
initial.epimarkers <- FindAllMarkers(b6.f420.combined.AEC_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
###remove contaminating proximal airway cells
b6.f420.combined.AEC_1=subset(b6.f420.combined.AEC_1, idents = "5", invert=T)
###cluster identification based on published data sets-PMIDs 32750316,32661339,32678092 
b6.f420.combined.AEC_1<-RenameIdents(b6.f420.combined.AEC_1, `0`="pAEC2", `1`="DATP", `2`="AEC2", `3`="AEC1", `4`="pAEC2", `6`="AEC1", `7`="cAEC2")
b6.f420.combined.AEC_1$celltype.sample <- paste(Idents(b6.f420.combined.AEC_1), b6.f420.combined.AEC_1$sample, sep = "_") 
b6.f420.combined.AEC_1$celltype <- Idents(b6.f420.combined.AEC_1)
### add module for DATP
p53_list<-filter(all_gene_sets, gs_name=="HALLMARK_P53_PATHWAY")
p53_list<-list(c(p53_list$gene_symbol))
b6.f420.combined.AEC_1<-AddModuleScore(b6.f420.combined.AEC_1, features =p53_list, name = 'p53_features')
senescence_list<-filter(all_gene_sets, gs_name=="FRIDMAN_SENESCENCE_UP")
senescence_list<-list(c(senescence_list$gene_symbol))
b6.f420.combined.AEC_1<-AddModuleScore(b6.f420.combined.AEC_1, features =senescence_list, name = 'senescence_features')
Distalairway.markers<-c("Lcn2","Lrg1","Cxcl17","Glrx","Krt8","Sfn","Krt19","Cldn4","Trp53","Ctgf","Sftpa1","Sftpb","Sftpc","Lpcat1","Etv5", "Abca3", "Pdpn", "Aqp5", "Hopx", "Cav1","Ager", "Akap5", "Clic5", "Top2a", "Mki67", "Cdk1")
###scTransform macs 
f420.macs<- SCTransform(f420.macs, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE)
###dim reduction
f420.macs<- RunPCA(f420.macs, verbose = FALSE)
f420.macs<- RunUMAP(f420.macs, dims = 1:30, verbose = FALSE)
###cluster 
f420.macs<- FindNeighbors(f420.macs, dims = 1:30, verbose = FALSE)
f420.macs<- FindClusters(f420.macs,resolution=0.4, verbose = FALSE)
mac.markers<-c("Itgam","Itgax", "Cd14", "Fcgr3", "Fcgr1", "Cd68", "Adgre1", "Ccr2", "Cx3cr1", "Ly6c1", "Ly6c2", "H2-Ab1", "H2-DMb1")
###Create and merge BALBC and k7M2 seurat objects followed by scTransform
BALBC<-CreateSeuratObject(counts = raw.balbc, min.cells = 3, min.features =200 )
BALBC$sample <- "BALBC"
BALBC[["percent.mt"]] <- PercentageFeatureSet(BALBC, pattern = "^mt-")
K7M2<-CreateSeuratObject(counts = raw.k7m2, min.cells = 3, min.features =200 )
K7M2$sample <- "K7M2"
K7M2[["percent.mt"]] <- PercentageFeatureSet(K7M2, pattern = "^mt-")
###merge control and met datasets 
balbc.k7m2.combined=merge(BALBC, y=K7M2, add.cell.ids = c("BALBC","K7M2"))
###scTransform merged data sets to include all genes in scaled data 
balbc.k7m2.combined<- SCTransform(balbc.k7m2.combined, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE)
###dim reduction
balbc.k7m2.combined <- RunPCA(balbc.k7m2.combined, verbose = FALSE)
balbc.k7m2.combined <- RunUMAP(balbc.k7m2.combined, dims = 1:30, verbose = FALSE)
###cluster low resolution "birds eye view" of lung microenvironment
balbc.k7m2.combined<- FindNeighbors(balbc.k7m2.combined, dims = 1:30, verbose = FALSE)
balbc.k7m2.combined<- FindClusters(balbc.k7m2.combined,resolution=0.2, verbose = FALSE)
###Remove  dead cells cluster 9  noted by high expression of mitochondrial as well as small clusters 18-20
balbc.k7m2.combined.subset<-subset(balbc.k7m2.combined, idents = c("9", "18", "19", "20"), invert=T)
balbc.k7m2.combined.subset<-RenameIdents(balbc.k7m2.combined.subset, `0`="Alveolar macrophage", `1`="Macrophage",`2`="T cell",`3`="Endothelial cell", `4`="K7M2", `5`="Granulocyte", `6`="Dendritic cell", `7`="B cell", `8`="Distal airway cell", `10`="Fibroblast", `11`="Endothelial cell", `12`="Macrophage",`13`="T cell",  `14`=" Proximal airway cell", `15`="Alveolar macrophage", `16`="Endothelial cell", `17`="Adipocyte")
###add cluster identity to metadata
balbc.k7m2.combined.subset$celltype.sample <- paste(Idents(balbc.k7m2.combined.subset), balbc.k7m2.combined.subset$sample, sep = "_") 
balbc.k7m2.combined.subset$celltype <- Idents(balbc.k7m2.combined.subset)
###add module scores
balbc.k7m2.combined.subset<-AddModuleScore(balbc.k7m2.combined.subset, features =Inflammation_list, name = 'Inflammation_features')
balbc.k7m2.combined.subset<-AddModuleScore(balbc.k7m2.combined.subset, features =Lung_fibrosis_list, name = 'Lung_fibrosis_features')
balbc.k7m2.combined.subset<-AddModuleScore(balbc.k7m2.combined.subset, features =TNFA_list, name = 'TNFA_features')
balbc.k7m2.combined.subset<-AddModuleScore(balbc.k7m2.combined.subset, features =ECM_degrad_list, name = 'ECM_degrad_features')
balbc.k7m2.combined.subset<-AddModuleScore(balbc.k7m2.combined.subset, features =ECM_org_list, name = 'ECM_org_features')
balbc.k7m2.combined.subset<-AddModuleScore(balbc.k7m2.combined.subset, features =Wound_response_list, name = 'Wound_response_features')
###subset epithelial cells and re-normalize
balbc.k7m2.combined.AEC=subset(balbc.k7m2.combined.subset, idents = "Distal airway cell")
balbc.k7m2.combined.AEC<- SCTransform(balbc.k7m2.combined.AEC, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE)
###dim reduction
balbc.k7m2.combined.AEC<- RunPCA(balbc.k7m2.combined.AEC, verbose = FALSE)
balbc.k7m2.combined.AEC<- RunUMAP(balbc.k7m2.combined.AEC, dims = 1:30, verbose = FALSE)
###cluster
balbc.k7m2.combined.AEC<- FindNeighbors(balbc.k7m2.combined.AEC, dims = 1:30, verbose = FALSE)
balbc.k7m2.combined.AEC<- FindClusters(balbc.k7m2.combined.AEC,resolution=0.6, verbose = FALSE)
###remove contaminating non-epithelial cells-immune cells
balbc.k7m2.combined.AEC=subset(balbc.k7m2.combined.AEC, idents = c("2","3", "4", "5", "6", "7"))
###re-cluster with pure epithelial cell populations 
balbc.k7m2.combined.AEC_1<- SCTransform(balbc.k7m2.combined.AEC, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE)
balbc.k7m2.combined.AEC_1<- RunPCA(balbc.k7m2.combined.AEC_1, verbose = FALSE)
balbc.k7m2.combined.AEC_1<- RunUMAP(balbc.k7m2.combined.AEC_1, dims = 1:30, verbose = FALSE)
balbc.k7m2.combined.AEC_1<- FindNeighbors(balbc.k7m2.combined.AEC_1, dims = 1:30, verbose = FALSE)
balbc.k7m2.combined.AEC_1<- FindClusters(balbc.k7m2.combined.AEC_1,resolution=0.6, verbose = FALSE)
###cluster identification based on published data sets-PMIDs 32750316,32661339,32678092 
balbc.k7m2.combined.AEC_1<-RenameIdents(balbc.k7m2.combined.AEC_1, `0`="pAEC2", `1`="DATP", `2`="AEC2", `3`="AEC2", `4`="AEC1",`5`="cAEC2")
balbc.k7m2.combined.AEC_1$celltype.sample <- paste(Idents(balbc.k7m2.combined.AEC_1), balbc.k7m2.combined.AEC_1$sample, sep = "_") 
balbc.k7m2.combined.AEC_1$celltype <- Idents(balbc.k7m2.combined.AEC_1)
###subset k7m2 macrophages
Idents(balbc.k7m2.combined.subset)="celltype.sample"
###remove "K7M2" from BALBC
balbc.k7m2.combined.subset=subset(balbc.k7m2.combined.subset, idents = "K7M2_BALBC", invert=T)
###scTransform macs 
k7m2.macs<- SCTransform(k7m2.macs, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE)
###dim reduction
k7m2.macs<- RunPCA(k7m2.macs, verbose = FALSE)
k7m2.macs<- RunUMAP(k7m2.macs, dims = 1:30, verbose = FALSE)
###cluster 
k7m2.macs<- FindNeighbors(k7m2.macs, dims = 1:30, verbose = FALSE)
k7m2.macs<- FindClusters(k7m2.macs,resolution=0.4, verbose = FALSE)
k7m2.mac.markers<-FindAllMarkers(k7m2.macs, only.pos = T, logfc.threshold = 0.5, min.pct = 0.25)
###mac cluster annotation based on PMID: 35690521, 32849616, 35365629
M1up_list<-filter(all_gene_sets, gs_name=="GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP")
M1up_list<-list(c(M1up_list$gene_symbol))
f420.macs<-AddModuleScore(f420.macs, features =M1up_list, name = 'm1_features')
phago_list<-filter(all_gene_sets, gs_name=="GOBP_POSITIVE_REGULATION_OF_PHAGOCYTOSIS")
phago_list<-list(c(phago_list$gene_symbol))
f420.macs<-AddModuleScore(f420.macs, features =phago_list, name = 'phagocytosis_features')
###remove contaminating endothelial, T cells
f420.mac.subset<-subset(f420.macs, idents = c("7", "9"), invert=T)
###mac markers based on 35690521 LA-TAM, ANGIO-TAM, IFN-TAM, INFLAM-TAM, PROLIF-TAM, Classical monocyte, nonclassical
TAM.markers<-c("Spp1", "Vegfa", "Fabp5", "Gpnmb", "Lgals3","Apoe","Pf4", "Mrc1", "Cd206","Arg1", "Cxcl1", "Cxcl2", "Cxcl3", "Il1a", "Il1b", "Inhba", "Ccl20", "Cxcl9", "Cxcl10", "Ifit1", "Ifitm1", "Tnfsf10", "Fn1", "Hmox1", "C1qa", "C1qb","Mki67","Ace", "Itgal", "Ear2", "Spn", "Itgam", "Itgax", "Cd14", "Fcgr3", "Fcgr1", "Cd68","Cd74", "Cd80", "Cd86", "Ccr2", "Cx3cr1", "Ly6c1", "Ly6c2", "H2-Ab1", "H2-DMb1")
###cluster annotation based on PMID: 35690521
dafdsfafs
