#Loading libraries 
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(MAST)

setwd("C:/Documents and Settings/PC/Documents/RStudio/D125fetal_retina/")

#Loading data
D125Cfetal <- Read10X("C:/Documents and Settings/PC/Documents/RStudio/D125fetal_retina/D125Cfetal_filtered_gene_bc_matrices/GRCh38/")
dim(D125Cfetal)

#Filtering (genes and barcodes)
D125CfetalS <- CreateSeuratObject(D125Cfetal) #we need TMEM119 & P2RY12 microglia markers
dim(D125CfetalS)

#Create the function for creating and processing the Seurat objects
ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat) #could be replaced with SCTransform
  Seurat <- RunPCA(Seurat)
  Seurat <- FindNeighbors(Seurat, dims = 1:20, reduction = "pca", verbose = F)
  Seurat <- FindClusters(Seurat, resolution = 1, verbose = F)
  Seurat <- RunUMAP(Seurat, dims = 1:20)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:20)
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}  

#Import S/G2M cell cycle genes
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

#Quantify the rb/mt proportion
D125CfetalS[["percent.rb"]] <- PercentageFeatureSet(D125CfetalS, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
D125CfetalS[["percent.mt"]] <- PercentageFeatureSet(D125CfetalS, pattern = "^MT-")

#Visualize mt & rb content 
FeatureScatter(D125CfetalS, "nCount_RNA", "percent.mt") + scale_x_log10()
FeatureScatter(D125CfetalS, "nFeature_RNA", "percent.mt") + scale_x_log10()
FeatureScatter(D125CfetalS, "nFeature_RNA", "percent.rb") + scale_x_log10()

##Normalization using NormalizeData()
D125CfetalS <- NormalizeData(D125CfetalS)

#Add cell cycle scores in metadata
D125CfetalS <- CellCycleScoring(D125CfetalS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

#Subsetting 
D125CfetalS <- subset(D125CfetalS, subset = nCount_RNA > 300 & nCount_RNA < 5000 
                      & nFeature_RNA > 300 & nFeature_RNA < 3000 
                      & percent.mt < 15 & percent.rb < 30)

options(future.globals.maxSize = 8000 * 1024^2)
D125CfetalS <- ScaleData(D125CfetalS, vars.to.regress = c("nCount_RNA", "percent.mt", "percent.rb", 
                                                          "S.Score", "G2M.Score"), verbose = T)

#If receive the error for nbins during the ScaleData, run NormalizeData before it
D125CfetalS <- ProcessSeu(D125CfetalS)

#Basic visualisation
FeaturePlot(D125CfetalS, features = c("RBPMS", "P2RY12", "TMEM119", "SOX2"))

DimPlot(D125CfetalS, reduction = "umap", label = TRUE, repel = TRUE)
FeaturePlot(D125CfetalS, reduction = "umap", features = "POU4F2")

D125CfetalS <- RenameIdents(D125CfetalS, "14" = "RGC", "18" = "RGC")
AverageExpression(D125CfetalS, features = "POU4F2")

#_______________________________________________________________________________
##Processing the 2nd dataset (D125P_fetal_retina) in the same way

#Loading data
D125Pfetal <- Read10X("C:/Documents and Settings/PC/Documents/RStudio/D125fetal_retina/D125Pfetal_filtered_gene_bc_matrices/GRCh38/")
dim(D125Pfetal)

##Filtering (genes and barcodes)
D125PfetalS <- CreateSeuratObject(D125Pfetal) 
dim(D125PfetalS)

#Quantify the rb/mt proportion
D125PfetalS[["percent.rb"]] <- PercentageFeatureSet(D125PfetalS, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
D125PfetalS[["percent.mt"]] <- PercentageFeatureSet(D125PfetalS, pattern = "^MT-")

##Normalization using NormalizeData()
D125PfetalS <- NormalizeData(D125PfetalS)

#Add cell cycle scores in metadata
D125PfetalS <- CellCycleScoring(D125PfetalS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

#Subsetting 
D125PfetalS <- subset(D125PfetalS, subset = nCount_RNA > 300 & nCount_RNA < 5000 
                      & nFeature_RNA > 300 & nFeature_RNA < 3000 
                      & percent.mt < 15 & percent.rb < 30)
dim(D125PfetalS)

options(future.globals.maxSize = 8000 * 1024^2)
D125PfetalS <- ScaleData(D125PfetalS, vars.to.regress = c("nCount_RNA", "percent.mt", "percent.rb", 
                                                          "S.Score", "G2M.Score"), verbose = T)

D125PfetalS <- ProcessSeu(D125PfetalS)

#Basic visualisation
FeaturePlot(D125PfetalS, features = c("RBPMS", "P2RY12", "TMEM119", "SOX2"))
DimPlot(D125PfetalS, reduction = "umap", label = TRUE, repel = TRUE)
#_______________________________________________________________________________
library(SeuratData)
library(patchwork)

#Add metadata
D125CfetalS <- AddMetaData(D125CfetalS, metadata = "Central", col.name = "Origin")
D125PfetalS <- AddMetaData(D125PfetalS, metadata = "Peripheral", col.name = "Origin")

###Integrating 2 datasets
#Create list of 2 Seurat objects
integration.list <- list(D125CfetalS, D125PfetalS)

#Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = integration.list)
data.anchors <- FindIntegrationAnchors(object.list = integration.list, anchor.features = features, 
                                       reduction = "rpca", k.filter = NA)

#Create an "integrated" data assay
data.integrated <- IntegrateData(anchorset = data.anchors)

#Specify that we will perform downstream analysis on the corrected data note that the
#original unmodified data still resides in the "RNA" assay
DefaultAssay(data.integrated) <- "RNA"
DefaultAssay(data.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
data.integrated <- ScaleData(data.integrated, verbose = FALSE)
data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
data.integrated <- FindNeighbors(data.integrated, reduction = "pca", dims = 1:30)
data.integrated <- FindClusters(data.integrated, resolution = 0.5)
data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30)
data.integrated <- RunTSNE(data.integrated, dims.use = 1:30)

#Basic visualization
DimPlot(data.integrated, reduction = "umap", split.by = "Origin")
DimPlot(data.integrated, reduction = "umap", group.by = "seurat_clusters", repel = TRUE, label = TRUE)
FeaturePlot(data.integrated, features = c("SOX2", "RBPMS", "RLBP1", "P2RY12", "TMEM119", "OPN1SW")) 

##Annotation
##To find the DEGs
allMarkers <- FindAllMarkers(data.integrated, max.cells.per.ident = 100, test.use = "MAST", only.pos = T)

DefaultAssay(data.integrated) <- "integrated"

allMarkers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(data.integrated, features = top10$gene) + NoLegend()

allMarkers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(data.integrated, features = top5$gene) + NoLegend()

allMarkers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
#_______________________________________________________________________________
DimPlot(data.integrated, reduction = "umap", repel = TRUE, label = TRUE)

FeaturePlot(data.integrated, features = c("RBPMS", "GAP43", "SNCG", "ELAVL4", "NEFL", "NEFM", "TUBB2A"), label = TRUE)
data.integrated <- RenameIdents(data.integrated, "9" = "RGC")
AverageExpression(data.integrated, features = c("NEFL", "POU4F2"))


FeaturePlot(data.integrated, features = c("AIF1", "CD74", "C1QA", "C1QB", "HLA-DPA1", 
                                          "HLA-DPB1", "HLA-DRA", "TGFBR1", "TMEM119", "CX3CR1", "P2RY12"), 
            repel = TRUE, label = TRUE)
data.integrated <- RenameIdents(data.integrated, "16" = "Microglia")
AverageExpression(data.integrated, features = c("P2RY12", "TMEM119", "C1QA"))


FeaturePlot(data.integrated, features = c("RCVRN", "PDE6H", "OPN1SW", "GNAT2", "GNGT2", "ARR3",
                                          "GUCA1A", "GUCA1B", "GUCA1C", "PDE5H", "PDE6C"), repel = TRUE, label = TRUE)
data.integrated <- RenameIdents(data.integrated, "10" = "Cones")
AverageExpression(data.integrated, features = "PDE6H")


FeaturePlot(data.integrated, features = c("ONECUT1", "ONECUT2", "LHX1", "CNTNAP2", "CACNG3", 
                                          "TFAP2B", "SOSTDC1", "TPM3", "PCSK1N", "TAGLN3"), repel = TRUE, label = TRUE)
data.integrated <- RenameIdents(data.integrated, "7" = "Horizontal")
AverageExpression(data.integrated, features = c("ONECUT1", "ONECUT2"))

ProcessSeuSub <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat) #could be replaced with SCTransform
  Seurat <- RunPCA(Seurat)
  Seurat <- FindNeighbors(Seurat, dims = 1:20, reduction = "pca", verbose = F)
  Seurat <- FindClusters(Seurat, resolution = 0.5, verbose = F)
  Seurat <- RunUMAP(Seurat, dims = 1:20)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:20)
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}  
horizontal.subset <- subset(data.integrated, subset = seurat_clusters == 7)
horizontal.subset <- ProcessSeuSub(horizontal.subset)
DimPlot(horizontal.subset, reduction = "umap", label = TRUE)
FeaturePlot(horizontal.subset, features = c("ONECUT1", "ONECUT2", "ONECUT3"))
FeaturePlot(horizontal.subset, features = c("LHX1", "ISL1"), label = TRUE)
FeaturePlot(horizontal.subset, features = c("LHX1", "PCDH9", "C1QL1", "SOSTDC1",
                                            "FAM135A", "SYNPR", "ISL1"), label = TRUE)
AverageExpression(horizontal.subset, features = c("LHX1", "ISL1"))
horizontal.subset <- RenameIdents(horizontal.subset, "0" = "H1", "1" = "H1", "2" = "H1", "3" = "H1", "4" = "H1")
horizontal.subset <- RenameIdents(horizontal.subset, "5" = "H2")

FeaturePlot(data.integrated, features = c("GYPC", "FBLN2", "SPARCL1", "TRH", "WIF1", "ALDH1L1", 
                                          "SOD2", "GLUL", "CLU", "NDRG2", "GFAP"), repel = TRUE, label = TRUE)
data.integrated <- RenameIdents(data.integrated, "15" = "Astrocytes")
AverageExpression(data.integrated, features = "GFAP")


FeaturePlot(data.integrated, features = c("FABP5", "NPVF", "CA2", "RLBP1", "CRABP1"), repel = TRUE, label = TRUE)
data.integrated <- RenameIdents(data.integrated, "12" = "Muller glia")
AverageExpression(data.integrated, features = c("NPVF", "RLBP1"))


#Markers which are highly expressed in both of Astrocytes & Muller glia) - prove of consistency
FeaturePlot(data.integrated, features = c("CRYAB", "GPX3", "PTGDS", "AQP4", "SLC1A3"), repel = TRUE, label = TRUE)
AverageExpression(data.integrated, features = c("CRYAB", "GPX3", "PTGDS", "AQP4", "SLC1A3"))

#Look at cluster 14
#Top upregulated genes in cluster 14: ISG15 (IFN-stimulated gene 15)
#MX1 (myxovirus (influenza virus) resistance 1, 
#PLSCR1 (phospholipid scramblase 1), IFITM3 (IFN Induced Transmembrane Protein 3), B2M (beta-2-microglobulin)
#These genes are involved in antiviral defence (against influenza virus)



FeaturePlot(data.integrated, features = c("IGSF21", "GSG1", "CA10", "CHN2", "NETO1",
                                          "SLC38A1", "GRM6", "PCP2", "TMEM215", 
                                          "TRPM1", "VSX1", "VSX2"), label = TRUE)
data.integrated <- RenameIdents(data.integrated, "3" = "Bipolar")                                          
AverageExpression(data.integrated, features = c("PCP2"))
AverageExpression(data.integrated, features = c("SLC38A1")) #co-expression in clusters 3 (bipolar), 
#RGC, Cones (but there isn't in Rods) & 2 (Amacrine?)
AverageExpression(data.integrated, features = "CA10") #co-expression in clusters 3 (bipolar) & RGC


FeaturePlot(data.integrated, features = c("LHX9", "RALYL", "MIR181A1HG", "RND3", "PAX6",
                                          "MEIS2", "GAD1", "GAD2", "C1QL2", "TFAP2A", "TFAP2B"), label = TRUE)
                                          
AverageExpression(data.integrated, features = "ZNF385D") #co-expression in clusters 2 (Amacrine?), 8 & 13 (?),
#Horizontal, Muller glia & RGC
AverageExpression(data.integrated, features = "PAX6") #co-expression in clusters 2 (Amacrine?), 8 & 13 (?),
#Horizontal, Muller glia & RGC
AverageExpression(data.integrated, features = "RALYL") #co-expression in clusters 2 (Amacrine?), 13 (?) & RGC
AverageExpression(data.integrated, features = "MIR181A1HG") #co-expression in clusters 2 (Amacrine?), 8, 13 (?) & RGC
AverageExpression(data.integrated, features = "MEIS2") #co-expression in clusters 2 (Amacrine?), Cones, 5, 0, 6 (Rods?), RGC


AverageExpression(data.integrated, features = c("GAD1", "GAD2", "C1QL2", "LHX9")) #specific markers for Amacrine cells
FeaturePlot(data.integrated, features = c("GAD1", "GAD2", "C1QL2", "LHX9", "RND3"), label = TRUE, repel = TRUE)
c("ZNF385D", "LMO4", "")
data.integrated <- RenameIdents(data.integrated, "2" = "Amacrine")  


FeaturePlot(data.integrated, features = c("NRL", "NR2E2", "GNAT1", "PDE6A", "IMPG1",
                                          "RHO", "CNGA1", "GNGT1", "SAG"), label = TRUE)
AverageExpression(data.integrated, features = c("PPEF2", "SAG"))
data.integrated <- RenameIdents(data.integrated, "5" = "Rods")


FeaturePlot(data.integrated, features = c("MKI67", "PCNA", "ASCL1", "SOX2", "HES1", "ID1"), label = TRUE)
data.integrated <- RenameIdents(data.integrated, "1" = "Progenitors", "4" = "Progenitors", "11" = "Progenitors")

#_______________________________________________________________________________

#Subsetting based on the clustering

levels(data.integrated)

horizontal.subset <- subset(data.integrated, idents = "Horizontal")
data.subset <- subset(data.integrated, idents = c("Progenitors", "Rods", "Cones", "Amacrine", "Horizontal", 
                                                  "Bipolar", "Muller glia", "Microglia", "Astrocytes", "RGC"))
data.integrated$labels <- data.integrated@active.ident
data.subset@active.ident <- data.subset$labels

# Run the standard workflow for visualization and clustering
data.subset <- ProcessSeuSub(data.subset)
DimPlot(data.subset, reduction = "umap", label.box = TRUE)
DefaultAssay(data.subset) <- "RNA"

FeaturePlot(data.subset, features = c("NEFL", "POU4F2"), label = TRUE)
data.subset <- RenameIdents(data.subset, "8" = "RGC")

FeaturePlot(data.subset, features = c("P2RY12", "TMEM119", "C1QA"), label = TRUE)
AverageExpression(data.subset, features = c("P2RY12", "TMEM119", "C1QA"))
data.subset <- RenameIdents(data.subset, "13" = "Microglia")

AverageExpression(data.subset, features = "PDE6H")
data.subset <- RenameIdents(data.subset, "9" = "Cones")

AverageExpression(data.subset, features = c("PPEF2", "SAG"))
data.subset <- RenameIdents(data.subset, "5" = "Rods")

FeaturePlot(data.subset, features = c("ONECUT1", "ONECUT2"), label = T)
data.subset <- RenameIdents(data.subset, "7" = "Horizontal")
FeaturePlot(data.subset, features = c("LHX1", "ISL1"))

AverageExpression(data.subset, features = "GFAP")
data.subset <- RenameIdents(data.subset, "12" = "Astrocytes")

AverageExpression(data.subset, features = c("NPVF", "RLBP1"))
data.subset <- RenameIdents(data.subset, "11" = "Muller glia")

AverageExpression(data.subset, features = "PCP2")
data.subset <- RenameIdents(data.subset, "2" = "Bipolar")                                

FeaturePlot(data.subset, features = c("GAD1", "GAD2", "C1QL2", "LHX9"), label = T)
AverageExpression(data.subset, features = c("GAD1", "GAD2", "C1QL2", "LHX9"))
data.subset <- RenameIdents(data.subset, "3" = "Amacrine", "6" = "Amacrine")

FeaturePlot(data.subset, features = c("MKI67", "PCNA", "ASCL1", "SOX2", "HES1", "ID1"), label = TRUE)
data.subset <- RenameIdents(data.subset, "0" = "Progenitors", "1" = "Progenitors", "4" = "Progenitors", "10" = "Progenitors")

levels(data.subset)

DimPlot(data.subset, reduction = "umap", label.box = TRUE)
