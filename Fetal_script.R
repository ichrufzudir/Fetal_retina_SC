library(DoubletFinder)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(Matrix)
library(MAST)

install.packages('remotes')
remotes::install_version("Seurat", "4.3.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))


RDoublet <- function(tmp){
  sweep.res.list <- paramSweep(tmp, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
  pKopt <- pKopt[1]
  homotypic.prop <- modelHomotypic(tmp$seurat_clusters)
  nExp_poi <- round(0.05*length(colnames(tmp)))  ## Assuming 10% doublet formation rate
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
  tmp <- doubletFinder(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
  return (tmp)
}
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 3000)
  Seurat <- ScaleData(Seurat) #could be replaced with SCTransform
  Seurat <- RunPCA(Seurat, npcs = 100)
  Seurat <- FindNeighbors(Seurat, dims = 1:100)
  Seurat <- FindClusters(Seurat, resolution = 1)
  Seurat <- RunUMAP(Seurat, dims = 1:100)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:100)
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}

#D80C
data_dir <- 'C://Documents and Settings/Nikita Bagaev/Documents/Datasets/Fetal retina/D80Cfetal_filtered_gene_bc_matrices/'
list.files(data_dir)
D80Cfetal <- Read10X(data.dir = data_dir)
D80CfetalS = CreateSeuratObject(counts = D80Cfetal)
D80CfetalS[["percent.rb"]] <- PercentageFeatureSet(D80CfetalS, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
D80CfetalS[["percent.mt"]] <- PercentageFeatureSet(D80CfetalS, pattern = "^MT-")
D80CfetalS <- CellCycleScoring(D80CfetalS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, nbin = 12 )
VlnPlot(D80CfetalS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
D80CfetalS <- subset(D80CfetalS, subset = nCount_RNA > 300 & nCount_RNA < 15000 & nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 15 & percent.rb < 40)
D80CfetalS <- ScaleData(D80CfetalS, verbose = T, vars.to.regress = c('percent.mt', "percent.rb","S.Score","G2M.Score"))
D80CfetalS <- ProcessSeu(D80CfetalS)
D80CfetalS <- RDoublet(D80CfetalS)
D80CfetalS <- subset(D80CfetalS, cells = colnames(D80CfetalS )[which(D80CfetalS [[]][12] == 'Singlet')])
D80CfetalS <- subset(D80CfetalS , cells = colnames(D80CfetalS )[which(D80CfetalS [[]][13] == 'Singlet')])
D80CfetalS <- ProcessSeu(D80CfetalS)
D80CfetalS$origin <- 'Thomas Reh Lab'
D80CfetalS$timepoint <- 'Day 80'
D80CfetalS$tissue <-  'Central retina'
D80CfetalS$donor <- '4'
SaveH5Seurat(D80CfetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/D80C.h5Seurat', overwrite = TRUE)
saveRDS(D80CfetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/D80C.rds')
gc()

#D82C
data_dir <- 'C://Documents and Settings/Nikita Bagaev/Documents/Datasets/Fetal retina/D82Cfetal_filtered_gene_bc_matrices/'
list.files(data_dir)
D80Cfetal <- Read10X(data.dir = data_dir)
D80CfetalS = CreateSeuratObject(counts = D80Cfetal)
D80CfetalS[["percent.rb"]] <- PercentageFeatureSet(D80CfetalS, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
D80CfetalS[["percent.mt"]] <- PercentageFeatureSet(D80CfetalS, pattern = "^MT-")
D80CfetalS <- CellCycleScoring(D80CfetalS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
VlnPlot(D80CfetalS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
D80CfetalS <- subset(D80CfetalS, subset = nCount_RNA > 300 & nCount_RNA < 15000 & nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 15 & percent.rb < 40)
D80CfetalS <- ScaleData(D80CfetalS, verbose = T, vars.to.regress = c('percent.mt', "percent.rb","S.Score","G2M.Score"))
D80CfetalS <- ProcessSeu(D80CfetalS)
D80CfetalS <- RDoublet(D80CfetalS)
D80CfetalS <- subset(D80CfetalS, cells = colnames(D80CfetalS )[which(D80CfetalS [[]][12] == 'Singlet')])
D80CfetalS <- subset(D80CfetalS , cells = colnames(D80CfetalS )[which(D80CfetalS [[]][13] == 'Singlet')])
D80CfetalS <- ProcessSeu(D80CfetalS)
D80CfetalS$origin <- 'Thomas Reh Lab'
D80CfetalS$timepoint <- 'Day 82'
D80CfetalS$tissue <-  'Central retina'
D80CfetalS$donor <- '5'
SaveH5Seurat(D80CfetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/D82C.h5Seurat', overwrite = TRUE)
saveRDS(D80CfetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/D82C.rds')
gc()

#D82P
data_dir <- 'C://Documents and Settings/Nikita Bagaev/Documents/Datasets/Fetal retina/D82Pfetal_filtered_gene_bc_matrices/'
list.files(data_dir)
D80Cfetal <- Read10X(data.dir = data_dir)
D80CfetalS = CreateSeuratObject(counts = D80Cfetal)
D80CfetalS[["percent.rb"]] <- PercentageFeatureSet(D80CfetalS, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
D80CfetalS[["percent.mt"]] <- PercentageFeatureSet(D80CfetalS, pattern = "^MT-")
D80CfetalS <- CellCycleScoring(D80CfetalS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
VlnPlot(D80CfetalS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
D80CfetalS <- subset(D80CfetalS, subset = nCount_RNA > 300 & nCount_RNA < 15000 & nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 15 & percent.rb < 40)
D80CfetalS <- ScaleData(D80CfetalS, verbose = T, vars.to.regress = c('percent.mt', "percent.rb","S.Score","G2M.Score"))
D80CfetalS <- ProcessSeu(D80CfetalS)
D80CfetalS <- RDoublet(D80CfetalS)
D80CfetalS <- subset(D80CfetalS, cells = colnames(D80CfetalS )[which(D80CfetalS [[]][12] == 'Singlet')])
D80CfetalS <- subset(D80CfetalS , cells = colnames(D80CfetalS )[which(D80CfetalS [[]][13] == 'Singlet')])
D80CfetalS <- ProcessSeu(D80CfetalS)
D80CfetalS$origin <- 'Thomas Reh Lab'
D80CfetalS$timepoint <- 'Day 82'
D80CfetalS$tissue <-  'Periferal retina'
D80CfetalS$donor <- '6'
SaveH5Seurat(D80CfetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/D82P.h5Seurat', overwrite = TRUE)
saveRDS(D80CfetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/D82P.rds')
gc()

#D105
data_dir <- 'C://Documents and Settings/Nikita Bagaev/Documents/Datasets/Fetal retina/D105fetal_filtered_gene_bc_matrices/'
list.files(data_dir)
D80Cfetal <- Read10X(data.dir = data_dir)
D80CfetalS = CreateSeuratObject(counts = D80Cfetal)
D80CfetalS[["percent.rb"]] <- PercentageFeatureSet(D80CfetalS, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
D80CfetalS[["percent.mt"]] <- PercentageFeatureSet(D80CfetalS, pattern = "^MT-")
D80CfetalS <- CellCycleScoring(D80CfetalS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
VlnPlot(D80CfetalS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
D80CfetalS <- subset(D80CfetalS, subset = nCount_RNA > 300 & nCount_RNA < 15000 & nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 15 & percent.rb < 40)
D80CfetalS <- ScaleData(D80CfetalS, verbose = T, vars.to.regress = c('percent.mt', "percent.rb","S.Score","G2M.Score"))
D80CfetalS <- ProcessSeu(D80CfetalS)
D80CfetalS <- RDoublet(D80CfetalS)
D80CfetalS <- subset(D80CfetalS, cells = colnames(D80CfetalS )[which(D80CfetalS [[]][12] == 'Singlet')])
D80CfetalS <- subset(D80CfetalS , cells = colnames(D80CfetalS )[which(D80CfetalS [[]][13] == 'Singlet')])
D80CfetalS <- ProcessSeu(D80CfetalS)
D80CfetalS$origin <- 'Thomas Reh Lab'
D80CfetalS$timepoint <- 'Day 105'
D80CfetalS$tissue <-  'Total retina'
D80CfetalS$donor <- '7'
SaveH5Seurat(D80CfetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/D105.h5Seurat', overwrite = TRUE)
saveRDS(D80CfetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/D105.rds')
gc()

#D125C
data_dir <- 'C://Documents and Settings/Nikita Bagaev/Documents/Datasets/Fetal retina/D125Cfetal_filtered_gene_bc_matrices/'
list.files(data_dir)
D80Cfetal <- Read10X(data.dir = data_dir)
D80CfetalS = CreateSeuratObject(counts = D80Cfetal)
D80CfetalS[["percent.rb"]] <- PercentageFeatureSet(D80CfetalS, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
D80CfetalS[["percent.mt"]] <- PercentageFeatureSet(D80CfetalS, pattern = "^MT-")
D80CfetalS <- CellCycleScoring(D80CfetalS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
VlnPlot(D80CfetalS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
D80CfetalS <- subset(D80CfetalS, subset = nCount_RNA > 300 & nCount_RNA < 15000 & nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 15 & percent.rb < 40)
D80CfetalS <- ScaleData(D80CfetalS, verbose = T, vars.to.regress = c('percent.mt', "percent.rb","S.Score","G2M.Score"))
D80CfetalS <- ProcessSeu(D80CfetalS)
D80CfetalS <- RDoublet(D80CfetalS)
D80CfetalS <- subset(D80CfetalS, cells = colnames(D80CfetalS )[which(D80CfetalS [[]][12] == 'Singlet')])
D80CfetalS <- subset(D80CfetalS , cells = colnames(D80CfetalS )[which(D80CfetalS [[]][13] == 'Singlet')])
D80CfetalS <- ProcessSeu(D80CfetalS)
D80CfetalS$origin <- 'Thomas Reh Lab'
D80CfetalS$timepoint <- 'Day 125'
D80CfetalS$tissue <-  'Central retina'
D80CfetalS$donor <- '8'
SaveH5Seurat(D80CfetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/D125C.h5Seurat', overwrite = TRUE)
saveRDS(D80CfetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/D125C.rds')
gc()

#D125P
data_dir <- 'C://Documents and Settings/Nikita Bagaev/Documents/Datasets/Fetal retina/D125Pfetal_filtered_gene_bc_matrices/'
list.files(data_dir)
D80Cfetal <- Read10X(data.dir = data_dir)
D80CfetalS = CreateSeuratObject(counts = D80Cfetal)
D80CfetalS[["percent.rb"]] <- PercentageFeatureSet(D80CfetalS, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
D80CfetalS[["percent.mt"]] <- PercentageFeatureSet(D80CfetalS, pattern = "^MT-")
D80CfetalS <- CellCycleScoring(D80CfetalS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, nbin = 12)
VlnPlot(D80CfetalS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
D80CfetalS <- subset(D80CfetalS, subset = nCount_RNA > 300 & nCount_RNA < 15000 & nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 15 & percent.rb < 40)
D80CfetalS <- ScaleData(D80CfetalS, verbose = T, vars.to.regress = c('percent.mt', "percent.rb","S.Score","G2M.Score"))
D80CfetalS <- ProcessSeu(D80CfetalS)
D80CfetalS <- RDoublet(D80CfetalS)
D80CfetalS <- subset(D80CfetalS, cells = colnames(D80CfetalS )[which(D80CfetalS [[]][12] == 'Singlet')])
D80CfetalS <- subset(D80CfetalS , cells = colnames(D80CfetalS )[which(D80CfetalS [[]][13] == 'Singlet')])
D80CfetalS <- ProcessSeu(D80CfetalS)
D80CfetalS$origin <- 'Thomas Reh Lab'
D80CfetalS$timepoint <- 'Day 125'
D80CfetalS$tissue <-  'Periferal retina'
D80CfetalS$donor <- '9'
SaveH5Seurat(D80CfetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/D125P.h5Seurat', overwrite = TRUE)
saveRDS(D80CfetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/D125P.rds')
gc()
