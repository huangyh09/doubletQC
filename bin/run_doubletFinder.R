## Command line for running doubletFinder
## https://github.com/chris-mcginnis-ucsf/DoubletFinder


## input arguments:
## input_directory, out_file_path, expected_doublet_rate
args <- commandArgs(TRUE)
if (length(args) < 1) {
    stop("require at least 1 parameter: input directory")
}
dat_dir <- normalizePath(args[1])
print(dat_dir)

if (length(args) == 1) {
    out_file = paste0(dat_dir, "/doubletFinder_table.tsv")
} else {
    out_file = args[2]
}

if (length(args) <= 2) {
    doublet_rate = NULL
} else {
    doublet_rate = numerical(args[3])
}


## load packages
library(ggpubr)
library(Matrix)
library(Seurat)

# dat_dir <- "/hps/nobackup/stegle/users/huangh/msclerosis/"
setwd(dat_dir)

dat10x <- Read10X(data.dir = dat_dir)
seuObj <- CreateSeuratObject(counts = dat10x, project = "MySample", min.cells = 3)
seuObj

## Preprocess
seuObj <- NormalizeData(seuObj, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", 
                               nfeatures = 2000)
top10 <- head(VariableFeatures(seuObj), 10)
seuObj <- ScaleData(seuObj)

## dimention reduction and clustering
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj))
seuObj <- FindNeighbors(seuObj, dims = 1:10)
seuObj <- FindClusters(seuObj, resolution = 0.5)

# saveRDS(pbmc, paste0(dat_dir, "/seurat_Object.precessed.rds"))

## run douletFinder
library(DoubletFinder)

## pK Identification (no ground-truth)
sweep.res <- paramSweep_v3(seuObj, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn_val <- find.pK(sweep.stats)


## Run DoubletFinder with varying classification stringencies 
if (is.null(doublet_rate)) {
    doublet_rate = ncol(seuObj) / 100000
}
doublet_rate = min(doublet_rate, 0.5)
nExp_poi <- round(doublet_rate * ncol(seuObj))
print(paste("Expected doublet rate:", doublet_rate))

pN = 0.25
pK = 0.09
name_part <- paste(pN, pK, nExp_poi, sep="_")

seuObj <- doubletFinder_v3(seuObj, PCs = 1:10, pN = pN, pK = pK, 
                           nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# saveRDS(seuObj, paste0(dat_dir, "/seurat_Object.precessed.rds"))

mean(seuObj[[paste0("DF.classifications_", name_part)]] == 'Doublet')

dbFinder_table = cbind(seuObj[[paste0("pANN_", name_part)]], 
                       seuObj[[paste0("DF.classifications_", name_part)]])


write.table(format(dbFinder_table, digits = 3), file = out_file, 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names=FALSE)
