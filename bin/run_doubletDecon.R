## Command line for running douletDecon
## https://github.com/EDePasquale/DoubletDecon


## input arguments:
## input_directory, out_file_path, expected_doublet_rate
args <- commandArgs(TRUE)
if (length(args) < 1) {
    stop("require at least 1 parameter: input directory")
}
dat_dir <- normalizePath(args[1])

if (length(args) == 1) {
    out_file = paste0(dat_dir, "/doubletDecon_table.tsv")
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
seuObj <- FindClusters(seuObj, resolution = 1.0)
print(paste(length(unique(seuObj@active.ident)), "clusters"))

# saveRDS(seuObj, paste0(dat_dir, "/seurat_Object.precessed2.rds"))
# seuObj <- readRDS(paste0(dat_dir, "/seurat_Object.precessed.rds"))

## run douletDecon
library(DoubletDecon)

newFiles = Improved_Seurat_Pre_Process(seuObj, num_genes=50, 
                                       write_files=FALSE)

# saveRDS(newFiles, paste0(dat_dir, "/doubletDecon.newFiles.rds"))

filename = "MySample"
results=Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                           groupsFile=newFiles$newGroupsFile, 
                           filename=filename, 
                           location=paste0(dirname(out_file), "/"),
                           fullDataFile=NULL, 
                           removeCC=FALSE, 
                           species="hsa", 
                           rhop=1.1, 
                           write=FALSE, 
                           PMF=TRUE, 
                           useFull=FALSE, 
                           heatmap=FALSE,
                           centroids=TRUE,
                           num_doubs=100, 
                           only50=FALSE,
                           min_uniq=4)

saveRDS(results, paste0(dat_dir, "/doubletDecon.results.rds"))
