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
    # out_file = args[2]
    out_file = paste0(normalizePath(dirname(args[2])), 
                      "/", basename(args[2]))
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

## Run DoubletFinder with varying classification stringencies 
if (is.null(doublet_rate)) {
    doublet_rate = ncol(seuObj) / 100000
}
doublet_rate = min(doublet_rate, 0.5)
nExp_poi <- round(doublet_rate * ncol(seuObj))
homotypic.prop <- modelHomotypic(seuObj$seurat_clusters)
print(paste("Doublet rate:", doublet_rate, "homo prop:", homotypic.prop))


## run 1: fixed parameters 
pN1 = 0.67
pK1 = round(0.5 / sqrt(ncol(seuObj)) * (1 - pN1), 4)
nExp_hete1 <- round(nExp_poi * (1 - homotypic.prop))
name_part1 <- paste(pN1, pK1, nExp_hete1, sep="_")
print(name_part1)

seuObj <- doubletFinder_v3(seuObj, PCs = 1:10, pN = pN1, pK = pK1, 
                           nExp = nExp_hete1, reuse.pANN = FALSE, sct = FALSE)

mean(seuObj[[paste0("DF.classifications_", name_part1)]] == 'Doublet')

dbFinder_table1 = cbind(seuObj[[paste0("pANN_", name_part1)]], 
                        seuObj[[paste0("DF.classifications_", name_part1)]])


## run 2: sweep the pK parameters
pN2 = 0.25 # 2 times doublets than input cells

sweep.res <- paramSweep_v3(seuObj, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn_val <- find.pK(sweep.stats)
bcmvn_val$pK <- as.numeric(as.character(bcmvn_val$pK))
# max nK=100; min nK=5
bcmvn_val = bcmvn_val[bcmvn_val$pK <= (150 / (ncol(seuObj) * (1 - pN2))), ]
bcmvn_val = bcmvn_val[bcmvn_val$pK >= (5 / (ncol(seuObj) * (1 - pN2))), ]
pK2 = bcmvn_val$pK[which.max(bcmvn_val$BCmetric)]


nExp_hete2 <- round(nExp_poi * (1 - homotypic.prop))
name_part2 <- paste(pN2, pK2, nExp_hete2, sep="_")
print(name_part2)

seuObj <- doubletFinder_v3(seuObj, PCs = 1:10, pN = pN2, pK = pK2, 
                           nExp = nExp_hete2, reuse.pANN = FALSE, sct = FALSE)

mean(seuObj[[paste0("DF.classifications_", name_part2)]] == 'Doublet')

dbFinder_table2 = cbind(seuObj[[paste0("pANN_", name_part2)]], 
                        seuObj[[paste0("DF.classifications_", name_part2)]])

## save data
out_table = cbind(rownames(dbFinder_table1), dbFinder_table1, dbFinder_table2)
colnames(out_table) <- c("cellID", "pANN_fix", "label_fix", "pANN_sweep", "label_sweep")

write.table(format(out_table, digits = 3), file = out_file, 
            quote = FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)


cat(paste0(
    paste("method", "pN", "pK", "nDoublet", sep="\t"), "\n",
    paste("fix", pN1, pK1, nExp_hete1, sep="\t"), "\n",
    paste("sweep", pN2, pK2, nExp_hete2, sep="\t"), "\n"),
    file=paste0(out_file, ".log"))
