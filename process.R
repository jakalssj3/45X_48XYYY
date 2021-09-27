
################################
###     load R libraries     ###
################################

require("Seurat")
require("SingleR")
require("AnnotationHub")
require("SingleCellExperiment")
require("data.table")
require("Matrix")

################################
###  environmental variables ###
################################

seed <- 123
set.seed(seed)
options(width=220)
options(scipen=5)

################################
###       set up paths       ###
################################

wpath <- getwd()
setwd(wpath)
datapath <- "input"
refdata <- "refs"

################################
###    prepare references    ###
################################

## gene annotations from 10X Genomics
gene.annotation.GFF <- readRDS(paste(wpath, refdata, "gene.annotation.GFF.rds", sep="/"))

## SingleR Monaco dataset
monaco <- readRDS(paste(wpath, refdata, "singleR.MonacoImmuneData.rds", sep="/"))

################################
###     custom functions     ###
################################

source("helper_functions.R")

################################
###   load CellRanger data   ###
################################

pbmc.raw <- Read10X(paste(datapath, sep="/"))
# dim(pbmc.raw) # 36601  3337
pbmc.raw <- pbmc.raw[rowSums(pbmc.raw>0)>0, ]
# dim(pbmc.raw) # 22933  3337

################################
###    output some stats     ###
################################

pbmc.raw.stats <- get_stats_with_dens(pbmc.raw)

mito.genes <- grep("^MT-", gene.annotation.GFF$symbol, value=TRUE)
mito.genes <- intersect(rownames(pbmc.raw), mito.genes)
pbmc.raw.stats$percent.mito <- (Matrix::colSums(pbmc.raw[mito.genes, ])*100)/Matrix::colSums(pbmc.raw)
#apply(pbmc.raw.stats[,c(2:3,6)], 2, summary)


################################
###   create Seurat object   ###
################################

maxdims <- 20
mingenes <- 500; maxgenes <- 5000
minumi <- 1000; maxumi <- 20000

keep_genes <- names(which(rowSums(pbmc.raw>0)>=10)) ## list of genes present in at least 10 cells
#length(keep_genes) # 15779

if (file.exists(paste(wpath, "pbmc.rds", sep="/"))) {
	pbmc <- readRDS(paste(wpath, "pbmc.rds", sep="/"))
} else {
	pbmc <- CreateSeuratObject(counts = pbmc.raw, project="45,X/48,XYYY")
	pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
	pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")
	pbmc[["percent.y"]] <- PercentageFeatureSet(pbmc, features = intersect(rownames(pbmc), subset(gene.annotation.GFF, V1=="chrY")$symbol))
	pbmc[["percent.x"]] <- PercentageFeatureSet(pbmc, features = intersect(rownames(pbmc), subset(gene.annotation.GFF, V1=="chrX")$symbol))
	#dim(pbmc) # 22933  3337

	pbmc <- subset(pbmc, features = keep_genes)
	pbmc <- subset(pbmc, subset = percent.mt < 10) # filter by mitochondrial transcript conten
	pbmc <- subset(pbmc, subset = nFeature_RNA >= mingenes & nFeature_RNA <= maxgenes) # filter by gene count
	pbmc <- subset(pbmc, subset = nCount_RNA >= minumi & nCount_RNA <= maxumi) # filter by UMI
	#dim(pbmc) # 15779  2936
  
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes
	pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) # assign cell-cycle scores
	pbmc$CC_Difference <- pbmc$S.Score - pbmc$G2M.Score

	pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
	pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(pbmc)
	pbmc <- ScaleData(pbmc, features = all.genes, vars.to.regress = c("percent.mt", "nCount_RNA"))

	pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
	pbmc <- JackStraw(pbmc, num.replicate = 100, dims=50)
	pbmc <- ScoreJackStraw(pbmc, dims = 1:50)
	pbmc <- RunUMAP(pbmc, dims = 1:maxdims)
	pbmc <- FindNeighbors(pbmc, dims = 1:maxdims)
	pbmc <- FindClusters(pbmc, resolution=c(1.0, 0.8, 0.6, 0.4, 0.2))

  saveRDS(pbmc, paste(wpath, "pbmc.rds", sep="/"))
}

pbmc$seurat_clusters <- pbmc$RNA_snn_res.0.6
#apply(pbmc@meta.data[,c(2:7,12)], 2, summary)

pbmc.cluster_stats <- get_cluster_stats(pbmc)


################################
###    assign cell labels    ###
################################

Idents(pbmc) <- pbmc$seurat_clusters
sce.pbmc <- as.SingleCellExperiment(pbmc)

pred.monaco <- SingleR(test=sce.pbmc, ref=monaco, assay.type.test=1, labels=monaco$label.main)
pred.monaco2 <- SingleR(test=sce.pbmc, ref=monaco, assay.type.test=1, labels=monaco$label.fine)
celllabels <- cbind("cell_type.monaco"=pred.monaco$labels, "cell_type.monacoF"=pred.monaco2$labels)
rownames(celllabels) <- rownames(pred.monaco)

pbmc <- AddMetaData(pbmc, as.data.frame(celllabels))
#table(pbmc$cell_type.monaco, pbmc$seurat_clusters)
#table(pbmc$cell_type.monacoF, pbmc$seurat_clusters)
#table(pbmc$cell_type.monacoF, pbmc$cell_type.monaco)

new_labels <- list("0"="Monocytes", "1"="CD4+ T cells", "2"="CD4+ T cells", "3"="CD8+ T cells", "4"="B cells", "5"="TFH/TREG",
	 	   "6"="Natural killer", "7"="Vδ2 T cells", "8"="B cells", "9"="Monocytes", "10"="Megakaryocyte progenitors", "11"="NA")
pbmc$new_labels <- unlist(new_labels[pbmc$seurat_clusters], use.names=F)

