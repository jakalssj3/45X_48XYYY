
message("\n## Setting up working environment...\n")

################################
###     load R libraries     ###
################################

library("Seurat", quietly=T)
library("SingleR", quietly=T)
library("AnnotationHub", quietly=T)
library("SingleCellExperiment", quietly=T)
library("data.table", quietly=T)
library("Matrix", quietly=T)
library("dplyr", quietly=T)

################################
###  environmental variables ###
################################

seed <- 42
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
###      get references      ###
################################

## gene annotations from 10X Genomics
gene.annotation.GFF <- readRDS(paste(wpath, refdata, "gene.annotation.GFF.rds", sep="/"))
par_genes <- read.csv(paste(wpath, refdata, "PAR_genes.txt", sep="/"), sep="\t", stringsAsFactors=F)

## SingleR Monaco dataset
monaco <- readRDS(paste(wpath, refdata, "singleR.MonacoImmuneData.rds", sep="/"))

################################
###  load custom functions   ###
################################

source("helper_functions.R")

################################
###   load CellRanger data   ###
################################

message("\n## Processing raw data...\n")

if (file.exists(paste(wpath, "pbmc.raw.rds", sep="/"))) {
	pbmc.raw <- readRDS(paste(wpath, "pbmc.raw.rds", sep="/"))
} else {
	pbmc.raw <- Read10X(paste(datapath, sep="/"))
	# dim(pbmc.raw) # 36601  3337
	pbmc.raw <- pbmc.raw[rowSums(pbmc.raw>0)>0, ] # remove genes with 0 counts only
	# dim(pbmc.raw) # 22933  3337
	saveRDS(pbmc.raw, paste(wpath, "pbmc.raw.rds", sep="/"))
}

################################
###    output some stats     ###
################################

message("\n## Calculating some basic statistics...\n")

pbmc.raw.stats <- get_stats_with_dens(pbmc.raw)

mito.genes <- subset(gene.annotation.GFF, V1 == "chrM")$symbol
mito.genes <- intersect(rownames(pbmc.raw), mito.genes)
pbmc.raw.stats$percent.mito <- (Matrix::colSums(pbmc.raw[mito.genes, ])*100)/Matrix::colSums(pbmc.raw)
print(apply(pbmc.raw.stats[,c(2:3,6)], 2, summary))

################################
###     quality control      ###
################################

message("\n## Performing quality control...\n")

maxdims <- 20
mingenes <- 500; maxgenes <- 5000
minumi <- 1000; maxumi <- 20000
maxmito <- 10

filteroutcells <- data.frame("bad_cells"=rbind(
		sum(pbmc.raw.stats$geneCount < mingenes),
		sum(pbmc.raw.stats$geneCount > maxgenes),
		sum(pbmc.raw.stats$umiCount < minumi),
		sum(pbmc.raw.stats$umiCount > maxumi),
		sum(pbmc.raw.stats$percent.mito > maxmito)))
rownames(filteroutcells) <- c(paste0("<",mingenes," genes"),paste0(">",maxgenes," genes"),
		paste0("<",minumi," UMI"),paste0(">",maxumi," UMI"), paste0(">",maxmito, "% mito"))
print(filteroutcells)

keep_genes <- names(which(rowSums(pbmc.raw>0)>=10)) ## list of genes present in at least 10 cells
filteroutgenes <- data.frame("genes"=rbind(nrow(pbmc.raw), length(keep_genes)))
rownames(filteroutgenes) <- c("all","in >=10 cells")
print(filteroutgenes)

################################
###   create Seurat object   ###
################################

message("\n## Generating Seurat object...\n")

if (file.exists(paste(wpath, "pbmc.rds", sep="/"))) {
	pbmc <- readRDS(paste(wpath, "pbmc.rds", sep="/"))
} else {
	pbmc <- CreateSeuratObject(counts = pbmc.raw, project = "45,X/48,XYYY")
	pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
	pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")
	pbmc[["percent.y"]] <- PercentageFeatureSet(pbmc, features = intersect(rownames(pbmc), subset(gene.annotation.GFF, V1=="chrY")$symbol))
	pbmc[["percent.x"]] <- PercentageFeatureSet(pbmc, features = intersect(rownames(pbmc), subset(gene.annotation.GFF, V1=="chrX")$symbol))
	#dim(pbmc) # 22933  3337

	pbmc <- subset(pbmc, features = keep_genes)
	pbmc <- subset(pbmc, subset = percent.mt < maxmito) # filter by mitochondrial transcript conten
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

	pbmc$seurat_clusters <- pbmc$RNA_snn_res.0.6

	saveRDS(pbmc, paste(wpath, "pbmc.rds", sep="/"))
}

apply(pbmc@meta.data[,c(2:7)], 2, summary)

pbmc.cluster_stats <- get_cluster_stats(pbmc)
print(pbmc.cluster_stats)

################################
###    assign cell labels    ###
################################

message("\n## Assigning cell labels...\n")

if (!"cell_type.monaco" %in% colnames(pbmc)) {
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

	new_labels <- list("0"="Monocytes", "1"="CD4+ T cells", "2"="CD4+ T cells", "3"="CD8+ T cells",
			"4"="B cells", "5"="TFH/TREG", "6"="Natural killer", "7"="Vδ2 T cells", "8"="B cells",
			"9"="Monocytes", "10"="Megakaryocyte progenitors", "11"="NA")

	pbmc$new_labels <- unlist(new_labels[pbmc$seurat_clusters], use.names=F)
	pbmc$new_labels2 <- pbmc$new_labels
	pbmc$new_labels2[pbmc$seurat_clusters=="4"] <- "B cells (clust_4)"
	pbmc$new_labels2[pbmc$seurat_clusters=="8"] <- "B cells (clust_8)"
	pbmc$new_labels2[pbmc$seurat_clusters=="0"] <- "Monocytes (clust_0)"
	pbmc$new_labels2[pbmc$seurat_clusters=="9"] <- "Monocytes (clust_9)"
}

print(table(pbmc$new_labels, pbmc$seurat_clusters, dnn=c("cell.type:","cluster:")))

################################
### sex chromosomes content  ###
################################

message("\n## Calculating sex chromosome statistics and defining cell karyotypes...\n")

y_genes <- intersect(rownames(pbmc), subset(gene.annotation.GFF, V1=="chrY")$symbol)
y_genes.raw <- intersect(rownames(pbmc.raw), subset(gene.annotation.GFF, V1=="chrY")$symbol)

y_stats <- get_y_stats(y_genes, pbmc)
y_stats.raw <- get_y_stats(y_genes.raw, pbmc.raw)

pbmc$karyo <- ifelse(colSums(pbmc@assays$RNA@counts[y_genes,])==0, "45X", "48XYYY")
pbmc$karyo.new_labels <- paste(pbmc$karyo, pbmc$new_labels, sep=": ")
pbmc$karyo.new_labels2 <- paste(pbmc$karyo, pbmc$new_labels2, sep=": ")
print(table(pbmc$karyo, pbmc$new_labels, dnn=c("karyotype:","cell.type:")))

x_genes <- intersect(rownames(pbmc), subset(gene.annotation.GFF, V1=="chrX")$symbol)
auto_genes <- intersect(rownames(pbmc), subset(gene.annotation.GFF, !V1 %in% c("chrY","chrX"))$symbol)
pbmc$YtoX <- colSums(pbmc[y_genes,])/colSums(pbmc[x_genes,])
pbmc$YtoAuto <- colSums(pbmc[y_genes,])/colSums(pbmc[auto_genes,])
#sapply(unique(pbmc$karyo), function(x) summary(subset(pbmc, karyo==x)$YtoX))
#sapply(unique(pbmc$karyo), function(x) summary(subset(pbmc, karyo==x)$YtoAuto))

#mos.stats <- rbind("ALL"=table(colSums(pbmc[y_genes,])==0), table(pbmc$seurat_clusters, colSums(pbmc[y_genes,])==0))
#mos.stats <- rbind("ALL"=table(colSums(pbmc[y_genes,])==0), table(pbmc$new_labels, colSums(pbmc[y_genes,])==0))
mos.stats <- rbind("ALL"=table(colSums(subset(pbmc, new_labels!="NA")[y_genes,])==0),
	table(subset(pbmc, new_labels!="NA")$new_labels, colSums(subset(pbmc, new_labels!="NA")[y_genes,])==0))
mos.stats <- cbind("ALL_cells"=rowSums(mos.stats), "%_of_ALL_cells"=round(rowSums(mos.stats)*100/sum(mos.stats[1,]),2), mos.stats)
mos.stats <- cbind(mos.stats, round(mos.stats[,c(3:4)]*100/mos.stats[1,1],2)) # regarding all cells
mos.stats <- cbind(mos.stats, round(mos.stats[,3]*100/mos.stats[1,3],2)) # regarding the karyotype
mos.stats <- cbind(mos.stats, round(mos.stats[,4]*100/mos.stats[1,4],2))
mos.stats <- mos.stats[,c(1:3,5,7,4,6,8)]
colnames(mos.stats) <- c("ALL_cells","%_of_ALL_cells","48XYYY_cells","48XYYY-%_of_ALL","%_of_48XYYY","45X_cells","45X-%_of_ALL","%_of_45X")
print(mos.stats)

################################
### differential expression  ###
################################

message("\n## Running differential expression analyses...\n")

#pbmc.de <- get_DEGs_for_karyo(pbmc, "new_labels", lfc=0.2, alpha=0.05)
pbmc.de <- get_DEGs_for_karyo(pbmc, "new_labels2", lfc=0.2, alpha=0.05) ## final DEG set
print(pbmc.de)

################################
### save the final R object  ###
################################

saveRDS(pbmc, paste(wpath, "pbmc.rds", sep="/"))

message("\n## All done!\n")

