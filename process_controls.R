
message("\n## Setting up working environment...\n")

################################
###     load R libraries     ###
################################

library("Seurat", quietly=T)
library("data.table", quietly=T)
library("Matrix", quietly=T)
library("SingleR", quietly=T)

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

## SingleR Monaco dataset
monaco <- readRDS(paste(wpath, refdata, "singleR.MonacoImmuneData.rds", sep="/"))

################################
###  retrieve controls data  ###
################################

message("\n## Processing controls' data...\n")

### Controls 1
### A single-cell atlas of the peripheral immune response to severe COVID-19
### healthy controls: 4 males (H2: 49 yo, H4: 49 yo, H5: 48 yo, H6: 37 yo)
### raw data available under: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150728

if (file.exists(paste(wpath, refdata, "controls1.rds", sep="/"))) {
	controls1 <- readRDS(paste(wpath, refdata, "controls1.rds", sep="/"))
} else {
	url1 <- "https://hosted-matrices-prod.s3-us-west-2.amazonaws.com/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection-25/blish_covid.seu.rds"
	download.file(url1, paste(wpath, refdata, "controls1.rds", sep="/"))
	controls1 <- readRDS(paste(wpath, refdata, "controls1.rds", sep="/"))
}
#dim(controls1) # 26361 44721
controls1$orig.cell_type1 <- controls1$singler
controls1$orig.cell_type2 <- controls1$cell.type
controls1$orig.cell_type3 <- controls1$cell.type.fine
controls1$orig.cell_type4 <- controls1$cell.type.coarse
#table(controls1$Status, controls1$Donor)
#table(controls1$Sex, controls1$Donor)
#table(controls1$singler, controls1$cell.type.coarse)

### Controls 2
### Sampling time-dependent artifacts in single-cell genomics studies
### Selecting PBMC samples that were preserved immediately (0 hours)
### raw data available under https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132065

controls2 <- readRDS(paste(wpath, refdata, "10X_pbmc_Seurat_donors_list_clustered.RDS", sep="/"))
controls2$orig.cell_type <- controls2$`cell_type`
# dim(controls2$male) # 9803 7471

### Controls 3
### Single cell transcriptomics reveals opioid usage evokes widespread suppression of antiviral gene program
### healthy controls: 5 males (M1: 24 yo, M2: 26 yo, M3: 33 yo, M4: 41 yo, M5: 45 yo)
### raw data available under: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128879

## manually download two files below-mentioned from the following URL containing the data related to the study
## Human_PBMC_Opiate_Metadata.txt and Human_PBMC_LPS_Processed_Log_UMI.txt
## https://singlecell.broadinstitute.org/single_cell/study/SCP589/ifnb-treatment-single-cell-transcriptomics-reveals-opioid-usage-evokes-widespread-suppression-of-antiviral-gene-program?cluster=tSNE%20projection%20of%20IFNb-treated%20PBMCs&spatialGroups=--&annotation=Condition--group--study&subsample=all#study-download
## and place these files under the refdata folder

controls3.meta <- fread(paste(wpath, refdata, "Human_PBMC_Opiate_Metadata.txt", sep="/"), sep="\t", header=T)
controls3.meta <- subset(controls3.meta, Condition=="Control" & Stim=="Naive")
controls3.meta <- controls3.meta[-grep("LPS", controls3.meta$`Cell Types`),]
controls3.meta$`Cell Types` <- gsub("Naive ", "", controls3.meta$`Cell Types`)

controls3 <- fread(paste(wpath, refdata, "Human_PBMC_LPS_Processed_Log_UMI.txt", sep="/"), sep="\t", header=T)
controls3 <- data.frame(controls3)
colnames(controls3) <- gsub("\\.","\\-",colnames(controls3))
controls3 <- controls3[,c("GENE",controls3.meta$NAME)]
rownames(controls3) <- controls3$GENE
controls3$GENE <- NULL

controls3 <- CreateSeuratObject(counts = controls3, project="pbmc.controls3")
controls3@meta.data <- cbind(controls3@meta.data, controls3.meta[match(colnames(controls3), controls3.meta$NAME),])
controls3$orig.cell_type <- controls3$`Cell Types`
#dim(controls3) #19012 23402

################################
###    process public data   ###
################################

message("\n## Creating Seurat objects...\n")

maxdims <- 20
mingenes <- 500; maxgenes <- 5000
minumi <- 1000; maxumi <- 20000
maxmito <- 10

controls1 <- subset(controls1, Status=="Healthy" & Sex=="M")
#dim(controls1) # 26361 12221

controls2 <- controls2$male
controls2 <- subset(controls2, condition=="0h")
#dim(controls2) #9803 3637

controls3 <- subset(controls3, Sex=="Male" & Stim=="Naive")
# dim(controls3) #19012 16376

controls3$percent.mt <- PercentageFeatureSet(controls3, pattern = "^MT-")
controls3 <- FindVariableFeatures(controls3, selection.method = "vst", nfeatures = 2000)
controls3 <- ScaleData(controls3, vars.to.regress = c("nCount_RNA", "percent.mt"))
controls3 <- RunPCA(controls3, features = VariableFeatures(object = controls3))
controls3 <- RunUMAP(controls3, dims = 1:maxdims)
controls3 <- FindNeighbors(controls3, dims = 1:maxdims)
controls3 <- FindClusters(controls3, resolution=c(1.0, 0.8, 0.6, 0.4, 0.2))
controls3$seurat_clusters <- controls3$RNA_snn_res.0.6

###

pbmc.public <- list("controls1"=controls1, "controls2"=controls2, "controls3"=controls3)

pbmc.public <- lapply(pbmc.public, function(obj) {
	obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
	obj <- subset(obj, percent.mt < maxmito & nCount_RNA > minumi & nCount_RNA < maxumi & nFeature_RNA > mingenes & nFeature_RNA < maxgenes)
	kg <- names(which(Matrix::rowSums(obj@assays$RNA@counts>0)>=10))
	obj <- obj[kg, ]
	obj@meta.data <- droplevels(obj@meta.data)
	return(obj)
})
# lapply(pbmc.public, dim)

################################
###   annotate public data   ###
################################

message("\n## Assigning cell labels...\n")

pbmc.public <- lapply(pbmc.public, function(obj) {
	Idents(obj) <- obj$seurat_clusters
	sce.obj <- as.SingleCellExperiment(obj)
	pred.monaco <- SingleR(test=sce.obj, ref=monaco, assay.type.test=1, labels=monaco$label.main)
	pred.monaco2 <- SingleR(test=sce.obj, ref=monaco, assay.type.test=1, labels=monaco$label.fine)
	celllabels <- cbind("cell_type.monaco"=pred.monaco$labels, "cell_type.monacoF"=pred.monaco2$labels)
	rownames(celllabels) <- rownames(pred.monaco)
	obj <- AddMetaData(obj, as.data.frame(celllabels))
	return(obj)
})

################################
### check Y/X in public data ###
################################

message("\n## Calculating sex chromosome statistics...\n")

pbmc.public <- lapply(pbmc.public, function(obj) {
	yg <- intersect(rownames(obj), subset(gene.annotation.GFF, V1=="chrY")$symbol)
	xg <- intersect(rownames(obj), subset(gene.annotation.GFF, V1=="chrX")$symbol)
	ag <- intersect(rownames(obj), subset(gene.annotation.GFF, !V1 %in% c("chrY","chrX"))$symbol)
	obj$hasY <- ifelse(colSums(obj[yg,])!=0, "Y","noY")
	obj$percent.y <- PercentageFeatureSet(obj, features = yg)
	obj$hasX <- ifelse(colSums(obj[xg,])!=0, "X","noX")
	obj$percent.x <- PercentageFeatureSet(obj, features = xg)
	# operate on normalized counts
	obj$YtoX <- colSums(obj@assays$RNA@data[yg,])/colSums(obj@assays$RNA[xg,])
	obj$YtoAuto <- colSums(obj@assays$RNA[yg,])/colSums(obj@assays$RNA[ag,])
	return(obj)
})

#################################
###   merge pbmc & controls   ###
#################################

message("\n## Merging with main PBMC dataset...\n")

pbmc <- readRDS(paste(wpath, "pbmc.rds", sep="/"))

pbmc.big <- merge(subset(pbmc, karyo=="45X"), y = c(subset(pbmc, karyo=="48XYYY"), pbmc.public[[1]], pbmc.public[[2]], pbmc.public[[3]]),
	add.cell.ids = c("45,X", "48,XYYY","controls1", "controls2", "controls3"), project = "45X/48XXX_vs_controls")
newids <- paste(pbmc.big$karyo, pbmc.big$batch, pbmc.big$Donor, pbmc.big$sample, pbmc.big$Age)
newids <- gsub("^H","MH",gsub("JULIA_","M",gsub("Control \\d ","M",gsub("NA | NA","",newids))))
pbmc.big$newids <- factor(newids)
levels(pbmc.big$newids) <- c("45X","48XYYY","M24","M26","M33","M41","M45","M03","M04","MH2","MH4","MH5","MH6")
pbmc.big$karyo[is.na(pbmc.big$karyo)] <- "46XY"
pbmc.big$karyo <- factor(pbmc.big$karyo, levels=c("45X","46XY","48XYYY"))

################################
### save the final R object  ###
################################

saveRDS(pbmc.big, paste(wpath, "pbmc.big.rds", sep="/"))

message("\n## All done!\n")

