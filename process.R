
################################
###     load R libraries     ###
################################

library("Seurat")
library("SingleR")
library("AnnotationHub")
library("SingleCellExperiment")
library("data.table")


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
gene.annotation.GFF <- readRDS(paste(wpath, refdata, "gene.annotation.GFF.rds",sep=""))

## SingleR Monaco dataset
monaco <- readRDS(paste(wpath, refdata, "singleR.MonacoImmuneData.rds",sep="/"))


################################
###     custom functions     ###
################################

get_stats_with_dens <- function(mat) {
	mat.genes_per_cell <- Matrix::colSums(mat > 0)
	mat.counts_per_cell <- Matrix::colSums(mat)
	mat.stats <- data.frame(cell=names(mat.genes_per_cell), geneCount=unname(mat.genes_per_cell),
							umiCount=unname(mat.counts_per_cell), stringsAsFactors=F)
	cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
	mat.dens <- densCols(mat.stats$geneCount, mat.stats$umiCount, colramp=colorRampPalette(c("black", "white")))
	mat.stats$dens <- col2rgb(mat.dens)[1,] + 1L
	mat.stats$col <- cols[mat.stats$dens]
	return(mat.stats)
}


################################
###   load CellRanger data   ###
################################

pbmc.raw <- Read10X(paste(datapath,"filtered_feature_bc_matrix",sep="/"))
# dim(pbmc.raw) # 36601  3337
pbmc.raw <- pbmc.raw[rowSums(pbmc.raw>0)>0, ]
# dim(pbmc.raw) # 22933  3337

pbmc.raw.stats <- get_stats_with_dens(pbmc.raw)
#apply(pbmc.raw.stats[,2:3], 2, summary)



