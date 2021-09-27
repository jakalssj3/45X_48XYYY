
get_stats_with_dens <- function(mat) {
	mat.genes_per_cell <- Matrix::colSums(mat > 0)
	mat.counts_per_cell <- Matrix::colSums(mat)
	mat.stats <- data.frame(cell=names(mat.genes_per_cell),
				geneCount=unname(mat.genes_per_cell),
				umiCount=unname(mat.counts_per_cell),
				stringsAsFactors=F)
	cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
	mat.dens <- densCols(mat.stats$geneCount, mat.stats$umiCount, colramp=colorRampPalette(c("black", "white")))
	mat.stats$dens <- col2rgb(mat.dens)[1,] + 1L
	mat.stats$col <- cols[mat.stats$dens]
	return(mat.stats)
}

get_cluster_stats <- function(object) {
	res <- as.data.frame(table(object$seurat_clusters))
	colnames(res) <- c("clusterID","cells")
	res$cells.perc <- round(res[,2]*100/ncol(object),2)
	res$genes.total <- lapply(res$clusterID, function(x) sum(rowSums(subset(object, seurat_clusters==x))!=0))
	res$genes.median <- lapply(res$clusterID, function(x) median(colSums(subset(object, seurat_clusters==x)@assays$RNA@counts!=0)))
	res$umi.median <- lapply(res$clusterID, function(x) median(colSums(subset(object, seurat_clusters==x)@assays$RNA@counts)))
	return(res)
}

