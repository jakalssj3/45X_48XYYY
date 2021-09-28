
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

get_y_stats <- function(genes, counts, annotations=gene.annotation.GFF) {
	y_stats <- subset(annotations, symbol %in% genes, select=c("symbol","description","gene_type"))
	y_stats$description <- gsub(" \\[.*", "", y_stats$description)
	if ("Seurat" %in% head(summary(counts),2)) {
		cwg=rowSums(counts@assays$RNA@counts[genes,]>0)
	} else {
		cwg=rowSums(counts[genes,]>0)
	}
	y_stats <- cbind(y_stats,
		"UMIsum"=round(rowSums(counts[y_stats$symbol,]),2),
		"UMImean"=round(rowMeans(counts[y_stats$symbol,]),2),
		"cellsWithGene"=cwg)
	y_stats$PctCellsWithGene <- paste(round(y_stats$cellsWithGene*100/ncol(counts),2),"%",sep="")
	y_stats <- y_stats[order(y_stats$UMIsum, decreasing=T),]
	return(y_stats)
}

get_DEGs_for_karyo <- function(object, labels, lfc=0.2, alpha=0.5, detest="wilcox") {
	degs <- lapply(unique(object@meta.data[[labels]]), function(x) {
		temp <- object[,object[[labels]]==x]
		Idents(temp) <- temp[[paste("karyo",labels,sep=".")]]
		f = setdiff(rownames(temp), c(y_genes, par_genes$Approved.symbol)) # we remove Y chromosome and PAR genes from testing
		res = Seurat::FindAllMarkers(temp, only.pos=T, verbose=T, logfc.threshold=lfc, min.pct=0.2, test.use=detest, features=f)
		if (nrow(res)!=0) {
			res = subset(res, p_val_adj < alpha)
			res = res %>% dplyr::group_by(cluster) %>% dplyr::arrange(cluster, desc(avg_logFC)) %>% as.data.frame()
		}
		return(res)
	})
	names(degs) <- unique(object@meta.data[[labels]])
	degs <- degs[sapply(degs, nrow)!=0] # remove cell types without DEGs
	return(degs)
}

