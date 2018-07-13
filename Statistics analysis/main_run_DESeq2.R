# set my own data files and parameters
#inputCts = "EC_table_counts.csv"
#inputAnno = "sample_annotation.csv"

#cutoff_p = 0.05	## threshold of adjusted p value
#cutoff_f = 2	## threshold of fold change (log2)

#cond1 = "Plin2-HF"
#cond2 = "WT-HF"

#outputECs = "DeSeq_results.csv"

main_run_DESeq2 <- function(inputCts, inputAnno, cond1, cond2, cutoff_p=0.05, cutoff_f=2, outputECs) {

	# load in the counts from a csv file
	countData <- as.matrix(read.csv(inputCts,sep=",",row.names=1)) 


	# load in the metadata from a csv file
	colData <- read.csv(inputAnno, row.names=1, sep=',') 

	# only pick the condition column of metadata
	colData <- colData[,c("condition"), drop=FALSE]

	# rearrange the columns of countData according to colData
	countData <- countData[, rownames(colData)] 


	# check if the condition names and the order are consistent
	all(rownames(colData) == colnames(countData))


	# set the input data for DESeq2
	library("DESeq2")
	dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition) 


	# perform a minimal pre-filtering to remove rows that have only 0 or 1 read. 
	dds <- dds[ rowSums(counts(dds)) > 1, ]


	#  pick the pair-wise conditions to compare with
	dds <- dds[, dds$condition %in% c(cond1,cond2) ]


	# remove those levels which do not have samples in the current DESeqDataSet:
	dds$condition <- droplevels(dds$condition)


	# run the DESeq pipeline 
	dds <- DESeq(dds)


	# generate results table with defined p value (defautly 0.1):
	res <- results(dds)
	#res <- results(dds, alpha=cutoff_p)


	# order our results table by the smallest adjusted p value:
	res <- res[order(res$padj),]


	## Merge with normalized count data
	resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
	names(resdata)[1] <- "enzyme"
	head(resdata)


	
	# How many adjusted p-values were less than a given cutoff?
	num_sigECs <- sum(res$padj < cutoff_p, na.rm=TRUE)


	# Plot dispersions
	plotDispEsts(dds, main="Dispersion plot")


	# Regularized log transformation for clustering/heatmaps, etc
	rld <- rlogTransformation(dds)
	head(assay(rld))
	hist(assay(rld))


	library(RColorBrewer)
	library(gplots)
	colData1 <- as.factor(dds$condition)
	(mycols <- brewer.pal(8, "Dark2")[1:length(unique(colData1))])

	# Sample distance heatmap
	sampleDists <- as.matrix(dist(t(assay(rld))))
	heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[colData1], RowSideColors=mycols[colData1],
          margin=c(10, 10), main="Sample Distance Matrix")


	## Examine plot of adjusted p-values
	hist(res$pvalue, breaks=50, col="grey")



	# select the gene with minimal adjusted pvalue to plot
	row_res <- rownames(res)
	plotCounts(dds, gene=row_res[which.min(res$padj)], intgroup="condition")


	## MA plot
	maplot(resdata, thresh=cutoff_p, main="MA Plot")


	## Volcano plot with "significant" enzymes labeled
	volcanoplot(resdata, lfcthresh=cutoff_f, sigthresh=cutoff_p, textcx=.8, xlim=c(-2.3, 2))



	# find the information about which variables and tests were used 
	summary <- mcols(res)$description


	# Exporting only the results which pass an adjusted p value threshold
	resSig <- subset(resdata, padj < cutoff_p & (log2FoldChange <= (-1*cutoff_f) | log2FoldChange >= cutoff_f))


	## Write results
	write.csv(resdata, file=outputECs, row.names=FALSE)


	return(list(dds = dds, res = res, summary = summary, resSig = resSig))
}