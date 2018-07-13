# set my own data files and parameters
#inputKEGG = "KEGG_pathway.csv"
#inputECs = "DeSeq_results.csv"

#cutoff_p = 0.05	## threshold of adjusted p value
#cutoff_f = 2	## threshold of fold change (log2)

#outputPath = "enrich_KEGG_pathway.csv"


main_enrichment_analysis <- function(inputKEGG, inputECs, cutoff_p, cutoff_f, outputPath) {

	# load in expressed enzymes from a csv file
	ECsData <- read.csv(inputECs, row.names=1, sep=',') 
	
	# select sig. enzymes
	resSig <- subset(ECsData, padj < cutoff_p & (log2FoldChange <= (-1*cutoff_f) | log2FoldChange >= cutoff_f))

	# get vectors of expressed and sig. enzymes
	exp_ECs <- rownames(ECsData)
	sig_ECs <- rownames(resSig)


	# load in KEGG pathway info. from a csv file
	pathData <- read.csv(inputKEGG, row.names=1, sep=',') 

	# number of KEGG pathways
	np <- nrow(pathData)


	# generate a vector of all enzymes in KEGG pathways, and
	# calculate the number of sig. enzymes appeared in each pathway

	enrich_path = pathData	
	enrich_path["num_exp_ECs"] <- 0
	enrich_path["num_sig_ECs"] <- 0

	ECs_all = character()

	for (i in 1:np){

		EC1 <- pathData[i,"ECs",drop=FALSE]
		EC2 <- unlist(strsplit(as.character(EC1$ECs),','))

		ECs_all = union(ECs_all,EC2)

		enrich_path[i,"num_exp_ECs"] <- length(EC2[EC2 %in% exp_ECs])
		enrich_path[i,"num_sig_ECs"] <- length(EC2[EC2 %in% sig_ECs])
	}


	# calculate p value for each pathway using enrichment hypergeomatric test

	M = length(ECs_all)
	K1 = length(exp_ECs)
	K = length(sig_ECs)

	pvalues_sig = vector(mode='double',length=np)
	pvalues_exp = vector(mode='double',length=np)

	for (i in 1:np){
		x <- enrich_path[i,"num_sig_ECs"]
		x1 <- enrich_path[i,"num_exp_ECs"]
		n <- enrich_path[i,"num_ECs"]
		pvalues_sig[i] = phyper(x-1,K,K1-K,x1,lower.tail=FALSE)
		pvalues_exp[i] = phyper(x1-1,K1,M-K1,n,lower.tail=FALSE)
	}

	# adjust p values based on Benjamini-Hochberg method
	padj_sig = p.adjust(pvalues_sig, method = 'BH', n = length(pvalues_sig))
	padj_exp = p.adjust(pvalues_exp, method = 'BH', n = length(pvalues_exp))


	## add more columns
	enrich_path["pvalues_exp"] <- pvalues_exp
	enrich_path["padj_exp"] <- padj_exp
	enrich_path["pvalues_sig"] <- pvalues_sig
	enrich_path["padj_sig"] <- padj_sig


	## Write results
	write.csv(enrich_path, file=outputPath, row.names=TRUE)

	return(enrich_path)
}