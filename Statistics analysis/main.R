#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#biocLite("GenomeInfoDb")
#biocLite("gplots")
#biocLite("calibrate")

graphics.off()
rm(list=ls())


source('my_functions.R')
source('main_run_DESeq2.R')
source('main_enrichment_analysis.R')


inputCts = "EC_table_counts.csv"
inputAnno = "sample_annotation.csv"

cutoff_p = 0.1	## threshold of adjusted p value
cutoff_f = 1	## threshold of fold change (log2)

cond1 = "Plin2-HF"
cond2 = "WT-HF"

outputECs = "DeSeq_results.csv"


inputKEGG = "KEGG_pathway.csv"

outputPath = "enrich_KEGG_pathway.csv"


R1 = main_run_DESeq2(inputCts,inputAnno,cond1,cond2,cutoff_p,cutoff_f,outputECs)


R2 = main_enrichment_analysis(inputKEGG, outputECs, cutoff_p, cutoff_f, outputPath)


