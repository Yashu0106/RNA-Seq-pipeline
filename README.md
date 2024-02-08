# RNA-Seq-pipeline
# RNA-Sequencing pipeline to determine differentially regulated genes.


# RNA-Seq general bioinformatics pipeline involved following (Song, et al. 2021) : 
# 1) RNA-Seq raw data quality control (data trimming). This step is to remove adaptors (arbitrarily unnecessary as some may agree)and low-quality reads. 

# 2) Transcriptome assembly / Alignment of sequence reads. As RNA-Seq generates millions of fragmented reads (e.g. Illumina sequencing platform generates reads length between 75 - 150bp), these fragmented reads need to be mapped to their proper locations before further counting. Aligner tool such as STAR is very computational expensive. It requires on average 30GB of memory. Therefore, Transcriptome assembly step is preferably performed on cloud computer. I used SHARCNET, a supercomputer system funded and shared by the universities of Ontario, Canada.

# 3) Read counts / gene feature counts and normalization. Once the fragmented reads are mapped properly, the number of reads will be quantified and such quantification is considered as 'raw counts'. Raw counts cannot be used directly to compared expressions between genes as they do not account for within- or between-sample effects such as read length, sequencing depth, and biological conditions. Normalization step such as RPKM and FPKM (poorly reviewed) or Med by DESeq and TMM by EdgeR (better reviewed) will be performed.

# 4) Identification of differentially expression genes (DEGs) and enrichment analysis. Identified read sequence and normalized read counts can be used to located associated gene annotations and be used to compare the expressional difference of genes.



# References: 
# Song, Y.; Hanner, R.H.; Meng, B. Probing into the effects of grapevine leafroll-associated viruses on the physiology, fruit quality and gene expression of grapes. Viruses 2021, 13, 593.
