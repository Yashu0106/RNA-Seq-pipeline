
library("DESeq2")
library("ggplot2","oligo","affycoretools")
getwd() 
#change directory to location of htseq_count result files 
setwd("/Users/Dir_where_your_HTSeqResults.txt/")
directory<-"/Users/Users/Dir_where_your_HTSeqResults.txt/"

######################################
##########DEG analysis ###############
######################################
sampleFiles_leaf_el31<-grep("31",
                             list.files(directory),
                             value=T) #grep all files belongs to EL31 treatment group.

sampleNames_leaf_el31<- c("sample_1","sample_2","sample_3",
                        "sample_4","sample_5","sample_6") #set names, in the sample order as file names.

sampleCondition_leaf_el31<-c("control", "LR3","control",
                             "LR3", "LR3","control") #set conditions to each sample file, in the same order as files appeared in your DIR.

sampleTable_leaf_el31=data.frame(sampleName=sampleNames_leaf_el31,
                                  fileName=sampleFiles_leaf_el31,
                                  condition=sampleCondition_leaf_el31) #build a sample table
sampleTable_leaf_el31

#if got character error, convert character to factors
#it is the problem with default "stringsAsFactors = TRUE"
#Set this behaviour off permanently for your session with options(stringsAsFactors = FALSE)
options(stringsAsFactors = FALSE)
str(sampleTable_leaf_el31)

sampleTable_leaf_el31$condition<-as.factor(sampleTable_leaf_el31$condition)
sampleTable_leaf_el31$fileName<-as.factor(sampleTable_leaf_el31$fileName)
sampleTable_leaf_el31$sampleName<-as.factor(sampleTable_leaf_el31$sampleName)
sampleCondition_leaf_el31<-as.factor(sampleCondition_leaf_el31)
str(sampleTable_leaf_el31)
#####################################

ddsHTseq_leaf_el31<-DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_leaf_el31,
                                                directory = directory,
                                                design=~condition) #Build DESeqDataSet 
ddsHTseq_leaf_el31
colData(ddsHTseq_leaf_el31)$condition<-factor(colData(ddsHTseq_leaf_el31)$condition,
                                               levels=c('control','LR3'))
#set 'control' at first: the comparison will be the last level of this variable over the reference level 
dds_leaf_el31<-DESeq(ddsHTseq_leaf_el31) #DESeq2 analysis 
res_leaf_el31<-results(dds_leaf_el31) 
head(res_leaf_el31)
tail(res_leaf_el31)
res_leaf_el31<-res_leaf_el31[order(res_leaf_el31$padj),] #order results by padj value (the most significant to the least significant)
head(res_leaf_el31)
summary(res_leaf_el31) 


mcols(res_leaf_el31)$description
resdata_leaf_el31<-merge(as.data.frame(res_leaf_el31),
                          as.data.frame(counts(dds_leaf_el31,normalized=T)),
                          by='row.names',sort=F) #add normalized reads of each samples to table
names(resdata_leaf_el31)[1]<-'gene'
head(resdata_leaf_el31)


############cutoff value padj 0.05, order genes by fold change################
write.csv(resdata_leaf_el31,file="el31leaf_.csv",
          row.names = F) #save to directory
write.table(as.data.frame(counts(dds_leaf_el31),normalized=T), 
            file = 'el31leaf_normalized_counts.txt', sep = '\t') #save normalized data separately, can be used for downstream functional annotation analysis

sum(resdata_leaf_el31$padj<0.05, na.rm=TRUE) #summary to see how many genes passed padj < 0.05, which can be determined as significantly expressed (i.e., DEGs).
resdata_leaf_el31=subset(resdata_leaf_el31,padj<0.05) #subset DEG with padj<0.05
resdata_leaf_el31<-resdata_leaf_el31[order(resdata_leaf_el31$log2FoldChange),] #order the DEGs by log2FoldChange to see which one is up-/down-regulated
write.csv(as.data.frame(resdata_leaf_el31),file="leaf_el31_DEG analysis_padj0.05.csv",row.names = F)
##Separate up and down csv file manually
##Followed by functional annotation such as GO annotation.


###################################################
################Make plots and figures#############
###################################################
#MA plot
plotMA(dds_leaf_el31,ylim=c(-6,6),main="MA plot of EL-31 leaf LR3+ vs. control")
dev.off()

#rld and vsd for normalization raw counts 
rld_leaf_el31<-rlogTransformation(dds_leaf_el31, blind=T)
vsd_leaf_el31<-varianceStabilizingTransformation(dds_leaf_el31, blind=T)

#PCA#

plotPCA(rld_leaf_el31, intgroup = c("condition"))+
  geom_label(aes(label=sampleTable_leaf_el31$sampleName))
dev.off()





