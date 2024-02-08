install.packages("dplyr")
install.packages("tidyverse")
install.packages("stringr")
install.packages("readxl")
install.packages("xlsx")
library(tidyverse)
library(stringr)
library("tidyr")
library("dplyr")
library(readxl)
library(xlsx)

setwd("/Users/dir_where_you_want_it_to_be/")
#Fetch your DEGs' GO annotations via Ensembl Biomart first, generate annotated file as GO.csv.
#Again, EL31L is name of my sample.

#Separate up-regulated genes and down-regulated genes into two files, easier to make analysis and graphs.
EL31L.up.mydata<-read.csv("EL31L.up.mydata.GO.long.csv",header = TRUE,)
EL31L.down.mydata<-read.csv("EL31L.down.mydata.GO.long.csv",header = TRUE,)

head(EL31L.up.mydata)
head(EL31L.down.mydata)
#########transform dataset#########
#collapse Gene.stable.ID column according to GO term column 
##Up-regulated DEGs##
EL31L.up.mydata.wide <- EL31L.up.mydata %>%
  group_by(GO.term.accession) %>%
  summarise(text=str_c(Gene.stable.ID,collapse=","))
colnames(EL31L.up.mydata.wide)<- c("GO term accession","Ensembl Gene ID mydata")
head(EL31L.up.mydata.wide)

write.csv(EL31L.up.mydata.wide,file="dir_where_you_want_it_to_be/EL31L.up.mydata.wide.csv",
          row.names = F) 

##Down-regulated DEGs##
EL31L.down.mydata.wide <- EL31L.down.mydata %>%
  group_by(GO.term.accession) %>%
  summarise(text=str_c(Gene.stable.ID,collapse=","))
colnames(EL31L.down.mydata.wide)<- c("GO term accession","Ensembl Gene ID")
head(EL31L.down.mydata.wide)
write.csv(EL31L.down.mydata.wide,file="dir_where_you_want_it_to_be/EL31L.down.mydata.wide.csv",
          row.names = F) #save to directory

#Via Ensembl Biomart, use GO term accession of 'mydata.wide.csv' file, 
#Fetch Ensembl ID via BioMart, to generate 'background.long.csv'  file. 
#This step is to generate a background reference list, i.e., SUPPOSED number of gene ID numbers matched to each GO ID)
#Now collapse using same strategy 
EL31L.up.background<-read.csv("EL31L.up.GO background.long.csv",header = TRUE)
EL31L.down.background<-read.csv("EL31L.down.GO background.long.csv",header = TRUE)
head(EL31L.up.background)
nrow(EL31L.up.background)
head(EL31L.down.background)
nrow(EL31L.down.background)

#now collapse gene ID based on GO term accession
EL31L.up.background.wide <- EL31L.up.background %>%
  group_by(GO.term.accession) %>%
  summarise(text=str_c(Gene.stable.ID,collapse=","))
colnames(EL31L.up.background.wide)<- c("GO term accession","Ensembl Gene ID background")

EL31L.down.background.wide <- EL31L.down.background %>%
  group_by(GO.term.accession) %>%
  summarise(text=str_c(Gene.stable.ID,collapse=","))
colnames(EL31L.down.background.wide)<- c("GO term accession","Ensembl Gene ID background")

head(EL31L.up.background.wide)
head(EL31L.down.background.wide)
nrow(EL31L.up.background.wide)
nrow(EL31L.down.background.wide)


write.csv(EL31L.up.background.wide,file="dir_where_you_want_it_to_be/EL31L.down.GO background.wide.csv",
          row.names = F) 
write.csv(EL31L.down.background.wide,file="dir_where_you_want_it_to_be/EL31L.down.GO background.wide.csv",
          row.names = F) 





#Merge mydata.wide.csv and background.wide.csv together by 'GO term accession'.

head(EL31L.up.mydata.wide)
head(EL31L.up.background.wide)
EL31L.up.enriched<-merge(EL31L.up.mydata.wide,EL31L.up.background.wide,by="GO term accession")
EL31L.down.enriched<-merge(EL31L.down.mydata.wide,EL31L.down.background.wide,by="GO term accession")
#Count the number of comma separated gene IDs of both mydata_input and background_input in a single cell, then add these numbers of counts into another column.
#This is your numbers used to conduct Fisher test. 
EL31L.up.enriched$EL31L.up.mydata_input<-count.fields(textConnection(EL31L.up.enriched$`Ensembl Gene ID mydata`),sep = ",")
EL31L.up.enriched$EL31L.up.background_input<-count.fields(textConnection(EL31L.up.enriched$`Ensembl Gene ID background`),sep = ",")
EL31L.down.enriched$EL31L.down.mydata_input<-count.fields(textConnection(EL31L.down.enriched$`Ensembl Gene ID`),sep = ",")
EL31L.down.enriched$EL31L.down.background_input<-count.fields(textConnection(EL31L.down.enriched$`Ensembl Gene ID background`),sep = ",")

head(EL31L.up.enriched[,c(1,4,5)],10)#print only 1, 4, 5 column, and first 5 rows
head(EL31L.down.enriched[,c(1,4,5)],10)
EL31L.up.enriched.backgroundComplete<-EL31L.up.enriched[c(1,4,5)]
EL31L.down.enriched.backgroundComplete<-EL31L.down.enriched[c(1,4,5)]
head(EL31L.up.enriched.backgroundComplete)
head(EL31L.down.enriched.backgroundComplete)
write.xlsx(EL31L.up.enriched.backgroundComplete,file="EL31L.up.enriched.kegg.backgrounComplete.xlsx")
write.xlsx(EL31L.down.enriched.backgroundComplete,file="EL31L.down.enriched.kegg.backgrounComplete.xlsx")


write.csv(EL31L.up.enriched,file="dir_where_you_want_it_to_be/EL31L.up.enriched.csv",
          row.names = F) 
write.csv(EL31L.down.enriched,file="dir_where_you_want_it_to_be/EL31L.down.enriched.csv",
          row.names = F) 

#Count the comma separated gene ID in a single cell, do it via excel. 
#Put in sample size (DEGs number), and population size (REF gene numbers) columns
#Save the file again
#Then conduct Fisher exact test: number of gene ID attached to this GO ID vs.number of gene ID supposed to attached to this GO ID in Vitis Vinifera reference gene 

#In excel, input DEGs number column
#Input background gene column 
#Input number of DEGs column (same number all way
#Input Reference list gene number column (same number 30661 all way)
#conduct fisher exact test in next column, using the four column build above
#Rank the p-values in excel
#Then do FDR adjust. 
#Extract the GOs with FDR <0.05 (FDR =pvalue*total rank numbers/rank nubmer of this sample)

#now save the data and extract 
install.packages("readxl")
install.packages("xlsx")
library(readxl)
library(xlsx)

EL31L.up.FDR<-read_excel("EL31L.up.enriched.xlsx")
EL31L.down.FDR<-read_excel("EL31L.down.enriched.xlsx")
head(EL31L.up.FDR)
head(EL31L.down.FDR)
nrow(EL31L.up.FDR)
nrow(EL31L.down.FDR)
EL31L.up.FDR=subset(EL31L.up.FDR,FDR<0.05) #subset GO ID with FDR<0.05 (i.e., the enriched GOs)
EL31L.up.FDR<-EL31L.up.FDR[order(EL31L.up.FDR$FDR),]#order results by FDR value (most significant to least)

EL31L.down.FDR=subset(EL31L.down.FDR,FDR<0.05) 
EL31L.down.FDR<-EL31L.down.FDR[order(EL31L.down.FDR$FDR),]
write.xlsx(EL31L.up.FDR,file="EL31L.up.enriched.FDR<0.05.xlsx")
write.xlsx(EL31L.down.FDR,file="EL31L.down.enriched.FDR<0.05.xlsx")

#Now fetch the annotation of these enriched GO ID via Biomart, then merge into enriched.final.xlsx
#Remember to remove duplicates in file retrieved from Biomart, don't know why they have duplicates, but remove it.
EL31L.up.FDR.anno<-read.csv("EL31L.up.enriched.final(FDR<0.05)_annotation.csv",
                            header=TRUE)
EL31L.down.FDR.anno<-read.csv("EL31L.down.enriched.final(FDR<0.05)_annotation.csv",
                              header=TRUE)


#make column name unique in both dataset 
head(EL31L.up.FDR.anno)
head(EL31L.up.FDR)
colnames(EL31L.up.FDR.anno)[1]<-"GO term accession"
EL31L.up.FDR.final<-merge(EL31L.up.FDR.anno,EL31L.up.FDR,by="GO term accession")
head(EL31L.up.FDR.final)

head(EL31L.down.FDR.anno)
head(EL31L.down.FDR)
colnames(EL31L.down.FDR.anno)[1]<-"GO term accession"
EL31L.down.FDR.final<-merge(EL31L.down.FDR.anno,EL31L.down.FDR,by="GO term accession")
head(EL31L.down.FDR.final)

#now you have your final version of enriched GO terms, with FDR values, and counted gene id.
write.xlsx(EL31L.up.FDR.final,file="EL31L.up.GO.final.xlsx")
write.xlsx(EL31L.down.FDR.final,file="EL31L.down.GO.final.xlsx")
















############################Doing the same thing for KEGG pathway enrichment########
#Fetch DEG'S KEGG annotation via Biomart, save as mydata.kegg.long.csv

EL31L.up.mydata.kegg<-read.csv("KEGG pathway enrichment/EL31L.up.mydata.kegg.long.csv",
                               header = TRUE,colClasses = c(KEGG.Pathway.and.Enzyme.ID = "character"))
#need to tell R that kegg pathway ID is read as 'characters', otherwise it will see ID as numeric and drop the leading zeros
EL31L.down.mydata.kegg<-read.csv("KEGG pathway enrichment/EL31L.down.mydata.kegg.long.csv",
                                 header = TRUE,colClasses = c(KEGG.Pathway.and.Enzyme.ID = "character"))

head(EL31L.up.mydata.kegg)
head(EL31L.down.mydata.kegg)
nrow(EL31L.up.mydata.kegg)
nrow(EL31L.down.mydata.kegg)
#########transform dataset#########
#collapse Gene.stable.ID column according to 'KEGG.Pathway.and.Enzyme.ID', making table frame wide
EL31L.up.mydata.kegg.wide <- EL31L.up.mydata.kegg %>%
  group_by(KEGG.Pathway.and.Enzyme.ID) %>%
  summarise(text=str_c(Gene.stable.ID,collapse=","))

colnames(EL31L.up.mydata.kegg.wide)<- c("KEGG pathway","Ensembl Gene ID mydata")
tail(EL31L.up.mydata.kegg.wide)
write.xlsx(EL31L.up.mydata.kegg.wide,file="KEGG pathway enrichment/EL31L.up.mydata.kegg.wide.xlsx") 
##DOWN##
EL31L.down.mydata.kegg.wide <- EL31L.down.mydata.kegg %>%
  group_by(KEGG.Pathway.and.Enzyme.ID) %>%
  summarise(text=str_c(Gene.stable.ID,collapse=","))

colnames(EL31L.down.mydata.kegg.wide)<- c("KEGG pathway","Ensembl Gene ID mydata")
tail(EL31L.down.mydata.kegg.wide)
write.xlsx(EL31L.down.mydata.kegg.wide,file="KEGG pathway enrichment/EL31L.down.mydata.kegg.wide.xlsx") #save to directory




#Via biomart, use KEGG ids of your wide form, fetch gene IDs from backgrou, generate background.kegg.long.csv
#(This is to generate background ref list, i.e. SUPPOSED number of gene ID numbers matched to each GO ID)
#now collapse using same strategy 
EL31L.up.background.kegg.long<-read.csv("KEGG pathway enrichment/EL31L.up.background.kegg.long.csv",
                                        header = TRUE,colClasses = c(KEGG.Pathway.and.Enzyme.ID = "character"))
#again tell R that your KEGG ID is characters, not numeric
EL31L.down.background.kegg.long<-read.csv("KEGG pathway enrichment/EL31L.down.background.kegg.long.csv",
                                          header = TRUE,colClasses = c(KEGG.Pathway.and.Enzyme.ID = "character"))
head(EL31L.up.background.kegg.long)
head(EL31L.down.background.kegg.long)
nrow(EL31L.up.background.kegg.long)
nrow(EL31L.down.background.kegg.long)

#now collapse gene ID based on 'KEGG id'KEGG.Pathway.and.Enzyme.ID'
EL31L.up.background.kegg.wide <- EL31L.up.background.kegg.long %>%
  group_by(KEGG.Pathway.and.Enzyme.ID) %>%
  summarise(text=str_c(Gene.stable.ID,collapse=","))
colnames(EL31L.up.background.kegg.wide)<- c("KEGG pathway","Ensembl Gene ID background")

EL31L.down.background.kegg.wide <- EL31L.down.background.kegg.long %>%
  group_by(KEGG.Pathway.and.Enzyme.ID) %>%
  summarise(text=str_c(Gene.stable.ID,collapse=","))
colnames(EL31L.down.background.kegg.wide)<- c("KEGG pathway","Ensembl Gene ID background")

head(EL31L.up.background.kegg.wide)
head(EL31L.down.background.kegg.wide)
nrow(EL31L.up.background.kegg.wide)
nrow(EL31L.down.background.kegg.wide)
write.xlsx(EL31L.up.background.kegg.wide,file="KEGG pathway enrichment/EL31L.up.background.kegg.wide.xlsx") 
write.xlsx(EL31L.down.background.kegg.wide,file="KEGG pathway enrichment/EL31L.down.background.kegg.wide.xlsx") 





#Merge mydata.kegg.wide with background.kegg.wide together by 'KEGG pathway' 
head(EL31L.up.mydata.kegg.wide)
head(EL31L.down.mydata.kegg.wide)
head(EL31L.up.background.kegg.wide)
head(EL31L.down.background.kegg.wide)

EL31L.up.enriched.kegg<-merge(EL31L.up.mydata.kegg.wide,EL31L.up.background.kegg.wide,by="KEGG pathway")
EL31L.down.enriched.kegg<-merge(EL31L.down.mydata.kegg.wide,EL31L.down.background.kegg.wide,by="KEGG pathway")

head(EL31L.up.enriched.kegg)
head(EL31L.down.enriched.kegg)
nrow(EL31L.up.enriched.kegg)
nrow(EL31L.down.enriched.kegg)
#Count the number of comma separated gene IDs of both mydata_input and backgroun_input in a single cell, then add these number of counts into another column
#this is your numbers used to calculate fisher test. 
EL31L.up.enriched.kegg$EL31L.up.mydata_input<-count.fields(textConnection(EL31L.up.enriched.kegg$`Ensembl Gene ID mydata`),sep = ",")
EL31L.up.enriched.kegg$EL31L.up.background_input<-count.fields(textConnection(EL31L.up.enriched.kegg$`Ensembl Gene ID background`),sep = ",")
EL31L.down.enriched.kegg$EL31L.down.mydata_input<-count.fields(textConnection(EL31L.down.enriched.kegg$`Ensembl Gene ID mydata`),sep = ",")
EL31L.down.enriched.kegg$EL31L.down.background_input<-count.fields(textConnection(EL31L.down.enriched.kegg$`Ensembl Gene ID background`),sep = ",")
tail(EL31L.up.enriched.kegg[,c(1,4,5)],10)#print only 1, 4, 5 column, and first 5 rows
tail(EL31L.down.enriched.kegg[,c(1,4,5)],10)

#the ensemble complete reference gene list column WILL be truncated, because it is too long! it won't really impact anything
write.xlsx(EL31L.up.enriched.kegg,file="KEGG pathway enrichment/EL31L.up.enriched.kegg.xlsx")
write.xlsx(EL31L.down.enriched.kegg,file="KEGG pathway enrichment/EL31L.down.enriched.kegg.xlsx")


#Count the comma separated gene ID in a single cell, via excel. 
#put in sample size (DEGs number), and population size (REF gene numbers) columns
#same the file again
#then do fisher exact test: number of gene ID attached to this GO ID vs. 
#number of gene ID supposed to attached to this GO ID in vitis vinifera reference gene 

#In excel, input DEGs number column
#input backgroun gene column 
#input number of DEGs column (same number all way
#input Reference list gene number clumn (same number 30661 all way)
#conduct fisher exact test in next column, using the four column build above
#!!!rank the p-values in excel
#then do FDR adjust. 
#extract the GOs with FDR <0.05 (FDR =pvalue*total rank numbers/rank nubmer of this sample)

#extract data and subset only FDR <0.05 ones (i.e., the enriched KEGG pathways)
EL31L.up.FDR.kegg<-read_excel("KEGG pathway enrichment/EL31L.up.enriched.kegg.xlsx")
EL31L.down.FDR.kegg<-read_excel("KEGG pathway enrichment/EL31L.down.enriched.kegg.xlsx")
head(EL31L.up.FDR.kegg)
head(EL31L.down.FDR.kegg)

EL31L.up.FDR.kegg=subset(EL31L.up.FDR.kegg,FDR<0.05) #subset GO ID with FDR<0.05 (i.e., the enriched GOs)
EL31L.up.FDR.kegg<-EL31L.up.FDR.kegg[order(EL31L.up.FDR.kegg$FDR),]#order results by FDR value (most significant to least)
EL31L.down.FDR.kegg=subset(EL31L.down.FDR.kegg,FDR<0.05) 
EL31L.down.FDR.kegg<-EL31L.down.FDR.kegg[order(EL31L.down.FDR.kegg$FDR),]
nrow(EL31L.up.FDR.kegg)
nrow(EL31L.down.FDR.kegg)

write.xlsx(EL31L.up.FDR.kegg,file="KEGG pathway enrichment/EL31L.up.FDR<0.05.kegg.xlsx")
write.xlsx(EL31L.down.FDR.kegg,file="KEGG pathway enrichment/EL31L.down.FDR<0.05.kegg.xlsx")



#Now you have a list of enriched KEGG pathways for your sample.
#Retrieve KEGG pathway name based on this KEGG pathway ID you have, all done.
kegg_idandname<-read.csv("/Users/alex/Desktop/PhD courses/Ph.D 2018-2019/RNA Seq 2021/2021_4_Ensembl GO annotation and own enrichment analysis/!KEGG ID and name list_cleaned!!!.csv",
                         header=TRUE,colClasses = c(Kegg.ID = "character"))
head(kegg_idandname)
colnames(kegg_idandname)<-c("Pathway_Rank","KEGG pathway",
                            "KEGG name rank_C","KEGG name rank_B","KEGG name rank_A") #'KEGG pathway' col name needs to be consistent
head(EL31L.up.FDR.kegg)
head(EL31L.down.FDR.kegg)

#Merge dataframes 
EL31L.up.final.kegg<-merge(kegg_idandname,EL31L.up.FDR.kegg,by="KEGG pathway")
head(EL31L.up.final.kegg)
EL31L.down.final.kegg<-merge(kegg_idandname,EL31L.down.FDR.kegg,by="KEGG pathway")
head(EL31L.down.final.kegg)


#now you have your final version of enriched GO terms, with FDR values, and counted gene id
write.xlsx(EL31L.up.final.kegg,file="KEGG pathway enrichment/EL31L.up.final.kegg.xlsx")
write.xlsx(EL31L.down.final.kegg,file="KEGG pathway enrichment/EL31L.down.final.kegg.xlsx")
