library(varbvs)
library(pracma)
library(glassoFast)
library(nor1mix)
library(roxygen2)
library(scEChIA)
library(reshape2)

#Processing spleen data
setwd('/Users/akankshaarora/Desktop/')
getwd()
cerebellum=read.table('cerebellum.txt',sep="\t",header=F)
head(cerebellum)
peaks=read.table('peaks_list.bed')
head(peaks)
peaks1 = cbind.data.frame(c(1:nrow(peaks)),paste(peaks[,1],peaks[,2],peaks[,3],sep="_"))
colnames(peaks1) = c("V1","V2")
finalmat=merge(peaks1, cerebellum, by='V1',all=F)
#final_mat1=finalmat[,c('V2.y', 'V1', 'V3')]
data=acast(finalmat,finalmat$V2.x~finalmat$V3)
View(data)

#write.table(data, file='cerebellum_data.tsv', quote=FALSE, sep='\t')
#Making count matrix, convert chr1_start_end to tsv by making them tab separated - chr1 start end
input_cer<-read.table('cerebellum_data.tsv',sep="\t",header=F)
head(input_cer)
input_cer_1<-subset(input_cer,input_cer$V1=="chr1")
View(input_cer_1)

#input$V5<-as.numeric(input$V5)

input_cer_1[,2:58] <- sapply(input_cer_1[,2:58],as.numeric)
sapply(input_cer_1, class)
#setwd(system.file("data", package="scEChIA")) #Keep your input data format like example data
genomic_region1 = read.table('GM12878_chr1.txt', sep="\t",header=F)
View(genomic_region2)
View(genomic_region1)
genomic_region1[is.na(genomic_region1)] = 0
genomic_region2 = read.table('IMR90_chr1.txt', sep="\t",header=F)
genomic_region2[is.na(genomic_region2)] = 0

gap = 25000 #Bin size 25kb (User can change it as per his choice)
patternf = 1 #Chromosome number (User can change it according to chromosome number)
chrNo = patternf
#data1 = read.table('GM1278_HiC.txt')
chrinfo = input_cer_1[, 1:3] #Chromosome location (chr, start, end)
chrNo = patternf
startCell = 4 #start sample column
endCell = 58 #end sample column
chromSize = 248956422 #Size of chromosome 1 (Change it according to chromosome number. Follow heading Chromosome Size of hg19)

rhomatrix = rhomatAvg(genomic_region1, genomic_region2, gap, patternf, input_cer_1, chrinfo)
predicted_interaction = Interaction_Prediction_1(chrinfo, input17, rhomatrix, chrNo, startCell, endCell, chromSize)

View(predicted_interaction)