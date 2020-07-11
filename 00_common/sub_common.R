library(Biostrings)
library(ggplot2)
library(gridExtra)
library(scales)
library(plyr)
library(reshape2)
library(proxy)
library(stats)
library(gplots)
library(RColorBrewer)
#library(dndscv)
library(caret)  
library(Rtsne)
library(factoextra)
library(stringi)
library(stringr)
#library(tidyverse)
library(VennDiagram) # use for up to 3 sets
#library(UpSetR) # use for more than 3 sets
library("VariantAnnotation")
library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)
library(LaF)
library("signature.tools.lib")
SampleMutNum <- function(total_muts, h, w, outputname){
  mut_number <- data.frame(table(total_muts$Sample))
  names(mut_number) <- c("Sample","Freq")
  # Write it into a file
  write.table(mut_number,paste0(outputname,".txt"),sep = "\t", col.names = T, row.names = F, quote = F)
  # Plot it to pdf file
  pdf(file=paste0(outputname,".pdf"), onefile=TRUE,height=h,width=w, useDingbats=FALSE)
  q <- ggplot(mut_number,aes(x=Sample,y=Freq,fill=Sample))+geom_bar(stat="identity", position=position_dodge())
  #q <- q+scale_fill_manual(values=c("orange","forestgreen","pink","royalblue","yellow"))
  q <- q+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=10),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  print(q)
  dev.off()
}
SampleMutNum_jitter <- function(total_muts,clone_col, subclone_col, h, w, outputname){
  mut_number <- data.frame(table(total_muts[,subclone_col], total_muts[,clone_col]))
  names(mut_number) <- c("Sample","Treatment","Freq")
  mut_number <- mut_number[mut_number$Freq>0,]
  # Write it into a file
  write.table(mut_number,paste0(outputname,".txt"),sep = "\t", col.names = T, row.names = F, quote = F)
  # Plot it to pdf file
  pdf(file=paste0(outputname,".pdf"), onefile=TRUE,height=h,width=w, useDingbats=FALSE)
  q <- ggplot(mut_number,aes(x=Treatment,y=Freq))+geom_jitter(position=position_jitter(0.1))+ 
    stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), 
                 geom="pointrange", color = "red")
  #q <- q+scale_fill_manual(values=c("orange","forestgreen","pink","royalblue","yellow"))
  q <- q+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  print(q)
  dev.off()
}
SampleMutNum_barplot <- function(total_muts,clone_col, subclone_col, h, w, outputname){
  mut_number <- data.frame(table(total_muts[,subclone_col], total_muts[,clone_col]))
  names(mut_number) <- c("Sample","Treatment","Freq")
  mut_number <- mut_number[mut_number$Freq>0,]
  # Write it into a file
  write.table(mut_number,paste0(outputname,".txt"),sep = "\t", col.names = T, row.names = F, quote = F)
  muts_summary <- ddply(mut_number,c("Treatment"),summarise,NChild=length(Freq),mean=mean(Freq),sd=sd(Freq),se=sd/sqrt(NChild))
  muts_summary[is.na(muts_summary)] <- 0
  muts_summary <- muts_summary[order(muts_summary$mean, decreasing = T),]
  write.table(mut_number,paste0(outputname,"_ddply.txt"),sep = "\t", col.names = T, row.names = F, quote = F)
  
  # Plot it to pdf file
  pdf(file=paste0(outputname,".pdf"), onefile=TRUE,height=h,width=w, useDingbats=FALSE)
  p <- ggplot(data=muts_summary, aes(x=Treatment, y=mean))+ geom_bar(stat="identity",position="dodge",fill="royal blue", width=.8)+xlab("KO gene")+ylab("Count")
#  p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",position=position_dodge(.9),width=.1)+scale_y_continuous()
  p <- p+scale_x_discrete(limits = muts_summary$Treatment,labels = muts_summary$Treatment)
  p <- p+geom_point(data=mut_number,aes(y=Freq,x=Treatment),colour="orange")
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=8,colour = "black",hjust = 1),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  print(p)
  dev.off()
}


muttype_freq_template <- read.table("/nfs/cancer_archive04/xz3/b_1176/00_common/muttype_freq_template.txt", sep = "\t", header = T, as.is = T)
#knockoutlist<- read.table("knockoutlist.txt", sep = "\t", header = F, as.is = T)
gen_muttype_selected <- function(CTsubs,selectedID){
  CTsubs <- CTsubs[CTsubs$VariantID %in% selectedID[,1],]
  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  
  CTsubs$muttype <- paste0(CTsubs$Ref, CTsubs$Alt, CTsubs$pre_context, CTsubs$rear_context)
  muttype_freq <- data.frame(table(CTsubs$muttype))
  names(muttype_freq) <- c("mut_type", "freq")
  muttype_freq$ref <-substr((muttype_freq$mut_type),1,1)
  muttype_freq$alt <-substr((muttype_freq$mut_type),2,2)
  muttype_freq$pre_context <-substr((muttype_freq$mut_type),3,3)
  muttype_freq$rear_context <-substr((muttype_freq$mut_type),4,4)
  muttype_freq$mutation <- paste0(muttype_freq$ref,">",muttype_freq$alt)
  muttype_freq$context <- paste0(muttype_freq$pre_context,muttype_freq$ref,muttype_freq$rear_context)
  muttype_freq$relfreq <- muttype_freq$freq/sum(muttype_freq$freq)
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(muttype_freq)
  
}

quick_pointplot <- function(inputdata,xcol,ycol,outputfilename){
  
  pdf(file=outputfilename, onefile=TRUE,height=5,width=6, useDingbats=FALSE)
  p <- ggplot(inputdata,aes(x=xcol,y=ycol))+geom_point(size=5)
  p <- p+theme(panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  p
  dev.off()
  
}
quick_barplot <- function(inputdata,xcol,ycol,h,w,outputfilename){
  
  pdf(file=outputfilename, onefile=TRUE,height=h,width=w, useDingbats=FALSE)
  p <- ggplot(inputdata,aes(x=inputdata[,xcol],y=inputdata[,ycol]))+geom_bar(stat="identity",position="dodge", width=.8,fill="blue")
  p <- p+scale_x_discrete(limits = as.character(inputdata[,xcol]))+xlab(xcol)+ylab(ycol)
  p <- p+theme(axis.text.x=element_text(size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  print(p)
  dev.off()
  
}
quick_jitterplot <- function(inputdata,xcol,ycol,xlabel_order, h, w, outputname){
 
  pdf(file=outputname, onefile=TRUE,height=h,width=w, useDingbats=FALSE)
  q <- ggplot(inputdata,aes(x=inputdata[,xcol],y=inputdata[,ycol]))+ geom_point(size = 2,shape=4)+#geom_jitter(position=position_jitter(0.1))+ 
    stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), 
                 geom="pointrange", color = "red")
  q <- q+scale_x_discrete(limits = xlabel_order)
  q <- q+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=10,colour = "black",hjust = 0.9),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  print(q)
  dev.off()
  
}



muttype_bysamples_freq <- function(allsubs,colnum,rownum,h,w,outputname){
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  
  muttype_freq_template <- read.table("../00_sampletracking/muttype_freq_template.txt", sep = "\t", header = T, as.is = T)
  samplelist <- data.frame(table(allsubs$Sample))
  names(samplelist) <- c("Sample","Freq")
  allsample_muttype_freq <- NULL
  for(i in 1:dim(samplelist)[1]){
    print(i)
    sample_subs <- allsubs[allsubs$Sample==as.character(samplelist[i,"Sample"]),]
    muttype_freq <- gen_muttype(sample_subs)
    muttype_freq <- merge(muttype_freq_template,muttype_freq[,c("mut_type","freq")],by="mut_type",all=T)
    muttype_freq[is.na(muttype_freq)] <- 0
    muttype_freq$Sample <- samplelist[i,"Sample"]
    muttype_freq$relfreq <- muttype_freq$freq/sum(muttype_freq$freq)
    allsample_muttype_freq <- rbind(allsample_muttype_freq,muttype_freq)
  }
  write.table(allsample_muttype_freq,paste0(outputname, "_muttype_freqbysample.txt"), sep="\t", row.names=F, col.names=T, quote=F)
  print(max(allsample_muttype_freq$freq))
  filename <- paste0(outputname, "_plot_bysample.pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=allsample_muttype_freq, aes(x=mut_type, y=freq,fill=mutation,width=0.5))+ geom_bar(position="dodge", stat="identity")+xlab("mutation type")+ylab("relative frequency")
  #+scale_x_discrete(breaks= c(muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1]),labels = as.character(muttype_freq[,"mutation"]))
 # p <- p+coord_cartesian(ylim=c(0, 300))
  p <- p+scale_x_discrete(labels = as.character(allsample_muttype_freq[,"context"]))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x= element_blank(),
              #axis.text.x=element_text(angle=90, vjust=0.5,size=5),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  p <- p+facet_wrap(~Sample,ncol=colnum,scales = "free")
  #expand_limits(x = muttype_freq[1,1], y = 0) 
  #p <- p+annotate("text", x=muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1], y=muttype_freq[,"mut_type"]+ max(muttype_freq[,"mut_type"])/50, label = as.character(muttype_freq[,"mutation"]), size=3)
  print(p)
  dev.off()
  
}
gen_muttype <- function(CTsubs){
  CTsubs_copy <- CTsubs
  
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  
  CTsubs$muttype <- paste0(CTsubs$Ref, CTsubs$Alt, CTsubs$pre_context, CTsubs$rear_context)
  muttype_freq <- data.frame(table(CTsubs$muttype))
  names(muttype_freq) <- c("mut_type", "freq")
  muttype_freq$ref <-substr((muttype_freq$mut_type),1,1)
  muttype_freq$alt <-substr((muttype_freq$mut_type),2,2)
  muttype_freq$pre_context <-substr((muttype_freq$mut_type),3,3)
  muttype_freq$rear_context <-substr((muttype_freq$mut_type),4,4)
  muttype_freq$mutation <- paste0(muttype_freq$ref,">",muttype_freq$alt)
  muttype_freq$context <- paste0(muttype_freq$pre_context,muttype_freq$ref,muttype_freq$rear_context)
  muttype_freq$relfreq <- muttype_freq$freq/sum(muttype_freq$freq)
  
  muttype_freq_template <- read.table("../00_sampletracking/muttype_freq_template.txt", sep = "\t", header = T, as.is = T)
  muttype_freq <- merge(muttype_freq_template,muttype_freq[,c("mut_type","freq")],by="mut_type",all=T)
  muttype_freq[is.na(muttype_freq)] <- 0
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(muttype_freq)
  
}
gen_6muttype <- function(CTsubs){
  CTsubs_copy <- CTsubs
  

  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  
  CTsubs$muttype <- paste0(CTsubs$Ref,">",CTsubs$Alt)
  muttype_freq <- data.frame(table(CTsubs$muttype))
  names(muttype_freq) <- c("mut_type", "freq")
  muttype_freq$ref <-substr((muttype_freq$mut_type),1,1)
  muttype_freq$alt <-substr((muttype_freq$mut_type),3,3)
  muttype_freq$mutation <- paste0(muttype_freq$ref,">",muttype_freq$alt)
  muttype_freq$relfreq <- muttype_freq$freq/sum(muttype_freq$freq)
  

  muttype_freq_template <- read.table("../00_common/muttype_freq_template.txt", sep = "\t", header = T, as.is = T)
  muttype_freq_template <- unique(muttype_freq_template[,c("ref","alt","mutation")])
  names(muttype_freq_template) <- c("ref","alt","mut_type")
  muttype_freq <- merge(muttype_freq_template,muttype_freq[,c("mut_type","freq")],by="mut_type",all=T)
  muttype_freq[is.na(muttype_freq)] <- 0
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(muttype_freq)
  
}

gen_muttype_new <- function(CTsubs){
  muttype_freq_template <- read.table("/nfs/cancer_archive04/xz3/b_1176/14_subs_0913/MutationType_template.txt", sep = "\t", header = T, as.is = T)

  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$MutationType <- paste0(CTsubs$pre_context,"[",CTsubs$Ref,">",CTsubs$Alt,"]",CTsubs$rear_context)
  sigfile_freq <- data.frame(table(CTsubs$Sample,CTsubs$MutationType))
  names(sigfile_freq) <- c("Sample","MutationType","Freq")
  control_sigset <-dcast(sigfile_freq,MutationType~Sample,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="MutationType",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(control_sigset)
  
}
gen_6muttype_new <- function(CTsubs){
  CTsubs_copy <- CTsubs
  
  muttype_freq_template <- read.table("/nfs/cancer_archive04/xz3/b_1176/14_subs_0913/MutationType_template.txt", sep = "\t", header = T, as.is = T)
  muttype_freq_template <- as.data.frame(unique(muttype_freq_template[,"Mutation"]))
  names(muttype_freq_template) <- "Mutation"
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$Mutation <- paste0(CTsubs$Ref,">",CTsubs$Alt)
  sigfile_freq <- data.frame(table(CTsubs$Sample,CTsubs$Mutation))
  names(sigfile_freq) <- c("Sample","Mutation","Freq")
  control_sigset <-dcast(sigfile_freq,Mutation~Sample,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="Mutation",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(control_sigset)
  
  
}

# choose the sample column
gen_muttype_sample <- function(CTsubs,sample){
  muttype_freq_template <- read.table("/nfs/cancer_archive04/xz3/b_1176/14_subs_0913/MutationType_template.txt", sep = "\t", header = T, as.is = T)
  
  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$MutationType <- paste0(CTsubs$pre_context,"[",CTsubs$Ref,">",CTsubs$Alt,"]",CTsubs$rear_context)
  sigfile_freq <- data.frame(table(CTsubs[,sample],CTsubs$MutationType))
  names(sigfile_freq) <- c("Sample","MutationType","Freq")
  control_sigset <-dcast(sigfile_freq,MutationType~Sample,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="MutationType",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(control_sigset)
  
}

# Generate 96 channel catalogue
GenCatalogue <- function(CTsubs, SampleCol){
  
  CTsubs[CTsubs$Chrom=="23","Chrom"]="X"
  CTsubs[CTsubs$Chrom=="24","Chrom"]="Y"
  
  # add 5' and 3' base information 
  CTsubs$pre_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos-1, CTsubs$Pos-1))
  CTsubs$rear_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos+1, CTsubs$Pos+1))
  
  mutation <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  base <- c("A","C","G","T")
  mutationtype <- NULL
  for(i in 1:length(mutation)){
    for(j in 1:length(base)){
      for(k in 1:length(base)){
        
        mutationtype <- c(mutationtype,paste0(base[j],"[",mutation[i],"]",base[k]))
      }
    }
    
  }
  
  muttype_freq_template <- data.frame("MutationType"=mutationtype)
  muttype_freq_template$Mutation <- substr(muttype_freq_template$MutationType,3,5)
    #read.table("./MutationType_template.txt", sep = "\t", header = T, as.is = T)
  
  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$MutationType <- paste0(CTsubs$pre_context,"[",CTsubs$Ref,">",CTsubs$Alt,"]",CTsubs$rear_context)
  sigfile_freq <- data.frame(table(CTsubs[,SampleCol],CTsubs$MutationType))
  names(sigfile_freq) <- c("Sample","MutationType","Freq")
  control_sigset <-dcast(sigfile_freq,MutationType~Sample,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="MutationType",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(control_sigset)
  
}

# Generate 32 trinucleotide context channel catalogue, pyrimidine 
Gen32Catalogue <- function(CTsubs, SampleCol){
  
  CTsubs[CTsubs$Chrom=="23","Chrom"]="X"
  CTsubs[CTsubs$Chrom=="24","Chrom"]="Y"
  
  # add 5' and 3' base information 
  CTsubs$pre_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos-1, CTsubs$Pos-1))
  CTsubs$rear_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos+1, CTsubs$Pos+1))
  
  pyrimidine <- c("C","T")
  base <- c("A","C","G","T")
  tri_pyrimidine <- NULL
  for(i in 1:length(pyrimidine)){
    for(j in 1:length(base)){
      for(k in 1:length(base)){
        
        tri_pyrimidine <- c(tri_pyrimidine,paste0(base[j],pyrimidine[i],base[k]))
      }
    }
    
  }
  
  muttype_freq_template <- data.frame("trinuc_pyrimidine"=tri_pyrimidine)

  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$trinuc_pyrimidine <- paste0(CTsubs$pre_context,CTsubs$Ref,CTsubs$rear_context)
  sigfile_freq <- data.frame(table(CTsubs[,SampleCol],CTsubs$trinuc_pyrimidine))
  names(sigfile_freq) <- c("Sample","trinuc_pyrimidine","Freq")
  control_sigset <-dcast(sigfile_freq,trinuc_pyrimidine~Sample,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="trinuc_pyrimidine",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(control_sigset)
  
}

# Generate 1536 channel catalogue
Gen1536Catalogue <- function(CTsubs, SampleCol){
  
  CTsubs[CTsubs$Chrom=="23","Chrom"]="X"
  CTsubs[CTsubs$Chrom=="24","Chrom"]="Y"
  
  # add 5' and 3' base information 
  CTsubs$pre_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos-2, CTsubs$Pos-1))
  CTsubs$rear_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos+1, CTsubs$Pos+2))
  
  mutation <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  base <- c("A","C","G","T")
  
  mutationtype <- NULL
  for(i in 1:length(mutation)){
    for(j in 1:length(base)){
      for(k in 1:length(base)){
        for(l in 1:length(base)){
          for(m in 1:length(base)){
            mutationtype <- c(mutationtype,paste0(base[k],base[j],"[",mutation[i],"]",base[l],base[m]))
          }
        }
        
      }
    }
    
  }
  
  muttype_freq_template <- data.frame("MutationType"=mutationtype)
  muttype_freq_template$Mutation <- substr(muttype_freq_template$MutationType,4,6)
  #read.table("./MutationType_template.txt", sep = "\t", header = T, as.is = T)
  
  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$MutationType <- paste0(CTsubs$pre_context,"[",CTsubs$Ref,">",CTsubs$Alt,"]",CTsubs$rear_context)
  sigfile_freq <- data.frame(table(CTsubs[,SampleCol],CTsubs$MutationType))
  names(sigfile_freq) <- c("Sample","MutationType","Freq")
  control_sigset <-dcast(sigfile_freq,MutationType~Sample,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="MutationType",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(control_sigset)
  
}
plot1536Countbasis <- function(muts_basis,colnum,h,w,outputname){
  
  muts_basis_melt <- melt(muts_basis,"MutationType")
  names(muts_basis_melt) <- c("MutationType","sample","count")
  muts_basis_melt$mutation <- substr(muts_basis_melt$MutationType,4,6)
  
  mutation_order <- muts_basis[,c("MutationType","MutationType")]
  names(mutation_order) <- c("MutationType","mutation")
  mutation_order$mutation <- substr(mutation_order$MutationType,4,6)
  mutation_order$pre1 <- substr(mutation_order$MutationType,2,2)
  mutation_order$rear1 <- substr(mutation_order$MutationType,8,8)
  mutation_order <- mutation_order[order(mutation_order$mutation,mutation_order$pre1,mutation_order$rear1),]
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  pdf(file=outputname, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=MutationType, y=count,fill=mutation))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Mutation Types")+ylab("Count")
  p <- p+scale_x_discrete(limits = as.character(mutation_order$MutationType))+ggtitle(outputname)
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(#axis.text.x=element_blank(),
    axis.text.x=element_text(size=3,angle=90,colour = "black"),
    axis.text.y=element_text(size=10,colour = "black"),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    plot.title = element_text(size=10),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill=NA))
  p <- p+facet_wrap(~sample,ncol=colnum,scales = "free")
  
  print(p)
  dev.off()
}

# Generate 5_96 channel catalogue
Gen5_96Catalogue <- function(CTsubs, SampleCol){
  
  CTsubs[CTsubs$Chrom=="23","Chrom"]="X"
  CTsubs[CTsubs$Chrom=="24","Chrom"]="Y"
  
  # add 5' and 3' base information 
  CTsubs$pre_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos-2, CTsubs$Pos-1))
  CTsubs$rear_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos+1, CTsubs$Pos+2))
  
  mutation <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  base <- c("A","C","G","T")
  
  mutationtype <- NULL
  for(i in 1:length(mutation)){
    for(j in 1:length(base)){
      for(k in 1:length(base)){
        
            mutationtype <- c(mutationtype,paste0(base[k],base[j],"[",mutation[i],"]"))
         
        
      }
    }
    
  }
  
  muttype_freq_template <- data.frame("MutationType"=mutationtype)
  muttype_freq_template$Mutation <- substr(muttype_freq_template$MutationType,4,6)
  #read.table("./MutationType_template.txt", sep = "\t", header = T, as.is = T)
  
  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$MutationType <- paste0(CTsubs$pre_context,"[",CTsubs$Ref,">",CTsubs$Alt,"]")
  sigfile_freq <- data.frame(table(CTsubs[,SampleCol],CTsubs$MutationType))
  names(sigfile_freq) <- c("Sample","MutationType","Freq")
  control_sigset <-dcast(sigfile_freq,MutationType~Sample,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="MutationType",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(control_sigset)
  
}
plot5_96Countbasis <- function(muts_basis,colnum,h,w,outputname){
  
  muts_basis_melt <- melt(muts_basis[,-2],"MutationType")
  names(muts_basis_melt) <- c("MutationType","sample","count")
  muts_basis_melt$mutation <- substr(muts_basis_melt$MutationType,4,6)
  mutation_order <- unique(muts_basis[,c("MutationType","MutationType")])
  names(mutation_order) <- c("MutationType","mutation")
  mutation_order$mutation <- substr(mutation_order$MutationType,4,6)
  mutation_order <- mutation_order[order(mutation_order$mutation),]
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  pdf(file=outputname, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=MutationType, y=count,fill=mutation))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Mutation Types")+ylab("Count")
  p <- p+scale_x_discrete(limits = as.character(mutation_order$MutationType))+ggtitle(outputname)
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(#axis.text.x=element_blank(),
    axis.text.x=element_text(size=5,angle=90,colour = "black"),
    axis.text.y=element_text(size=10,colour = "black"),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    plot.title = element_text(size=10),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill=NA))
  p <- p+facet_wrap(~sample,ncol=colnum,scales = "free")
  
  print(p)
  dev.off()
}

# Generate 3_96 channel catalogue
Gen3_96Catalogue <- function(CTsubs, SampleCol){
  
  CTsubs[CTsubs$Chrom=="23","Chrom"]="X"
  CTsubs[CTsubs$Chrom=="24","Chrom"]="Y"
  
  # add 5' and 3' base information 
  CTsubs$pre_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos-2, CTsubs$Pos-1))
  CTsubs$rear_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos+1, CTsubs$Pos+2))
  
  mutation <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  base <- c("A","C","G","T")
  
  mutationtype <- NULL
  for(i in 1:length(mutation)){
    for(l in 1:length(base)){
      for(m in 1:length(base)){
        
        mutationtype <- c(mutationtype,paste0("[",mutation[i],"]",base[l],base[m]))
        
        
      }
    }
    
  }
  
  muttype_freq_template <- data.frame("MutationType"=mutationtype)
  muttype_freq_template$Mutation <- substr(muttype_freq_template$MutationType,2,4)
  #read.table("./MutationType_template.txt", sep = "\t", header = T, as.is = T)
  
  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$MutationType <- paste0("[",CTsubs$Ref,">",CTsubs$Alt,"]",CTsubs$rear_context)
  sigfile_freq <- data.frame(table(CTsubs[,SampleCol],CTsubs$MutationType))
  names(sigfile_freq) <- c("Sample","MutationType","Freq")
  control_sigset <-dcast(sigfile_freq,MutationType~Sample,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="MutationType",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(control_sigset)
  
}
plot3_96Countbasis <- function(muts_basis,colnum,h,w,outputname){
  
  muts_basis_melt <- melt(muts_basis[,-2],"MutationType")
  names(muts_basis_melt) <- c("MutationType","sample","count")
  muts_basis_melt$mutation <- substr(muts_basis_melt$MutationType,2,4)
  mutation_order <- unique(muts_basis[,c("MutationType","MutationType")])
  names(mutation_order) <- c("MutationType","mutation")
  mutation_order$mutation <- substr(mutation_order$MutationType,2,4)
  mutation_order <- mutation_order[order(mutation_order$mutation),]
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  pdf(file=outputname, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=MutationType, y=count,fill=mutation))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Mutation Types")+ylab("Count")
  p <- p+scale_x_discrete(limits = as.character(mutation_order$MutationType))+ggtitle(outputname)
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(#axis.text.x=element_blank(),
    axis.text.x=element_text(size=5,angle=90,colour = "black"),
    axis.text.y=element_text(size=10,colour = "black"),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    plot.title = element_text(size=10),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill=NA))
  p <- p+facet_wrap(~sample,ncol=colnum,scales = "free")
  
  print(p)
  dev.off()
}

# Summarize number of -10 and +10 base around the mutation Unfinished
Count_SeqContext <- function(CTsubs,Contextlength,w,h,outputname){
  CTsubs[CTsubs$Chrom=="23","Chrom"]="X"
  CTsubs[CTsubs$Chrom=="24","Chrom"]="Y"
  
  # add 5' and 3' base information 
  CTsubs$pre_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos-Contextlength, CTsubs$Pos-1))
  CTsubs$rear_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos+1, CTsubs$Pos+Contextlength))
  
  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$Mutation <- paste0(CTsubs$Ref,">",CTsubs$Alt)
  mutation <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  for(i in 1:6){
    CTsubs_i <- CTsubs[CTsubs$Mutation==mutation[i],]
    sc5_all <- data.frame(matrix(unlist(strsplit(CTsubs_i$pre_context,split="")), ncol = Contextlength, byrow = T))
    sc3_all <- data.frame(matrix(unlist(strsplit(CTsubs_i$rear_context,split="")), ncol = Contextlength, byrow = T))
    
    sc5_summary <- NULL
    sc3_summary <- NULL
    
    for(j in 1:Contextlength){
      sc5_j <- table(sc5_all[,j])
      sc5_summary <- rbind(sc5_summary,sc5_j)
      
      sc3_j <- table(sc3_all[,j])
      sc3_summary <- rbind(sc3_summary,sc3_j)
    }
    sc5_summary <- data.frame(sc5_summary)
    sc5_summary$pos <- paste0("-",seq(from=Contextlength,to=1))
    
    sc3_summary <- data.frame(sc3_summary)
    sc3_summary$pos <- paste0("+",seq(from=1,to=Contextlength))
    
    sc_all <- rbind(sc5_summary,sc3_summary)
    if(substr(mutation[i],1,1)=="C"){
      sc_all <- rbind(sc_all,c(0,dim(CTsubs_i)[1],0,0,0))
      
    }else {
      sc_all <- rbind(sc_all,c(0,0,0,dim(CTsubs_i)[1]))
      
    }
    sc_all_melt <- melt(sc_all,c("pos"))
    names(sc_all_melt) <- c("pos","nucleotide","freq")
    pdf(file=paste0(outputname,"_",substr(mutation[i],1,1),"_",substr(mutation[i],3,3),"_",Contextlength,".pdf"), onefile=TRUE,width=w,height=h)
    p <- ggplot(data=sc_all_melt, aes(x=pos, y=freq,fill=nucleotide,width=1))+ geom_bar(stat="identity")+xlab("Position")+ylab("")
    #  p <- p+coord_cartesian(ylim=c(0, max(muttype_freq[,"freq"])))
    p <- p+scale_x_discrete(limits = c(seq(from=-Contextlength,to=-1),0,paste0("+",seq(from=1,to=Contextlength))))
#    p <- p+scale_fill_manual(values=mypalette)
    p <- p+theme(axis.text.x=element_text(size=10,colour = "black"),
                 axis.text.y=element_text(size=10,colour = "black"),
                 axis.title.x = element_text(size=15),
                 axis.title.y = element_text(size=15),
                 plot.title = element_text(size=10),
                 panel.grid.minor.x=element_blank(),
                 panel.grid.major.x=element_blank(),
                 panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank(),
                 panel.background = element_rect(fill = "white"),
                 panel.border = element_rect(colour = "black", fill=NA))
    
    #expand_limits(x = muttype_freq[1,1], y = 0) 
    #p <- p+annotate("text", x=muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1], y=muttype_freq[,"mut_type"]+ max(muttype_freq[,"mut_type"])/50, label = as.character(muttype_freq[,"mutation"]), size=3)
    print(p)
    dev.off()
    
  }
  
}

# Find -10 and +10 base around the mutation
Find_SeqContext <- function(CTsubs,Contextlength){
  CTsubs[CTsubs$Chrom=="23","Chrom"]="X"
  CTsubs[CTsubs$Chrom=="24","Chrom"]="Y"
  
  # add 5' and 3' base information 
  CTsubs$pre_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos-Contextlength, CTsubs$Pos-1))
  CTsubs$rear_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos+1, CTsubs$Pos+Contextlength))
  
  CTsubs_copy <- CTsubs
  
  CTsubs$Alt_py <- CTsubs$Alt
  CTsubs$pre_context_py <- CTsubs$pre_context
  CTsubs$rear_context_py <- CTsubs$rear_context
  CTsubs$Ref_py <- CTsubs$Ref
  
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt_py <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context_py <- as.character(reverseComplement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context_py <- as.character(reverseComplement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref_py <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$Mutation <- paste0(CTsubs$Ref_py,">",CTsubs$Alt_py)
  
  return(CTsubs)
}


# Assign MutationType for each mutation
AssignMutType <- function(CTsubs){
  
  CTsubs[CTsubs$Chrom=="23","Chrom"]="X"
  CTsubs[CTsubs$Chrom=="24","Chrom"]="Y"
  
  # add 5' and 3' base information 
  CTsubs$pre_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos-1, CTsubs$Pos-1))
  CTsubs$rear_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos+1, CTsubs$Pos+1))
  
  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  CTsubs$MutationType <- paste0(CTsubs$pre_context,"[",CTsubs$Ref,">",CTsubs$Alt,"]",CTsubs$rear_context)
  return(CTsubs)
  
}


plot_muttype <- function(muttype_freq,h,w,outputname){
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  muttype_freq_order <- muttype_freq[order(muttype_freq$mutation),]
  
  pdf(file=outputname, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muttype_freq_order, aes(x=MutationType, y=freq,fill=mutation,width=0.5))+ geom_bar(position="dodge", stat="identity")+xlab("mutation type")+ylab("")
 #  p <- p+coord_cartesian(ylim=c(0, max(muttype_freq[,"freq"])))
  p <- p+scale_x_discrete(limits = as.character(muttype_freq_order$MutationType))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  #expand_limits(x = muttype_freq[1,1], y = 0) 
  #p <- p+annotate("text", x=muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1], y=muttype_freq[,"mut_type"]+ max(muttype_freq[,"mut_type"])/50, label = as.character(muttype_freq[,"mutation"]), size=3)
  print(p)
  dev.off()
  
}

plot_mutationprofile <- function(CTsubs,SampleCol,colnum,h,w,outputname){
  
  a <- GenCatalogue(CTsubs,SampleCol)
  plotCountbasis(a[,-2],colnum,h,w,outputname)
}

count_subnum <- function(sublist,allsamplelist,outputname){
  subs_clean_sample_subnum <- data.frame(table(sublist$Sample_ko))
  names(subs_clean_sample_subnum) <- c("Sample_ko","Freq")
  subs_clean_sample_subnum <- merge(allsamplelist,subs_clean_sample_subnum,all.x=TRUE)
  subs_clean_sample_subnum[is.na(subs_clean_sample_subnum)] <- 0
  subs_clean_sample_subnum$Sample_ko=as.character(subs_clean_sample_subnum$Sample_ko)
  subs_clean_sample_subnum <- subs_clean_sample_subnum[order(subs_clean_sample_subnum$Sample_ko),]
  write.table(subs_clean_sample_subnum,paste0(gsub(".pdf","",outputname),".txt"),sep="\t",col.names = T, row.names = F, quote = F)
  
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
                  "#CC79A7","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  pdf(file=outputname, onefile=TRUE,height=4,width=15, useDingbats=FALSE)
  q <- ggplot(subs_clean_sample_subnum,aes(x=Sample_ko,y=Freq,fill=ko_gene))+geom_bar(stat="identity")+ scale_fill_manual(values=cbbPalette)
  q <- q+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=6),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA),
               legend.position="none")
  print(q)
  dev.off()
  
  
}

count_subnum_doublings <- function(sublist,allsamplelist,outputname){
  subs_clean_sample_subnum <- data.frame(table(sublist$Sample_ko))
  names(subs_clean_sample_subnum) <- c("Sample_ko","Freq")
  subs_clean_sample_subnum <- merge(allsamplelist,subs_clean_sample_subnum,all.x=TRUE)
  subs_clean_sample_subnum[is.na(subs_clean_sample_subnum)] <- 0
  subs_clean_sample_subnum$Sample_ko=as.character(subs_clean_sample_subnum$Sample_ko)
  subs_clean_sample_subnum <- subs_clean_sample_subnum[order(subs_clean_sample_subnum$Sample_ko),]
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
                  "#CC79A7","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  pdf(file=outputname, onefile=TRUE,height=4,width=15, useDingbats=FALSE)
  q <- ggplot(subs_clean_sample_subnum,aes(x=Sample_ko,y=Freq,fill=ko_gene))+geom_bar(stat="identity")+ scale_fill_manual(values=cbbPalette)
  q <- q+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA),
               legend.position="none")
  q <- q+facet_grid(.~doublings,scales="free",space="free")
  print(q)
  dev.off()
  
  
}

count_avgsubnum_KO <- function(sublist,kodoubling,outputname){
  subs_clean_sample_subnum <- data.frame(table(sublist$ko_gene))
  names(subs_clean_sample_subnum) <- c("ko_gene","Freq")
  subs_clean_sample_subnum <- merge(kodoubling,subs_clean_sample_subnum,all.x=TRUE)
  subs_clean_sample_subnum[is.na(subs_clean_sample_subnum)] <- 0
  subs_clean_sample_subnum$avgnum <- subs_clean_sample_subnum$Freq/subs_clean_sample_subnum[,2]
  subs_clean_sample_subnum <- subs_clean_sample_subnum[order(subs_clean_sample_subnum$avgnum, decreasing=TRUE),]
  write.table(subs_clean_sample_subnum,paste0(gsub(".pdf","",outputname),".txt"),sep="\t",col.names = T, row.names = F, quote = F)
  
  pdf(file=outputname, onefile=TRUE,height=3,width=9, useDingbats=FALSE)
  q <- ggplot(subs_clean_sample_subnum,aes(x=ko_gene,y=avgnum))+geom_bar(stat="identity",fill="blue")+ scale_fill_manual(values=cbbPalette)
  q <- q+scale_x_discrete(limits = subs_clean_sample_subnum$ko_gene)
  q <- q+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=10),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA),
               legend.position="none")
  print(q)
  dev.off()
  
  
}
count_ddplysubnum_KO <- function(sublist,allsamplelist,kodoubling,outputname){
  subs_clean_sample_subnum <- data.frame(table(sublist$Sample_ko))
  names(subs_clean_sample_subnum) <- c("Sample_ko","Freq")
  subs_clean_sample_subnum <- merge(allsamplelist,subs_clean_sample_subnum,all.x=TRUE)
  subs_clean_sample_subnum[is.na(subs_clean_sample_subnum)] <- 0
  subs_clean_sample_subnum$Sample_ko=as.character(subs_clean_sample_subnum$Sample_ko)
  
  subs_clean_sample_subnum_ddply <- ddply(subs_clean_sample_subnum,c("ko_gene"),summarise,NChild=length(Freq),total=sum(Freq),mean=mean(Freq),sd=sd(Freq),se=sd/sqrt(NChild))

  subs_clean_sample_subnum_ddply <- subs_clean_sample_subnum_ddply[order(subs_clean_sample_subnum_ddply$mean, decreasing=TRUE),]
  write.table(subs_clean_sample_subnum_ddply,paste0(gsub(".pdf","",outputname),".txt"),sep="\t",col.names = T, row.names = F, quote = F)
  
  hline_y <- subs_clean_sample_subnum_ddply[subs_clean_sample_subnum_ddply$ko_gene=="ATP284","mean"]
  
  #MAD <- mad(subs_clean_sample_subnum$Freq,center = hline_y) # equals to 1.4826*cMedian(abs(x - center)) MAD: Median Absolute Deviation
  #uplevel <- hline_y+2.5*MAD
  #lowlevel <- hline_y-2.5*MAD
  #print(paste0("ATP284_mean=", hline_y, "; MAD=", MAD, "; uplevel=", uplevel, "; lowlevel=",lowlevel))
  uplevel <- hline_y + (hline_y-min(subs_clean_sample_subnum_ddply$mean))
  lowlevel<-hline_y -(hline_y-min(subs_clean_sample_subnum_ddply$mean))
  print(paste0("ATP284_mean=", hline_y, "; Min=", min(subs_clean_sample_subnum_ddply$mean), "; uplevel=", uplevel, "; lowlevel=",lowlevel))
  pdf(file=outputname, onefile=TRUE,height=7,width=8, useDingbats=FALSE)
  q <- ggplot(subs_clean_sample_subnum_ddply,aes(x=ko_gene,y=mean))+geom_bar(stat="identity",position="dodge",fill="lightgreen")
  q <- q+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",position=position_dodge(.9),width=.2)
  q <- q+scale_x_discrete(limits = subs_clean_sample_subnum_ddply$ko_gene)
  q <- q+geom_hline(yintercept=hline_y, linetype="dashed", color = "blue")
  q <- q+geom_hline(yintercept=floor(lowlevel), linetype="dashed", color = "red")
  q <- q+geom_hline(yintercept=floor(uplevel), linetype="dashed", color = "red")
  q <- q+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=10),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA),
               legend.position="none")
  print(q)
  dev.off()
  
  
}


cal_culturesig <- function(muts,cutoff,outputname){
  knockoutlist <- data.frame(table(muts$knockout))
  muts$mut_ID <- paste0(muts$Chrom,"_",muts$Pos)
  culture_muts <- NULL
  for(i in 1:(dim(knockoutlist)[1]-1)){
    muts_knockout1 <- muts[muts$knockout==as.character(knockoutlist[i,1]),]
    for(j in (i+1):dim(knockoutlist)[1]){
      muts_knockout2 <- muts[muts$knockout==as.character(knockoutlist[j,1]),]
      
      overlap_muts <- muts_knockout2[muts_knockout2$mut_ID %in% muts_knockout1$mut_ID,]
      
      if(dim(overlap_muts)[1]>=cutoff){
        print(paste0(knockoutlist[i,1],"_",knockoutlist[j,1],":",dim(overlap_muts)[1]))
        
        culture_muts <- rbind(culture_muts,overlap_muts)
      }
    }
  }
  a=unique(culture_muts[,c("Chrom","Pos","Ref","Alt","pre_context","rear_context")])
  b=gen_muttype(a)
  plot_muttype(b,outputname)
}


average_muttype_bysamples_freq <- function(allsubs,h,w,outputname){
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  
  muttype_freq_template <- read.table("../00_common/muttype_freq_template.txt", sep = "\t", header = T, as.is = T)
  samplelist <- data.frame(table(allsubs$Sample))
  names(samplelist) <- c("Sample","Freq")
  allsample_muttype_freq <- NULL
  for(i in 1:dim(samplelist)[1]){
    print(i)
    sample_subs <- allsubs[allsubs$Sample==as.character(samplelist[i,"Sample"]),]
    muttype_freq <- gen_muttype(sample_subs)
    muttype_freq <- merge(muttype_freq_template,muttype_freq[,c("mut_type","freq")],by="mut_type",all=T)
    muttype_freq[is.na(muttype_freq)] <- 0
    muttype_freq$Sample <- samplelist[i,"Sample"]
    muttype_freq$relfreq <- muttype_freq$freq/sum(muttype_freq$freq)
    allsample_muttype_freq <- rbind(allsample_muttype_freq,muttype_freq)
  }
  allsample_muttype_freq <- allsample_muttype_freq[,c("mut_type","mutation","context","freq","Sample")]
  allsample_muttype_freq_cast <- dcast(allsample_muttype_freq,mut_type+mutation+context~Sample,value.var = "freq")
  allsample_muttype_freq_cast$total <- rowSums(allsample_muttype_freq_cast[,4:dim(allsample_muttype_freq_cast)[2]])
  write.table(allsample_muttype_freq_cast[,c(1:3,dim(allsample_muttype_freq_cast)[2])],sub(".pdf",".txt",outputname), sep="\t", row.names=F, col.names=T, quote=F)
  pdf(file=outputname, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=allsample_muttype_freq_cast, aes(x=mut_type, y=total,fill=mutation,width=0.5))+ geom_bar(position="dodge", stat="identity")+xlab("mutation type")+ylab("frequency")
  #+scale_x_discrete(breaks= c(muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1]),labels = as.character(muttype_freq[,"mutation"]))
  # p <- p+coord_cartesian(ylim=c(0, 300))
  p <- p+scale_x_discrete(labels = as.character(allsample_muttype_freq[,"context"]))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(#axis.text.x= element_blank(),
               axis.text.x=element_text(angle=90, vjust=0.5,size=5),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  #p <- p+facet_wrap(~Sample,ncol=colnum,scales = "free")
  #expand_limits(x = muttype_freq[1,1], y = 0) 
  #p <- p+annotate("text", x=muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1], y=muttype_freq[,"mut_type"]+ max(muttype_freq[,"mut_type"])/50, label = as.character(muttype_freq[,"mutation"]), size=3)
  print(p)
  dev.off()
  
}

average_6muttype_bysamples_freq <- function(allsubs,h,w,outputname){
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  
  muttype_freq_template <- read.table("../00_common/muttype_freq_template.txt", sep = "\t", header = T, as.is = T)
  muttype_freq_template <- unique(muttype_freq_template[,c("ref","alt","mutation")])
  names(muttype_freq_template) <- c("ref","alt","mut_type")
  
  samplelist <- data.frame(table(allsubs$Sample))
  names(samplelist) <- c("Sample","Freq")
  allsample_muttype_freq <- NULL
  for(i in 1:dim(samplelist)[1]){
    print(i)
    sample_subs <- allsubs[allsubs$Sample==as.character(samplelist[i,"Sample"]),]
    muttype_freq <- gen_6muttype(sample_subs)
    muttype_freq <- merge(muttype_freq_template,muttype_freq[,c("mut_type","freq")],by="mut_type",all=T)
    muttype_freq[is.na(muttype_freq)] <- 0
    muttype_freq$Sample <- samplelist[i,"Sample"]
    muttype_freq$relfreq <- muttype_freq$freq/sum(muttype_freq$freq)
    allsample_muttype_freq <- rbind(allsample_muttype_freq,muttype_freq)
  }
  allsample_muttype_freq <- allsample_muttype_freq[,c("mut_type","freq","Sample")]
  allsample_muttype_freq_cast <- dcast(allsample_muttype_freq,mut_type~Sample,value.var = "freq")
  allsample_muttype_freq_cast$total <- rowSums(allsample_muttype_freq_cast[,4:dim(allsample_muttype_freq_cast)[2]])
  allsample_muttype_freq_cast$percentage <- allsample_muttype_freq_cast$total/sum(allsample_muttype_freq_cast$total)
  
  write.table(allsample_muttype_freq_cast[,c(1:3,dim(allsample_muttype_freq_cast)[2])],sub(".pdf",".txt",outputname), sep="\t", row.names=F, col.names=T, quote=F)
  pdf(file=outputname, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=allsample_muttype_freq_cast, aes(x=mut_type, y=total,fill=mut_type,width=0.5))+ geom_bar(position="dodge", stat="identity")+xlab("mutation type")+ylab("frequency")
  #+scale_x_discrete(breaks= c(muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1]),labels = as.character(muttype_freq[,"mutation"]))
  # p <- p+coord_cartesian(ylim=c(0, 300))
  p <- p+scale_x_discrete(labels = as.character(allsample_muttype_freq[,"mut_type"]))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(#axis.text.x= element_blank(),
    axis.text.x=element_text(angle=90, vjust=0.5,size=10,colour = "black"),
    axis.text.y=element_text(size=10,colour = "black"),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    plot.title = element_text(size=10),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill=NA))
  #p <- p+facet_wrap(~Sample,ncol=colnum,scales = "free")
  #expand_limits(x = muttype_freq[1,1], y = 0) 
  #p <- p+annotate("text", x=muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1], y=muttype_freq[,"mut_type"]+ max(muttype_freq[,"mut_type"])/50, label = as.character(muttype_freq[,"mutation"]), size=3)
  print(p)
  dev.off()
  
}

plotbasis_average_se <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,2:dim(muts_basis)[2]])/(dim(muts_basis)[2]-1)
  
  muts_basis_melt <- melt(muts_basis,"MutationType")
  names(muts_basis_melt) <- c("MutationType","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("MutationType"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N),snr=abs(mean/sd))
  muts_basis_melt_summary_percentage <- muts_basis_melt_summary
  muts_basis_melt_summary_percentage$mean <- muts_basis_melt_summary$mean/mean_parentmuts
  muts_basis_melt_summary_percentage$sd <- muts_basis_melt_summary$sd/mean_parentmuts
  muts_basis_melt_summary_percentage$se <- muts_basis_melt_summary$se/mean_parentmuts
  muts_basis_melt_summary_percentage$snr <- muts_basis_melt_summary_percentage$mean/muts_basis_melt_summary_percentage$sd
  muts_basis_melt_summary_percentage$mutation <- substr(muts_basis_melt_summary_percentage$MutationType,3,5)
  muts_basis_melt_summary_percentage <- muts_basis_melt_summary_percentage[order(muts_basis_melt_summary_percentage$mutation),]
  
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary_percentage, aes(x=MutationType, y=mean,fill=mutation))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Mutation Types")+ylab("Percentage")
  p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="blue",position=position_dodge(.9),size=.2,width=0.5)+scale_y_continuous(limits=c(0, 0.2),breaks=seq(0, 0.2, 0.05),labels=percent)
  p <- p+scale_x_discrete(limits = as.character(muts_basis_melt_summary_percentage$MutationType))+ggtitle(outputname)
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  print(p)
  dev.off()
}
plotCountbasis_average_se <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,2:dim(muts_basis)[2]])/(dim(muts_basis)[2]-1)
  
  muts_basis_melt <- melt(muts_basis,"MutationType")
  names(muts_basis_melt) <- c("MutationType","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("MutationType"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N),snr=abs(mean/sd))
  muts_basis_melt_summary$mutation <- substr(muts_basis_melt_summary$MutationType,3,5)
  muts_basis_melt_summary <- muts_basis_melt_summary[order(muts_basis_melt_summary$mutation),]
  
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary, aes(x=MutationType, y=mean,fill=mutation))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Mutation Types")+ylab("Count")
  p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="blue",position=position_dodge(.9),size=.2,width=0.5) #+scale_y_continuous(limits=c(0, 350),breaks=seq(0, 350, 100))
  p <- p+scale_x_discrete(limits = as.character(muts_basis_melt_summary$MutationType))+ggtitle(outputname)
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=4,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  print(p)
  dev.off()
#  return(muts_basis_melt_summary)
  write.table(muts_basis_melt_summary,paste0(outputname, ".txt"), sep = "\t", col.names = T, row.names = F, quote = F)
}
plotCountbasis <- function(muts_basis,colnum,h,w,outputname){
  
  muts_basis_melt <- melt(muts_basis,"MutationType")
  names(muts_basis_melt) <- c("MutationType","sample","count")
  muts_basis_melt$mutation <- substr(muts_basis_melt$MutationType,3,5)
  
  mutation_order <- muts_basis[,c("MutationType","MutationType")]
  names(mutation_order) <- c("MutationType","mutation")
  mutation_order$mutation <- substr(mutation_order$MutationType,3,5)
  mutation_order <- mutation_order[order(mutation_order$mutation),]
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  pdf(file=outputname, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=MutationType, y=count,fill=mutation))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Mutation Types")+ylab("Count")
  p <- p+scale_x_discrete(limits = as.character(mutation_order$MutationType))+ggtitle(outputname)
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(#axis.text.x=element_blank(),
               axis.text.x=element_text(size=5,angle=90,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  p <- p+facet_wrap(~sample,ncol=colnum,scales = "free")
  
  print(p)
  dev.off()
}
plotPercentagebasis <- function(muts_basis,colnum,h,w,outputname){
  
  muts_basis[,-which(colnames(muts_basis)=="MutationType")] <- muts_basis[,-which(colnames(muts_basis)=="MutationType")]/colSums(muts_basis[,-which(colnames(muts_basis)=="MutationType")])[col(muts_basis[,-which(colnames(muts_basis)=="MutationType")])]
  muts_basis_melt <- melt(muts_basis,"MutationType")
  names(muts_basis_melt) <- c("MutationType","sample","count")
  muts_basis_melt$mutation <- substr(muts_basis_melt$MutationType,3,5)
  
  mutation_order <- muts_basis[,c("MutationType","MutationType")]
  names(mutation_order) <- c("MutationType","mutation")
  mutation_order$mutation <- substr(mutation_order$MutationType,3,5)
  mutation_order <- mutation_order[order(mutation_order$mutation),]
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  pdf(file=outputname, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=MutationType, y=count,fill=mutation))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Mutation Types")+ylab("Percentage")
#  p <- p+scale_x_discrete(limits = as.character(mutation_order$MutationType))+ggtitle(outputname)+scale_y_continuous(limits=c(0, 0.2),breaks=seq(0, 0.2, 0.05),labels=percent)
  p <- p+scale_x_discrete(limits = as.character(mutation_order$MutationType))+ggtitle(outputname)+scale_y_continuous(limits=c(0, 0.5),breaks=seq(0, 0.5, 0.1),labels=percent)
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5,colour = "black"),
               axis.text.y=element_text(size=15,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  p <- p+facet_wrap(~sample,ncol=colnum,scales = "free")
  
  print(p)
  dev.off()
}

sd_highD <- function(target_matrix, direction){
  
  centroid <- rowMeans(target_matrix)
  target_matrix_diff <- target_matrix-centroid
  
  #### Project target_matrix_diff to given direction
  # project_vector <- apply(target_matrix_diff,2,function(x) project(x, direction, type='length'))
  cos_vector_square <- apply(target_matrix_diff,2,function(x) (norm(as.matrix(x),"f")*cos_similarity(x,direction))^2)
  sd_matrix <- sqrt(sum(cos_vector_square)/(dim(target_matrix_diff)[2]-1))
  #sd_matrix <- 
  return(sd_matrix)
  
}
bootstrapGenomesfun <- function(genomes){
  
  return(apply(genomes, 2, function(x) rmultinom(1, sum(x), x)))
}
cos_similarity <- function(v1,v2){
  v1v2 <- sum(v1*v2)
  v1_length <- sqrt(sum(v1*v1))
  v2_length <- sqrt(sum(v2*v2))
  return(v1v2/v1_length/v2_length)
}

ExtractSig <- function(control_profile_set, compound_profile_set, sampling_number,boundary,h,w,outputname) {
  
  centroid_control <- rowMeans(control_profile_set)
  sd_control <- apply(control_profile_set,1,sd)
  boundary_control <- centroid_control+boundary*sd_control
  
  centroid_compound <- rowMeans(compound_profile_set)
  
  #RepControl <- matrix(rep(centroid_control,sampling_number),ncol = sampling_number)
  #bootstrapControl <- bootstrapGenomesfun(RepControl)
  diff_mean_all <- NULL
  for(bt_num in 1:sampling_number){
    RepCompound <- matrix(rep(centroid_compound,dim(compound_profile_set)[2]),ncol = dim(compound_profile_set)[2])
    bootstrapCompound <- bootstrapGenomesfun(RepCompound)
    
    diff_all_boundary <- bootstrapCompound-boundary_control
    diff_all <- bootstrapCompound-centroid_control
    
    diff_all[which(diff_all_boundary<0)] <- 0
    
    diff_mean_all <- cbind(diff_mean_all,rowMeans(diff_all))
  }
  diff_mean_all <- as.data.frame(diff_mean_all)
  names(diff_mean_all) <- paste0("s",c(1:dim(diff_mean_all)[2]))
  rownames(diff_mean_all) <- rownames(control_profile_set)
  diff_mean_all$MutationType <- rownames(control_profile_set)
  
  
  #diff_all <- bootstrapCompound-bootstrapControl
  #diff_all_positive <- as.data.frame(diff_all[,which(apply(diff_all,2,min)>=0)])
  # diff_all_percentage <- diff_all/centroid_control
  #  diff_all_positive <- as.data.frame(diff_all[,which(apply(diff_all_percentage,2,min)>=0)])
  
  muts_basis_melt <- melt(diff_mean_all,"MutationType")
  names(muts_basis_melt) <- c("MutationType","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("MutationType"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N))
  muttype_freq_template <- read.table("/nfs/cancer_archive04/xz3/b_1176/14_subs_0913/MutationType_template.txt", sep = "\t", header = T, as.is = T)
  muts_basis_melt_summary <- merge(muttype_freq_template, muts_basis_melt_summary,by="MutationType",all.x=T)
  muts_basis_melt_summary[is.na(muts_basis_melt_summary)] <- 0
  muts_basis_melt_summary <- muts_basis_melt_summary[order(muts_basis_melt_summary$Mutation),]
  
  write.table(muts_basis_melt_summary,paste0(outputname, "_exposure.txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  
  diff_all_percentage <- diff_mean_all[,1:(dim(diff_mean_all)[2]-1)]/colSums(diff_mean_all[,1:(dim(diff_mean_all)[2]-1)])[col(diff_mean_all[,1:(dim(diff_mean_all)[2]-1)])]
  rownames(diff_all_percentage) <- rownames(control_profile_set)
  diff_all_percentage$MutationType <- rownames(control_profile_set)
  
  muts_basis_melt_percentage <- melt(diff_all_percentage,"MutationType")
  names(muts_basis_melt_percentage) <- c("MutationType","sample","count")
  muts_basis_melt_summary_percentage <- ddply(muts_basis_melt_percentage,c("MutationType"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N))
  muttype_freq_template <- read.table("/nfs/cancer_archive04/xz3/b_1176/14_subs_0913/MutationType_template.txt", sep = "\t", header = T, as.is = T)
  muts_basis_melt_summary_percentage <- merge(muttype_freq_template, muts_basis_melt_summary_percentage,by="MutationType",all.x=T)
  muts_basis_melt_summary_percentage[is.na(muts_basis_melt_summary_percentage)] <- 0
  muts_basis_melt_summary_percentage <- muts_basis_melt_summary_percentage[order(muts_basis_melt_summary_percentage$Mutation),]
  write.table(muts_basis_melt_summary_percentage,paste0(outputname, "_percentage.txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  
  filename <- paste0(outputname, "_Signature.pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary, aes(x=MutationType, y=mean,fill=Mutation))+ geom_bar(stat="identity",position="dodge")+xlab("Substitution Types")+ylab("Count")
  p <- p+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),color="blue",position=position_dodge(.9),size=.2,width=0.5)#+scale_y_continuous(limits=c(-1,1),breaks=(seq(-1,1,0.2)),labels = scales::percent)
  p <- p+scale_x_discrete(limits = as.character(muts_basis_melt_summary$MutationType))+ggtitle(paste0(outputname,"_exposure"))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  
  o <- ggplot(data=muts_basis_melt_summary_percentage, aes(x=MutationType, y=mean,fill=Mutation))+ geom_bar(stat="identity",position="dodge")+xlab("Substitution Types")+ylab("Percentage")
  o <- o+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),color="blue",position=position_dodge(.9),size=.2,width=0.5)+scale_y_continuous(labels = scales::percent) #limits=c(0,0.3),breaks=(seq(0,1,0.3)),
  o <- o+scale_x_discrete(limits = as.character(muts_basis_melt_summary$MutationType))+ggtitle(paste0(outputname,"_signature"))
  o <- o+scale_fill_manual(values=mypalette)
  o <- o+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  grid.arrange(p,o,ncol=1)
  dev.off()
  
}
ExtractSig_centroid_subs <- function(control_profile_set, compound_profile_set, sampling_number,boundary,muttype_template,h,w,outputname) {
  
  # Test convergence of sig
  if(TRUE){
    result_profile <- NULL
    for(j in 1:dim(compound_profile_set)[2]){
      a <- RemoveBackground_single(control_profile_set,compound_profile_set[,j],1.65)
      result_profile <- cbind(result_profile,a)
    }
    
    min_b <- 1
    max_b <- 0
    for(j in 1:dim(compound_profile_set)[2]){
      if((j+1)<=dim(compound_profile_set)[2]){
        for(k in (j+1):dim(compound_profile_set)[2]){
          b <- cos_similarity(result_profile[,j],result_profile[,k])
          min_b <- ifelse(b<min_b,b,min_b)
          max_b <- ifelse(b>max_b,b,max_b)
          print(paste0("min_b:",min_b,"; max_b:",max_b))
        }
      }
    }
    stability_sig <- c(min_b,max_b)
    
  }
  
  
  compound_sig <- RemoveBackground_centroid(control_profile_set,compound_profile_set,1000,1.65)
  compound_sig$MutationType <- rownames(compound_sig)
  muts_basis_melt_summary <- merge(muttype_template, compound_sig,by="MutationType",all.x=T)
  muts_basis_melt_summary[is.na(muts_basis_melt_summary)] <- 0
  muts_basis_melt_summary <- muts_basis_melt_summary[order(muts_basis_melt_summary$Mutation),]
  
  muts_basis_melt_summary$percentage <- muts_basis_melt_summary[,"centroid"]/sum(muts_basis_melt_summary[,"centroid"])
  muts_basis_melt_summary$percentage_sd <- muts_basis_melt_summary[,"sd"]/sum(muts_basis_melt_summary[,"centroid"])
  
  write.table(muts_basis_melt_summary,paste0(outputname, "_exposure.txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  
  filename <- paste0(outputname, "_Signature.pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary, aes(x=MutationType, y=centroid,fill=Mutation))+ geom_bar(stat="identity",position="dodge")+xlab("Substitution Types")+ylab("Count")
  p <- p+geom_errorbar(aes(ymin=centroid-sd,ymax=centroid+sd),color="blue",position=position_dodge(.9),size=.2,width=0.5)#+scale_y_continuous(limits=c(-1,1),breaks=(seq(-1,1,0.2)),labels = scales::percent)
  p <- p+scale_x_discrete(limits = as.character(muts_basis_melt_summary$MutationType))+ggtitle(paste0(outputname,"_exposure"))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=4.5,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  
  o <- ggplot(data=muts_basis_melt_summary, aes(x=MutationType, y=percentage,fill=Mutation))+ geom_bar(stat="identity",position="dodge")+xlab("Substitution Types")+ylab("Percentage")
  o <- o+geom_errorbar(aes(ymin=percentage-percentage_sd,ymax=percentage+percentage_sd),color="blue",position=position_dodge(.9),size=.2,width=0.5)+scale_y_continuous(labels = scales::percent) #limits=c(0,0.3),breaks=(seq(0,1,0.3)),
  o <- o+scale_x_discrete(limits = as.character(muts_basis_melt_summary$MutationType))+ggtitle(paste0(outputname,"_signature","(stability=",round(max_b,2),")"))
  o <- o+scale_fill_manual(values=mypalette)
  o <- o+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=4.5,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  grid.arrange(p,o,ncol=1)
  dev.off()
  
  return(stability_sig)
  
}
RemoveBackground_single <- function(background_profile_set, sig_profile,boundary){
  
  # Range of background
  centroid_background <- rowMeans(background_profile_set)
  sd_background <- apply(background_profile_set,1,sd)
  boundary_background <- centroid_background+boundary*sd_background
  
  
  diff_all_boundary <- sig_profile-boundary_background
  diff_all <- sig_profile-centroid_background
  
  diff_all[which(diff_all_boundary<0)] <- 0
  return(diff_all)
  
}
RemoveBackground_centroid <- function(background_profile_set, sig_profile_set, sampling_number,boundary){
  
  # Range of background
  centroid_background <- rowMeans(background_profile_set)
  sd_background <- apply(background_profile_set,1,sd)
  boundary_background <- centroid_background+boundary*sd_background
  
  # Bootstrap subclones 
  diff_mean_all <- NULL
  for(bt_num in 1:sampling_number){
    
    RepCompound <- sig_profile_set[,sample(1:dim(sig_profile_set)[2],dim(sig_profile_set)[2],replace = T)]
    bootstrapCompound <- bootstrapGenomesfun(RepCompound)
    
    diff_all_boundary <- rowMeans(bootstrapCompound)-boundary_background
    diff_all <- rowMeans(bootstrapCompound)-centroid_background
    
    diff_all[which(diff_all_boundary<0)] <- 0
    diff_mean_all <- cbind(diff_mean_all,diff_all)
    
  }
  
  diff_centroid_all <- rowMeans(diff_mean_all)
  diff_quantile_all_sd <- apply(diff_mean_all,1,sd)
  
  pure_sig <- data.frame(cbind(diff_centroid_all,diff_quantile_all_sd))
  names(pure_sig) <- c("centroid","sd")
  
  return(pure_sig)
  
}
FindDinucleotides <- function(allsubs) {
  dinuclist <- NULL
  samplelist <- data.frame(table(allsubs$Sample))
  for(i in 1:dim(samplelist)[1]){
    print(i)
    samplesubs <- allsubs[allsubs$Sample==samplelist[i,1],]
    a <- samplesubs$Pos
    b <- samplesubs$Ref
    c <- samplesubs$Alt
    samplesubs$Pos_neighbor <- c(a[-1],0)
    samplesubs$neigbor_dist <- samplesubs$Pos-samplesubs$Pos_neighbor
    samplesubs$Ref_neighbor <- c(b[-1],"N")
    samplesubs$Alt_neighbor <- c(c[-1],"N")
    
    dinuc_index <- which(samplesubs$neigbor_dist==-1)
    print(length(dinuc_index))
    if(length(dinuc_index)>0){
      dinuclist <- rbind(dinuclist,samplesubs[dinuc_index,])
      dinuclist <- rbind(dinuclist,samplesubs[dinuc_index+1,])
    }
  }
  return(dinuclist)
}
bootstrapGenomesfun2 <- function(genomes,n){
  
  return(apply(genomes, 2, function(x) rmultinom(1, n, x)))
}
# For a given profile, calculate the similarity of bootstraped samples and the given profile with different mutation numbers
Bootstrap_profile_similarity <- function(profile,bsnum=20,count_max,gap=30,output){
  profile_replicates <- matrix(rep(profile,bsnum),ncol = bsnum)
  test_num <- c(seq(10,count_max,by=gap),count_max)
  result <- NULL
  for(i in 1:length(test_num)){
    profile_bootstraps <- bootstrapGenomesfun2(profile_replicates,test_num[i])
    cossim <- apply(profile_bootstraps,2,function(x) cos_similarity(x,profile))
    result <- rbind(result,c(test_num[i],mean(cossim),sd(cossim)))
  }
  result <- as.data.frame(result)
  names(result) <- c("count","simi_mean","simi_sd")
  pdf(file=paste0(output,".pdf"), onefile=TRUE,height=5,width=5, useDingbats=FALSE)
  q <- ggplot(result, aes(x=count, y=simi_mean)) + geom_point()
  q <- q+geom_errorbar(aes(ymin=simi_mean-simi_sd,ymax=simi_mean+simi_sd),color="black",position=position_dodge(.9),width=.1)+scale_y_continuous()
  q <- q+theme(axis.text.x=element_text(size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  print(q)
  dev.off()
  return(result)
}

# Return single vector
RemoveBackground_vector_single <- function(background_profile, sig_profile,sampling_number, start_num,boundary){
  
  # Remove weak mutation types in sig_profile
  removeWeakMutationTypes <- 0.01
  genomesOriginal <- as.data.frame(sig_profile)
  Totalmutations <- sum(sig_profile)
  removeMutations_max <- removeWeakMutationTypes * Totalmutations
  removeIdx <- which(cumsum(rowSums(genomesOriginal[,1:dim(genomesOriginal)[2],drop = FALSE])[order(rowSums(genomesOriginal[,1:dim(genomesOriginal)[2],drop = FALSE]))])<= min(removeMutations_max,4))
  mutationTypesToRemoveSet <- order(rowSums(genomesOriginal[,1:dim(genomesOriginal)[2],drop = FALSE]))[removeIdx]
  genomesReducted <- as.data.frame(genomesOriginal[-mutationTypesToRemoveSet,])
  reducedMutationtypes <- dim(genomesReducted)[1]
  
  # bootstrap sig profile
  centroid_sig <- rowMeans(genomesReducted[,1:dim(genomesReducted)[2],drop = FALSE])
  RepSig <- matrix(rep(centroid_sig,sampling_number),ncol = sampling_number)
  Sig_bootstraps <- bootstrapGenomesfun2(RepSig,sum(centroid_sig))
  
  i=start_num
  reachLimit <- FALSE
  
  while(i<sum(centroid_sig) & !reachLimit){
    # Remove the same weak mutation types in background profile
    backgroundReducted <- background_profile[-mutationTypesToRemoveSet]
    RepControl <- matrix(rep(backgroundReducted,sampling_number),ncol = sampling_number)
    bg_bootstraps <- bootstrapGenomesfun2(RepControl,i)
    
    # Range of background
    centroid_background <- rowMeans(bg_bootstraps)
    sd_background <- apply(bg_bootstraps,1,sd)
    boundary_background <- centroid_background+boundary*sd_background
    
    # Range of sig
     centroid_sig <- rowMeans(Sig_bootstraps)
    sd_sig <- apply(Sig_bootstraps,1,sd)
    boundary_sig <- centroid_sig+boundary*sd_sig
    
    
    diff_all_boundary <- boundary_sig-centroid_background
    diff_all <- centroid_sig-centroid_background
    
    if(length(which(diff_all_boundary<0))>0){
      reachLimit <- TRUE
    }
    
    if(length(which(diff_all_boundary<0))==0){
      diff_all_boundary_save <- diff_all_boundary
      diff_all_save <- diff_all
      diff_all_save[which(diff_all_save<0)] <- 0
    }
    
    
    
    i = i+1
    
  }
  # Add Weak mutations for KO expoure
  exposure <- rep(0,96)
  origArrayIndex <- 1
  for(i in 1:96){
    if(! i %in% mutationTypesToRemoveSet){
      exposure[i] <- diff_all_save[origArrayIndex]
      origArrayIndex=origArrayIndex+1
    }
  }
  
  # Add Weak mutations for background
  background_exposure<- rep(0,96)
  origArrayIndex <- 1
  for(i in 1:96){
    if(! i %in% mutationTypesToRemoveSet){
      background_exposure[i] <- centroid_background[origArrayIndex]
      origArrayIndex=origArrayIndex+1
    }
  }
  
  
  return(data.frame("KO_exposure"=exposure,"background_exposure"=background_exposure))
  
  
}


# Wrap for subtracting the background from the target signal to estimate the KO-associated signature
Wrap_KOSig <- function(MutCatalogue,bg_column,ko_column,sampling_number, start_num,boundary,outputname){
  
  KOSig <- RemoveBackground_vector_single(MutCatalogue[,bg_column], MutCatalogue[,ko_column],sampling_number, start_num,boundary)
  KOSig$MutationType <- MutCatalogue[,"MutationType"]
  write.table(KOSig,paste0(outputname,".txt"),sep = "\t",col.names = T, row.names = F, quote = F)
  plotCountbasis(KOSig,1,6,9,paste0(outputname,".pdf"))
  plotPercentagebasis(KOSig,1,6,9,paste0(outputname,"_percentage.pdf"))
  
}



#########################################
#
# Analysis share mutations
#
#########################################
# Remove mutations shown in parental clone
Remove_parentMuts <- function(sublist,outputname){
  samplelist <- data.frame(table(sublist$Sample))
  sublist$IsParentMut <- "DenovoMuts"
  for(i in 1:dim(samplelist)[1]){
    current_parent <- unique(sublist[sublist$Sample==as.character(samplelist[i,1]),"parent"])
    sublist[sublist$Sample==samplelist[i,1] & sublist$mut_ID %in% sublist[sublist$Sample==current_parent,"mut_ID"],"IsParentMut"] <- "ParentSharedMuts"
  }
  
  a <- data.frame(table(sublist$Sample,sublist$IsParentMut))
  names(a) <- c("Sample","IsParentMut","count")
  a_dcast <- dcast(a,Sample~IsParentMut,value.var = "count")
  write.table(a_dcast,outputname,sep="\t",col.names = T, row.names = F, quote = F)
  
  b <- sublist[sublist$IsParentMut=="DenovoMuts" | sublist$parent == sublist$Sample, ]
  return(b)
}
# Remove child mutations shown from parental clone
Remove_childMuts <- function(sublist,outputname){
  parentlist <- data.frame(table(sublist$parent))
  sublist$IschildMuts <- "No"
  for(i in 1:dim(parentlist)[1]){
    current_parent <- sublist[sublist$Sample==as.character(parentlist[i,1]),]
    parents_childmut <- sublist[sublist$parent==as.character(parentlist[i,1]) & sublist$Sample!=as.character(parentlist[i,1]),]
    sublist[sublist$Sample==samplelist[i,1] & sublist$mut_ID %in% sublist[sublist$Sample==current_parent,"mut_ID"],"IsParentMut"] <- "Yes"
  }
  
  a <- data.frame(table(sublist$Sample,sublist$IsParentMut))
  write.table(a,"IsParentMuts.txt",sep="\t",col.names = T, row.names = F, quote = F)
  
  b <- sublist[sublist$IsParentMut=="No" | sublist$parent == sublist$Sample, ]
  return(b)
}

Summary_SharedMutNumber <- function(sublist){
  samplelist <- data.frame(table(sublist$Sample))
  sublist$IsParentMut <- "No"
  for(i in 1:dim(samplelist)[1]){
    current_parent <- unique(sublist[sublist$Sample==as.character(samplelist[i,1]),"parent"])
    sublist[sublist$Sample==samplelist[i,1] & sublist$mut_ID %in% sublist[sublist$Sample==current_parent,"mut_ID"],"IsParentMut"] <- "Yes"
  }
  
  a <- data.frame(table(sublist$Sample,sublist$IsParentMut))
  names(a) <- c("Sample","ParentMut","count")
  write.table(a,"SharedMutNumber.txt",sep="\t",col.names = T, row.names = F, quote = F)
  
  # Number of childs containing the same mutation with the same parent
  ChildSubs <- sublist[sublist$clonelevel=="child",]
  ParentSubs <- sublist[sublist$clonelevel=="parent",]
  
  parentlist <- data.frame(table(sublist$parent))
  ParentSubsChildinfoAll <- NULL
  ChildSubsinfoAll <- NULL
  
  for(i in 1:dim(parentlist)[1]){
    print(paste0(i,"/",dim(parentlist)[1]))
    ParentSubs_p <- ParentSubs[ParentSubs$parent==parentlist[i,1],]
    ChildSubs_p <- ChildSubs[ChildSubs$parent==parentlist[i,1],]
    
    ChildSubs_p_freq <- data.frame(table(ChildSubs_p$mut_ID))
    names(ChildSubs_p_freq) <- c("mut_ID","freq")
    #ChildSubsAvg <- ddply(ChildSubs_p, c("mut_ID"),summarise, N=length(PM.Tum),mean=mean(PM.Tum),sd=sd(PM.Tum),se=sd/sqrt(N))
    ParentSubs_childinfo <- merge(ParentSubs_p,ChildSubs_p_freq,by="mut_ID",all.x=T)
    ParentSubs_childinfo[is.na(ParentSubs_childinfo)] <- 0
    ParentSubsChildinfoAll <- rbind(ParentSubsChildinfoAll,ParentSubs_childinfo)
    
    
    # Number of child containing the same parent mutation for each parent
    
    # ChildSubs_p_freq <- data.frame(table(ChildSubs_p$mut_ID))
    
    ChildSubs_info <- merge(ChildSubs_p,ChildSubs_p_freq,by="mut_ID",all.x=T)
    ChildSubs_info[is.na(ChildSubs_info)] <- 0
    ChildSubs_info$parent <- 0
    ChildSubs_info[ChildSubs_info$mut_ID%in%ParentSubs_p$mut_ID,"parent"] <- 1
    ChildSubsinfoAll <- rbind(ChildSubsinfoAll,ChildSubs_info)
    
  }
  ParentSubsChildinfoAll_summary <- data.frame(table(ParentSubsChildinfoAll$Sample,ParentSubsChildinfoAll$freq))
  names(ParentSubsChildinfoAll_summary) <- c("Sample","ShownInchildnum","freq")
  write.table(ParentSubsChildinfoAll_summary,"ParentSubsChildinfoAll_summary.txt",sep="\t",col.names = T, row.names = F, quote = F)
  write.table(ParentSubsChildinfoAll,"ParentSubsChildinfoAll.txt",sep="\t",col.names = T, row.names = F, quote = F)
  
  ChildSubsinfoAll_summary <- data.frame(table(ChildSubsinfoAll$Sample,ChildSubsinfoAll$freq))
  names(ChildSubsinfoAll_summary) <- c("Sample","ShownInchildnum","freq")
  write.table(ChildSubsinfoAll_summary,"ChildSubsinfoAll_summary.txt",sep="\t",col.names = T, row.names = F, quote = F)
  write.table(ChildSubsinfoAll,"ChildSubsinfoAll.txt",sep="\t",col.names = T, row.names = F, quote = F)
  
  
}

# Remove mutations list A from mutation list B
Remove_Muts <- function(sublistA, sublistB){
  muts_sublistA$mut_pos <- paste0(muts_sublistA$Chrom,"_",muts_sublistA$Pos)
  muts_sublistB$mut_pos <- paste0(muts_sublistB$Chrom,"_",muts_sublistB$Pos)
  
  muts_denovo <- muts_sublistB[!muts_sublistB$mut_pos %in% muts_sublistA$mut_pos,]
  return(muts_denovo)
}

# group is usually the same clone (with the same treatment)
Remove_sharedmuts <- function(sublist,group){
  sublist$mut_ID <- paste(sublist$Chrom,sublist$Pos,sublist[,group],sep = "_")
  mutslist <- data.frame(table(sublist$mut_ID))
  names(mutslist) <- c("mut_ID","Freq")
  sublist_denovo <- sublist[sublist$mut_ID %in% mutslist[mutslist$Freq==1,"mut_ID"],]
  
  return(sublist_denovo) 
}

# shared mutation by all subclones in a treatment (subclone_col usually use "Sample" column)
Sharedmuts <- function(sublist,subclone_col){
  sublist$mut_ID <- paste(sublist$Chrom,sublist$Pos,sep = "_")
  mutslist <- data.frame(table(sublist$mut_ID))
  names(mutslist) <- c("mut_ID","Freq")
  subclone_num <- length(table(sublist[,subclone_col]))
  sublist_shared <- sublist[sublist$mut_ID %in% mutslist[mutslist$Freq==subclone_num,"mut_ID"],]
  
  return(sublist_shared) 
}

# Find shared mutation number between two samples
# subclone_col is usually the Sample column
# clone_col is usually the treatment column 
ShareMuts_summary <- function(mutlist,clone_col, subclone_col){
  mutlist$mut_pos <- paste0(mutlist$Chrom,"_",mutlist$Pos)
  subclone_list <- data.frame(table(mutlist[,subclone_col]))
  names(subclone_list) <- c("Subclone","all_num")
  share_mut_list <- NULL
  share_mut_num <- NULL
  for(i in 1:(dim(subclone_list)[1]-1)){
    print(i)
    i_muts <- mutlist[mutlist[,subclone_col]==as.character(subclone_list[i,"Subclone"]),]
    for(j in (i+1):dim(subclone_list)[1]){
      j_muts <- mutlist[mutlist[,subclone_col]==as.character(subclone_list[j,"Subclone"]),]
      shared <- length(intersect(i_muts$mut_pos,j_muts$mut_pos))
      share_mut_list <- rbind(share_mut_list,c(as.character(subclone_list[i,"Subclone"]),as.character(subclone_list[j,"Subclone"])))
      share_mut_num <- c(share_mut_num,shared)
    }
  }
  share_mut_list <- data.frame(share_mut_list)
  share_mut_list$shared_num <- share_mut_num
  names(share_mut_list) <- c("SubcloneA","SubcloneB","shared_num")
  share_mut_list_copy <- share_mut_list
  subclone_clone <- data.frame(unique(mutlist[,c(clone_col, subclone_col)]))
  share_mut_list$cloneA <- apply(share_mut_list,1,function(x) subclone_clone[subclone_clone[,subclone_col]==x[1],clone_col])
  share_mut_list$cloneB <- apply(share_mut_list,1,function(x) subclone_clone[subclone_clone[,subclone_col]==x[2],clone_col])
  share_mut_list$sameCondition <- "No"
  share_mut_list[share_mut_list$cloneA==share_mut_list$cloneB,]$sameCondition <- "Yes"
  
  write.table(share_mut_list,"share_mut_list.txt",sep = "\t",col.names = T, row.names = F, quote = F)
  
}
ShareMuts_UpSetR <- function(mutlist,h,w,outputname){
  mutlist$mut_pos <- paste0(mutlist$Chrom,"_",mutlist$Pos)
  listInput <- list()
  samplelist <- data.frame(table(mutlist$Sample))
  for(i in 1:dim(samplelist)[1]){
    name <- paste0(samplelist[i,1])
    listInput[[name]] <- mutlist[mutlist$Sample==samplelist[i,1],"mut_pos"]
  }
  pdf(file=paste0(outputname,".pdf"), onefile=F,width=w,height=h)
  upset(fromList(listInput), sets=as.character(samplelist[,1]), nsets = dim(samplelist)[1],nintersects = 100,order.by = "freq",keep.order = TRUE)
  dev.off()
  
}

#########################################
#
# Topography analysis
#
#########################################
regufea_length <- read.table("/nfs/cancer_archive04/xz3/a_1242/00_common/RegulatoryFeatures_length.txt", sep = "\t", header = T, as.is = T)
replitime_length <- read.table("/nfs/cancer_archive04/xz3/a_1242/00_common/MCF7_RepliTime_length.txt", sep = "\t", header = T, as.is = T)
transtrand_length <- read.table("/nfs/cancer_archive04/xz3/a_1242/00_common/TranscriptionalStrand_length.txt", sep = "\t", header = T, as.is = T)
replitrand_length <- read.table("/nfs/cancer_archive04/xz3/a_1242/00_common/ReplicativeStrand_length.txt", sep = "\t", header = T, as.is = T)
NucleoPos_length <- read.table("/nfs/cancer_archive04/xz3/a_1242/00_common/NucleosomePos_length41.txt", sep = "\t", header = T, as.is = T)

ATGC <- read.table("/nfs/cgp_signatures/Analysis/xz3/00_Resource/GenomeFeature/00_common/ATGC_counts_genome.txt", sep = "\t", header = T, as.is = T)




# bed file is 0-based, half-closed-half-open 
Tab2Bed_generic <- function(inputfile,ChromCol,Pos1Col,Pos2Col,VID,outputname){ # bed file is 0-based, half-closed-half-open 
  muts_bed <- inputfile[,c(ChromCol,Pos1Col,Pos2Col,VID)]
  muts_bed[,ChromCol] <- paste0("chr",as.character(muts_bed[,ChromCol]))
  muts_bed[muts_bed[,ChromCol]=="chr23",ChromCol] <- "chrX"
  muts_bed[muts_bed[,ChromCol]=="chr24",ChromCol] <- "chrY"
  muts_bed[,2] <- muts_bed[,2]-1
  write.table(muts_bed,outputname,sep="\t",col.names = F, row.names = F, quote = F)
}

Tab2Bed <- function(muts, outputname){ # bed file is 0-based, half-closed-half-open 
  muts_bed <- muts[,c("Chrom","Pos","Pos","VariantID")]
  muts_bed$Chrom <- paste0("chr",muts_bed$Chrom)
  muts_bed[muts_bed$Chrom=="chr23","Chrom"] <- "chrX"
  muts_bed[muts_bed$Chrom=="chr24","Chrom"] <- "chrY"
  muts_bed[,2] <- muts_bed[,2]-1
  write.table(muts_bed,outputname,sep="\t",col.names = F, row.names = F, quote = F)
}
AddStrandInfo <- function(mutfile1, mutfile2,muts_context,outputname){
  mut_strand1 <- read.table(mutfile1,sep = "\t",header = F,as.is = T)
  mut_strand2 <- read.table(mutfile2,sep = "\t",header = F,as.is = T)
  mut_strand <- rbind(mut_strand1,mut_strand2)
  mut_strand_short <- mut_strand[,c(4,8)]
  names(mut_strand_short) <- c("VariantID","Strand_original")
  
  muts_withstrandinfo <- merge(muts_context,mut_strand_short,by="VariantID",all.x=T)
  muts_withstrandinfo[is.na(muts_withstrandinfo)] <- "others"
  muts_withstrandinfo[(muts_withstrandinfo$Strand_original == "Leading" | muts_withstrandinfo$Strand_original == 1),]$Strand_original <- "leading_uts"
  muts_withstrandinfo[(muts_withstrandinfo$Strand_original == "Lagging" | muts_withstrandinfo$Strand_original == -1),]$Strand_original <- "lagging_ts"
  
  muts_withstrandinfo$Strand <- muts_withstrandinfo$Strand_original
  
  CTsubs_copy <- muts_withstrandinfo
  CTsubs <- muts_withstrandinfo
  CTsubs$Alt2 <- CTsubs$Alt
  CTsubs$Ref2 <- CTsubs$Ref
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt2 <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref2 <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  
  CTsubs[(CTsubs$Ref2!=CTsubs$Ref & CTsubs$Strand_original=="leading_uts"),]$Strand <- "lagging_ts"
  CTsubs[(CTsubs$Ref2!=CTsubs$Ref & CTsubs$Strand_original=="lagging_ts"),]$Strand <- "leading_uts"
  
  
  write.table(CTsubs,outputname,sep="\t",col.names = T, row.names = F, quote = F)
  return(CTsubs)
}

AddRepliTimeInfo <- function(mutlist, mutBedfile,featureBedfile,intersectResultfile,outputfilename){
  Tab2Bed(mutlist,mutBedfile)
  intersectBed_command <- paste0("intersectBed -a ",mutBedfile," -b ",featureBedfile," -wo > ", intersectResultfile)
  
  #"nohup intersectBed -a denovo_muts_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/MCF7.all.compact.bed -wo > subs_ReplicatingTime.txt &"
  try(system(intersectBed_command))
  
  subs_ReplicatingTime <- read.table(intersectResultfile,sep = "\t",header = F,as.is = T)
  subs_ReplicatingTime_short <- subs_ReplicatingTime[,c(4,9)]
  names(subs_ReplicatingTime_short) <- c("VariantID","ReplicatingTime")
  subs_ReplicatingTime_short$ReplicatingTime=abs(subs_ReplicatingTime_short$ReplicatingTime-11)
  denovo_muts_rt <- merge(mutlist,subs_ReplicatingTime_short,by="VariantID",all.x=T)
  denovo_muts_rt[is.na(denovo_muts_rt)] <- "others"
  write.table(denovo_muts_rt,outputfilename,sep = "\t", col.names = T, row.names = F, quote = F)
  
}
AddStrandInfo_intersect <- function(mutlist, mutBedfile,featureBedfile_strand1,featureBedfile_strand2,intersectResultfile,outputfilename){
  Tab2Bed(mutlist,mutBedfile)
  intersectBed_command_1 <- paste0("intersectBed -a ",mutBedfile," -b ",featureBedfile_strand1," -wo > ", paste0(intersectResultfile,"_1.txt"))
  intersectBed_command_2 <- paste0("intersectBed -a ",mutBedfile," -b ",featureBedfile_strand2," -wo > ", paste0(intersectResultfile,"_2.txt"))
  
  # nohup intersectBed -a denovo_muts_control_mutagen_bed.txt -b ./00_data/TranscribStrand.uts.txt -wo > subs_control_mutagen_uts.txt &
  # nohup intersectBed -a denovo_muts_control_mutagen_bed.txt -b ./00_data/TranscribStrand.ts.txt -wo > subs_control_mutagen_ts.txt &
  try(system(intersectBed_command_1))
  try(system(intersectBed_command_2))
  
  AddStrandInfo(paste0(intersectResultfile,"_1.txt"),paste0(intersectResultfile,"_2.txt"),mutlist,outputfilename)
}
AddFeatureInfo <- function(mutlist, mutBedfile,featureBedfile,intersectResultfile,outputfilename){
  Tab2Bed(mutlist,mutBedfile)
  intersectBed_command <- paste0("intersectBed -a ",mutBedfile," -b ",featureBedfile," -wo > ", intersectResultfile)
  
  #"nohup intersectBed -a denovo_muts_bed.txt -b /nfs/cancer_archive04/xz3/a_1242/00_common/MCF7.all.compact.bed -wo > subs_ReplicatingTime.txt &"
  try(system(intersectBed_command))
  
  subs_Feature <- read.table(intersectResultfile,sep = "\t",header = F,as.is = T)
  subs_Feature_short <- subs_Feature[,c(4,8)]
  names(subs_Feature_short) <- c("VariantID","Feature")
  denovo_muts_rt <- merge(mutlist,subs_Feature_short,by="VariantID",all.x=T)
  denovo_muts_rt[is.na(denovo_muts_rt)] <- "NonRegulatory"
  write.table(denovo_muts_rt,outputfilename,sep = "\t", col.names = T, row.names = F, quote = F)
  
}

Rintersect <- function(mutlist, mutBedfile,featureBedfile,intersectResultfile){
  Tab2Bed(mutlist,mutBedfile)
  intersectBed_command <- paste0("intersectBed -a ",mutBedfile," -b ",featureBedfile," -wo > ", paste0(intersectResultfile,".txt"))
  
  try(system(intersectBed_command))
}

Rintersect_generic <- function(mutlist,ChromCol,Pos1Col,Pos2Col,VID, mutBedfile,featureBedfile,intersectResultfile){
  Tab2Bed_generic(mutlist,ChromCol,Pos1Col,Pos2Col,VID, mutBedfile)
  intersectBed_command <- paste0("intersectBed -a ",mutBedfile," -b ",featureBedfile," -wo > ", paste0(intersectResultfile,".txt"))
  
  try(system(intersectBed_command))
}

# Includding simulation on replication timing distribtution of mutSig
ReplicationTiming_Sample <- function(muts_replitime,SampleCol,ReptimeCol, DoSimulation="TRUE", outputname){
  
  replitime_observed <- data.frame(table(muts_replitime[,SampleCol], muts_replitime[,ReptimeCol]))
  names(replitime_observed) <- c(SampleCol,"ReplicatingTime","observed_Freq")
  replitime_observed <- replitime_observed[replitime_observed[,ReptimeCol] != "others",]
  # simulate the expected distribution 
  if(DoSimulation){
    bootstrap_num <- 100
    muts_replitime[muts_replitime$Chrom=="23","Chrom"]="X"
    muts_replitime[muts_replitime$Chrom=="24","Chrom"]="Y"
    
    trinuc_catalogue <- Gen32Catalogue(muts_replitime[muts_replitime$ReplicatingTime !="others",],SampleCol)
    trinuc_replitime <- read.table("/nfs/cgp_signatures/Analysis/xz3/00_Resource/GenomeFeature/00_common/trinuc_replitime_1early.txt", sep = "\t", header = T, as.is = T)
    
    # According to distribution of trinuc on replication timing regions (trinuc_replitime), bootstrap this distribution
    # bootstrapping method to evaluate the difference between expected distribution and observed one
    # bootstrapping is used to construct a population of expected distributions
    
    
    # loop for sample
    for(i in 2:dim(trinuc_catalogue)[2]){
      current_sample <- trinuc_catalogue[,c(1,i)]
      replitime_expected <- data.frame("ReplicatingTime"=c("1","2","3","4","5","6","7","8","9","10"))
      
      # loop for bootstrap numbers
      for(j in 1:bootstrap_num){
        bootstrap_replitime_all <- NULL
        
        # loop for trinucleotides
        for(k in 1:dim(current_sample)[1]){
          current_context <- current_sample[k,1]
          trinuc_replitime_current <- trinuc_replitime[trinuc_replitime$trinuc_pyrimidine==current_context,]
          bootstrap_replitime=sample(trinuc_replitime_current$ReplicatingTime, current_sample[k,2],replace = T, prob = trinuc_replitime_current$sum/sum(trinuc_replitime_current$sum))
          bootstrap_replitime_all <- c(bootstrap_replitime_all,bootstrap_replitime)
        } # k
        current_reptime <- data.frame(table(bootstrap_replitime_all))
        names(current_reptime) <- c("ReplicatingTime",paste0("bs",j))
        replitime_expected <- merge(replitime_expected,current_reptime,by="ReplicatingTime",all.x=T)
      } # j
      
      replitime_expected$expected_mean <- rowMeans(replitime_expected[,2:dim(replitime_expected)[2]])
      replitime_expected$expected_sd <- apply(replitime_expected[,2:dim(replitime_expected)[2]], 1, sd)
      replitime_observed_sample <- replitime_observed[replitime_observed[,SampleCol]==colnames(current_sample)[2],]
      
      # combine simulation (expected) and observed distributions together for each sample
      replitime_observed_expected <- merge(replitime_observed_sample,replitime_expected[,c("ReplicatingTime","expected_mean","expected_sd")],by="ReplicatingTime")
      replitime_observed_expected$ReplicatingTime <- as.numeric(as.character(replitime_observed_expected$ReplicatingTime))
      replitime_observed_expected <- replitime_observed_expected[order(replitime_observed_expected$ReplicatingTime),]
      # plot
      filename <- paste0(outputname,"_",colnames(current_sample)[2],"_simu.pdf")
      pdf(file=filename, onefile=TRUE,width=6,height=6, useDingbats=FALSE)
      p <- ggplot(data=replitime_observed_expected, aes(x=ReplicatingTime, y=observed_Freq))+ geom_bar(stat="identity",position="dodge",fill="#3cb44b")+xlab("Replication timing region")+ylab("Count")
      p <- p+geom_point(data=replitime_observed_expected,aes(x=ReplicatingTime,y=expected_mean),color="blue")+geom_line(data=replitime_observed_expected,aes(x=ReplicatingTime,y=expected_mean, group=1),color="blue")
      p <- p+geom_errorbar(data=replitime_observed_expected,aes(ymin=expected_mean-expected_sd,ymax=expected_mean+expected_sd),color="blue",position=position_dodge(.9),size=.2,width=0.5) #+scale_y_continuous(limits=c(0, 350),breaks=seq(0, 350, 100))
      p <- p+ggtitle(paste0(colnames(current_sample)[2]))
      p <- p+theme(axis.text.x=element_text(size=15,colour = "black"),
                   axis.text.y=element_text(size=15,colour = "black"),
                   axis.title.x = element_text(size=15),
                   axis.title.y = element_text(size=15),
                   plot.title = element_text(size=10),
                   panel.grid.minor.x=element_blank(),
                   panel.grid.major.x=element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   panel.background = element_rect(fill = "white"),
                   panel.border = element_rect(colour = "black", fill=NA))
      print(p)
      dev.off()
      write.table(replitime_observed_expected,paste0(outputname,"_",colnames(current_sample)[2],"_simu.txt"),sep = "\t",col.names = T, row.names = F, quote = F)
      
    } # i
    
  }else{
    
    # loop for sample
    for(i in 2:dim(trinuc_catalogue)[2]){
      current_sample <- trinuc_catalogue[,c(1,i)]
      replitime_observed_sample <- replitime_observed[replitime_observed[,SampleCol]==colnames(current_sample)[2],]
      replitime_observed_sample$ReplicatingTime <- as.numeric(as.character(replitime_observed_sample$ReplicatingTime))
      replitime_observed_sample <- replitime_observed_sample[order(replitime_observed_sample$ReplicatingTime),]
      
      # plot
      filename <- paste0(outputname,"_",colnames(current_sample)[2],"_observed.pdf")
      pdf(file=filename, onefile=TRUE,width=6,height=6, useDingbats=FALSE)
      p <- ggplot(data=replitime_observed_sample, aes(x=ReplicatingTime, y=observed_Freq))+ geom_bar(stat="identity",position="dodge",fill="#3cb44b")+xlab("Replication timing region")+ylab("Count")
      p <- p+ggtitle(paste0(colnames(current_sample)[2]))
      p <- p+theme(axis.text.x=element_text(size=15,colour = "black"),
                   axis.text.y=element_text(size=15,colour = "black"),
                   axis.title.x = element_text(size=15),
                   axis.title.y = element_text(size=15),
                   plot.title = element_text(size=10),
                   panel.grid.minor.x=element_blank(),
                   panel.grid.major.x=element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   panel.background = element_rect(fill = "white"),
                   panel.border = element_rect(colour = "black", fill=NA))
      print(p)
      dev.off()
      write.table(replitime_observed_sample,paste0(outputname,"_",colnames(current_sample)[2],"_observed.txt"),sep = "\t",col.names = T, row.names = F, quote = F)
      
      
    }
  }
  
}

# Compare feature with others; 
Feature_Sample <- function(muts_feature,SampleCol,FeatureCol,trinuc_feature, FeatureName,outputname){
  
  feature_observed <- data.frame(table(muts_feature[,SampleCol], muts_feature[,FeatureCol]))
  names(feature_observed) <- c(SampleCol,"Feature","observed_Freq")
  # simulate the expected distribution 

    bootstrap_num <- 100
    muts_feature[muts_feature$Chrom=="23","Chrom"]="X"
    muts_feature[muts_feature$Chrom=="24","Chrom"]="Y"
    
    trinuc_catalogue <- Gen32Catalogue(muts_feature,SampleCol)
   # trinuc_feature <- read.table(TrinucFeatureFile, sep = "\t", header = F, as.is = T)
   # names(TrinucFeature) <- c("trinuc_pyrimidine","Feature","sum")
    # According to distribution of trinuc on replication timing regions (trinuc_replitime), bootstrap this distribution
    # bootstrapping method to evaluate the difference between expected distribution and observed one
    # bootstrapping is used to construct a population of expected distributions
    
    # loop for sample
    for(i in 2:dim(trinuc_catalogue)[2]){
      current_sample <- trinuc_catalogue[,c(1,i)]
      feature_expected <- data.frame(table(muts_feature[,FeatureCol]))
      names(feature_expected) <- c("Feature","Observed")
      # loop for bootstrap numbers
      for(j in 1:bootstrap_num){
        bootstrap_Feature_all <- NULL
        
        # loop for trinucleotides
        for(k in 1:dim(current_sample)[1]){
          current_context <- as.character(current_sample[k,1])
          trinuc_feature_current <- trinuc_feature[trinuc_feature$trinuc_pyrimidine==current_context,]
          bootstrap_Feature=sample(trinuc_feature_current$Feature, current_sample[k,2],replace = T, prob = trinuc_feature_current$sum/sum(trinuc_feature_current$sum))
          bootstrap_Feature_all <- c(bootstrap_Feature_all,bootstrap_Feature)
        } # k
        current_feature <- data.frame(table(bootstrap_Feature_all))
        names(current_feature) <- c("Feature",paste0("bs",j))
        feature_expected <- merge(feature_expected,current_feature,by="Feature",all.x=T)
      } # j
      feature_expected <- feature_expected[,-2]
      feature_expected$expected_mean <- rowMeans(feature_expected[,2:dim(feature_expected)[2]])
      feature_expected$expected_sd <- apply(feature_expected[,2:dim(feature_expected)[2]], 1, sd)
      feature_expected_sample <- feature_observed[feature_observed[,SampleCol]==colnames(current_sample)[2],]
      
      # combine simulation (expected) and observed distributions together for each sample
      feature_observed_expected <- merge(feature_expected_sample,feature_expected[,c("Feature","expected_mean","expected_sd")],by="Feature")
       # plot
      filename <- paste0(outputname,"_",colnames(current_sample)[2],"_simu.pdf")
      pdf(file=filename, onefile=TRUE,width=6,height=6, useDingbats=FALSE)
      p <- ggplot(data=feature_observed_expected, aes(x=Feature, y=observed_Freq))+ geom_bar(stat="identity",position="dodge",fill="#3cb44b")+xlab("Feature")+ylab("Count")
      p <- p+geom_point(data=feature_observed_expected,aes(x=Feature,y=expected_mean),color="blue")
      p <- p+geom_errorbar(data=feature_observed_expected,aes(ymin=expected_mean-expected_sd,ymax=expected_mean+expected_sd),color="blue",position=position_dodge(.9),size=.2,width=0.5) #+scale_y_continuous(limits=c(0, 350),breaks=seq(0, 350, 100))
      p <- p+ggtitle(paste0(colnames(current_sample)[2]))
      p <- p+theme(axis.text.x=element_text(size=15,colour = "black"),
                   axis.text.y=element_text(size=15,colour = "black"),
                   axis.title.x = element_text(size=15),
                   axis.title.y = element_text(size=15),
                   plot.title = element_text(size=10),
                   panel.grid.minor.x=element_blank(),
                   panel.grid.major.x=element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   panel.background = element_rect(fill = "white"),
                   panel.border = element_rect(colour = "black", fill=NA))
      print(p)
      dev.off()
      
      feature_observed_expected$chisq_pvalue <- chisq.test(c(feature_observed_expected[1,"observed_Freq"], feature_observed_expected[2,"observed_Freq"]), p=c(feature_observed_expected[1,"expected_mean"], feature_observed_expected[2,"expected_mean"]),rescale.p = TRUE)$p.value
      feature_observed_expected$binom_pvalue <- binom.test(feature_observed_expected[feature_observed_expected$Feature==FeatureName,"observed_Freq"],sum(feature_observed_expected$observed_Freq),p=feature_observed_expected[feature_observed_expected$Feature==FeatureName,"expected_mean"]/feature_observed_expected[feature_observed_expected$Feature=="NonRegulatory","expected_mean"])$p.value
      
      write.table(feature_observed_expected,paste0(outputname,"_",colnames(current_sample)[2],"_simu.txt"),sep = "\t",col.names = T, row.names = F, quote = F)
      
    } # i
    
   
    
      

      
    
       
    
    
}

Feature_Sample_chisq <- function(muts_feature,SampleCol,FeatureCol,trinuc_feature,outputname){
  
  feature_observed <- data.frame(table(muts_feature[,SampleCol], muts_feature[,FeatureCol]))
  names(feature_observed) <- c(SampleCol,"Feature","observed_Freq")
  # simulate the expected distribution 
  
  bootstrap_num <- 100
  muts_feature[muts_feature$Chrom=="23","Chrom"]="X"
  muts_feature[muts_feature$Chrom=="24","Chrom"]="Y"
  
  trinuc_catalogue <- Gen32Catalogue(muts_feature,SampleCol)
  # trinuc_feature <- read.table(TrinucFeatureFile, sep = "\t", header = F, as.is = T)
  # names(TrinucFeature) <- c("trinuc_pyrimidine","Feature","sum")
  # According to distribution of trinuc on replication timing regions (trinuc_replitime), bootstrap this distribution
  # bootstrapping method to evaluate the difference between expected distribution and observed one
  # bootstrapping is used to construct a population of expected distributions
  
  # loop for sample
  for(i in 2:dim(trinuc_catalogue)[2]){
    current_sample <- trinuc_catalogue[,c(1,i)]
    feature_expected <- data.frame(table(muts_feature[,FeatureCol]))
    names(feature_expected) <- c("Feature","Observed")
    # loop for bootstrap numbers
    for(j in 1:bootstrap_num){
      bootstrap_Feature_all <- NULL
      
      # loop for trinucleotides
      for(k in 1:dim(current_sample)[1]){
        current_context <- as.character(current_sample[k,1])
        trinuc_feature_current <- trinuc_feature[trinuc_feature$trinuc_pyrimidine==current_context,]
        bootstrap_Feature=sample(trinuc_feature_current$Feature, current_sample[k,2],replace = T, prob = trinuc_feature_current$sum/sum(trinuc_feature_current$sum))
        bootstrap_Feature_all <- c(bootstrap_Feature_all,bootstrap_Feature)
      } # k
      current_feature <- data.frame(table(bootstrap_Feature_all))
      names(current_feature) <- c("Feature",paste0("bs",j))
      feature_expected <- merge(feature_expected,current_feature,by="Feature",all.x=T)
    } # j
    feature_expected <- feature_expected[,-2]
    feature_expected$expected_mean <- rowMeans(feature_expected[,2:dim(feature_expected)[2]])
    feature_expected$expected_sd <- apply(feature_expected[,2:dim(feature_expected)[2]], 1, sd)
    feature_expected_sample <- feature_observed[feature_observed[,SampleCol]==colnames(current_sample)[2],]
    
    # combine simulation (expected) and observed distributions together for each sample
    feature_observed_expected <- merge(feature_expected_sample,feature_expected[,c("Feature","expected_mean","expected_sd")],by="Feature")
    # plot
    filename <- paste0(outputname,"_",colnames(current_sample)[2],"_simu.pdf")
    pdf(file=filename, onefile=TRUE,width=6,height=6, useDingbats=FALSE)
    p <- ggplot(data=feature_observed_expected, aes(x=Feature, y=observed_Freq))+ geom_bar(stat="identity",position="dodge",fill="#3cb44b")+xlab("Feature")+ylab("Count")
    p <- p+geom_point(data=feature_observed_expected,aes(x=Feature,y=expected_mean),color="blue")
    p <- p+geom_errorbar(data=feature_observed_expected,aes(ymin=expected_mean-expected_sd,ymax=expected_mean+expected_sd),color="blue",position=position_dodge(.9),size=.2,width=0.5) #+scale_y_continuous(limits=c(0, 350),breaks=seq(0, 350, 100))
    p <- p+ggtitle(paste0(colnames(current_sample)[2]))
    p <- p+theme(axis.text.x=element_text(size=15,colour = "black"),
                 axis.text.y=element_text(size=15,colour = "black"),
                 axis.title.x = element_text(size=15),
                 axis.title.y = element_text(size=15),
                 plot.title = element_text(size=10),
                 panel.grid.minor.x=element_blank(),
                 panel.grid.major.x=element_blank(),
                 panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank(),
                 panel.background = element_rect(fill = "white"),
                 panel.border = element_rect(colour = "black", fill=NA))
    print(p)
    dev.off()
    
    feature_observed_expected$chisq_pvalue <- chisq.test(c(feature_observed_expected[1,"observed_Freq"], feature_observed_expected[2,"observed_Freq"]), p=c(feature_observed_expected[1,"expected_mean"], feature_observed_expected[2,"expected_mean"]),rescale.p = TRUE)$p.value
 #   feature_observed_expected$binom_pvalue <- binom.test(feature_observed_expected[feature_observed_expected$Feature==FeatureName,"observed_Freq"],sum(feature_observed_expected$observed_Freq),p=feature_observed_expected[feature_observed_expected$Feature==FeatureName,"expected_mean"]/feature_observed_expected[feature_observed_expected$Feature=="NonRegulatory","expected_mean"])$p.value
    
    write.table(feature_observed_expected,paste0(outputname,"_",colnames(current_sample)[2],"_simu.txt"),sep = "\t",col.names = T, row.names = F, quote = F)
    
  } # i
  
  

  
}

plotCountbasis_average_sd_feature <- function(muts_basis,muttype_template,selectcolor,h,w,outputname){
  
  muts_basis[is.na(muts_basis)] <- 0
  mean_parentmuts <- sum(muts_basis[,2:dim(muts_basis)[2]])/(dim(muts_basis)[2]-1)
  
  muts_basis_melt <- melt(muts_basis,"targetfeature")
  names(muts_basis_melt) <- c("targetfeature","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("targetfeature"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N),snr=abs(mean/sd))
  
  muts_basis_melt_summary$mean_perc <- muts_basis_melt_summary$mean/mean_parentmuts
  muts_basis_melt_summary$sd_perc <- muts_basis_melt_summary$sd/mean_parentmuts
  muts_basis_melt_summary$se_perc <- muts_basis_melt_summary$se/mean_parentmuts
  muts_basis_melt_summary$snr_perc <- muts_basis_melt_summary$mean/muts_basis_melt_summary$sd
  muts_basis_melt_summary <- muts_basis_melt_summary[order(muts_basis_melt_summary$targetfeature),]
  
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h, useDingbats=FALSE)
  p <- ggplot(data=muts_basis_melt_summary, aes(x=targetfeature, y=mean))+ geom_bar(stat="identity",position="dodge", width=.8,fill=selectcolor)+xlab("Mutation Types")+ylab("Count")
  p <- p+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),color="blue",position=position_dodge(.9),size=.2,width=0.5) #+scale_y_continuous(limits=c(0, 350),breaks=seq(0, 350, 100))
  p <- p+scale_x_discrete(limits = as.character(muts_basis_melt_summary$targetfeature))+ggtitle(outputname)
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  q <- ggplot(data=muts_basis_melt_summary, aes(x=targetfeature, y=mean_perc))+ geom_bar(stat="identity",position="dodge", width=.8,fill=selectcolor)+xlab("Mutation Types")+ylab("Percentage")
  q <- q+geom_errorbar(aes(ymin=mean_perc-sd_perc,ymax=mean_perc+sd_perc),color="blue",position=position_dodge(.9),size=.2,width=0.5) #+scale_y_continuous(limits=c(0, 350),breaks=seq(0, 350, 100))
  q <- q+scale_x_discrete(limits = as.character(muts_basis_melt_summary$targetfeature))+ggtitle(outputname)
  q <- q+scale_fill_manual(values=mypalette)
  q <- q+theme(axis.text.x=element_text(size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  grid.arrange(p,q,ncol=1)
  dev.off()
  write.table(muts_basis_melt_summary,paste0("targetfeature","_",outputname, "_basis.txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  
  return(muts_basis_melt_summary)
  
  
  
  
}
# Plot Strand bias
StrandBias_6muttype <- function(denovo_muts_feature, SampleCol,feature,featurelength,h,w,colnum, outputname){
  
  #featurelength$all <-rowSums(featurelength[,2:5])
  
  # 1: uts/leading
  # -1: ts/lagging
  #featurelength_strand <- data.frame("Strand"=c("1","-1"),"C"=c(sum(featurelength[,"C"]), sum(featurelength[,"G"])),"T"=c(sum(featurelength[,"T"]), sum(featurelength[,"A"])))
  featurelength_strand <- data.frame("Ref"=c("C","T"),"uts_leading_wg"=c(sum(featurelength[,"C"]), sum(featurelength[,"T"])), "ts_lagging_wg"=c(sum(featurelength[,"G"]), sum(featurelength[,"A"])))
  
  denovo_muts_feature <- denovo_muts_feature[denovo_muts_feature[,feature] != "others",]
  denovo_muts_feature$Mutation <- paste0(denovo_muts_feature$Ref2,">",denovo_muts_feature$Alt2)
  denovo_muts_dis <- data.frame(table(denovo_muts_feature[,SampleCol],denovo_muts_feature[,feature], denovo_muts_feature[,"Mutation"]))
  names(denovo_muts_dis) <- c("Sample",feature,"Mutation", "Freq")
  
  gtc <- denovo_muts_dis
  gtc$targetfeature <- gtc[,feature]
  
  filename=paste0(outputname,".pdf")
  pdf(file=filename, onefile=TRUE,height=h,width = w) 
  d1 <- ggplot(gtc,aes(x=Mutation,y=Freq,fill=targetfeature))+geom_bar(stat="identity",position="dodge")
  d1 <- d1+ylab("Count")+theme(axis.title.x = element_text(size=15),
                               axis.title.y = element_text(size=15),
                               plot.title = element_text(size=10),
                               axis.text.x=element_text(angle=90, vjust=0.5,colour = "black"),
                               axis.text.y=element_text(colour = "black"),
                               panel.grid.minor.x=element_blank(),
                               panel.grid.major.x=element_blank(),
                               panel.grid.major.y = element_blank(),
                               panel.grid.minor.y = element_blank(),
                               panel.background = element_rect(fill = "white"),
                               panel.border = element_rect(colour = "black", fill=NA))
  d1 <- d1+facet_wrap(~Sample,ncol=colnum,scales="free")
  print(d1)
  dev.off()
  
  gtc_dcast <- dcast(gtc,Sample+Mutation~Strand, value.var="Freq")
  gtc_dcast$Ref <- substr(gtc_dcast$Mutation,1,1)
  gtc_dcast <- merge(gtc_dcast, featurelength_strand, by="Ref", all.x=T)
  gtc_dcast$chisq_pvalue <- 1
  for(i in 1:dim(gtc_dcast)[1]){
    if(gtc_dcast[i,"leading_uts"]+gtc_dcast[i,"lagging_ts"]>0){
      gtc_dcast[i,"chisq_pvalue"] <- chisq.test(c(gtc_dcast[i,"leading_uts"], gtc_dcast[i,"lagging_ts"]), p=c(gtc_dcast[i,"uts_leading_wg"], gtc_dcast[i,"ts_lagging_wg"]),rescale.p = TRUE)$p.value
      
    }
    
  }
  gtc_dcast$P_adjust <- p.adjust(gtc_dcast$chisq_pvalue,method = "BH")
  gtc_dcast$flag <- ""
  gtc_dcast[gtc_dcast$P_adjust<=0.05,"flag"] <- "*"
  gtc_dcast[gtc_dcast$P_adjust<=0.01,"flag"] <- "**"
  gtc_dcast[gtc_dcast$P_adjust<=0.001,"flag"] <- "***"
  
  write.table(gtc_dcast,paste0(outputname, "_chisq_adjust.txt"),sep = "\t",col.names = T, row.names = F, quote = F)
  
  
}
StrandBias_96muttype <- function(denovo_muts_feature, SampleCol,feature,featurelength,h,w,colnum, outputname){
  
  #featurelength$all <-rowSums(featurelength[,2:5])
  
  # 1: uts/leading
  # -1: ts/lagging
  #featurelength_strand <- data.frame("Strand"=c("1","-1"),"C"=c(sum(featurelength[,"C"]), sum(featurelength[,"G"])),"T"=c(sum(featurelength[,"T"]), sum(featurelength[,"A"])))
  featurelength_strand <- data.frame("Ref"=c("C","T"),"uts_leading_wg"=c(sum(featurelength[,"C"]), sum(featurelength[,"T"])), "ts_lagging_wg"=c(sum(featurelength[,"G"]), sum(featurelength[,"A"])))
  
  denovo_muts_feature <- denovo_muts_feature[denovo_muts_feature[,feature] != "others",]
  denovo_muts_feature$MutationType <- paste0(denovo_muts_feature$pre_context,"[",denovo_muts_feature$Ref2,">",denovo_muts_feature$Alt2,"]",denovo_muts_feature$rear_context)

  denovo_muts_dis <- data.frame(table(denovo_muts_feature[,SampleCol],denovo_muts_feature[,feature], denovo_muts_feature[,"MutationType"]))
  names(denovo_muts_dis) <- c("Sample",feature,"MutationType", "Freq")
  denovo_muts_dis$Mutation <- substr(denovo_muts_dis$MutationType,3,5)
  mutation_order <- unique(denovo_muts_dis[,c("MutationType","Mutation")])
  mutation_order <- mutation_order[order(mutation_order$Mutation),]
  gtc <- denovo_muts_dis
  gtc$targetfeature <- gtc[,feature]
  
  filename=paste0(outputname,".pdf")
  pdf(file=filename, onefile=TRUE,height=h,width = w) 
  d1 <- ggplot(gtc,aes(x=MutationType,y=Freq,fill=targetfeature))+geom_bar(stat="identity",position="dodge")
  d1 <- d1+scale_x_discrete(limits = as.character(mutation_order$MutationType))+ggtitle(outputname)
  d1 <- d1+ylab("Count")+theme(axis.title.x = element_text(size=15),
                               axis.title.y = element_text(size=15),
                               plot.title = element_text(size=10),
                               axis.text.x=element_text(angle=90, vjust=0.5,colour = "black",size = 5),
                               axis.text.y=element_text(colour = "black"),
                               panel.grid.minor.x=element_blank(),
                               panel.grid.major.x=element_blank(),
                               panel.grid.major.y = element_blank(),
                               panel.grid.minor.y = element_blank(),
                               panel.background = element_rect(fill = "white"),
                               panel.border = element_rect(colour = "black", fill=NA))
  d1 <- d1+facet_wrap(~Sample,ncol=colnum,scales="free")
  print(d1)
  dev.off()
  
  gtc_dcast <- dcast(gtc,Sample+MutationType~Strand, value.var="Freq")
  gtc_dcast$Ref <- substr(gtc_dcast$MutationType,1,1)
  gtc_dcast <- merge(gtc_dcast, featurelength_strand, by="Ref", all.x=T)
  gtc_dcast$chisq_pvalue <- 1
  for(i in 1:dim(gtc_dcast)[1]){
    if(gtc_dcast[i,"leading_uts"]+gtc_dcast[i,"lagging_ts"]>0){
      gtc_dcast[i,"chisq_pvalue"] <- chisq.test(c(gtc_dcast[i,"leading_uts"], gtc_dcast[i,"lagging_ts"]), p=c(gtc_dcast[i,"uts_leading_wg"], gtc_dcast[i,"ts_lagging_wg"]),rescale.p = TRUE)$p.value
      
    }
    
  }
  gtc_dcast$P_adjust <- p.adjust(gtc_dcast$chisq_pvalue,method = "BH")
  gtc_dcast$flag <- ""
  gtc_dcast[gtc_dcast$P_adjust<=0.05,"flag"] <- "*"
  gtc_dcast[gtc_dcast$P_adjust<=0.01,"flag"] <- "**"
  gtc_dcast[gtc_dcast$P_adjust<=0.001,"flag"] <- "***"
  
  write.table(gtc_dcast,paste0(outputname, "_chisq_adjust.txt"),sep = "\t",col.names = T, row.names = F, quote = F)
  
  
}

# odds ratio
StrandBias_oddsratio <- function(denovo_muts_feature, SampleCol,feature,featurelength,h,w,colnum, outputname){
  
  #featurelength$all <-rowSums(featurelength[,2:5])
  
  # 1: uts/leading
  # -1: ts/lagging
  #featurelength_strand <- data.frame("Strand"=c("1","-1"),"C"=c(sum(featurelength[,"C"]), sum(featurelength[,"G"])),"T"=c(sum(featurelength[,"T"]), sum(featurelength[,"A"])))
  featurelength_strand <- data.frame("Ref"=c("C","T"),"uts_leading_wg"=c(sum(featurelength[,"C"]), sum(featurelength[,"T"])), "ts_lagging_wg"=c(sum(featurelength[,"G"]), sum(featurelength[,"A"])))
  
  denovo_muts_feature <- denovo_muts_feature[denovo_muts_feature[,feature] != "others",]
  denovo_muts_feature$Mutation <- paste0(denovo_muts_feature$Ref2,">",denovo_muts_feature$Alt2)
  denovo_muts_dis <- data.frame(table(denovo_muts_feature[,SampleCol],denovo_muts_feature[,feature], denovo_muts_feature[,"Mutation"]))
  names(denovo_muts_dis) <- c("Sample",feature,"Mutation", "Freq")
  
  gtc <- denovo_muts_dis
  gtc$targetfeature <- gtc[,feature]
  
  gtc_dcast <- dcast(gtc,Sample+Mutation~Strand, value.var="Freq")
  gtc_dcast$Ref <- substr(gtc_dcast$Mutation,1,1)
  gtc_dcast <- merge(gtc_dcast, featurelength_strand, by="Ref", all.x=T)
  gtc_dcast$oddsratio <- 1
  gtc_dcast$lowerconfint <- 1
  gtc_dcast$higherconfint <- 1
  for(i in 1:dim(gtc_dcast)[1]){
    if(gtc_dcast[i,"leading_uts"]*gtc_dcast[i,"lagging_ts"]>0){
      
      M <- matrix(c(gtc_dcast[i,"lagging_ts"], gtc_dcast[i,"leading_uts"], (gtc_dcast[i,"ts_lagging_wg"]-gtc_dcast[i,"lagging_ts"]), (gtc_dcast[i,"uts_leading_wg"]-gtc_dcast[i,"leading_uts"])), ncol = 2)
      b <- data.frame(t(calcOddsRatio(M)))
      names(b) <- c("oddsratio","lowerconfint","higherconfint")
      gtc_dcast[i,"oddsratio"] <- b$oddsratio
      gtc_dcast[i,"lowerconfint"] <- b$lowerconfint
      gtc_dcast[i,"higherconfint"] <- b$higherconfint
      #    gtc_dcast[i,"chisq_pvalue"] <- chisq.test(c(gtc_dcast[i,"leading_uts"], gtc_dcast[i,"lagging_ts"]), p=c(gtc_dcast[i,"uts_leading_wg"], gtc_dcast[i,"ts_lagging_wg"]),rescale.p = TRUE)$p.value
      
    }
    
  }
  gtc_dcast$flag <- ""
  gtc_dcast[gtc_dcast$lowerconfint>1 | gtc_dcast$higherconfint<1,"flag"] <- "*"
  write.table(gtc_dcast,paste0(outputname, "_oddsratio.txt"),sep = "\t",col.names = T, row.names = F, quote = F)
  
  filename=paste0(outputname,"_oddsratio.pdf")
  pdf(file=filename, onefile=TRUE,height=h,width = w) 
  d1 <- ggplot(gtc_dcast,aes(x=oddsratio,y=Mutation))+geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") 
  d1 <- d1 + geom_errorbarh(aes(xmax = higherconfint, xmin = lowerconfint), size = .5, height = .2, color = "gray50")
  d1 <- d1 + geom_point(size = 3.5, color = "orange") +ylab("") +xlab("Odds Ratio")
  #  d1 <- d1 + coord_trans(x = scales:::exp_trans(10)) + scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1),
  #                   limits = log10(c(0.09,2.5)))
  d1 <- d1+theme(axis.title.x = element_text(size=15),
                 axis.title.y = element_text(size=15),
                 plot.title = element_text(size=10),
                 axis.text.x=element_text(angle=90, vjust=0.5,colour = "black"),
                 axis.text.y=element_text(colour = "black"),
                 panel.grid.minor.x=element_blank(),
                 panel.grid.major.x=element_blank(),
                 panel.grid.major.y = element_blank(),
                 panel.background = element_rect(fill = "white"),
                 panel.border = element_rect(colour = "black", fill=NA))
  d1 <- d1+facet_wrap(~Sample,ncol=colnum,scales="free")
  print(d1)
  dev.off()
  
  
  
}

calcOddsRatio <- function(mymatrix,alpha=0.05,referencerow=2,quiet=FALSE)
{
  numrow <- nrow(mymatrix)
  myrownames <- rownames(mymatrix)
  
  for (i in 1:numrow)
  {
    rowname <- myrownames[i]
    DiseaseUnexposed <- mymatrix[referencerow,1]
    ControlUnexposed <- mymatrix[referencerow,2]
    if (i != referencerow)
    {
      DiseaseExposed <- mymatrix[i,1]
      ControlExposed <- mymatrix[i,2]
      
      totExposed <- DiseaseExposed + ControlExposed
      totUnexposed <- DiseaseUnexposed + ControlUnexposed
      
      probDiseaseGivenExposed <- DiseaseExposed/totExposed
      probDiseaseGivenUnexposed <- DiseaseUnexposed/totUnexposed
      probControlGivenExposed <- ControlExposed/totExposed
      probControlGivenUnexposed <- ControlUnexposed/totUnexposed
      
      # calculate the odds ratio
      oddsRatio <- (probDiseaseGivenExposed*probControlGivenUnexposed)/
        (probControlGivenExposed*probDiseaseGivenUnexposed)
      if (quiet == FALSE)
      {
        print(paste("category =", rowname, ", odds ratio = ",oddsRatio))
      }
      
      # calculate a confidence interval
      confidenceLevel <- (1 - alpha)*100
      sigma <- sqrt((1/DiseaseExposed)+(1/ControlExposed)+
                      (1/DiseaseUnexposed)+(1/ControlUnexposed))
      # sigma is the standard error of our estimate of the log of the odds ratio
      z <- qnorm(1-(alpha/2))
      lowervalue <- oddsRatio * exp(-z * sigma)
      uppervalue <- oddsRatio * exp( z * sigma)
      if (quiet == FALSE)
      {
        print(paste("category =", rowname, ", ", confidenceLevel,
                    "% confidence interval = [",lowervalue,",",uppervalue,"]"))
        return(c(oddsRatio,lowervalue,uppervalue))
      }
    }
  }
  if (quiet == TRUE && numrow == 2) # If there are just two treatments (exposed/nonexposed)
  {
    return(oddsRatio)
  }
}


RemoveBackground <- function(background_profile, sig_profile,chan_num=96,sampling_number, start_num,boundary){
  
  # Remove weak mutation types in sig_profile
  removeWeakMutationTypes <- 0.01
  genomesOriginal <- as.data.frame(sig_profile)
  Totalmutations <- sum(sig_profile)
  removeMutations_max <- removeWeakMutationTypes * Totalmutations
  removeIdx <- which(cumsum(rowSums(genomesOriginal[,1:dim(genomesOriginal)[2],drop = FALSE])[order(rowSums(genomesOriginal[,1:dim(genomesOriginal)[2],drop = FALSE]))])<= min(removeMutations_max,4))
  mutationTypesToRemoveSet <- order(rowSums(genomesOriginal[,1:dim(genomesOriginal)[2],drop = FALSE]))[removeIdx]
  genomesReducted <- as.data.frame(genomesOriginal[-mutationTypesToRemoveSet,])
  reducedMutationtypes <- dim(genomesReducted)[1]
  
  # bootstrap sig profile
  centroid_sig <- rowMeans(genomesReducted[,1:dim(genomesReducted)[2],drop = FALSE])
  RepSig <- matrix(rep(centroid_sig,sampling_number),ncol = sampling_number)
  Sig_bootstraps <- bootstrapGenomesfun2(RepSig,sum(centroid_sig))
  
  i=start_num
  reachLimit <- FALSE
  
  while(i<sum(centroid_sig) & !reachLimit){
    # Remove the same weak mutation types in background profile
    backgroundReducted <- background_profile[-mutationTypesToRemoveSet]
    RepControl <- matrix(rep(backgroundReducted,sampling_number),ncol = sampling_number)
    bg_bootstraps <- bootstrapGenomesfun2(RepControl,i)
    
    # Range of background
    centroid_background <- rowMeans(bg_bootstraps)
    sd_background <- apply(bg_bootstraps,1,sd)
    boundary_background <- centroid_background+boundary*sd_background
    
    # Range of sig
    centroid_sig <- rowMeans(Sig_bootstraps)
    sd_sig <- apply(Sig_bootstraps,1,sd)
    boundary_sig <- centroid_sig+boundary*sd_sig
    
    
    diff_all_boundary <- boundary_sig-centroid_background
    diff_all <- centroid_sig-centroid_background
    
    if(length(which(diff_all_boundary<0))>0){
      reachLimit <- TRUE
    }
    
    if(length(which(diff_all_boundary<0))==0){
      diff_all_boundary_save <- diff_all_boundary
      diff_all_save <- diff_all
      diff_all_save[which(diff_all_save<0)] <- 0
    }
    
    
    
    i = i+1
    
  }
  # Add Weak mutations for KO expoure
  exposure <- rep(0,chan_num)
  origArrayIndex <- 1
  for(i in 1:chan_num){
    if(! i %in% mutationTypesToRemoveSet){
      exposure[i] <- diff_all_save[origArrayIndex]
      origArrayIndex=origArrayIndex+1
    }
  }
  
  # Add Weak mutations for background
  background_exposure<- rep(0,chan_num)
  origArrayIndex <- 1
  for(i in 1:chan_num){
    if(! i %in% mutationTypesToRemoveSet){
      background_exposure[i] <- centroid_background[origArrayIndex]
      origArrayIndex=origArrayIndex+1
    }
  }
  
  
  return(data.frame("KO_exposure"=exposure,"background_exposure"=background_exposure))
  
  
}

AssessEnrichment <- function(mutlist, featurefile, rangelength,binsize,h,w,outputname){
  
  featurelist <- read.table(featurefile,sep = "\t",header = F, as.is = T)
  featurelist[,2] <- featurelist[,2]-rangelength
  featurelist[,3] <- featurelist[,3]+rangelength
  write.table(featurelist,"newfeaturefile.txt",sep = "\t", col.names = F, row.names = F, quote = F)
  Rintersect(mutlist,"muts_bed.txt","newfeaturefile.txt","muts_feature_flank")
  muts_feature_flank1000 <- read.table("muts_feature_flank.txt", sep = "\t", header = F, as.is = T)
  
  # The distance between mutation to the center (dc)
  muts_feature_flank1000$dc <- muts_feature_flank1000[,2]-floor((muts_feature_flank1000[,6]+muts_feature_flank1000[,7])/2)
  filename=paste0(outputname,".pdf")
  pdf(file=filename, onefile=TRUE,width=8,height=3) 
  dd <- ggplot(muts_feature_flank1000,aes(x=dc))+geom_histogram(binwidth = binsize)
  print(dd)
  dev.off()
  
}

#########################################
#
# Signature comparison
#
#########################################


Cossimi_CompareSig <- function(target_sig,cosmic_sig,h,w,text_size,x_text,outputname){
  
  muttype <- colnames(cosmic_sig)[1]
  sig_all <- merge(cosmic_sig, target_sig, by=muttype)
  cosmic_sig_new <- sig_all[,2:dim(cosmic_sig)[2]]
  target_sig_new <- sig_all[,(dim(cosmic_sig)[2]+1):dim(sig_all)[2]]
  simi_matrix <- NULL
  for(i in 1:dim(cosmic_sig_new)[2]){
    print(i)
    simi_m <- apply(target_sig_new,2,function(x) abs(cos_similarity(cosmic_sig_new[,i],x)))
    simi_matrix <- rbind(simi_matrix, simi_m)
  }
  simi_matrix <- data.frame(simi_matrix)
  simi_matrix$CosmicSig <- colnames(cosmic_sig_new)
  simi_matrix_melt <- melt(simi_matrix)
  names(simi_matrix_melt) <- c("CosmicSig","Sample","similarity")
  simi_matrix_melt$similarity <- round(simi_matrix_melt$similarity,2)
  cosmig_pos <- simi_matrix$CosmicSig
  cosmig_label <- c(1:30)
  pdf(file=paste0(outputname,".pdf"), onefile=TRUE,width = w,height = h)
  g <-ggplot(simi_matrix_melt, aes(y=factor(Sample), x=factor(CosmicSig))) + geom_tile(aes(fill=similarity),colour="white")
  g <- g+xlab(x_text)+ylab("Samples")+scale_fill_gradient2(low ="blue",mid = "white",high="red",space="Lab")
  g <- g+scale_x_discrete(limits = cosmig_pos, labels = cosmig_label)
  g <- g+geom_text(aes(label=paste(similarity)),size=text_size)
  g <- g+theme(axis.text.x=element_text(size=10, vjust=0.5,colour="black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),                                                               
               axis.title.y = element_text(size=15),                                                               
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA)
  )
  
  print(g)
  dev.off()
  write.table(simi_matrix,paste0(outputname,".txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  
}

# any two signature sets
Cossimi_CompareSig2 <- function(target_sig,cosmic_sig,h,w,text_size,x_text,outputname){
  
  muttype <- colnames(cosmic_sig)[1]
  sig_all <- merge(cosmic_sig, target_sig, by=muttype)
  cosmic_sig_new <- data.frame(sig_all[,2:dim(cosmic_sig)[2]])
  
  target_sig_new <- sig_all[,(dim(cosmic_sig)[2]+1):dim(sig_all)[2]]
  simi_matrix <- NULL
  for(i in 1:dim(cosmic_sig_new)[2]){
    print(i)
    simi_m <- apply(target_sig_new,2,function(x) abs(cos_similarity(cosmic_sig_new[,i],x)))
    simi_matrix <- rbind(simi_matrix, simi_m)
  }
  simi_matrix <- data.frame(simi_matrix)
  simi_matrix$CosmicSig <- colnames(cosmic_sig)[-1]
  simi_matrix_melt <- melt(simi_matrix,c("CosmicSig"))
  names(simi_matrix_melt) <- c("CosmicSig","Sample","similarity")
  simi_matrix_melt$similarity <- round(simi_matrix_melt$similarity,2)
  cosmig_pos <- simi_matrix$CosmicSig
  pdf(file=paste0(outputname,".pdf"), onefile=TRUE,width = w,height = h)
  g <-ggplot(simi_matrix_melt, aes(y=factor(Sample), x=factor(CosmicSig))) + geom_tile(aes(fill=similarity),colour="white")
  g <- g+xlab(x_text)+ylab("Samples")+scale_fill_gradient2(high="red", low="white",limits=c(0.5, 1),midpoint=0.5) #scale_fill_gradient2(low ="blue",mid = "white",high="red",space="Lab")
  g <- g+scale_x_discrete(limits = cosmig_pos)
  g <- g+geom_text(aes(label=paste(similarity)),size=text_size)
  g <- g+theme(axis.text.x=element_text(size=10, vjust=0.5,colour="black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),                                                               
               axis.title.y = element_text(size=15),                                                               
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA)
  )
  
  print(g)
  dev.off()
  write.table(simi_matrix,paste0(outputname,".txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  return(simi_matrix_melt)
}
Cossimi_CompareSig2_noplot <- function(target_sig,cosmic_sig,outputname){
  
  muttype <- colnames(cosmic_sig)[1]
  sig_all <- merge(cosmic_sig, target_sig, by=muttype)
  cosmic_sig_new <- data.frame(sig_all[,2:dim(cosmic_sig)[2]])
  
  target_sig_new <- sig_all[,(dim(cosmic_sig)[2]+1):dim(sig_all)[2]]
  simi_matrix <- NULL
  for(i in 1:dim(cosmic_sig_new)[2]){
    print(i)
    simi_m <- apply(target_sig_new,2,function(x) abs(cos_similarity(cosmic_sig_new[,i],x)))
    simi_matrix <- rbind(simi_matrix, simi_m)
  }
  simi_matrix <- data.frame(simi_matrix)
  simi_matrix$CosmicSig <- colnames(cosmic_sig)[-1]
  simi_matrix_melt <- melt(simi_matrix,c("CosmicSig"))
  names(simi_matrix_melt) <- c("CosmicSig","Sample","similarity")
  simi_matrix_melt$similarity <- round(simi_matrix_melt$similarity,2)
  cosmig_pos <- simi_matrix$CosmicSig
  simi_matrix_dcast <- dcast(simi_matrix_melt,Sample~CosmicSig)
  write.table(simi_matrix_dcast,paste0(outputname,".txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  return(simi_matrix_melt)
}


# compare signature within a data set
Cossimi_CompareSig3 <- function(target_sig,h,w,text_size,x_text,outputname){
  
  cossimil <- as.matrix(proxy::simil(as.matrix(target_sig), diag=FALSE, upper=FALSE, method="cosine",auto_convert_data_frames = FALSE))
  cossimil[lower.tri(cossimil,diag=TRUE)]=NA
  cossimil <- as.data.frame(as.table(cossimil))
  cossimil=na.omit(cossimil)
  names(cossimil) <- c("sample1","sample2","simil")
  cossimil$simil <- round(cossimil$simil,2)
  write.table(cossimil,paste0(outputname,".txt"),sep="\t",col.names = T, row.names = F, quote = F)
  filename=paste0(outputname,".pdf")
  pdf(file=filename, onefile=TRUE,width = w,height =h)
  g1 <-ggplot(cossimil, aes(x=sample1, y=sample2)) + geom_tile(aes(fill=simil),colour="white")+geom_text(aes(label=paste(simil)),size=text_size)
  g1 <-g1 +scale_fill_gradient2(high="red", low="white",space="Lab",limits=c(0, 1))
  #g1 <-g1 +scale_x_discrete(limits = unique(as.character(cossimil$sample1)))
  #g1 <-g1 +scale_y_discrete(limits = unique(as.character(cossimil$sample2[order(as.character(cossimil$sample2))])))
  
  g1 <-g1 +theme(axis.text.x=element_text(size=x_text,colour = "black"),
                 axis.text.y=element_text(size=10,colour = "black"),
                 plot.title = element_text(size=10),
                 panel.grid.minor.x=element_blank(),
                 panel.grid.major.x=element_blank(),
                 panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank(),
                 panel.background = element_rect(fill = "white"),
                 panel.border = element_rect(colour = "black", fill=NA))
  print(g1)
  dev.off()
  return(cossimil)
}





# compare two signatures,return cosine similarity
Cossimi_CompareSig4 <- function(target_sig,cosmic_sig){
  
  muttype <- colnames(cosmic_sig)[1]
  sig_all <- merge(cosmic_sig, target_sig, by=muttype)
  return(abs(cos_similarity(sig_all[,2],sig_all[,3])))
}

# compare one signature with another signature set,return cosine similarity
Cossimi_CompareSig5 <- function(target_sig,cosmic_sig){
  
  muttype <- colnames(cosmic_sig)[1]
  sig_all <- merge(cosmic_sig, target_sig, by=muttype)
  cosmic_sig_new <- data.frame(sig_all[,2:dim(cosmic_sig)[2]])
  
  target_sig_new <- sig_all[,(dim(cosmic_sig)[2]+1):dim(sig_all)[2]]
  
  simi_m <- apply(cosmic_sig_new,2,function(x) abs(cos_similarity(target_sig_new,x)))
  
  
  return(simi_m)
}

#########################################
#
# Double substitutions
#
#########################################

FindDinucleotides <- function(allsubs) {
  dinuclist <- NULL
  samplelist <- data.frame(table(allsubs$Sample))
  for(i in 1:dim(samplelist)[1]){
    print(i)
    samplesubs <- allsubs[allsubs$Sample==samplelist[i,1],]
    a <- samplesubs$Pos
    b <- samplesubs$Ref
    c <- samplesubs$Alt
    samplesubs$Pos_neighbor <- c(a[-1],0)
    samplesubs$neigbor_dist <- samplesubs$Pos-samplesubs$Pos_neighbor
    samplesubs$Ref_neighbor <- c(b[-1],"N")
    samplesubs$Alt_neighbor <- c(c[-1],"N")
    
    dinuc_index <- which(samplesubs$neigbor_dist==-1)
    print(length(dinuc_index))
    if(length(dinuc_index)>0){
      dinuclist <- rbind(dinuclist,samplesubs[dinuc_index,])
      dinuclist <- rbind(dinuclist,samplesubs[dinuc_index+1,])
    }
  }
  return(dinuclist)
}

# for a given list of subs of a sample,
FindAdjacentMut_sample <- function(allsubs) {
  dinuclist <- NULL   
  samplesubs <- allsubs[order(allsubs$Chrom, allsubs$Pos),]
  a <- samplesubs$Pos
  
  samplesubs$Pos_neighbor <- c(a[-1],0)
  samplesubs$neigbor_dist <- samplesubs$Pos-samplesubs$Pos_neighbor
  
  dinuc_index <- which(samplesubs$neigbor_dist==-1)
  print(length(dinuc_index))
  if(length(dinuc_index)>0){
    dinuclist <- rbind(dinuclist,samplesubs[dinuc_index,])
    dinuclist <- rbind(dinuclist,samplesubs[dinuc_index+1,])
  }
  return(dinuclist)
}


gen_muttype_Dinucleotides <- function(CTsubs,SelectCol){
  muttype_freq_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/12_subs/Dinucleiotides/Dinucleotides_template.txt", sep = "\t", header = T, as.is = T)
  
  sigfile_freq <- data.frame(table(CTsubs[,SelectCol],CTsubs$dinuc_mutation_final))
  names(sigfile_freq) <- c("SelectPara","MutationType","Freq")
  control_sigset <-dcast(sigfile_freq,MutationType~SelectPara,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="MutationType",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  return(control_sigset)
  
}
plot_allsample_Dinucleotides_profile <- function(muttype_catalogue,colnum,h,w,outputname){
  
  muts_basis_melt <- melt(muttype_catalogue,id=c("MutationType","Ref"))
  names(muts_basis_melt) <- c("MutationType","Ref","sample","count")
  
  mypalette <- c("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4","#46f0f0","#f032e6","#d2f53c","#fabebe")
  
  pdf(file=outputname, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=MutationType, y=count,fill=Ref,width=0.5))+ geom_bar(position="dodge", stat="identity")+xlab("Dinucleotide mutation type")+ylab("")
  #  p <- p+coord_cartesian(ylim=c(0, max(muttype_freq[,"freq"])))
  p <- p+scale_x_discrete(limits = as.character(muttype_catalogue$MutationType))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  p <- p+facet_wrap(~sample,ncol=colnum,scales = "free")
  print(p)
  dev.off()
} 

gen_10Ref_Dinucleotides <- function(CTsubs,SelectCol){
  muttype_freq_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/12_subs/Dinucleiotides/Dinucleotides_10Ref.txt", sep = "\t", header = T, as.is = T)
  
  sigfile_freq <- data.frame(table(CTsubs[,SelectCol],CTsubs$dinuc_Ref_final))
  names(sigfile_freq) <- c("SelectPara","Ref","Freq")
  control_sigset <-dcast(sigfile_freq,Ref~SelectPara,value.var = "Freq")
  control_sigset <- merge(muttype_freq_template,control_sigset,by="Ref",all.x=T)
  control_sigset[is.na(control_sigset)] <- 0
  return(control_sigset)
  
}
plot_allsample_Dinucleotides_10Ref <- function(muttype_catalogue,colnum,h,w,outputname){
  
  muts_basis_melt <- melt(muttype_catalogue,id=c("Ref","Ref2"))
  names(muts_basis_melt) <- c("Ref","Ref2","sample","count")
  
  mypalette <- c("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4","#46f0f0","#f032e6","#d2f53c","#fabebe")
  
  pdf(file=outputname, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=Ref, y=count,fill=Ref,width=0.5))+ geom_bar(position="dodge", stat="identity")+xlab("Dinucleotide mutation type")+ylab("")
  #  p <- p+coord_cartesian(ylim=c(0, max(muttype_freq[,"freq"])))
  p <- p+scale_x_discrete(limits = as.character(muttype_catalogue$Ref))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  p <- p+facet_wrap(~sample,ncol=colnum,scales = "free")
  print(p)
  dev.off()
} 




