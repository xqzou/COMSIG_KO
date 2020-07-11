library(reshape2)
library(ggplot2)
library(plyr)
library("VariantAnnotation")
library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)
library(stringi)
#library("tidyverse")

bootstrapGenomesfun2 <- function(genomes,n){
  
  return(apply(genomes, 2, function(x) rmultinom(1, n, x)))
}
cos_similarity <- function(v1,v2){
  v1v2 <- sum(v1*v2)
  v1_length <- sqrt(sum(v1*v1))
  v2_length <- sqrt(sum(v2*v2))
  return(v1v2/v1_length/v2_length)
}
# compare one signature set (target_sig) with another signature set (cosmic_sig),return cosine similarity
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
  simi_matrix_dcast <- dcast(simi_matrix_melt,Sample~CosmicSig,value.var="similarity")
  write.table(simi_matrix_dcast,paste0(outputname,".txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  return(simi_matrix_melt)
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
  p <- p+scale_x_discrete(limits = as.character(mutation_order$MutationType))+ggtitle(outputname)+scale_y_continuous(limits=c(0, 0.21),breaks=seq(0, 0.2, 0.05),labels=percent)
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
Wrap_KOSig <- function(MutCatalogue,bg_column,ko_column,sampling_number, start_num,boundary,outputname){
  
  KOSig <- RemoveBackground_vector_single(MutCatalogue[,bg_column], MutCatalogue[,ko_column],sampling_number, start_num,boundary)
  KOSig$MutationType <- MutCatalogue[,"MutationType"]
  write.table(KOSig,paste0(outputname,".txt"),sep = "\t",col.names = T, row.names = F, quote = F)
  plotCountbasis(KOSig,1,6,9,paste0(outputname,".pdf"))
  plotPercentagebasis(KOSig,1,6,9,paste0(outputname,"_percentage.pdf"))
  
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
  p <- p+theme(axis.text.x=element_blank(),
               #axis.text.x=element_text(size=5,angle=90,colour = "black"),
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
  CTsubs$pre_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos-1, CTsubs$Pos-1))
  CTsubs$rear_context <- as.character(getSeq(Hsapiens, paste0('chr',CTsubs$Chrom), CTsubs$Pos+1, CTsubs$Pos+1))
  CTsubs_copy <- CTsubs
  
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt2 <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref2 <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  
  CTsubs[(CTsubs$Ref2!=CTsubs$Ref & CTsubs$Strand_original=="leading_uts"),]$Strand <- "lagging_ts"
  CTsubs[(CTsubs$Ref2!=CTsubs$Ref & CTsubs$Strand_original=="lagging_ts"),]$Strand <- "leading_uts"
  CTsubs$Mutation <- paste0(CTsubs$Ref2,">",CTsubs$Alt2)
  CTsubs$MutationType <- paste0(CTsubs$pre_context,"[",CTsubs$Ref2,">",CTsubs$Alt2,"]",CTsubs$rear_context)
  
  
  write.table(CTsubs,outputname,sep="\t",col.names = T, row.names = F, quote = F)
  return(CTsubs)
}
StrandBias_oddsratio <- function(denovo_muts_feature, SampleCol,feature,featurelength,h,w,colnum, outputname){
  
  #featurelength$all <-rowSums(featurelength[,2:5])
  
  # 1: uts/leading
  # -1: ts/lagging
  #featurelength_strand <- data.frame("Strand"=c("1","-1"),"C"=c(sum(featurelength[,"C"]), sum(featurelength[,"G"])),"T"=c(sum(featurelength[,"T"]), sum(featurelength[,"A"])))
  # featurelength_strand <- data.frame("Ref"=c("C","T"),"uts_leading_wg"=c(sum(featurelength[,"C"]), sum(featurelength[,"T"])), "ts_lagging_wg"=c(sum(featurelength[,"G"]), sum(featurelength[,"A"])))
  featurelength_strand <- data.frame("Ref"=c("C","T"),
                                     "uts_leading_wg"=c((featurelength[featurelength$Strand==1,"C"]+featurelength[featurelength$Strand==-1,"G"]), (featurelength[featurelength$Strand==1,"T"]+featurelength[featurelength$Strand==-1,"A"])),
                                     "ts_lagging_wg"=c((featurelength[featurelength$Strand==1,"G"]+featurelength[featurelength$Strand==-1,"C"]), (featurelength[featurelength$Strand==1,"A"]+featurelength[featurelength$Strand==-1,"T"])))
  
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
StrandBias_6muttype <- function(denovo_muts_feature, SampleCol,feature,featurelength,h,w,colnum, outputname){
  
  #featurelength$all <-rowSums(featurelength[,2:5])
  
  # 1: uts/leading/G4
  # -1: ts/lagging/NonG4
  #featurelength_strand <- data.frame("Strand"=c("1","-1"),"C"=c(sum(featurelength[,"C"]), sum(featurelength[,"G"])),"T"=c(sum(featurelength[,"T"]), sum(featurelength[,"A"])))
  #  featurelength_strand <- data.frame("Ref"=c("C","T"),"uts_leading_wg"=c(sum(featurelength[,"C"]), sum(featurelength[,"T"])), "ts_lagging_wg"=c(sum(featurelength[,"G"]), sum(featurelength[,"A"])))
  featurelength_strand <- data.frame("Ref"=c("C","T"),
                                     "uts_leading_wg"=c((featurelength[featurelength$Strand==1,"C"]+featurelength[featurelength$Strand==-1,"G"]), (featurelength[featurelength$Strand==1,"T"]+featurelength[featurelength$Strand==-1,"A"])),
                                     "ts_lagging_wg"=c((featurelength[featurelength$Strand==1,"G"]+featurelength[featurelength$Strand==-1,"C"]), (featurelength[featurelength$Strand==1,"A"]+featurelength[featurelength$Strand==-1,"T"])))
  
  denovo_muts_feature <- denovo_muts_feature[denovo_muts_feature[,feature] != "others",]
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
  #  featurelength_strand <- data.frame("Ref"=c("C","T"),"uts_leading_wg"=c(sum(featurelength[,"C"]), sum(featurelength[,"T"])), "ts_lagging_wg"=c(sum(featurelength[,"G"]), sum(featurelength[,"A"])))
  # 
  
  denovo_muts_feature <- denovo_muts_feature[denovo_muts_feature[,feature] != "others",]
  
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
  gtc_dcast$Mutation <- paste0(substr(gtc_dcast$MutationType,3,3), substr(gtc_dcast$MutationType,4,4),substr(gtc_dcast$MutationType,5,5))
  write.table(gtc_dcast,paste0(outputname, "_freq.txt"),sep = "\t",col.names = T, row.names = F, quote = F)
  
  #names(featurelength) <- c("seq_pyrimidine","plus_strand","minus_strand")
  gtc_dcast$seq_pyrimidine <- paste0(substr(as.character(gtc_dcast$MutationType),1,1),substr(as.character(gtc_dcast$MutationType),3,3),substr(as.character(gtc_dcast$MutationType),7,7))
  
  gtc_dcast <- merge(gtc_dcast, featurelength, by="seq_pyrimidine", all.x=T)
  gtc_dcast$chisq_pvalue <- 1
  for(i in 1:dim(gtc_dcast)[1]){
    if(gtc_dcast[i,"leading_uts"]+gtc_dcast[i,"lagging_ts"]>0){
      gtc_dcast[i,"chisq_pvalue"] <- chisq.test(c(gtc_dcast[i,"leading_uts"], gtc_dcast[i,"lagging_ts"]), p=c(gtc_dcast[i,"leading_uts_wg"], gtc_dcast[i,"lagging_ts_wg"]),rescale.p = TRUE)$p.value
      
    }
    
  }
  gtc_dcast$P_adjust <- p.adjust(gtc_dcast$chisq_pvalue,method = "BH")
  gtc_dcast$flag <- ""
  gtc_dcast[gtc_dcast$P_adjust<=0.05,"flag"] <- "*"
  gtc_dcast[gtc_dcast$P_adjust<=0.01,"flag"] <- "**"
  gtc_dcast[gtc_dcast$P_adjust<=0.001,"flag"] <- "***"
  
  write.table(gtc_dcast,paste0(outputname, "_chisq_adjust.txt"),sep = "\t",col.names = T, row.names = F, quote = F)
  
  
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
Find_RepeatSeq <- function(CTsubs){
  
  #  CTsubs <- Find_SeqContext(CTsubs,10)
  CTsubs$pre_pos <- 10
  CTsubs$pre_repeat <- ""
  CTsubs$rear_repeat <- ""
  CTsubs$pre_pos <- 10
  CTsubs$rear_pos <- 1
  for(i in 1:dim(CTsubs)[1]){
    print(i)
    # 5' pre - max
    base_5 <- max(which(strsplit(CTsubs[i,]$pre_context_py, "")[[1]]!=CTsubs[i,]$Ref_py))
    
    # 3' rear - max
    base_3 <- min(which(strsplit(CTsubs[i,]$rear_context_py, "")[[1]]!=CTsubs[i,]$Ref_py))
    
    # 5' base of the repeat
    CTsubs[i,]$pre_pos <- base_5
    CTsubs[i,]$rear_pos <- base_3
    CTsubs[i,]$pre_repeat <- substr(CTsubs[i,]$pre_context_py,max(which(strsplit(CTsubs[i,]$pre_context_py, "")[[1]]!=CTsubs[i,]$Ref_py)),max(which(strsplit(CTsubs[i,]$pre_context_py, "")[[1]]!=CTsubs[i,]$Ref_py)))
    CTsubs[i,]$rear_repeat <- substr(CTsubs[i,]$rear_context_py,min(which(strsplit(CTsubs[i,]$rear_context_py, "")[[1]]!=CTsubs[i,]$Ref_py)),min(which(strsplit(CTsubs[i,]$rear_context_py, "")[[1]]!=CTsubs[i,]$Ref_py)))
    CTsubs[i,]$repeat_length <- base_3-base_5+10
    
  }
  write.table(CTsubs,"CTsubs_repeat.txt", sep = "\t", col.names = T, row.names = F, quote = F)
  # return(CTsubs)
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
  p <- ggplot(data=muts_basis_melt, aes(x=MutationType, y=count,fill=mutation))+ geom_bar(stat="identity",position="dodge", width=.7)+xlab("Mutation Types")+ylab("Count")
  p <- p+scale_x_discrete(limits = as.character(mutation_order$MutationType))+ggtitle(outputname)
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(#axis.text.x=element_blank(),
    axis.text.x=element_text(size=3,angle=90,colour = "black",vjust=0.5,hjust=0.9),
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


###########################
#  MMRDetect
##########################

MMRDetect.compute.variables <- function(subs, indels, tissue_type,MMR_sig_indel, tissue_subsig96, conversion_mtx){
  
  
  sub_summary <- data.frame(table(subs$Sample))
  names(sub_summary) <- c("Sample","sub_num")
  indel_classied <- indel_classifier(indels)
  
  indel_classied_rep <- indel_classied[!indel_classied$indeltype_short%in%c("[-]Mh","[->1]Other","[+]Mh","[+>1]Other","Complex","[-C]Rep=1","[-T]Rep=1","[+C]Rep=0","[+T]Rep=0", "[+C]Rep=1","[+T]Rep=1"),]
  indel_classied_rep_summary <- data.frame(table(indel_classied_rep$Sample))
  names(indel_classied_rep_summary) <- c("Sample","RepIndel_num")
  
  muts_summary <- merge(sub_summary,indel_classied_rep_summary,by="Sample")
  
  write.table(muts_summary,paste0("muts_summary_",tissue_type,".txt"),sep = "\t", col.names = T, row.names = F, quote = F)
  #sample_highburden <- muts_summary[muts_summary$sub_num>=190 & muts_summary$RepIndel_num>=80,"Sample"]
  
  
  # Generate catalouge for subs_highburden
  sub_catalouge <- GenCatalogue(subs,"Sample")
  sub_catalouge <- sub_catalouge[,-2]
  
  # Step 2:  Signature fitting for subs using Andrea's tissue-specific signatures 
  selected_tissueSig <- tissue_subsig96[,c(1,which(sub("_[^_]+$","",names(tissue_subsig96))==tissue_type))]
  
  mut_sig <- merge(selected_tissueSig,sub_catalouge,by="MutationType")
  mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
  mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
  sig_cat <- mut_sig[,2:dim(selected_tissueSig)[2]]
  mut_cat <- mut_sig[,(dim(selected_tissueSig)[2]+1):dim(mut_sig)[2]]
  
  a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
  
  tissue_exposure <- t(a$E_median_filtered)
  # conversion matrix
  conversion_mtx_tissue <- conversion_mtx[rownames(conversion_mtx)%in% names(selected_tissueSig),] 
  conversion_mtx_tissue <- conversion_mtx_tissue[colnames(tissue_exposure),]
  
  # convert tissue-specific signatures to reference sigs
  refsig_exposure <- as.matrix(tissue_exposure) %*% as.matrix(conversion_mtx_tissue)
  
  write.table(refsig_exposure,paste0("exposure_",tissue_type,".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
  
  # calculate cosine similarity of indel profile and  MLH1/MSH2/MSH6 PMS2 signatures 
  Sample_MMR <- refsig_exposure[,c("RefSig.MMR1","RefSig.MMR2")]
  
  cossim_allsample <- NULL
  cossim_MMR <- NULL
  cossim_sample <- NULL
  for(i in 2:dim(MMR_sig_indel)[2]){
    for(j in 1:dim(Sample_MMR)[1]){
      current_sample_indel <- indel_classied[indel_classied$Sample==rownames(Sample_MMR)[j],]
      if(dim(current_sample_indel)[1]>1){
        cossim <- Generate_CossimVector_SingleSample_RepIndel2(current_sample_indel,MMR_sig_indel[,c(1,i)])    
        cossim_allsample <- rbind(cossim_allsample,cossim)
        cossim_MMR <- c(cossim_MMR,names(MMR_sig_indel)[i])
        cossim_sample <- c(cossim_sample,rownames(Sample_MMR)[j])
        
      }
      
    }
    
  }
  cossim_allsample <- data.frame(cossim_allsample)
  names(cossim_allsample) <- c("Del_rep","Ins_rep")
  cossim_allsample$Sample <- cossim_sample
  cossim_allsample$MMRgene <- cossim_MMR
  
  Del_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Del_rep_mean=mean(Del_rep),Del_rep_sd=sd(Del_rep))
  Ins_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Ins_rep_mean=mean(Ins_rep),Ins_rep_sd=sd(Ins_rep))
  
  
  MMRsig_2 <- merge(Del_rep_mean[,-2],Ins_rep_mean[,-2],by="Sample")
  MMRsig_2 <- merge(MMRsig_2,indel_classied_rep_summary,by="Sample")
  
  Sample_MMR <- as.data.frame(Sample_MMR)
  Sample_MMR$Sample <- rownames(Sample_MMR)
  Sample_MMR$MMR_sum <- Sample_MMR$RefSig.MMR1+Sample_MMR$RefSig.MMR2
  MMRsig_2 <- merge(MMRsig_2,Sample_MMR,by="Sample")
  write.table(MMRsig_2,paste0("sample_feature_summary_",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)
  
  
  return(MMRsig_2) 
}
