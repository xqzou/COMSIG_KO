library(ggplot2)
library(plyr) 
library(reshape2)
library(gridExtra)
library(scales)
library(VennDiagram) # use for up to 3 sets
library(Rtsne)
library(factoextra)
library("VariantAnnotation")
library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)
SampleMutNum <- function(total_muts, h, w, outputname){
  mut_number <- data.frame(table(total_muts$Sample))
  names(mut_number) <- c("Sample","Freq")
  # Write it into a file
  write.table(mut_number,paste0(outputname,".txt"),sep = "\t", col.names = T, row.names = F, quote = F)
  # Plot it to pdf file
  pdf(file=paste0(outputname,".pdf"), onefile=TRUE,height=h,width=w, useDingbats=FALSE)
  q <- ggplot(mut_number,aes(x=Sample,y=Freq,fill=Sample))+geom_bar(stat="identity", position=position_dodge())
  #q <- q+scale_fill_manual(values=c("orange","forestgreen","pink","royalblue","yellow"))
  q <- q+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=6,colour = "black",hjust = 1),
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
  q <- q+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=6,colour = "black",hjust = 1),
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

quick_barplot <- function(inputdata,xcol,ycol,h,w,outputfilename){
  
  pdf(file=outputfilename, onefile=TRUE,height=h,width=w, useDingbats=FALSE)
  p <- ggplot(inputdata,aes(x=inputdata[,xcol],y=inputdata[,ycol]))+geom_bar(stat="identity",position="dodge", width=.6,fill="blue")
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



##########################
#
#  Indel classification
#
##########################
prepare.indel.df.tab <- function(indel.data) {
  
  if (nrow(indel.data)>0) {
    
    
    indel.data$ref.length <- nchar(indel.data$Ref)
    indel.data$alt.length <- nchar(indel.data$Alt)
    indel.data$indel.length <- abs(indel.data$ref.length - indel.data$alt.length)
    
    #
    indel.data$Type <- NULL
    indel.data[indel.data$ref.length==1,"Type"] <- "Ins"
    indel.data[indel.data$alt.length==1,"Type"] <- "Del"
    indel.data[indel.data$alt.length!=1 & indel.data$ref.length!=1,"Type"] <- "Complex"
    
    # sequence of change
    indel.data$change <- NULL
    indel.data[indel.data$Type=="Complex","change"] <- substr( as.character(indel.data[indel.data$Type=='Complex',"Ref"]),1,1e5)
    indel.data[indel.data$Type=="Ins","change"] <- substr( as.character(indel.data[indel.data$Type=='Ins',"Alt"]),2,1e5)
    indel.data[indel.data$Type=="Del","change"] <- substr( as.character(indel.data[indel.data$Type=='Del',"Ref"]),2,1e5)
    
    # 5'-neighbor base
    indel.data$pre_context <- substr( as.character(indel.data[,"Ref"]),1,1)
    
    # 27bp before and after change                     
    indel.data$extend5 = indel.data$Pos-indel.data$indel.length-25;
    indel.data$extend3 = indel.data$Pos+indel.data$indel.length + indel.data$indel.length+25;
    
    
    indel.data$slice5 <- as.character(getSeq(Hsapiens, paste0('chr',indel.data$Chrom), indel.data$extend5, indel.data$Pos))
    indel.data$slice3 <- NULL
    indel.data[indel.data$Type=="Del","slice3"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Del","Chrom"]), indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1, indel.data[indel.data$Type=="Del","extend3"]))
    indel.data[indel.data$Type=="Ins","slice3"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Ins","Chrom"]), indel.data[indel.data$Type=="Ins","Pos"]+1, indel.data[indel.data$Type=="Ins","extend3"]))
    
    # 1bp before and after change                     
    indel.data$slice5_1bp <- as.character(getSeq(Hsapiens, paste0('chr',indel.data$Chrom), indel.data$Pos, indel.data$Pos))
    indel.data$slice3_1bp <- NULL
    indel.data[indel.data$Type=="Del","slice3_1bp"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Del","Chrom"]), indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1, indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1))
    indel.data[indel.data$Type=="Ins","slice3_1bp"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Ins","Chrom"]), indel.data[indel.data$Type=="Ins","Pos"]+1, indel.data[indel.data$Type=="Ins","Pos"]+1))
    
    # 2bp before and after change                     
    
    #    indel.data$slice5_2bp <- as.character(getSeq(Hsapiens, paste0('chr',indel.data$Chrom), indel.data$Pos-1, indel.data$Pos))
    #    indel.data$slice3_2bp <- NULL
    #    indel.data[indel.data$Type=="Del","slice3_2bp"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Del","Chrom"]), indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1, indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+2))
    #    indel.data[indel.data$Type=="Ins","slice3_2bp"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Ins","Chrom"]), indel.data[indel.data$Type=="Ins","Pos"]+1, indel.data[indel.data$Type=="Ins","Pos"]+2))
    
    
    # 3bp before and after change                     
    
    #    indel.data$slice5_3bp <- as.character(getSeq(Hsapiens, paste0('chr',indel.data$Chrom), indel.data$Pos-2, indel.data$Pos))
    #    indel.data$slice3_3bp <- NULL
    #    indel.data[indel.data$Type=="Del","slice3_3bp"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Del","Chrom"]), indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1, indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+3))
    #    indel.data[indel.data$Type=="Ins","slice3_3bp"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Ins","Chrom"]), indel.data[indel.data$Type=="Ins","Pos"]+1, indel.data[indel.data$Type=="Ins","Pos"]+3))
    
    # Pyrimidine represnetation for 1bp indels
    indel.data$change.pyr <- indel.data$change
    indel.data$slice3_1bp_pyr <- indel.data$slice3_1bp
    indel.data$slice5_1bp_pyr <- indel.data$slice5_1bp
    
    #    indel.data$slice3_2bp_pyr <- indel.data$slice3_2bp
    #    indel.data$slice5_2bp_pyr <- indel.data$slice5_2bp
    
    indel.data[indel.data$change=="A","change.pyr"] <- "T"
    indel.data[indel.data$change=="A","slice5_1bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="A","slice3_1bp"])))
    indel.data[indel.data$change=="A","slice3_1bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="A","slice5_1bp"])))
    
    #    indel.data[indel.data$change=="A","slice5_2bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="A","slice3_2bp"])))
    #    indel.data[indel.data$change=="A","slice3_2bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="A","slice5_2bp"])))
    
    
    indel.data[indel.data$change=="G","change.pyr"] <- "C"
    indel.data[indel.data$change=="G","slice5_1bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="G","slice3_1bp"])))
    indel.data[indel.data$change=="G","slice3_1bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="G","slice5_1bp"])))
    
    #    indel.data[indel.data$change=="G","slice5_2bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="G","slice3_2bp"])))
    #    indel.data[indel.data$change=="G","slice3_2bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="G","slice5_2bp"])))
    
    
    # indel.df needs following columns:
    # indel.type
    # change
    # slice3
    # slice5
    # indel.length
    
    return(indel.data)
    
    
  }
}
# micro-homology between the deleted sequence and 3' sequence (number of bases)
mhcaller <- function(d, prime3) {
  # d is the change of a deletion
  # prime 3 is the sequence 3' of the deleted segment
  countmh = 0
  seq = '-'
  # Dealing with microhomology or lack of microhmomology in the first position
  
  
  for (i in 1:(nchar(d))) { # for all substrings
    
    if (substr(d, 1, i)==substr(prime3,1,i)) {
      countmh <- countmh + 1
      seq <- substr(d, 1, i)
    }
  }
  result <- list()
  result$countmh <- countmh # number of basepairs matching between the deleted fragment and the 3' segment
  result$seq <- seq
  result
}


# counts how many times the pattern is repeated within the string
tandemcount <- function(pat,string) {
  sum = 0
  
  start.pos <- seq(from=1, to=nchar(string)-nchar(pat)+1, by=nchar(pat))
  
  for (s in start.pos) {
    if (substr(string, s, s+nchar(pat)-1) == pat) {
      sum <- sum + 1
    }
  }
  sum
  
}

# tandemcount_v2, by Xueqing, only counts continuous repeated patterns
tandemcount_v2 <- function(pat,string) {
  sum = 0
  
  start.pos <- seq(from=1, to=nchar(string)-nchar(pat)+1, by=nchar(pat))
  s = 1
  
  while (substr(string, s, s+nchar(pat)-1) == pat) {
    sum <- sum + 1
    s = s + nchar(pat)
  }
  
  sum
  
}


# get all factors of an integer
get_all_factors <- function(n)
{
  prime_factor_tables <- lapply(
    setNames(n, n), 
    function(i)
    {
      if(i == 1) return(data.frame(x = 1L, freq = 1L))
      plyr::count(as.integer(gmp::factorize(i)))
    }
  )
  lapply(
    prime_factor_tables, 
    function(pft)
    {
      powers <- plyr::alply(pft, 1, function(row) row$x ^ seq.int(0L, row$freq))
      power_grid <- do.call(expand.grid, powers)
      sort(unique(apply(power_grid, 1, prod)))
    }
  )
}


# ----------------------------------------------------------
# This finds the smallest repeating subunit of the deletion
# ----------------------------------------------------------
findsmallestrep <- function(d) {
  d.factors <- as.numeric(sort(unlist(get_all_factors(nchar(d))), decreasing=TRUE))
  rep.unit <- ''
  
  for (f in d.factors) {
    no.repeats <- nchar(d)/f
    rep.string <- paste0(rep(substring(d,1,f), no.repeats), collapse='')
    if (d==rep.string) {
      rep.unit <- substring(d,1,f)
    }
  }
  rep.unit
}


# ------------------------------------------------ 
# This is the Repeat caller
# ------------------------------------------------      

# Taken in the deletion sequence, the 3' context and the length of the deletion
# I have joined deletion and 3' context. To remove this, the rules for microhom/rep calling need to change
# 
repcaller <- function(d, prime3, prime5, l) {
  # d : deletion
  # prime3 : 3' context
  # prime5: 5' context
  # l : length of change
  
  result <- list() 
  result$countrep = 0
  result$rep<-''
  # This is for counting single base repeats
  # if the length of change is 1
  if (l==1) {
    countrep <- 0
    i <- 1
    while (substr(d,1,1)==substr(prime3,i,i)) {
      countrep <- countrep +1
      i <- i + 1
    }
    result$countrep <- countrep
    result$rep <- d
    return(result)
  } else if (d==substr(prime3, 1, nchar(d))) { # This is for counting whole deletion/DI repeats that are present in the flanking region  
    
    countrep = tandemcount_v2(d,prime3)
    rep = findsmallestrep(d)
    countrep = max(countrep,tandemcount_v2(rep,prime3))
    result$countrep <- countrep
    result$rep <- rep
    return(result)
  } else {   # This is for counting anything in between 1bp and whole del repetition                                     
    rep = '-'
    
    for (t in seq(from=(nchar(d)-1), to=2)) {  # Look for repeats of 2bp to n-1 bp of the length of the indel
      if (grepl(substring(d,1,t), prime3)) {
        countrep = tandemcount_v2(substr(d, 1, t),prime3)
        rep = findsmallestrep(substr(d, 1, t))
        unit = tandemcount_v2(rep,d)*nchar(rep)
        # The false calls arise in examples such as these : del = AACCCCATCTCTACT; 3' = AAAATTACAAACAAAT; rep = 'A'; repcount in 3' = 4 which is greater than MH = 2; Therefore, it is called Repeat-mediated
        # In fact, it should be repeat count = 0; Therefore, call should be microhomology mediated. 
        # To do this, compare check how far the repeat stretched into the indel itself. Eg: 'A' is counted twice in the deletion. So compare del[:2] to del[2:4]. If they are the same,then keep it, else false
        if (substr(d,1,unit) == substr(d, unit+1,unit*2)) {
          #print countrep, rep, tandemcount(rep,d), tandemcount(rep,prime3), "Repeat", unit, d[:unit], d[unit:unit*2]
          countrep = max(countrep, tandemcount_v2(rep,prime3))         
          result$countrep <- countrep
          result$rep <- rep
          return(result)
        } else {
          result$countrep <- 0
          result$rep <- '-'
          return(result)
        }
      } # endif
    } # end for
    result$countrep <- 0 # in case no repeat 2bp or longer is fount
    result$rep <- '-'
    return(result)
  } # end else
} # end repcaller
repcaller_del <- function(d, prime3, prime5, l){
  # d : deletion
  # prime3 : 3' context
  # prime5: 5' context
  # l : length of change
  
  prime3 <- paste0(d,prime3)
  result <- list() 
  result$countrep = 0
  result$rep<-''
  # This is for counting single base repeats
  # if the length of change is 1
  if (l==1) {
    countrep <- 0
    i <- 1
    while (substr(d,1,1)==substr(prime3,i,i)) {
      countrep <- countrep +1
      i <- i + 1
    }
    result$countrep <- countrep
    result$rep <- d
    return(result)
  } else if (d==substr(prime3, 1, nchar(d))) { # This is for counting whole deletion/DI repeats that are present in the flanking region  
    
    countrep = tandemcount_v2(d,prime3)
    rep = findsmallestrep(d)
    countrep = max(countrep,tandemcount_v2(rep,prime3))
    result$countrep <- countrep
    result$rep <- rep
    return(result)
  } else {   # This is for counting anything in between 1bp and whole del repetition                                     
    rep = '-'
    
    for (t in seq(from=(nchar(d)-1), to=2)) {  # Look for repeats of 2bp to n-1 bp of the length of the indel
      if (grepl(substring(d,1,t), prime3)) {
        countrep = tandemcount_v2(substr(d, 1, t),prime3)
        rep = findsmallestrep(substr(d, 1, t))
        unit = tandemcount_v2(rep,d)*nchar(rep)
        # The false calls arise in examples such as these : del = AACCCCATCTCTACT; 3' = AAAATTACAAACAAAT; rep = 'A'; repcount in 3' = 4 which is greater than MH = 2; Therefore, it is called Repeat-mediated
        # In fact, it should be repeat count = 0; Therefore, call should be microhomology mediated. 
        # To do this, compare check how far the repeat stretched into the indel itself. Eg: 'A' is counted twice in the deletion. So compare del[:2] to del[2:4]. If they are the same,then keep it, else false
        if (substr(d,1,unit) == substr(d, unit+1,unit*2)) {
          #print countrep, rep, tandemcount(rep,d), tandemcount(rep,prime3), "Repeat", unit, d[:unit], d[unit:unit*2]
          countrep = max(countrep, tandemcount_v2(rep,prime3))         
          result$countrep <- countrep
          result$rep <- rep
          return(result)
        } else {
          result$countrep <- 0
          result$rep <- '-'
          return(result)
        }
      } # endif
    } # end for
    result$countrep <- 0 # in case no repeat 2bp or longer is fount
    result$rep <- '-'
    return(result)
  } # end else
} # end repcaller

# Given microhomology is in bases, repeat should be in bases too
# If there is a single repeat of the indel 3' of it, then it should be labelled as Repeat-mediated not MH. Eg : TTTA	TTTATTATTAAGATTTTTAAATTTTAATT has 4bp MH and 1 repeat of TTTA. Counting itself, this is a repeat of 2, so it is repeat-mediated.15.05.14
# Except if it is a single base from a longer indel that is repeating, then it is treated as MH

finalcaller <- function(mhcount, replength, rept) {
  # mhcount
  # replength
  # rept
  
  
  if (replength >= mhcount) {
    if ((replength/nchar(rept)) >= 1) {
      return("Repeat-mediated")
    } else if (mhcount > 0) {
      return("Microhomology-mediated")
    } else {
      return("None")
    }
  } else {
    if (mhcount > 0) {
      return("Microhomology-mediated")
    } else { 
      return("None")
    }
  }
  
}

mh <- function(indel.df) {
  
  # indel.df needs following columns:
  # indel.type
  # change
  # slice3
  # slice5
  # indel.length
  
  
  classification <- rep (NA, nrow(indel.df))
  repcount_all <- rep (NA, nrow(indel.df))
  if (nrow(indel.df)>0) {
    for (i in 1:nrow(indel.df)) {
      print(i)
      if (indel.df$Type[i]=='Del') { # the classification is only for deletions
        
        as = as.character(indel.df$change[i]) # The actual deletion
        bs = as.character(indel.df$slice3[i])            
        cs = as.character(indel.df$slice5[i]) # Sequence 5' to deletion            
        
        mhcount = mhcaller(as,bs)$countmh  #  number of microhomology bases            
        
        #Look for Microhomology first and then for repeats - tandem/normal
        r = repcaller(as,bs,cs,indel.df$indel.length[i])
        repcount <- r$countrep # number of times there is a repeat
        rept <- r$rep # the repeat sequence
        
        finalcall = finalcaller(mhcount, repcount*nchar(rept), rept)
        classification[i] <- finalcall
        repcount_all[i] <- repcount
        mhcount_all[i] <- mhcount
        #print a, b, mhcount, mh, repcount*len(repeat), repeat, finalcall
      } # end if
    } # end for
    
    classification[is.na(classification)] <- '-'
    repcount_all[is.na(repcount_all)] <- '-'
    indel.df$repcount <- repcount_all
    indel.df$mhcount <- mhcount_all
    indel.df$classification <- classification
  } else {
  }
  
  indel.df
}  # end the main mh function  

# add repeat insertion by Xueqing
mh_indel_v1 <- function(indel.df) {
  
  # indel.df needs following columns:
  # indel.type
  # change
  # slice3
  # slice5
  # indel.length
  
  
  classification <- rep (NA, nrow(indel.df))
  repcount_all <- rep (0, nrow(indel.df))
  mhcount_all <- rep (0, nrow(indel.df))
  
  if (nrow(indel.df)>0) {
    for (i in 1:nrow(indel.df)) {
      print(i)
      if (indel.df$Type[i]=='Del') { # the classification is only for deletions
        
        as = as.character(indel.df$change[i]) # The actual deletion
        bs = as.character(indel.df$slice3[i])            
        cs = as.character(indel.df$slice5[i]) # Sequence 5' to deletion            
        
        mhcount = mhcaller(as,bs)$countmh  #  number of microhomology bases            
        
        #Look for Microhomology first and then for repeats - tandem/normal
        r = repcaller(as,bs,cs,indel.df$indel.length[i])
        repcount <- r$countrep # number of times there is a repeat
        rept <- r$rep # the repeat sequence
        
        finalcall = finalcaller(mhcount, repcount*nchar(rept), rept)
        classification[i] <- finalcall
        repcount_all[i] <- repcount
        mhcount_all[i] <- mhcount
        
        #print a, b, mhcount, mh, repcount*len(repeat), repeat, finalcall
      } # end if
      
      if (indel.df$Type[i]=='Ins') { # the classification is only for insertions
        
        as = as.character(indel.df$change[i]) # The actual deletion
        bs = as.character(indel.df$slice3[i])            
        cs = as.character(indel.df$slice5[i]) # Sequence 5' to deletion            
        
        mhcount = mhcaller(as,bs)$countmh  #  number of microhomology bases            
        
        r = repcaller(as,bs,cs,indel.df$indel.length[i])
        repcount <- r$countrep # number of times there is a repeat
        rept <- r$rep # the repeat sequence
        finalcall = finalcaller(mhcount, repcount*nchar(rept), rept)
        
        classification[i] <- finalcall
        repcount_all[i] <- repcount
        mhcount_all[i] <- mhcount
        
      }
    } # end for
    
    classification[is.na(classification)] <- '-'
    #  repcount_all[is.na(repcount_all)] <- '-'
    #  mhcount_all[is.na(mhcount_all)] <- '-'
    indel.df$repcount <- repcount_all
    indel.df$mhcount <- mhcount_all
    indel.df$classification <- classification
  } else {
  }
  
  indel.df
}  # end the main mh function  

# Joined deletion and 3' context for deletions, insertion use the original one
# add repeat insertion by Xueqing Oct 12, 2018, repcount includes change
mh_indel_v2 <- function(indel.df) {
  
  # indel.df needs following columns:
  # indel.type
  # change
  # slice3
  # slice5
  # indel.length
  
  
  classification <- rep (NA, nrow(indel.df))
  repcount_all <- rep (0, nrow(indel.df))
  mhcount_all <- rep (0, nrow(indel.df))
  slice3_nonrep_all <- rep (NA, nrow(indel.df))
  if (nrow(indel.df)>0) {
    for (i in 1:nrow(indel.df)) {
      print(i)
      if (indel.df$Type[i]=='Del') { # the classification is only for deletions
        
        as = as.character(indel.df$change[i]) # The actual deletion
        bs = as.character(indel.df$slice3[i])            
        cs = as.character(indel.df$slice5[i]) # Sequence 5' to deletion            
        
        mhcount = mhcaller(as,bs)$countmh  #  number of microhomology bases            
        
        #Look for Microhomology first and then for repeats - tandem/normal
        r_ref = repcaller_del(as,bs,cs,indel.df$indel.length[i])
        repcount_ref <- r_ref$countrep # number of times there is a repeat
        
        r = repcaller(as,bs,cs,indel.df$indel.length[i])
        repcount <- r$countrep # number of times there is a repeat
        rept <- r$rep # the repeat sequence
        
        finalcall = finalcaller(mhcount, repcount*nchar(rept), rept)
        classification[i] <- finalcall
        repcount_all[i] <- repcount_ref
        mhcount_all[i] <- mhcount
        slice3_nonrep_all[i] <- substr(bs,repcount*nchar(rept)+1,repcount*nchar(rept)+1)
        #print a, b, mhcount, mh, repcount*len(repeat), repeat, finalcall
      } # end if
      
      if (indel.df$Type[i]=='Ins') { # the classification is only for insertions
        
        as = as.character(indel.df$change[i]) # The actual deletion
        bs = as.character(indel.df$slice3[i])            
        cs = as.character(indel.df$slice5[i]) # Sequence 5' to deletion            
        
        mhcount = mhcaller(as,bs)$countmh  #  number of microhomology bases            
        
        r = repcaller(as,bs,cs,indel.df$indel.length[i])
        repcount <- r$countrep # number of times there is a repeat
        rept <- r$rep # the repeat sequence
        finalcall = finalcaller(mhcount, repcount*nchar(rept), rept)
        
        classification[i] <- finalcall
        repcount_all[i] <- repcount
        mhcount_all[i] <- mhcount
        slice3_nonrep_all[i] <- substr(bs,repcount*nchar(rept)+1,repcount*nchar(rept)+1)
        
      }
    } # end for
    
    classification[is.na(classification)] <- '-'
    slice3_nonrep_all[is.na(slice3_nonrep_all)] <- '-'
    #  mhcount_all[is.na(mhcount_all)] <- '-'
    indel.df$repcount <- repcount_all
    indel.df$mhcount <- mhcount_all
    indel.df$slice3_nonrep <- slice3_nonrep_all
    indel.df$slice3_nonrep_pyr <- indel.df$slice3_nonrep
    indel.df$slice5_nonrep_pyr <- indel.df$slice5_1bp_pyr
    indel.df[indel.df$change=="A","slice5_nonrep_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.df[indel.df$change=="A","slice3_nonrep"])))
    indel.df[indel.df$change=="A","slice3_nonrep_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.df[indel.df$change=="A","slice5_1bp"])))
    indel.df[indel.df$change=="G","slice5_nonrep_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.df[indel.df$change=="G","slice3_nonrep"])))
    indel.df[indel.df$change=="G","slice3_nonrep_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.df[indel.df$change=="G","slice5_1bp"])))
    
    
    indel.df$classification <- classification
  } else {
  }
  
  indel.df
}  # end the main mh function  


subtype_indel <- function(indel.df) {
  
  indel.df$Subtype <- NULL
  indel.df$TypeS <- ""
  indel.df[indel.df$Type=="Ins","TypeS"] <- "+"
  indel.df[indel.df$Type=="Del","TypeS"] <- "-"
  indel.df[indel.df$Type=="Complex","TypeS"] <- "Complex"
  
  indel.new <- NULL
  #[+C], [+T], [-C], [-T]
  indel.a1 <- subset(indel.df,Type=="Ins" & indel.length==1 & repcount==0)
  indel.a1$Subtype <- paste0(indel.a1$slice5_1bp_pyr,"|[",indel.a1$TypeS,indel.a1$change.pyr,"]|",indel.a1$slice3_1bp_pyr)
  indel.a1$classification <- paste0("[",indel.a1$TypeS,indel.a1$change.pyr,"]")
  indel.new <- rbind(indel.new,indel.a1)
  
  indel.a2 <- subset(indel.df,Type=="Del" & indel.length==1 & repcount==1)
  indel.a2$Subtype <- paste0(indel.a2$slice5_1bp_pyr,"|[",indel.a2$TypeS,indel.a2$change.pyr,"]|",indel.a2$slice3_1bp_pyr)
  indel.a2$classification <- paste0("[",indel.a2$TypeS,indel.a2$change.pyr,"]")
  indel.new <- rbind(indel.new,indel.a2)
  
  indel.b1 <- subset(indel.df,Type=="Ins" & indel.length==1 & repcount>0 & repcount<5)
  indel.b1$Subtype <- paste0(indel.b1$slice5_nonrep_pyr,"|[",indel.b1$TypeS,indel.b1$change.pyr,"]Rep<5|",indel.b1$slice3_nonrep_pyr)
  indel.b1$classification <- paste0("[",indel.b1$TypeS,indel.b1$change.pyr,"]")
  indel.new <- rbind(indel.new,indel.b1)
  
  indel.b2 <- subset(indel.df,Type=="Del" & indel.length==1 & repcount>1 & repcount<5)
  indel.b2$Subtype <- paste0(indel.b2$slice5_nonrep_pyr,"|[",indel.b2$TypeS,indel.b2$change.pyr,"]Rep<5|",indel.b2$slice3_nonrep_pyr)
  indel.b2$classification <- paste0("[",indel.b2$TypeS,indel.b2$change.pyr,"]")
  indel.new <- rbind(indel.new,indel.b2)
  
  
  indel.c <- subset(indel.df,Type%in%c("Ins","Del") & indel.length==1 & repcount>=5 & repcount<8)
  indel.c$Subtype <- paste0("[",indel.c$TypeS,indel.c$change.pyr,"]Rep=",indel.c$repcount)
  indel.c$classification <- paste0("[",indel.c$TypeS,indel.c$change.pyr,"]")
  indel.new <- rbind(indel.new,indel.c)
  
  indel.d <- subset(indel.df,Type%in%c("Ins","Del") & indel.length==1 & repcount>=8)
  indel.d$Subtype <- paste0("[",indel.d$TypeS,indel.d$change.pyr,"]Rep>=8")
  indel.d$classification <- paste0("[",indel.d$TypeS,indel.d$change.pyr,"]")
  indel.new <- rbind(indel.new,indel.d)
  
  # [+>1]Mh, [->1]Mh
  indel.e1 <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Microhomology-mediated" & mhcount==1)
  indel.e1$Subtype <-paste0("[",indel.e1$TypeS,"]Mh=",indel.e1$mhcount)
  indel.e1$classification <- paste0("[",indel.e1$TypeS,"]Mh")
  indel.new <- rbind(indel.new,indel.e1)
  
  indel.e2 <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Microhomology-mediated" & mhcount==2)
  indel.e2$Subtype <-paste0("[",indel.e2$TypeS,"]Mh=",indel.e2$mhcount)
  indel.e2$classification <- paste0("[",indel.e2$TypeS,"]Mh")
  indel.new <- rbind(indel.new,indel.e2)
  
  indel.e3 <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Microhomology-mediated" & mhcount==3)
  indel.e3$Subtype <-paste0("[",indel.e3$TypeS,"]Mh=",indel.e3$mhcount)
  indel.e3$classification <- paste0("[",indel.e3$TypeS,"]Mh")
  indel.new <- rbind(indel.new,indel.e3)
  
  
  indel.f <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Microhomology-mediated" & mhcount>=4 & mhcount<20)
  indel.f$Subtype <-paste0("[",indel.f$TypeS,"]4<=Mh<20")
  indel.f$classification <- paste0("[",indel.f$TypeS,"]Mh")
  indel.new <- rbind(indel.new,indel.f)
  
  indel.g <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Microhomology-mediated" & mhcount>=20)
  indel.g$Subtype <-paste0("[",indel.g$TypeS,"]Mh>=20")
  indel.g$classification <- paste0("[",indel.g$TypeS,"]Mh")
  indel.new <- rbind(indel.new,indel.g)
  
  # [+>1]Rep, [->1]Rep
  indel.h <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Repeat-mediated" & repcount<5)
  indel.h$Subtype <-paste0("[",indel.h$TypeS,">1]Rep<=4")
  indel.h$classification <- paste0("[",indel.h$TypeS,">1]Rep")
  indel.new <- rbind(indel.new,indel.h)
  
  indel.i <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Repeat-mediated" & repcount>=5)
  indel.i$Subtype <-paste0("[",indel.i$TypeS,">1]Rep>=5")
  indel.i$classification <- paste0("[",indel.i$TypeS,">1]Rep")
  indel.new <- rbind(indel.new,indel.i)
  
  # [+>1]Other, [->1]Other
  indel.j <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="None")
  indel.j$Subtype <-paste0("[",indel.j$TypeS,">1]Other")
  indel.j$classification <- paste0("[",indel.j$TypeS,">1]Other")
  indel.new <- rbind(indel.new,indel.j)
  
  # Complex
  indel.k <- subset(indel.df,Type=="Complex")
  indel.k$Subtype <-"Complex"
  indel.k$classification <- "Complex"
  indel.new <- rbind(indel.new,indel.k)
  
  #a=data.frame(table(indel.new$VariantID))
  #a=a[order(a$Freq,decreasing=T),]
  return(indel.new)
}
subtype_indel69 <- function(indel.df) {
  
  indel.df$Subtype <- NULL
  indel.df$TypeS <- ""
  indel.df[indel.df$Type=="Ins","TypeS"] <- "+"
  indel.df[indel.df$Type=="Del","TypeS"] <- "-"
  indel.df[indel.df$Type=="Complex","TypeS"] <- "Complex"
  
  indel.new <- NULL
  #[+C], [+T], [-C], [-T]
  
  indel.b1 <- subset(indel.df,Type=="Ins" & indel.length==1 & repcount<5)
  indel.b1$Subtype <- paste0(indel.b1$slice5_nonrep_pyr,"|[",indel.b1$TypeS,indel.b1$change.pyr,"]Rep<5|",indel.b1$slice3_nonrep_pyr)
  indel.b1$classification <- paste0("[",indel.b1$TypeS,indel.b1$change.pyr,"]")
  indel.new <- rbind(indel.new,indel.b1)
  
  indel.b2 <- subset(indel.df,Type=="Del" & indel.length==1 & repcount<5)
  indel.b2$Subtype <- paste0(indel.b2$slice5_nonrep_pyr,"|[",indel.b2$TypeS,indel.b2$change.pyr,"]Rep<5|",indel.b2$slice3_nonrep_pyr)
  indel.b2$classification <- paste0("[",indel.b2$TypeS,indel.b2$change.pyr,"]")
  indel.new <- rbind(indel.new,indel.b2)
  
  
  indel.c <- subset(indel.df,Type%in%c("Ins","Del") & indel.length==1 & repcount>=5 & repcount<8)
  indel.c$Subtype <- paste0("[",indel.c$TypeS,indel.c$change.pyr,"]Rep=",indel.c$repcount)
  indel.c$classification <- paste0("[",indel.c$TypeS,indel.c$change.pyr,"]")
  indel.new <- rbind(indel.new,indel.c)
  
  indel.d <- subset(indel.df,Type%in%c("Ins","Del") & indel.length==1 & repcount>=8)
  indel.d$Subtype <- paste0("[",indel.d$TypeS,indel.d$change.pyr,"]Rep>=8")
  indel.d$classification <- paste0("[",indel.d$TypeS,indel.d$change.pyr,"]")
  indel.new <- rbind(indel.new,indel.d)
  
  # [+>1]Mh, [->1]Mh
  indel.e1 <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Microhomology-mediated" & mhcount==1)
  indel.e1$Subtype <-paste0("[",indel.e1$TypeS,"]Mh=",indel.e1$mhcount)
  indel.e1$classification <- paste0("[",indel.e1$TypeS,"]Mh")
  indel.new <- rbind(indel.new,indel.e1)
  
  indel.e2 <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Microhomology-mediated" & mhcount==2)
  indel.e2$Subtype <-paste0("[",indel.e2$TypeS,"]Mh=",indel.e2$mhcount)
  indel.e2$classification <- paste0("[",indel.e2$TypeS,"]Mh")
  indel.new <- rbind(indel.new,indel.e2)
  
  indel.e3 <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Microhomology-mediated" & mhcount==3)
  indel.e3$Subtype <-paste0("[",indel.e3$TypeS,"]Mh=",indel.e3$mhcount)
  indel.e3$classification <- paste0("[",indel.e3$TypeS,"]Mh")
  indel.new <- rbind(indel.new,indel.e3)
  
  
  indel.f <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Microhomology-mediated" & mhcount>=4 & mhcount<20)
  indel.f$Subtype <-paste0("[",indel.f$TypeS,"]4<=Mh<20")
  indel.f$classification <- paste0("[",indel.f$TypeS,"]Mh")
  indel.new <- rbind(indel.new,indel.f)
  
  indel.g <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Microhomology-mediated" & mhcount>=20)
  indel.g$Subtype <-paste0("[",indel.g$TypeS,"]Mh>=20")
  indel.g$classification <- paste0("[",indel.g$TypeS,"]Mh")
  indel.new <- rbind(indel.new,indel.g)
  
  # [+>1]Rep, [->1]Rep
  indel.h <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Repeat-mediated" & repcount<5)
  indel.h$Subtype <-paste0("[",indel.h$TypeS,">1]Rep<=4")
  indel.h$classification <- paste0("[",indel.h$TypeS,">1]Rep")
  indel.new <- rbind(indel.new,indel.h)
  
  indel.i <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Repeat-mediated" & repcount>=5)
  indel.i$Subtype <-paste0("[",indel.i$TypeS,">1]Rep>=5")
  indel.i$classification <- paste0("[",indel.i$TypeS,">1]Rep")
  indel.new <- rbind(indel.new,indel.i)
  
  # [+>1]Other, [->1]Other
  indel.j <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="None")
  indel.j$Subtype <-paste0("[",indel.j$TypeS,">1]Other")
  indel.j$classification <- paste0("[",indel.j$TypeS,">1]Other")
  indel.new <- rbind(indel.new,indel.j)
  
  # Complex
  indel.k <- subset(indel.df,Type=="Complex")
  indel.k$Subtype <-"Complex"
  indel.k$classification <- "Complex"
  indel.new <- rbind(indel.new,indel.k)
  
  #a=data.frame(table(indel.new$VariantID))
  #a=a[order(a$Freq,decreasing=T),]
  return(indel.new)
}
subtype_indel_mmrd <- function(indel.df) {
  
  indel.df$Subtype <- NULL
  indel.df$TypeS <- ""
  indel.df[indel.df$Type=="Ins","TypeS"] <- "+"
  indel.df[indel.df$Type=="Del","TypeS"] <- "-"
  indel.df[indel.df$Type=="Complex","TypeS"] <- "Complex"
  
  indel.new <- NULL
  #[+C], [+T], [-C], [-T]
  
  #  1 rep
  indel.c <- subset(indel.df,Type%in%c("Ins","Del") & indel.length==1)
  indel.c$Subtype <- paste0(indel.c$slice5_nonrep_pyr,"|[",indel.c$TypeS,indel.c$change.pyr,"]Rep=",indel.c$repcount,"|",indel.c$slice3_nonrep_pyr)
  
  indel.c$classification <- paste0("[",indel.c$TypeS,indel.c$change.pyr,"]")
  indel.new <- rbind(indel.new,indel.c)
  
  # [+>1]Mh, [->1]Mh
  indel.e1 <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Microhomology-mediated")
  indel.e1$Subtype <-paste0("[",indel.e1$TypeS,"]Mh")
  
  indel.e1$classification <- paste0("[",indel.e1$TypeS,"]Mh")
  indel.new <- rbind(indel.new,indel.e1)
  
  
  # [+>1]Rep, [->1]Rep
  indel.h <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Repeat-mediated" & repcount<5)
  indel.h$Subtype <-paste0("[",indel.h$TypeS,">1]Rep<=4")
  indel.h$classification <- paste0("[",indel.h$TypeS,">1]Rep")
  indel.new <- rbind(indel.new,indel.h)
  
  indel.i <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="Repeat-mediated" & repcount>=5)
  indel.i$Subtype <-paste0("[",indel.i$TypeS,">1]Rep>=5")
  indel.i$classification <- paste0("[",indel.i$TypeS,">1]Rep")
  indel.new <- rbind(indel.new,indel.i)
  
  # [+>1]Other, [->1]Other
  indel.j <- subset(indel.df,Type%in%c("Ins","Del") & indel.length>1 & classification=="None")
  if(dim(indel.j)[1]>0){
    indel.j$Subtype <-paste0("[",indel.j$TypeS,">1]Other")
    indel.j$classification <- paste0("[",indel.j$TypeS,">1]Other")
    indel.new <- rbind(indel.new,indel.j)
  }

  
  # Complex
  indel.k <- subset(indel.df,Type=="Complex")
  if(dim(indel.k)[1]>0){
    indel.k$Subtype <-"Complex"
    indel.k$classification <- "Complex"
    indel.new <- rbind(indel.new,indel.k)
    
  }
  

  #a=data.frame(table(indel.new$VariantID))
  #a=a[order(a$Freq,decreasing=T),]
  return(indel.new)
}


indel_classifier <- function(indels){
  
  indel.data <- indels
  indel.data[indel.data$Chrom=="23","Chrom"]="X"
  indel.data[indel.data$Chrom=="24","Chrom"]="Y"
  
  # convert formats, and find context of the indels
  indel.df <- prepare.indel.df.tab(indel.data)
  indel.df.max100 <- indel.df[indel.df$indel.length<=100,]
  # indel classification
  indel.classified.df <- mh_indel_v2(indel.df.max100)
  
  # remove the insertions with repcount =10
  indel.classified.df <- indel.classified.df[indel.classified.df$repcount<10,]
  
  indel.classified_details <- subtype_indel_mmrd(indel.classified.df)
  
  indel.classified_details[indel.classified_details$Chrom=="X","Chrom"]=23
  indel.classified_details[indel.classified_details$Chrom=="Y","Chrom"]=24
  indel_templateMMR <- read.table("../00_common/indel_templateMMRD.txt",sep = "\t",header = T, as.is = T)
  indel_template2 <- indel_templateMMR
  names(indel_template2)[1] <- c("Subtype")
  indel.classified_details <- merge(indel.classified_details,indel_template2,by="Subtype")
  
  return(indel.classified_details)
  
}

indel_classifier15 <- function(indels){
  
  indel.data <- indels
  indel.data[indel.data$Chrom=="23","Chrom"]="X"
  indel.data[indel.data$Chrom=="24","Chrom"]="Y"
  
  # convert formats, and find context of the indels
  indel.df <- prepare.indel.df.tab(indel.data)
  indel.df.max100 <- indel.df[indel.df$indel.length<=100,]
  # indel classification
  indel.classified.df <- mh_indel_v2(indel.df.max100)
  
  # remove the insertions with repcount =10
  indel.classified.df <- indel.classified.df[indel.classified.df$repcount<10,]
  

  indel.classified_details <- subtype_indel69(indel.classified.df)
  
  indel.classified_details[indel.classified_details$Chrom=="X","Chrom"]=23
  indel.classified_details[indel.classified_details$Chrom=="Y","Chrom"]=24
  indel_template2 <- read.table("../00_common/indel_template69.txt",sep = "\t",header = T, as.is = T)

  names(indel_template2)[1] <- c("Subtype")
  indel.classified_details <- merge(indel.classified_details,indel_template2,by="Subtype")
  indel.classified_details$Subtype <- indel.classified_details$indeltype_short
  
  
  return(indel.classified_details)
  
}
##########################
#
#  Indel profiles
#
##########################

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


plotbasis_average_se_indel_10types <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  
  muts_basis_melt <- melt(muts_basis,"indelsubtype")
  names(muts_basis_melt) <- c("indelsubtype","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("indelsubtype"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N))
  muts_basis_melt_summary_percentage <- muts_basis_melt_summary
  muts_basis_melt_summary_percentage$mean <- muts_basis_melt_summary$mean/mean_parentmuts
  muts_basis_melt_summary_percentage$sd <- muts_basis_melt_summary$sd/mean_parentmuts
  muts_basis_melt_summary_percentage$se <- muts_basis_melt_summary$se/mean_parentmuts
  muts_basis_melt_summary_percentage$indeltype <- c("Mh-mediated_Del","Other_Del","Other_Del","Other_Del","Rep-mediated_Del","Rep-mediated_Del","Rep-mediated_Del","Insertion","Insertion","Insertion")
  
  
  indel_mypalette <- c("#D55E00","#F0E442","#CC79A7","#009E73")
  indel_positions <- c("Ins_-_1_C","Ins_-_1_T","Ins_-_>1bp","Del_Microhomology-mediated_>1bp","Del_Repeat-mediated_1_C","Del_Repeat-mediated_1_T","Del_Repeat-mediated_>1bp","Del_None_1_C","Del_None_1_T","Del_None_>1bp") 
  indel_labels <- c("Insertion=C","Insertion=T","Insertion>=2bp","Mh-mediated Del","Rep-mediated Del.=C","Rep-mediated Del.=T","Rep-mediated Del.>=2bp","Other Del.=C","Other Del.=T","Other Del.>=2bp") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary_percentage, aes(x=indelsubtype, y=mean,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Percentage")
  p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="blue",position=position_dodge(.9),width=.2)+scale_y_continuous(limits=c(0,1),breaks=(seq(0,1,0.2)),labels = scales::percent)
  p <- p+scale_x_discrete(limits = indel_positions,labels=indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=7,colour = "black"),
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
plotCountbasis_average_se_indel_10types <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  
  muts_basis_melt <- melt(muts_basis,"indelsubtype")
  names(muts_basis_melt) <- c("indelsubtype","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("indelsubtype"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N))
  muts_basis_melt_summary$indeltype <- c("Mh-mediated_Del","Other_Del","Other_Del","Other_Del","Rep-mediated_Del","Rep-mediated_Del","Rep-mediated_Del","Insertion","Insertion","Insertion")
  
  
  indel_mypalette <- c("#D55E00","#F0E442","#CC79A7","#009E73")
  indel_positions <- c("Ins_-_1_C","Ins_-_1_T","Ins_-_>1bp","Del_Microhomology-mediated_>1bp","Del_Repeat-mediated_1_C","Del_Repeat-mediated_1_T","Del_Repeat-mediated_>1bp","Del_None_1_C","Del_None_1_T","Del_None_>1bp") 
  indel_labels <- c("Insertion=C","Insertion=T","Insertion>=2bp","Mh-mediated Del","Rep-mediated Del.=C","Rep-mediated Del.=T","Rep-mediated Del.>=2bp","Other Del.=C","Other Del.=T","Other Del.>=2bp") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary, aes(x=indelsubtype, y=mean,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Count")
  p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="blue",position=position_dodge(.9),width=.2)  +scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels=indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=7,colour = "black"),
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
  return(muts_basis_melt_summary)
}

plotbasis_average_se_indel_25types <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  
  muts_basis_melt <- melt(muts_basis,"indelsubtype")
  names(muts_basis_melt) <- c("indelsubtype","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("indelsubtype"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N))
  muts_basis_melt_summary_percentage <- muts_basis_melt_summary
  muts_basis_melt_summary_percentage$mean <- muts_basis_melt_summary$mean/mean_parentmuts
  muts_basis_melt_summary_percentage$sd <- muts_basis_melt_summary$sd/mean_parentmuts
  muts_basis_melt_summary_percentage$se <- muts_basis_melt_summary$se/mean_parentmuts
  #muts_basis_melt_summary_percentage$indeltype <- sub(".*_","",muts_basis_melt_summary_percentage$indelsubtype)
  
  muts_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/21_indels/23_results_25channels/indel_25channels.txt",sep = "\t",header = T, as.is = T)
  muts_basis_melt_summary_percentage <- merge(muts_template, muts_basis_melt_summary_percentage,by="indelsubtype",all.x=T)
  muts_basis_melt_summary_percentage[is.na(muts_basis_melt_summary_percentage)] <- 0
  
  # D55E00(orange, OtherIns); F0E442(yellow, MhDel); CC79A7(pink, OtherDel); 009E73(green, RepDel); 0072B2(blue, RepIns)
  
  indel_mypalette <- c("#F0E442","#CC79A7","#D55E00","#009E73","#0072B2")
  indel_positions <- c("[+C]_RepIns","[+T]_RepIns","[+>1]_RepIns","[+C]A_OtherIns", "[+C]C_OtherIns", "[+C]G_OtherIns", "[+C]T_OtherIns","[+T]A_OtherIns", "[+T]C_OtherIns", "[+T]G_OtherIns","[+T]T_OtherIns","[+>1]_OtherIns",
                       "[-]_MhDel", "[-C]_RepDel", "[-T]_RepDel", "[->1]_RepDel", "[-C]A_OtherDel", "[-C]C_OtherDel", "[-C]G_OtherDel", "[-C]T_OtherDel", "[-T]A_OtherDel", "[-T]C_OtherDel", "[-T]G_OtherDel","[-T]T_OtherDel", "[->1]_OtherDel") 
  #  indel_labels <- c("Insertion=C","Insertion=T","Insertion>=2bp","Mh-mediated Del","Rep-mediated Del.=C","Rep-mediated Del.=T","Rep-mediated Del.>=2bp","Other Del.=C","Other Del.=T","Other Del.>=2bp") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary_percentage, aes(x=indelsubtype, y=mean,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Percentage")
  p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="blue",position=position_dodge(.9),width=.2)+scale_y_continuous(limits=c(0,1),breaks=(seq(0,1,0.2)),labels = scales::percent)
  p <- p+scale_x_discrete(limits = indel_positions)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=7,colour = "black"),
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
plotCountbasis_average_se_indel_25types <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  
  muts_basis_melt <- melt(muts_basis,"indelsubtype")
  names(muts_basis_melt) <- c("indelsubtype","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("indelsubtype"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N))
  #muts_basis_melt_summary$indeltype <- sub(".*_","",muts_basis_melt_summary$indelsubtype)
  muts_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/21_indels/23_results_25channels/indel_25channels.txt",sep = "\t",header = T, as.is = T)
  muts_basis_melt_summary <- merge(muts_template, muts_basis_melt_summary,by="indelsubtype",all.x=T)
  muts_basis_melt_summary[is.na(muts_basis_melt_summary)] <- 0
  
  
  # D55E00(orange, OtherIns); F0E442(yellow, MhDel); CC79A7(pink, OtherDel); 009E73(green, RepDel); 0072B2(blue, RepIns)
  
  indel_mypalette <- c("#F0E442","#CC79A7","#D55E00","#009E73","#0072B2")
  indel_positions <- c("[+C]_RepIns","[+T]_RepIns","[+>1]_RepIns","[+C]A_OtherIns", "[+C]C_OtherIns", "[+C]G_OtherIns", "[+C]T_OtherIns","[+T]A_OtherIns", "[+T]C_OtherIns", "[+T]G_OtherIns","[+T]T_OtherIns","[+>1]_OtherIns",
                       "[-]_MhDel", "[-C]_RepDel", "[-T]_RepDel", "[->1]_RepDel", "[-C]A_OtherDel", "[-C]C_OtherDel", "[-C]G_OtherDel", "[-C]T_OtherDel", "[-T]A_OtherDel", "[-T]C_OtherDel", "[-T]G_OtherDel","[-T]T_OtherDel", "[->1]_OtherDel") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary, aes(x=indelsubtype, y=mean,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Count")
  p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="blue",position=position_dodge(.9),width=.2)  +scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=7,colour = "black"),
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
  return(muts_basis_melt_summary)
}

plotbasis_average_se_indel_13types <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  
  muts_basis_melt <- melt(muts_basis,"indelsubtype")
  names(muts_basis_melt) <- c("indelsubtype","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("indelsubtype"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N))
  muts_basis_melt_summary_percentage <- muts_basis_melt_summary
  muts_basis_melt_summary_percentage$mean <- muts_basis_melt_summary$mean/mean_parentmuts
  muts_basis_melt_summary_percentage$sd <- muts_basis_melt_summary$sd/mean_parentmuts
  muts_basis_melt_summary_percentage$se <- muts_basis_melt_summary$se/mean_parentmuts
  #muts_basis_melt_summary_percentage$indeltype <- sub(".*_","",muts_basis_melt_summary_percentage$indelsubtype)
  
  muts_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/21_indels/24_results_13channels/indel_13channels.txt",sep = "\t",header = T, as.is = T)
  muts_basis_melt_summary_percentage <- merge(muts_template, muts_basis_melt_summary_percentage,by="indelsubtype",all.x=T)
  muts_basis_melt_summary_percentage[is.na(muts_basis_melt_summary_percentage)] <- 0
  
  # D55E00(orange, OtherIns); F0E442(yellow, MhDel); CC79A7(pink, OtherDel); 009E73(green, RepDel); 0072B2(blue, RepIns)
  
  indel_mypalette <- c("#F0E442","#CC79A7","#D55E00","#009E73","#0072B2")
  indel_positions <- c("[+C]_OtherIns","[+T]_OtherIns","[+>1]_OtherIns","[+C]_RepIns","[+T]_RepIns","[+>1]_RepIns",
                       "[-C]_OtherDel", "[-T]_OtherDel", "[->1]_OtherDel","[-C]_RepDel", "[-T]_RepDel", "[->1]_RepDel","[-]_MhDel") 
  #  indel_labels <- c("Insertion=C","Insertion=T","Insertion>=2bp","Mh-mediated Del","Rep-mediated Del.=C","Rep-mediated Del.=T","Rep-mediated Del.>=2bp","Other Del.=C","Other Del.=T","Other Del.>=2bp") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary_percentage, aes(x=indelsubtype, y=mean,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Percentage")
  p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="blue",position=position_dodge(.9),width=.2)+scale_y_continuous(limits=c(0,1),breaks=(seq(0,1,0.2)),labels = scales::percent)
  p <- p+scale_x_discrete(limits = indel_positions)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=7,colour = "black"),
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
plotCountbasis_average_se_indel_13types <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  
  muts_basis_melt <- melt(muts_basis,"indelsubtype")
  names(muts_basis_melt) <- c("indelsubtype","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("indelsubtype"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N))
  #muts_basis_melt_summary$indeltype <- sub(".*_","",muts_basis_melt_summary$indelsubtype)
  muts_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/21_indels/24_results_13channels/indel_13channels.txt",sep = "\t",header = T, as.is = T)
  muts_basis_melt_summary <- merge(muts_template, muts_basis_melt_summary,by="indelsubtype",all.x=T)
  muts_basis_melt_summary[is.na(muts_basis_melt_summary)] <- 0
  
  
  # D55E00(orange, OtherIns); F0E442(yellow, MhDel); CC79A7(pink, OtherDel); 009E73(green, RepDel); 0072B2(blue, RepIns)
  
  indel_mypalette <- c("#F0E442","#CC79A7","#D55E00","#009E73","#0072B2")
  indel_positions <- c("[+C]_OtherIns","[+T]_OtherIns","[+>1]_OtherIns","[+C]_RepIns","[+T]_RepIns","[+>1]_RepIns",
                       "[-C]_OtherDel", "[-T]_OtherDel", "[->1]_OtherDel","[-C]_RepDel", "[-T]_RepDel", "[->1]_RepDel","[-]_MhDel") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary, aes(x=indelsubtype, y=mean,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Count")
  p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="blue",position=position_dodge(.9),width=.2)  +scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=7,colour = "black"),
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
  return(muts_basis_melt_summary)
}

plotbasis_average_se_indel_13types_6 <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  
  muts_basis_melt <- melt(muts_basis,"indelsubtype")
  names(muts_basis_melt) <- c("indelsubtype","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("indelsubtype"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N))
  muts_basis_melt_summary_percentage <- muts_basis_melt_summary
  muts_basis_melt_summary_percentage$mean <- muts_basis_melt_summary$mean/mean_parentmuts
  muts_basis_melt_summary_percentage$sd <- muts_basis_melt_summary$sd/mean_parentmuts
  muts_basis_melt_summary_percentage$se <- muts_basis_melt_summary$se/mean_parentmuts
  #muts_basis_melt_summary_percentage$indeltype <- sub("\\_.*","",muts_basis_melt_summary_percentage$indelsubtype)
  
  muts_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/21_indels/25_results_13channels_6types/indel_13channels_6.txt",sep = "\t",header = T, as.is = T)
  muts_basis_melt_summary_percentage <- merge(muts_template, muts_basis_melt_summary_percentage,by="indelsubtype",all.x=T)
  muts_basis_melt_summary_percentage[is.na(muts_basis_melt_summary_percentage)] <- 0
  
  # E69F00(light orange, +C); D55E00(dark orange, +T); F0E442(yellow, +>1); 56B4E9(light blue, -C); 0072B2(dark blue, -T); 009E73(green, ->1); CC79A7(pink, MhDel)
  #indel_mypalette <- c("#009E73","#56B4E9","#0072B2","#F0E442","#E69F00","#D55E00")
  
  # 009E73 green, ->1; 56B4E9 light blue -C; E69F00 light orange -T; CC79A7 pink +>1; 0072B2 dark blue +C; D55E00 dark orange +T
  indel_mypalette <- c("#009E73","#56B4E9","#E69F00","#CC79A7","#0072B2","#D55E00")
  
  
  indel_positions <- c("[+C]_OtherIns","[+C]_RepIns","[+T]_OtherIns","[+T]_RepIns","[+>1]_OtherIns","[+>1]_RepIns",
                       "[-C]_OtherDel", "[-C]_RepDel", "[-T]_OtherDel","[-T]_RepDel", "[->1]_OtherDel","[->1]_RepDel","[->1]_MhDel") 
  #  indel_labels <- c("Insertion=C","Insertion=T","Insertion>=2bp","Mh-mediated Del","Rep-mediated Del.=C","Rep-mediated Del.=T","Rep-mediated Del.>=2bp","Other Del.=C","Other Del.=T","Other Del.>=2bp") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary_percentage, aes(x=indelsubtype, y=mean,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Percentage")
  p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="blue",position=position_dodge(.9),width=.2)+scale_y_continuous(limits=c(0,1),breaks=(seq(0,1,0.2)),labels = scales::percent)
  p <- p+scale_x_discrete(limits = indel_positions)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=7,colour = "black"),
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
plotCountbasis_average_se_indel_13types_6 <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  
  muts_basis_melt <- melt(muts_basis,"indelsubtype")
  names(muts_basis_melt) <- c("indelsubtype","sample","count")
  muts_basis_melt_summary <- ddply(muts_basis_melt,c("indelsubtype"),summarise, N=length(count),mean=mean(count),sd=sd(count),se=sd/sqrt(N))
  #muts_basis_melt_summary$indeltype <- sub(".*_","",muts_basis_melt_summary$indelsubtype)
  muts_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/21_indels/25_results_13channels_6types/indel_13channels_6.txt",sep = "\t",header = T, as.is = T)
  muts_basis_melt_summary <- merge(muts_template, muts_basis_melt_summary,by="indelsubtype",all.x=T)
  muts_basis_melt_summary[is.na(muts_basis_melt_summary)] <- 0
  
  # E69F00(light orange, +C); D55E00(dark orange, +T); F0E442(yellow, +>1); 56B4E9(light blue, -C); 0072B2(dark blue, -T); 009E73(green, ->1); CC79A7(pink, MhDel)
  #indel_mypalette <- c("#009E73","#56B4E9","#0072B2","#F0E442","#E69F00","#D55E00")
  
  # 009E73 green, ->1; 56B4E9 light blue -C; E69F00 light orange -T; CC79A7 pink +>1; 0072B2 dark blue +C; D55E00 dark orange +T
  indel_mypalette <- c("#009E73","#56B4E9","#E69F00","#CC79A7","#0072B2","#D55E00")
  
  
  indel_positions <- c("[+C]_OtherIns","[+C]_RepIns","[+T]_OtherIns","[+T]_RepIns","[+>1]_OtherIns","[+>1]_RepIns",
                       "[-C]_OtherDel", "[-C]_RepDel", "[-T]_OtherDel","[-T]_RepDel", "[->1]_OtherDel","[->1]_RepDel","[->1]_MhDel") 
  #  indel_labels <- c("Insertion=C","Insertion=T","Insertion>=2bp","Mh-mediated Del","Rep-mediated Del.=C","Rep-mediated Del.=T","Rep-mediated Del.>=2bp","Other Del.=C","Other Del.=T","Other Del.>=2bp") 
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt_summary, aes(x=indelsubtype, y=mean,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Count")
  p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="blue",position=position_dodge(.9),width=.2)  +scale_y_continuous(limits=c(0,110),breaks=(seq(0,110,20)))
  p <- p+scale_x_discrete(limits = indel_positions)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=7,colour = "black"),
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
  return(muts_basis_melt_summary)
}

plotCountbasis_aggragated_indel_105types_6 <- function(muts_basis,h,w,outputname){
  #mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_basis$aggragate <- rowSums(muts_basis[,-1])
  
  muts_template <-  read.table("/nfs/cancer_archive04/xz3/b_1176/25_indels_20180907/indel_template105.txt",sep = "\t",header = T, as.is = T)

  muts_basis <- merge(muts_template, muts_basis,by="indelsubtype",all.x=T)
  muts_basis[is.na(muts_basis)] <- 0
  
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  indel_mypalette_fill <- c("pink", "deeppink","lightcyan","skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","green","purple")
  indel_mypalette_outline <- c("grey", "white","black")
  # Outline: Insertion: black; Deletion: grey
  # 
  # Fill: [+C]: lightcyan; [+C]Rep<5: skyblue; [+C]Rep>=5: royalblue; [+T]: lightgoldenrod; [+T]Rep<5: goldenrod1; [+T]Rep>=5: darkorgane; 
  #        [+]Mh: deeppink; [+>1]Rep: hotpink; [+>1]Other: pink
  #        Complex:	green
  
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_positions <- c("A|[+C]|A","A|[+C]|G","A|[+C]|T","G|[+C]|A","G|[+C]|G","G|[+C]|T","T|[+C]|A","T|[+C]|G","T|[+C]|T",
                       "A|[+C]Rep<5|A","A|[+C]Rep<5|G","A|[+C]Rep<5|T","G|[+C]Rep<5|A","G|[+C]Rep<5|G","G|[+C]Rep<5|T","T|[+C]Rep<5|A","T|[+C]Rep<5|G","T|[+C]Rep<5|T",
                       "[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep>=8",
                       "A|[+T]|A","A|[+T]|C","A|[+T]|G","C|[+T]|A","C|[+T]|C","C|[+T]|G","G|[+T]|A","G|[+T]|C","G|[+T]|G",
                       "A|[+T]Rep<5|A","A|[+T]Rep<5|C","A|[+T]Rep<5|G","C|[+T]Rep<5|A","C|[+T]Rep<5|C","C|[+T]Rep<5|G","G|[+T]Rep<5|A","G|[+T]Rep<5|C","G|[+T]Rep<5|G",
                       "[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep>=8",
                       "[+]Mh=1","[+]Mh=2","[+]Mh=3","[+]4<=Mh<20","[+]Mh>=20","[+>1]Rep<=4","[+>1]Rep>=5","[+>1]Other",
                       "A|[-C]|A","A|[-C]|G","A|[-C]|T","G|[-C]|A","G|[-C]|G","G|[-C]|T","T|[-C]|A","T|[-C]|G","T|[-C]|T",
                       "A|[-C]Rep<5|A","A|[-C]Rep<5|G","A|[-C]Rep<5|T","G|[-C]Rep<5|A","G|[-C]Rep<5|G","G|[-C]Rep<5|T","T|[-C]Rep<5|A","T|[-C]Rep<5|G","T|[-C]Rep<5|T",
                       "[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep>=8",
                       "A|[-T]|A","A|[-T]|C","A|[-T]|G","C|[-T]|A","C|[-T]|C","C|[-T]|G","G|[-T]|A","G|[-T]|C","G|[-T]|G",
                       "A|[-T]Rep<5|A","A|[-T]Rep<5|C","A|[-T]Rep<5|G","C|[-T]Rep<5|A","C|[-T]Rep<5|C","C|[-T]Rep<5|G","G|[-T]Rep<5|A","G|[-T]Rep<5|C","G|[-T]Rep<5|G",
                       "[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep>=8",
                       "[-]Mh=1","[-]Mh=2","[-]Mh=3","[-]4<=Mh<20","[-]Mh>=20","[->1]Rep<=4","[->1]Rep>=5","[->1]Other",
                       "Complex") 
  
  indel_labels <- c("A|[+C]|A","A|[+C]|G","A|[+C]|T","G|[+C]|A","G|[+C]|G","G|[+C]|T","T|[+C]|A","T|[+C]|G","T|[+C]|T",
                    "A|[+C]Rep<5|A","A|[+C]Rep<5|G","A|[+C]Rep<5|T","G|[+C]Rep<5|A","G|[+C]Rep<5|G","G|[+C]Rep<5|T","T|[+C]Rep<5|A","T|[+C]Rep<5|G","T|[+C]Rep<5|T",
                    "[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep>=8",
                    "A|[+T]|A","A|[+T]|C","A|[+T]|G","C|[+T]|A","C|[+T]|C","C|[+T]|G","G|[+T]|A","G|[+T]|C","G|[+T]|G",
                    "A|[+T]Rep<5|A","A|[+T]Rep<5|C","A|[+T]Rep<5|G","C|[+T]Rep<5|A","C|[+T]Rep<5|C","C|[+T]Rep<5|G","G|[+T]Rep<5|A","G|[+T]Rep<5|C","G|[+T]Rep<5|G",
                    "[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep>=8",
                    "[+]Mh=1","[+]Mh=2","[+]Mh=3","[+]4<=Mh<20","[+]Mh>=20","[+>1]Rep<=4","[+>1]Rep>=5","[+>1]Other",
                    "A|[-C]|A","A|[-C]|G","A|[-C]|T","G|[-C]|A","G|[-C]|G","G|[-C]|T","T|[-C]|A","T|[-C]|G","T|[-C]|T",
                    "A|[-C]Rep<5|A","A|[-C]Rep<5|G","A|[-C]Rep<5|T","G|[-C]Rep<5|A","G|[-C]Rep<5|G","G|[-C]Rep<5|T","T|[-C]Rep<5|A","T|[-C]Rep<5|G","T|[-C]Rep<5|T",
                    "[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep>=8",
                    "A|[-T]|A","A|[-T]|C","A|[-T]|G","C|[-T]|A","C|[-T]|C","C|[-T]|G","G|[-T]|A","G|[-T]|C","G|[-T]|G",
                    "A|[-T]Rep<5|A","A|[-T]Rep<5|C","A|[-T]Rep<5|G","C|[-T]Rep<5|A","C|[-T]Rep<5|C","C|[-T]Rep<5|G","G|[-T]Rep<5|A","G|[-T]Rep<5|C","G|[-T]Rep<5|G",
                    "[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep>=8",
                    "[-]Mh=1","[-]Mh=2","[-]Mh=3","[-]4<=Mh<20","[-]Mh>=20","[->1]Rep<=4","[->1]Rep>=5","[->1]Other",
                    "Complex") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis[,c("indelsubtype","indeltype","type","aggragate")], aes(x=indelsubtype, y=aggragate,fill=indeltype, colour=type))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Count")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  p <- p+scale_color_manual(values=indel_mypalette_outline)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10,colour = "black",hjust=1),
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
  
  return(muts_basis)
}
plotCountbasis_single_indel_105types_6 <- function(muts_basis,h,w,outputname){
  #mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_template <-  read.table("/nfs/cancer_archive04/xz3/b_1176/25_indels_20180907/indel_template105.txt",sep = "\t",header = T, as.is = T)
  
  muts_basis <- merge(muts_template, muts_basis,by="indelsubtype",all.x=T)
  muts_basis[is.na(muts_basis)] <- 0
  names(muts_basis) <- c("indelsubtype","indeltype","type","freq")
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  indel_mypalette_fill <- c("pink", "deeppink","lightcyan","skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","green","purple")
  indel_mypalette_outline <- c("grey", "white","black")
  # Outline: Insertion: black; Deletion: grey
  # 
  # Fill: [+C]: lightcyan; [+C]Rep<5: skyblue; [+C]Rep>=5: royalblue; [+T]: lightgoldenrod; [+T]Rep<5: goldenrod1; [+T]Rep>=5: darkorgane; 
  #        [+]Mh: deeppink; [+>1]Rep: hotpink; [+>1]Other: pink
  #        Complex:	green
  
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_positions <- c("A|[+C]|A","A|[+C]|G","A|[+C]|T","G|[+C]|A","G|[+C]|G","G|[+C]|T","T|[+C]|A","T|[+C]|G","T|[+C]|T",
                       "A|[+C]Rep<5|A","A|[+C]Rep<5|G","A|[+C]Rep<5|T","G|[+C]Rep<5|A","G|[+C]Rep<5|G","G|[+C]Rep<5|T","T|[+C]Rep<5|A","T|[+C]Rep<5|G","T|[+C]Rep<5|T",
                       "[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep>=8",
                       "A|[+T]|A","A|[+T]|C","A|[+T]|G","C|[+T]|A","C|[+T]|C","C|[+T]|G","G|[+T]|A","G|[+T]|C","G|[+T]|G",
                       "A|[+T]Rep<5|A","A|[+T]Rep<5|C","A|[+T]Rep<5|G","C|[+T]Rep<5|A","C|[+T]Rep<5|C","C|[+T]Rep<5|G","G|[+T]Rep<5|A","G|[+T]Rep<5|C","G|[+T]Rep<5|G",
                       "[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep>=8",
                       "[+]Mh=1","[+]Mh=2","[+]Mh=3","[+]4<=Mh<20","[+]Mh>=20","[+>1]Rep<=4","[+>1]Rep>=5","[+>1]Other",
                       "A|[-C]|A","A|[-C]|G","A|[-C]|T","G|[-C]|A","G|[-C]|G","G|[-C]|T","T|[-C]|A","T|[-C]|G","T|[-C]|T",
                       "A|[-C]Rep<5|A","A|[-C]Rep<5|G","A|[-C]Rep<5|T","G|[-C]Rep<5|A","G|[-C]Rep<5|G","G|[-C]Rep<5|T","T|[-C]Rep<5|A","T|[-C]Rep<5|G","T|[-C]Rep<5|T",
                       "[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep>=8",
                       "A|[-T]|A","A|[-T]|C","A|[-T]|G","C|[-T]|A","C|[-T]|C","C|[-T]|G","G|[-T]|A","G|[-T]|C","G|[-T]|G",
                       "A|[-T]Rep<5|A","A|[-T]Rep<5|C","A|[-T]Rep<5|G","C|[-T]Rep<5|A","C|[-T]Rep<5|C","C|[-T]Rep<5|G","G|[-T]Rep<5|A","G|[-T]Rep<5|C","G|[-T]Rep<5|G",
                       "[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep>=8",
                       "[-]Mh=1","[-]Mh=2","[-]Mh=3","[-]4<=Mh<20","[-]Mh>=20","[->1]Rep<=4","[->1]Rep>=5","[->1]Other",
                       "Complex") 
  
  indel_labels <- c("A|[+C]|A","A|[+C]|G","A|[+C]|T","G|[+C]|A","G|[+C]|G","G|[+C]|T","T|[+C]|A","T|[+C]|G","T|[+C]|T",
                    "A|[+C]Rep<5|A","A|[+C]Rep<5|G","A|[+C]Rep<5|T","G|[+C]Rep<5|A","G|[+C]Rep<5|G","G|[+C]Rep<5|T","T|[+C]Rep<5|A","T|[+C]Rep<5|G","T|[+C]Rep<5|T",
                    "[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep>=8",
                    "A|[+T]|A","A|[+T]|C","A|[+T]|G","C|[+T]|A","C|[+T]|C","C|[+T]|G","G|[+T]|A","G|[+T]|C","G|[+T]|G",
                    "A|[+T]Rep<5|A","A|[+T]Rep<5|C","A|[+T]Rep<5|G","C|[+T]Rep<5|A","C|[+T]Rep<5|C","C|[+T]Rep<5|G","G|[+T]Rep<5|A","G|[+T]Rep<5|C","G|[+T]Rep<5|G",
                    "[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep>=8",
                    "[+]Mh=1","[+]Mh=2","[+]Mh=3","[+]4<=Mh<20","[+]Mh>=20","[+>1]Rep<=4","[+>1]Rep>=5","[+>1]Other",
                    "A|[-C]|A","A|[-C]|G","A|[-C]|T","G|[-C]|A","G|[-C]|G","G|[-C]|T","T|[-C]|A","T|[-C]|G","T|[-C]|T",
                    "A|[-C]Rep<5|A","A|[-C]Rep<5|G","A|[-C]Rep<5|T","G|[-C]Rep<5|A","G|[-C]Rep<5|G","G|[-C]Rep<5|T","T|[-C]Rep<5|A","T|[-C]Rep<5|G","T|[-C]Rep<5|T",
                    "[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep>=8",
                    "A|[-T]|A","A|[-T]|C","A|[-T]|G","C|[-T]|A","C|[-T]|C","C|[-T]|G","G|[-T]|A","G|[-T]|C","G|[-T]|G",
                    "A|[-T]Rep<5|A","A|[-T]Rep<5|C","A|[-T]Rep<5|G","C|[-T]Rep<5|A","C|[-T]Rep<5|C","C|[-T]Rep<5|G","G|[-T]Rep<5|A","G|[-T]Rep<5|C","G|[-T]Rep<5|G",
                    "[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep>=8",
                    "[-]Mh=1","[-]Mh=2","[-]Mh=3","[-]4<=Mh<20","[-]Mh>=20","[->1]Rep<=4","[->1]Rep>=5","[->1]Other",
                    "Complex") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis, aes(x=indelsubtype, y=freq,fill=indeltype, colour=type))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Count")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  p <- p+scale_color_manual(values=indel_mypalette_outline)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10,colour = "black",hjust=1),
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
  
  return(muts_basis)
}
plotCountbasis_indel_105types_6 <- function(muts_basis,colnum,h,w,outputname){
  #mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_template <-  read.table("/nfs/cancer_archive04/xz3/b_1176/25_indels_20180907/indel_template105.txt",sep = "\t",header = T, as.is = T)
  
  muts_basis_melt <- melt(muts_basis,"indelsubtype")
  
  muts_basis_melt <- merge(muts_template, muts_basis_melt,by="indelsubtype",all.x=T)
  muts_basis_melt[is.na(muts_basis_melt)] <- 0
  names(muts_basis_melt) <- c("indelsubtype","indeltype","type","Sample","freq")
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  indel_mypalette_fill <- c("pink", "deeppink","lightcyan","skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","green","purple")
  indel_mypalette_outline <- c("grey", "white","black")
  # Outline: Insertion: black; Deletion: grey
  # 
  # Fill: [+C]: lightcyan; [+C]Rep<5: skyblue; [+C]Rep>=5: royalblue; [+T]: lightgoldenrod; [+T]Rep<5: goldenrod1; [+T]Rep>=5: darkorgane; 
  #        [+]Mh: deeppink; [+>1]Rep: hotpink; [+>1]Other: pink
  #        Complex:	green
  
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_positions <- c("A|[+C]|A","A|[+C]|G","A|[+C]|T","G|[+C]|A","G|[+C]|G","G|[+C]|T","T|[+C]|A","T|[+C]|G","T|[+C]|T",
                       "A|[+C]Rep<5|A","A|[+C]Rep<5|G","A|[+C]Rep<5|T","G|[+C]Rep<5|A","G|[+C]Rep<5|G","G|[+C]Rep<5|T","T|[+C]Rep<5|A","T|[+C]Rep<5|G","T|[+C]Rep<5|T",
                       "[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep>=8",
                       "A|[+T]|A","A|[+T]|C","A|[+T]|G","C|[+T]|A","C|[+T]|C","C|[+T]|G","G|[+T]|A","G|[+T]|C","G|[+T]|G",
                       "A|[+T]Rep<5|A","A|[+T]Rep<5|C","A|[+T]Rep<5|G","C|[+T]Rep<5|A","C|[+T]Rep<5|C","C|[+T]Rep<5|G","G|[+T]Rep<5|A","G|[+T]Rep<5|C","G|[+T]Rep<5|G",
                       "[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep>=8",
                       "[+]Mh=1","[+]Mh=2","[+]Mh=3","[+]4<=Mh<20","[+]Mh>=20","[+>1]Rep<=4","[+>1]Rep>=5","[+>1]Other",
                       "A|[-C]|A","A|[-C]|G","A|[-C]|T","G|[-C]|A","G|[-C]|G","G|[-C]|T","T|[-C]|A","T|[-C]|G","T|[-C]|T",
                       "A|[-C]Rep<5|A","A|[-C]Rep<5|G","A|[-C]Rep<5|T","G|[-C]Rep<5|A","G|[-C]Rep<5|G","G|[-C]Rep<5|T","T|[-C]Rep<5|A","T|[-C]Rep<5|G","T|[-C]Rep<5|T",
                       "[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep>=8",
                       "A|[-T]|A","A|[-T]|C","A|[-T]|G","C|[-T]|A","C|[-T]|C","C|[-T]|G","G|[-T]|A","G|[-T]|C","G|[-T]|G",
                       "A|[-T]Rep<5|A","A|[-T]Rep<5|C","A|[-T]Rep<5|G","C|[-T]Rep<5|A","C|[-T]Rep<5|C","C|[-T]Rep<5|G","G|[-T]Rep<5|A","G|[-T]Rep<5|C","G|[-T]Rep<5|G",
                       "[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep>=8",
                       "[-]Mh=1","[-]Mh=2","[-]Mh=3","[-]4<=Mh<20","[-]Mh>=20","[->1]Rep<=4","[->1]Rep>=5","[->1]Other",
                       "Complex") 
  
  indel_labels <- c("A|[+C]|A","A|[+C]|G","A|[+C]|T","G|[+C]|A","G|[+C]|G","G|[+C]|T","T|[+C]|A","T|[+C]|G","T|[+C]|T",
                    "A|[+C]Rep<5|A","A|[+C]Rep<5|G","A|[+C]Rep<5|T","G|[+C]Rep<5|A","G|[+C]Rep<5|G","G|[+C]Rep<5|T","T|[+C]Rep<5|A","T|[+C]Rep<5|G","T|[+C]Rep<5|T",
                    "[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep>=8",
                    "A|[+T]|A","A|[+T]|C","A|[+T]|G","C|[+T]|A","C|[+T]|C","C|[+T]|G","G|[+T]|A","G|[+T]|C","G|[+T]|G",
                    "A|[+T]Rep<5|A","A|[+T]Rep<5|C","A|[+T]Rep<5|G","C|[+T]Rep<5|A","C|[+T]Rep<5|C","C|[+T]Rep<5|G","G|[+T]Rep<5|A","G|[+T]Rep<5|C","G|[+T]Rep<5|G",
                    "[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep>=8",
                    "[+]Mh=1","[+]Mh=2","[+]Mh=3","[+]4<=Mh<20","[+]Mh>=20","[+>1]Rep<=4","[+>1]Rep>=5","[+>1]Other",
                    "A|[-C]|A","A|[-C]|G","A|[-C]|T","G|[-C]|A","G|[-C]|G","G|[-C]|T","T|[-C]|A","T|[-C]|G","T|[-C]|T",
                    "A|[-C]Rep<5|A","A|[-C]Rep<5|G","A|[-C]Rep<5|T","G|[-C]Rep<5|A","G|[-C]Rep<5|G","G|[-C]Rep<5|T","T|[-C]Rep<5|A","T|[-C]Rep<5|G","T|[-C]Rep<5|T",
                    "[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep>=8",
                    "A|[-T]|A","A|[-T]|C","A|[-T]|G","C|[-T]|A","C|[-T]|C","C|[-T]|G","G|[-T]|A","G|[-T]|C","G|[-T]|G",
                    "A|[-T]Rep<5|A","A|[-T]Rep<5|C","A|[-T]Rep<5|G","C|[-T]Rep<5|A","C|[-T]Rep<5|C","C|[-T]Rep<5|G","G|[-T]Rep<5|A","G|[-T]Rep<5|C","G|[-T]Rep<5|G",
                    "[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep>=8",
                    "[-]Mh=1","[-]Mh=2","[-]Mh=3","[-]4<=Mh<20","[-]Mh>=20","[->1]Rep<=4","[->1]Rep>=5","[->1]Other",
                    "Complex") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=indelsubtype, y=freq,fill=indeltype, colour=type))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Count")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  p <- p+scale_color_manual(values=indel_mypalette_outline)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10,colour = "black",hjust=1),
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
  
  
  print(p)
  dev.off()
  
 # return(muts_basis)
}

plotCountbasis_aggragated_indel_69types_6 <- function(muts_basis,h,w,outputname){
  #mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_basis$aggragate <- rowSums(muts_basis[,-1])
  
  muts_template <-  read.table("/nfs/cancer_archive04/xz3/b_1176/25_indels_20180907/indel_template105.txt",sep = "\t",header = T, as.is = T)
  
  muts_basis <- merge(muts_template, muts_basis,by="indelsubtype",all.x=T)
  muts_basis[is.na(muts_basis)] <- 0
  
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  indel_mypalette_fill <- c("pink", "deeppink","skyblue","royalblue","goldenrod1","darkorange","green","purple")
  indel_mypalette_outline <- c("grey", "white","black")
  # Outline: Insertion: black; Deletion: grey
  # 
  # Fill: [+C]: lightcyan; [+C]Rep<5: skyblue; [+C]Rep>=5: royalblue; [+T]: lightgoldenrod; [+T]Rep<5: goldenrod1; [+T]Rep>=5: darkorgane; 
  #        [+]Mh: deeppink; [+>1]Rep: hotpink; [+>1]Other: pink
  #        Complex:	green
  
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_positions <- c("A|[+C]Rep<5|A","A|[+C]Rep<5|G","A|[+C]Rep<5|T","G|[+C]Rep<5|A","G|[+C]Rep<5|G","G|[+C]Rep<5|T","T|[+C]Rep<5|A","T|[+C]Rep<5|G","T|[+C]Rep<5|T",
                       "[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep>=8",
                       "A|[+T]Rep<5|A","A|[+T]Rep<5|C","A|[+T]Rep<5|G","C|[+T]Rep<5|A","C|[+T]Rep<5|C","C|[+T]Rep<5|G","G|[+T]Rep<5|A","G|[+T]Rep<5|C","G|[+T]Rep<5|G",
                       "[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep>=8",
                       "[+]Mh=1","[+]Mh=2","[+]Mh=3","[+]4<=Mh<20","[+]Mh>=20","[+>1]Rep<=4","[+>1]Rep>=5","[+>1]Other",
                       "A|[-C]Rep<5|A","A|[-C]Rep<5|G","A|[-C]Rep<5|T","G|[-C]Rep<5|A","G|[-C]Rep<5|G","G|[-C]Rep<5|T","T|[-C]Rep<5|A","T|[-C]Rep<5|G","T|[-C]Rep<5|T",
                       "[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep>=8",
                       "A|[-T]Rep<5|A","A|[-T]Rep<5|C","A|[-T]Rep<5|G","C|[-T]Rep<5|A","C|[-T]Rep<5|C","C|[-T]Rep<5|G","G|[-T]Rep<5|A","G|[-T]Rep<5|C","G|[-T]Rep<5|G",
                       "[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep>=8",
                       "[-]Mh=1","[-]Mh=2","[-]Mh=3","[-]4<=Mh<20","[-]Mh>=20","[->1]Rep<=4","[->1]Rep>=5","[->1]Other",
                       "Complex") 
  
  indel_labels <- c("A|[+C]Rep<5|A","A|[+C]Rep<5|G","A|[+C]Rep<5|T","G|[+C]Rep<5|A","G|[+C]Rep<5|G","G|[+C]Rep<5|T","T|[+C]Rep<5|A","T|[+C]Rep<5|G","T|[+C]Rep<5|T",
                    "[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep>=8",
                    "A|[+T]Rep<5|A","A|[+T]Rep<5|C","A|[+T]Rep<5|G","C|[+T]Rep<5|A","C|[+T]Rep<5|C","C|[+T]Rep<5|G","G|[+T]Rep<5|A","G|[+T]Rep<5|C","G|[+T]Rep<5|G",
                    "[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep>=8",
                    "[+]Mh=1","[+]Mh=2","[+]Mh=3","[+]4<=Mh<20","[+]Mh>=20","[+>1]Rep<=4","[+>1]Rep>=5","[+>1]Other",
                    "A|[-C]Rep<5|A","A|[-C]Rep<5|G","A|[-C]Rep<5|T","G|[-C]Rep<5|A","G|[-C]Rep<5|G","G|[-C]Rep<5|T","T|[-C]Rep<5|A","T|[-C]Rep<5|G","T|[-C]Rep<5|T",
                    "[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep>=8",
                    "A|[-T]Rep<5|A","A|[-T]Rep<5|C","A|[-T]Rep<5|G","C|[-T]Rep<5|A","C|[-T]Rep<5|C","C|[-T]Rep<5|G","G|[-T]Rep<5|A","G|[-T]Rep<5|C","G|[-T]Rep<5|G",
                    "[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep>=8",
                    "[-]Mh=1","[-]Mh=2","[-]Mh=3","[-]4<=Mh<20","[-]Mh>=20","[->1]Rep<=4","[->1]Rep>=5","[->1]Other",
                    "Complex") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis[,c("indelsubtype","indeltype","type","aggragate")], aes(x=indelsubtype, y=aggragate,fill=indeltype, colour=type))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Count")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  p <- p+scale_color_manual(values=indel_mypalette_outline)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10,colour = "black",hjust=1),
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
  
  return(muts_basis)
}
plotCountbasis_single_indel_69types_6 <- function(muts_basis,h,w,outputname){
  #mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_template <-  read.table("/nfs/cancer_archive04/xz3/b_1176/25_indels_20180907/indel_template105.txt",sep = "\t",header = T, as.is = T)
  
  muts_basis <- merge(muts_template, muts_basis,by="indelsubtype",all.x=T)
  muts_basis[is.na(muts_basis)] <- 0
  names(muts_basis) <- c("indelsubtype","indeltype","type","freq")
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  indel_mypalette_fill <- c("pink", "deeppink","skyblue","royalblue","goldenrod1","darkorange","green","purple")
  indel_mypalette_outline <- c("grey", "white","black")
  # Outline: Insertion: black; Deletion: grey
  # 
  # Fill: [+C]: lightcyan; [+C]Rep<5: skyblue; [+C]Rep>=5: royalblue; [+T]: lightgoldenrod; [+T]Rep<5: goldenrod1; [+T]Rep>=5: darkorgane; 
  #        [+]Mh: deeppink; [+>1]Rep: hotpink; [+>1]Other: pink
  #        Complex:	green
  
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_positions <- c("A|[+C]Rep<5|A","A|[+C]Rep<5|G","A|[+C]Rep<5|T","G|[+C]Rep<5|A","G|[+C]Rep<5|G","G|[+C]Rep<5|T","T|[+C]Rep<5|A","T|[+C]Rep<5|G","T|[+C]Rep<5|T",
                       "[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep>=8",
                       "A|[+T]Rep<5|A","A|[+T]Rep<5|C","A|[+T]Rep<5|G","C|[+T]Rep<5|A","C|[+T]Rep<5|C","C|[+T]Rep<5|G","G|[+T]Rep<5|A","G|[+T]Rep<5|C","G|[+T]Rep<5|G",
                       "[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep>=8",
                       "[+]Mh=1","[+]Mh=2","[+]Mh=3","[+]4<=Mh<20","[+]Mh>=20","[+>1]Rep<=4","[+>1]Rep>=5","[+>1]Other",
                       "A|[-C]Rep<5|A","A|[-C]Rep<5|G","A|[-C]Rep<5|T","G|[-C]Rep<5|A","G|[-C]Rep<5|G","G|[-C]Rep<5|T","T|[-C]Rep<5|A","T|[-C]Rep<5|G","T|[-C]Rep<5|T",
                       "[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep>=8",
                       "A|[-T]Rep<5|A","A|[-T]Rep<5|C","A|[-T]Rep<5|G","C|[-T]Rep<5|A","C|[-T]Rep<5|C","C|[-T]Rep<5|G","G|[-T]Rep<5|A","G|[-T]Rep<5|C","G|[-T]Rep<5|G",
                       "[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep>=8",
                       "[-]Mh=1","[-]Mh=2","[-]Mh=3","[-]4<=Mh<20","[-]Mh>=20","[->1]Rep<=4","[->1]Rep>=5","[->1]Other",
                       "Complex") 
  
  indel_labels <- c("A|[+C]Rep<5|A","A|[+C]Rep<5|G","A|[+C]Rep<5|T","G|[+C]Rep<5|A","G|[+C]Rep<5|G","G|[+C]Rep<5|T","T|[+C]Rep<5|A","T|[+C]Rep<5|G","T|[+C]Rep<5|T",
                    "[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep>=8",
                    "A|[+T]Rep<5|A","A|[+T]Rep<5|C","A|[+T]Rep<5|G","C|[+T]Rep<5|A","C|[+T]Rep<5|C","C|[+T]Rep<5|G","G|[+T]Rep<5|A","G|[+T]Rep<5|C","G|[+T]Rep<5|G",
                    "[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep>=8",
                    "[+]Mh=1","[+]Mh=2","[+]Mh=3","[+]4<=Mh<20","[+]Mh>=20","[+>1]Rep<=4","[+>1]Rep>=5","[+>1]Other",
                    "A|[-C]Rep<5|A","A|[-C]Rep<5|G","A|[-C]Rep<5|T","G|[-C]Rep<5|A","G|[-C]Rep<5|G","G|[-C]Rep<5|T","T|[-C]Rep<5|A","T|[-C]Rep<5|G","T|[-C]Rep<5|T",
                    "[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep>=8",
                    "A|[-T]Rep<5|A","A|[-T]Rep<5|C","A|[-T]Rep<5|G","C|[-T]Rep<5|A","C|[-T]Rep<5|C","C|[-T]Rep<5|G","G|[-T]Rep<5|A","G|[-T]Rep<5|C","G|[-T]Rep<5|G",
                    "[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep>=8",
                    "[-]Mh=1","[-]Mh=2","[-]Mh=3","[-]4<=Mh<20","[-]Mh>=20","[->1]Rep<=4","[->1]Rep>=5","[->1]Other",
                    "Complex") 
  
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis, aes(x=indelsubtype, y=freq,fill=indeltype, colour=type))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Count")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  p <- p+scale_color_manual(values=indel_mypalette_outline)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10,colour = "black",hjust=1),
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
  
  return(muts_basis)
}
plotCountbasis_indel_69types_6 <- function(muts_basis,colnum,h,w,outputname){
  #mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_template <-  read.table("/nfs/cancer_archive04/xz3/b_1176/25_indels_20180907/indel_template105.txt",sep = "\t",header = T, as.is = T)
  
  muts_basis_melt <- melt(muts_basis,"indelsubtype")
  
  muts_basis_melt <- merge(muts_template, muts_basis_melt,by="indelsubtype",all.x=T)
  muts_basis_melt[is.na(muts_basis_melt)] <- 0
  names(muts_basis_melt) <- c("indelsubtype","indeltype","type","Sample","freq")
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  indel_mypalette_fill <- c("pink", "deeppink","skyblue","royalblue","goldenrod1","darkorange","green","purple")
  indel_mypalette_outline <- c("grey", "white","black")
  # Outline: Insertion: black; Deletion: grey
  # 
  # Fill: [+C]: lightcyan; [+C]Rep<5: skyblue; [+C]Rep>=5: royalblue; [+T]: lightgoldenrod; [+T]Rep<5: goldenrod1; [+T]Rep>=5: darkorgane; 
  #        [+]Mh: deeppink; [+>1]Rep: hotpink; [+>1]Other: pink
  #        Complex:	green
  
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_positions <- c("A|[+C]Rep<5|A","A|[+C]Rep<5|G","A|[+C]Rep<5|T","G|[+C]Rep<5|A","G|[+C]Rep<5|G","G|[+C]Rep<5|T","T|[+C]Rep<5|A","T|[+C]Rep<5|G","T|[+C]Rep<5|T",
                       "[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep>=8",
                       "A|[+T]Rep<5|A","A|[+T]Rep<5|C","A|[+T]Rep<5|G","C|[+T]Rep<5|A","C|[+T]Rep<5|C","C|[+T]Rep<5|G","G|[+T]Rep<5|A","G|[+T]Rep<5|C","G|[+T]Rep<5|G",
                       "[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep>=8",
                       "[+]Mh=1","[+]Mh=2","[+]Mh=3","[+]4<=Mh<20","[+]Mh>=20","[+>1]Rep<=4","[+>1]Rep>=5","[+>1]Other",
                       "A|[-C]Rep<5|A","A|[-C]Rep<5|G","A|[-C]Rep<5|T","G|[-C]Rep<5|A","G|[-C]Rep<5|G","G|[-C]Rep<5|T","T|[-C]Rep<5|A","T|[-C]Rep<5|G","T|[-C]Rep<5|T",
                       "[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep>=8",
                       "A|[-T]Rep<5|A","A|[-T]Rep<5|C","A|[-T]Rep<5|G","C|[-T]Rep<5|A","C|[-T]Rep<5|C","C|[-T]Rep<5|G","G|[-T]Rep<5|A","G|[-T]Rep<5|C","G|[-T]Rep<5|G",
                       "[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep>=8",
                       "[-]Mh=1","[-]Mh=2","[-]Mh=3","[-]4<=Mh<20","[-]Mh>=20","[->1]Rep<=4","[->1]Rep>=5","[->1]Other",
                       "Complex") 
  
  indel_labels <- c("A|[+C]Rep<5|A","A|[+C]Rep<5|G","A|[+C]Rep<5|T","G|[+C]Rep<5|A","G|[+C]Rep<5|G","G|[+C]Rep<5|T","T|[+C]Rep<5|A","T|[+C]Rep<5|G","T|[+C]Rep<5|T",
                    "[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep>=8",
                    "A|[+T]Rep<5|A","A|[+T]Rep<5|C","A|[+T]Rep<5|G","C|[+T]Rep<5|A","C|[+T]Rep<5|C","C|[+T]Rep<5|G","G|[+T]Rep<5|A","G|[+T]Rep<5|C","G|[+T]Rep<5|G",
                    "[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep>=8",
                    "[+]Mh=1","[+]Mh=2","[+]Mh=3","[+]4<=Mh<20","[+]Mh>=20","[+>1]Rep<=4","[+>1]Rep>=5","[+>1]Other",
                    "A|[-C]Rep<5|A","A|[-C]Rep<5|G","A|[-C]Rep<5|T","G|[-C]Rep<5|A","G|[-C]Rep<5|G","G|[-C]Rep<5|T","T|[-C]Rep<5|A","T|[-C]Rep<5|G","T|[-C]Rep<5|T",
                    "[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep>=8",
                    "A|[-T]Rep<5|A","A|[-T]Rep<5|C","A|[-T]Rep<5|G","C|[-T]Rep<5|A","C|[-T]Rep<5|C","C|[-T]Rep<5|G","G|[-T]Rep<5|A","G|[-T]Rep<5|C","G|[-T]Rep<5|G",
                    "[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep>=8",
                    "[-]Mh=1","[-]Mh=2","[-]Mh=3","[-]4<=Mh<20","[-]Mh>=20","[->1]Rep<=4","[->1]Rep>=5","[->1]Other",
                    "Complex") 
  
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=indelsubtype, y=freq,fill=indeltype, colour=type))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Count")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  p <- p+scale_color_manual(values=indel_mypalette_outline)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10,colour = "black",hjust=1),
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
  
  
  print(p)
  dev.off()
  
  # return(muts_basis)
}

plotCountbasis_aggragated_indel_15types_6 <- function(muts_basis,h,w,outputname){
  #mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_basis$aggragate <- rowSums(muts_basis[,-1])
  
  muts_template <-  read.table("/nfs/cancer_archive04/xz3/b_1176/26_indels_fullscale/indel_template15.txt",sep = "\t",header = T, as.is = T)
  
  muts_basis <- merge(muts_template, muts_basis,by="indelsubtype",all.x=T)
  muts_basis[is.na(muts_basis)] <- 0
  
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  #indel_mypalette_fill <- c("pink", "deeppink","skyblue","royalblue","goldenrod1","darkorange","green","purple")
  # Outline: Insertion: black; Deletion: grey
  # 
  # Fill: [+C]: lightcyan; [+C]Rep<5: skyblue; [+C]Rep>=5: royalblue; [+T]: lightgoldenrod; [+T]Rep<5: goldenrod1; [+T]Rep>=5: darkorgane; 
  #        [+]Mh: deeppink; [+>1]Rep: hotpink; [+>1]Other: pink
  #        Complex:	green
  
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_mypalette_fill <- c("orchid","pink", "lightgreen","skyblue","skyblue","orange","orange",
                            "purple","deeppink","darkgreen","blue","blue","tomato","tomato",
                            "grey")
  indel_positions <- c("[+C]Rep<5","[+C]Rep>=5","[+T]Rep<5","[+T]Rep>=5","[+>1]Rep","[+]Mh","[+>1]Other",
                       "[-C]Rep<5","[-C]Rep>=5","[-T]Rep<5","[-T]Rep>=5","[->1]Rep","[-]Mh","[->1]Other",
                       "Complex") 
  
  indel_labels <- c("[+C]Rep<5","[+C]Rep>=5","[+T]Rep<5","[+T]Rep>=5","[+>1]Rep","[+]Mh","[+>1]Other",
                    "[-C]Rep<5","[-C]Rep>=5","[-T]Rep<5","[-T]Rep>=5","[->1]Rep","[-]Mh","[->1]Other",
                    "Complex") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis[,c("indelsubtype","indeltype","type","aggragate")], aes(x=indelsubtype, y=aggragate,fill=indelsubtype))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Count")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10,colour = "black",hjust=1),
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
  
  return(muts_basis)
}
plotCountbasis_single_indel_15types_6 <- function(muts_basis,h,w,outputname){
  #mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_template <-  read.table("/nfs/cancer_archive04/xz3/b_1176/26_indels_fullscale/indel_template15.txt",sep = "\t",header = T, as.is = T)
  
  muts_basis <- merge(muts_template, muts_basis,by="indelsubtype",all.x=T)
  muts_basis[is.na(muts_basis)] <- 0
  names(muts_basis) <- c("indelsubtype","indeltype","type","freq")
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  #indel_mypalette_fill <- c("pink", "deeppink","skyblue","royalblue","goldenrod1","darkorange","green","purple")
  # Outline: Insertion: black; Deletion: grey
  # 
  # Fill: [+C]: lightcyan; [+C]Rep<5: skyblue; [+C]Rep>=5: royalblue; [+T]: lightgoldenrod; [+T]Rep<5: goldenrod1; [+T]Rep>=5: darkorgane; 
  #        [+]Mh: deeppink; [+>1]Rep: hotpink; [+>1]Other: pink
  #        Complex:	green
  
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_mypalette_fill <- c("orchid","pink", "lightgreen","skyblue","skyblue","orange","orange",
                            "purple","deeppink","darkgreen","blue","blue","tomato","tomato",
                            "grey")
  indel_positions <- c("[+C]Rep<5","[+C]Rep>=5","[+T]Rep<5","[+T]Rep>=5","[+>1]Rep","[+]Mh","[+>1]Other",
                       "[-C]Rep<5","[-C]Rep>=5","[-T]Rep<5","[-T]Rep>=5","[->1]Rep","[-]Mh","[->1]Other",
                       "Complex") 
  
  indel_labels <- c("[+C]Rep<5","[+C]Rep>=5","[+T]Rep<5","[+T]Rep>=5","[+>1]Rep","[+]Mh","[+>1]Other",
                    "[-C]Rep<5","[-C]Rep>=5","[-T]Rep<5","[-T]Rep>=5","[->1]Rep","[-]Mh","[->1]Other",
                    "Complex") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis, aes(x=indelsubtype, y=freq,fill=indelsubtype))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Count")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10,colour = "black",hjust=1),
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
  
  return(muts_basis)
}
plotCountbasis_indel_15types_6 <- function(muts_basis,colnum,h,w,outputname){
  #mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_template <-  read.table("/nfs/cancer_archive04/xz3/b_1176/26_indels_fullscale/indel_template15.txt",sep = "\t",header = T, as.is = T)
  
  muts_basis_melt <- melt(muts_basis,"indelsubtype")
  
  muts_basis_melt <- merge(muts_template, muts_basis_melt,by="indelsubtype",all.x=T)
  muts_basis_melt[is.na(muts_basis_melt)] <- 0
  names(muts_basis_melt) <- c("indelsubtype","indeltype","type","Sample","freq")
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  #indel_mypalette_fill <- c("pink", "deeppink","skyblue","royalblue","goldenrod1","darkorange","green","purple")
  muts_basis_melt$Sample <- as.character(muts_basis_melt$Sample)
  # Outline: Insertion: black; Deletion: grey
  # 
  # Fill: [+C]: lightcyan; [+C]Rep<5: skyblue; [+C]Rep>=5: royalblue; [+T]: lightgoldenrod; [+T]Rep<5: goldenrod1; [+T]Rep>=5: darkorgane; 
  #        [+]Mh: deeppink; [+>1]Rep: hotpink; [+>1]Other: pink
  #        Complex:	green
  
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_mypalette_fill <- c("orchid","pink", "lightgreen","skyblue","skyblue","orange","orange",
                            "purple","deeppink","darkgreen","blue","blue","tomato","tomato",
                            "grey")
  indel_positions <- c("[+C]Rep<5","[+C]Rep>=5","[+T]Rep<5","[+T]Rep>=5","[+>1]Rep","[+]Mh","[+>1]Other",
                       "[-C]Rep<5","[-C]Rep>=5","[-T]Rep<5","[-T]Rep>=5","[->1]Rep","[-]Mh","[->1]Other",
                       "Complex") 
  
  indel_labels <- c("[+C]Rep<5","[+C]Rep>=5","[+T]Rep<5","[+T]Rep>=5","[+>1]Rep","[+]Mh","[+>1]Other",
                    "[-C]Rep<5","[-C]Rep>=5","[-T]Rep<5","[-T]Rep>=5","[->1]Rep","[-]Mh","[->1]Other",
                    "Complex") 
  
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=indelsubtype, y=freq,fill=indelsubtype))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Count")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=5,colour = "black",hjust=1),
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
  
  
  print(p)
  dev.off()
  
  # return(muts_basis)
}
plotbasis_indel_15types_6 <- function(muts_basis,colnum,h,w,outputname){

  muts_basis_percentage <- muts_basis
  muts_basis_percentage[,1:(dim(muts_basis_percentage)[2]-1)] <- muts_basis_percentage[,1:(dim(muts_basis_percentage)[2]-1)]/colSums(muts_basis_percentage[,1:(dim(muts_basis_percentage)[2]-1)])[col(muts_basis_percentage[,1:(dim(muts_basis_percentage)[2]-1)])]
  
  muts_template <-  read.table("/nfs/cancer_archive04/xz3/b_1176/26_indels_fullscale/indel_template15.txt",sep = "\t",header = T, as.is = T)
  
  muts_basis_melt <- melt(muts_basis_percentage,"indelsubtype")
  
  muts_basis_melt <- merge(muts_template, muts_basis_melt,by="indelsubtype",all.x=T)
  muts_basis_melt[is.na(muts_basis_melt)] <- 0
  names(muts_basis_melt) <- c("indelsubtype","indeltype","type","Sample","freq")
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  #indel_mypalette_fill <- c("pink", "deeppink","skyblue","royalblue","goldenrod1","darkorange","green","purple")
  muts_basis_melt$Sample <- as.character(muts_basis_melt$Sample)
  # Outline: Insertion: black; Deletion: grey
  # 
  # Fill: [+C]: lightcyan; [+C]Rep<5: skyblue; [+C]Rep>=5: royalblue; [+T]: lightgoldenrod; [+T]Rep<5: goldenrod1; [+T]Rep>=5: darkorgane; 
  #        [+]Mh: deeppink; [+>1]Rep: hotpink; [+>1]Other: pink
  #        Complex:	green
  
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_mypalette_fill <- c("orchid","pink", "lightgreen","skyblue","skyblue","orange","orange",
                            "purple","deeppink","darkgreen","blue","blue","tomato","tomato",
                            "grey")
  indel_positions <- c("[+C]Rep<5","[+C]Rep>=5","[+T]Rep<5","[+T]Rep>=5","[+>1]Rep","[+]Mh","[+>1]Other",
                       "[-C]Rep<5","[-C]Rep>=5","[-T]Rep<5","[-T]Rep>=5","[->1]Rep","[-]Mh","[->1]Other",
                       "Complex") 
  
  indel_labels <- c("[+C]Rep<5","[+C]Rep>=5","[+T]Rep<5","[+T]Rep>=5","[+>1]Rep","[+]Mh","[+>1]Other",
                    "[-C]Rep<5","[-C]Rep>=5","[-T]Rep<5","[-T]Rep>=5","[->1]Rep","[-]Mh","[->1]Other",
                    "Complex") 
  
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=indelsubtype, y=freq,fill=indelsubtype))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Percentage")
  p <- p+scale_y_continuous(limits=c(0, 0.8),breaks=seq(0, 0.8, 0.2),labels=percent)
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  p <- p+theme(axis.text.x=element_blank(),
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
  
  
  print(p)
  dev.off()
  
  # return(muts_basis)
}
plotbasis_single_indel_15types_6 <- function(muts_basis,h,w,outputname){
  #mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_template <-  read.table("/nfs/cancer_archive04/xz3/b_1176/26_indels_fullscale/indel_template15.txt",sep = "\t",header = T, as.is = T)
  
  muts_basis <- merge(muts_template, muts_basis,by="indelsubtype",all.x=T)
  muts_basis[is.na(muts_basis)] <- 0
  names(muts_basis) <- c("indelsubtype","indeltype","type","freq")
  muts_basis$percentage <- muts_basis$freq/sum(muts_basis$freq)
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  #indel_mypalette_fill <- c("pink", "deeppink","skyblue","royalblue","goldenrod1","darkorange","green","purple")
  #indel_mypalette_outline <- c("grey", "red","black")
  # Outline: Insertion: black; Deletion: grey
  # 
  # Fill: [+C]: lightcyan; [+C]Rep<5: skyblue; [+C]Rep>=5: royalblue; [+T]: lightgoldenrod; [+T]Rep<5: goldenrod1; [+T]Rep>=5: darkorgane; 
  #        [+]Mh: deeppink; [+>1]Rep: hotpink; [+>1]Other: pink
  #        Complex:	green
  
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  # blue +C;  tomato +T;  darkgreen [+>1]Rep;deeppink [+>1]Other; purple [+]Mh
  # skyblue -C; orange -T; lightgreen [->1]Rep ; pink [->1]Other; orchid [-]Mh
  # grey Complex
  
  indel_mypalette_fill <- c("orchid","pink", "lightgreen","skyblue","skyblue","orange","orange",
                            "purple","deeppink","darkgreen","blue","blue","tomato","tomato",
                            "grey")
  indel_positions <- c("[+C]Rep<5","[+C]Rep>=5","[+T]Rep<5","[+T]Rep>=5","[+>1]Rep","[+]Mh","[+>1]Other",
                       "[-C]Rep<5","[-C]Rep>=5","[-T]Rep<5","[-T]Rep>=5","[->1]Rep","[-]Mh","[->1]Other",
                       "Complex") 
  
  indel_labels <- c("[+C]Rep<5","[+C]Rep>=5","[+T]Rep<5","[+T]Rep>=5","[+>1]Rep","[+]Mh","[+>1]Other",
                    "[-C]Rep<5","[-C]Rep>=5","[-T]Rep<5","[-T]Rep>=5","[->1]Rep","[-]Mh","[->1]Other",
                    "Complex") 
  
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis, aes(x=indelsubtype, y=percentage,fill=indelsubtype))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Percentage")
  p <- p+scale_y_continuous(limits=c(0, 0.8),breaks=seq(0, 0.8, 0.2),labels=percent)
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10,colour = "black",hjust=1),
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
  
  return(muts_basis)
}

gen_indelmuttype_new <- function(muts_list){
  indel_catalogue <- data.frame(table(muts_list$Sample,muts_list$Subtype))
  names(indel_catalogue) <- c("subclone","indelsubtype","freq")
  indel_catalogue <- dcast(indel_catalogue,indelsubtype~subclone,value.var="freq")
  indel_catalogue <- merge(indel_template,indel_catalogue,by="indelsubtype",all.x=T)
  indel_catalogue[is.na(indel_catalogue)] <- 0
  return(indel_catalogue)
}
gen_indelmuttype_15 <- function(muts_list){
  indel_template <- read.table("../00_common/indel_template15.txt",sep = "\t",header = T, as.is = T)
  indel_catalogue <- data.frame(table(muts_list$Sample,muts_list$Subtype))
  names(indel_catalogue) <- c("subclone","indelsubtype","freq")
  indel_catalogue <- dcast(indel_catalogue,indelsubtype~subclone,value.var="freq")
  indel_catalogue <- merge(indel_template,indel_catalogue,by="indelsubtype",all.x=T)
  indel_catalogue[is.na(indel_catalogue)] <- 0
  return(indel_catalogue)
}
gen_indelmuttype_sample_15 <- function(muts_list,sample){
  indel_template <- read.table("/nfs/cancer_archive04/xz3/b_1176/26_indels_fullscale/indel_template15.txt",sep = "\t",header = T, as.is = T)
  indel_catalogue <- data.frame(table(muts_list[,sample],muts_list$Subtype))
  names(indel_catalogue) <- c("subclone","indelsubtype","freq")
  indel_catalogue <- dcast(indel_catalogue,indelsubtype~subclone,value.var="freq")
  indel_catalogue <- merge(indel_template,indel_catalogue,by="indelsubtype",all.x=T)
  indel_catalogue[is.na(indel_catalogue)] <- 0
  return(indel_catalogue)
}

gen_indelmuttype_MMRD <- function(muts_list, Sample_col, muttype_col){
  indel_template <- read.table("../00_common/indel_templateMMRD.txt",sep = "\t",header = T, as.is = T)
  indel_template_uniq <- unique(indel_template[,c(muttype_col,"type")])
  names(indel_template_uniq) <- c("indelsubtype","type")
  indel_catalogue <- data.frame(table(muts_list[,Sample_col],muts_list[,muttype_col]))
  names(indel_catalogue) <- c("subclone","indelsubtype","freq")
  indel_catalogue <- dcast(indel_catalogue,indelsubtype~subclone,value.var="freq")
  indel_catalogue <- merge(indel_template_uniq,indel_catalogue,by="indelsubtype",all.x=T)
  indel_catalogue[is.na(indel_catalogue)] <- 0
  return(indel_catalogue)
}
plotCountbasis_indel_45types_6 <- function(muts_basis,colnum,h,w,outputname){
  indel_template <- read.table("../00_common/indel_templateMMRD.txt",sep = "\t",header = T, as.is = T)
  indel_template_uniq <- unique(indel_template[,c("indeltype_short","type")])
  names(indel_template_uniq) <- c("indelsubtype","type")
  
  muts_basis_melt <- melt(muts_basis,c("indelsubtype"))
  
  muts_basis_melt <- merge(indel_template_uniq, muts_basis_melt,by=c("indelsubtype"),all.x=T)
  muts_basis_melt[is.na(muts_basis_melt)] <- 0
  names(muts_basis_melt) <- c("indelsubtype","type","Sample","freq")
  muts_basis_melt$indeltype <- sub("\\=.*","",muts_basis_melt$indelsubtype)
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  muts_basis_melt$Sample <- as.character(muts_basis_melt$Sample)
  
  # 009E73 green, ->1; 56B4E9 light blue -C; E69F00 light orange -T; CC79A7 pink +>1; 0072B2 dark blue +C; D55E00 dark orange +T
  #indel_mypalette_fill <- c("#009E73","#56B4E9","#E69F00","#CC79A7","#0072B2","#D55E00")
  
  
  # blue +C;  tomato +T;  darkgreen [+>1]Rep;deeppink [+>1]Other; purple [+]Mh
  # skyblue -C; orange -T; lightgreen [->1]Rep ; pink [->1]Other; orchid [-]Mh
  # grey Complex
  # indel_mypalette_fill <-c("blue","tomato","darkgreen","deeppink","purple",
  #                        "skyblue","orange","lightgreen","pink","orchid",
  #                       "grey")
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  indel_mypalette_fill <-c("orchid","pink","lightgreen","skyblue","orange",
                           "purple","deeppink","darkgreen","blue","tomato",
                           "grey")
  
  indel_positions <- c("[+C]Rep=0","[+C]Rep=1","[+C]Rep=2","[+C]Rep=3","[+C]Rep=4","[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep=8","[+C]Rep=9",
                       "[+T]Rep=0","[+T]Rep=1","[+T]Rep=2","[+T]Rep=3","[+T]Rep=4","[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep=8","[+T]Rep=9",
                       "[+>1]Rep","[+>1]Other","[+]Mh",
                       "[-C]Rep=1","[-C]Rep=2","[-C]Rep=3","[-C]Rep=4","[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep=8","[-C]Rep=9",
                       "[-T]Rep=1","[-T]Rep=2","[-T]Rep=3","[-T]Rep=4","[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep=8","[-T]Rep=9",
                       "[->1]Rep","[->1]Other","[-]Mh","Complex") 
  
  indel_labels <- c("[+C]Rep=0","[+C]Rep=1","[+C]Rep=2","[+C]Rep=3","[+C]Rep=4","[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep=8","[+C]Rep=9",
                    "[+T]Rep=0","[+T]Rep=1","[+T]Rep=2","[+T]Rep=3","[+T]Rep=4","[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep=8","[+T]Rep=9",
                    "[+>1]Rep","[+>1]Other","[+]Mh",
                    "[-C]Rep=1","[-C]Rep=2","[-C]Rep=3","[-C]Rep=4","[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep=8","[-C]Rep=9",
                    "[-T]Rep=1","[-T]Rep=2","[-T]Rep=3","[-T]Rep=4","[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep=8","[-T]Rep=9",
                    "[->1]Rep","[->1]Other","[-]Mh","Complex") 
  
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=indelsubtype, y=freq,fill=indeltype))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Count")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  #  p <- p+scale_color_manual(values=indel_mypalette_outline)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5, hjust=0.9,colour = "black"),
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
  
  
  print(p)
  dev.off()
  
  # return(muts_basis)
}
plotbasis_indel_45types_6 <- function(muts_basis,colnum,h,w,outputname){
  
  muts_basis_percentage <- muts_basis
  muts_basis_percentage[,2:(dim(muts_basis_percentage)[2])] <- round(muts_basis_percentage[,2:(dim(muts_basis_percentage)[2])]/colSums(muts_basis_percentage[,2:(dim(muts_basis_percentage)[2])])[col(muts_basis_percentage[,2:(dim(muts_basis_percentage)[2])])],2)
  
  
  indel_template <- read.table("/nfs/cancer_archive04/xz3/b_1176/29_indels_final03262019/20_MMRD_sig/indel_templateMMRD.txt",sep = "\t",header = T, as.is = T)
  indel_template_uniq <- unique(indel_template[,c("indeltype_short","type")])
  names(indel_template_uniq) <- c("indelsubtype","type")
  
  muts_basis_melt <- melt(muts_basis_percentage,c("indelsubtype"))
  
  muts_basis_melt <- merge(indel_template_uniq, muts_basis_melt,by=c("indelsubtype"),all.x=T)
  muts_basis_melt[is.na(muts_basis_melt)] <- 0
  names(muts_basis_melt) <- c("indelsubtype","type","Sample","freq")
  muts_basis_melt$indeltype <- sub("\\=.*","",muts_basis_melt$indelsubtype)
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  muts_basis_melt$Sample <- as.character(muts_basis_melt$Sample)
  
  # 009E73 green, ->1; 56B4E9 light blue -C; E69F00 light orange -T; CC79A7 pink +>1; 0072B2 dark blue +C; D55E00 dark orange +T
  #indel_mypalette_fill <- c("#009E73","#56B4E9","#E69F00","#CC79A7","#0072B2","#D55E00")
  
  
  # blue +C;  tomato +T;  darkgreen [+>1]Rep;deeppink [+>1]Other; purple [+]Mh
  # skyblue -C; orange -T; lightgreen [->1]Rep ; pink [->1]Other; orchid [-]Mh
  # grey Complex
  # indel_mypalette_fill <-c("blue","tomato","darkgreen","deeppink","purple",
  #                        "skyblue","orange","lightgreen","pink","orchid",
  #                       "grey")
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  indel_mypalette_fill <-c("orchid","pink","lightgreen","skyblue","orange",
                           "purple","deeppink","darkgreen","blue","tomato",
                           "grey")
  
  indel_positions <- c("[+C]Rep=0","[+C]Rep=1","[+C]Rep=2","[+C]Rep=3","[+C]Rep=4","[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep=8","[+C]Rep=9",
                       "[+T]Rep=0","[+T]Rep=1","[+T]Rep=2","[+T]Rep=3","[+T]Rep=4","[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep=8","[+T]Rep=9",
                       "[+>1]Rep","[+>1]Other","[+]Mh",
                       "[-C]Rep=1","[-C]Rep=2","[-C]Rep=3","[-C]Rep=4","[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep=8","[-C]Rep=9",
                       "[-T]Rep=1","[-T]Rep=2","[-T]Rep=3","[-T]Rep=4","[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep=8","[-T]Rep=9",
                       "[->1]Rep","[->1]Other","[-]Mh","Complex") 
  
  indel_labels <- c("[+C]Rep=0","[+C]Rep=1","[+C]Rep=2","[+C]Rep=3","[+C]Rep=4","[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep=8","[+C]Rep=9",
                    "[+T]Rep=0","[+T]Rep=1","[+T]Rep=2","[+T]Rep=3","[+T]Rep=4","[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep=8","[+T]Rep=9",
                    "[+>1]Rep","[+>1]Other","[+]Mh",
                    "[-C]Rep=1","[-C]Rep=2","[-C]Rep=3","[-C]Rep=4","[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep=8","[-C]Rep=9",
                    "[-T]Rep=1","[-T]Rep=2","[-T]Rep=3","[-T]Rep=4","[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep=8","[-T]Rep=9",
                    "[->1]Rep","[->1]Other","[-]Mh","Complex") 
  
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=indelsubtype, y=freq,fill=indeltype))+ geom_bar(stat="identity",position="dodge", width=.7)+xlab("Indel Types")+ylab("Count")
  p <- p+scale_y_continuous(limits=c(0, 0.3),breaks=seq(0, 0.3, 0.1),labels = scales::percent_format(accuracy=1)) # 
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  #  p <- p+scale_color_manual(values=indel_mypalette_outline)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5, hjust=0.9,colour = "black"),
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
  
  
  print(p)
  dev.off()
  
  # return(muts_basis)
}


plot_sub_indel_profile_sample <- function(sub_catalouge, indel_catalouge,h,w,outputname){
  
  # subs
  muts_basis_melt <- melt(sub_catalouge,"MutationType")
  names(muts_basis_melt) <- c("MutationType","sample","count")
  muts_basis_melt$mutation <- substr(muts_basis_melt$MutationType,3,5)
  
  mutation_order <- sub_catalouge[,c("MutationType","MutationType")]
  names(mutation_order) <- c("MutationType","mutation")
  mutation_order$mutation <- substr(mutation_order$MutationType,3,5)
  mutation_order <- mutation_order[order(mutation_order$mutation),]
  
  # indels
  indel_template <- read.table("/nfs/cancer_archive04/xz3/b_1176/29_indels_final03262019/20_MMRD_sig/indel_templateMMRD.txt",sep = "\t",header = T, as.is = T)
  indel_template_uniq <- unique(indel_template[,c("indeltype_short","type")])
  names(indel_template_uniq) <- c("indelsubtype","type")
  
  indels_basis_melt <- melt(indel_catalouge,c("indelsubtype"))
  
  indels_basis_melt <- merge(indel_template_uniq, indels_basis_melt,by=c("indelsubtype"),all.x=T)
  indels_basis_melt[is.na(indels_basis_melt)] <- 0
  names(indels_basis_melt) <- c("indelsubtype","type","Sample","freq")
  indels_basis_melt$indeltype <- sub("\\=.*","",indels_basis_melt$indelsubtype)
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  indels_basis_melt$Sample <- as.character(indels_basis_melt$Sample)
  
  indel_mypalette_fill <-c("orchid","pink","lightgreen","skyblue","orange",
                           "purple","deeppink","darkgreen","blue","tomato",
                           "grey")
  
  indel_positions <- c("[+C]Rep=0","[+C]Rep=1","[+C]Rep=2","[+C]Rep=3","[+C]Rep=4","[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep=8","[+C]Rep=9",
                       "[+T]Rep=0","[+T]Rep=1","[+T]Rep=2","[+T]Rep=3","[+T]Rep=4","[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep=8","[+T]Rep=9",
                       "[+>1]Rep","[+>1]Other","[+]Mh",
                       "[-C]Rep=1","[-C]Rep=2","[-C]Rep=3","[-C]Rep=4","[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep=8","[-C]Rep=9",
                       "[-T]Rep=1","[-T]Rep=2","[-T]Rep=3","[-T]Rep=4","[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep=8","[-T]Rep=9",
                       "[->1]Rep","[->1]Other","[-]Mh","Complex") 
  
  indel_labels <- c("[+C]Rep=0","[+C]Rep=1","[+C]Rep=2","[+C]Rep=3","[+C]Rep=4","[+C]Rep=5","[+C]Rep=6","[+C]Rep=7","[+C]Rep=8","[+C]Rep=9",
                    "[+T]Rep=0","[+T]Rep=1","[+T]Rep=2","[+T]Rep=3","[+T]Rep=4","[+T]Rep=5","[+T]Rep=6","[+T]Rep=7","[+T]Rep=8","[+T]Rep=9",
                    "[+>1]Rep","[+>1]Other","[+]Mh",
                    "[-C]Rep=1","[-C]Rep=2","[-C]Rep=3","[-C]Rep=4","[-C]Rep=5","[-C]Rep=6","[-C]Rep=7","[-C]Rep=8","[-C]Rep=9",
                    "[-T]Rep=1","[-T]Rep=2","[-T]Rep=3","[-T]Rep=4","[-T]Rep=5","[-T]Rep=6","[-T]Rep=7","[-T]Rep=8","[-T]Rep=9",
                    "[->1]Rep","[->1]Other","[-]Mh","Complex") 
  
  
  filename <- paste0(outputname)
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  p <- ggplot(data=muts_basis_melt, aes(x=MutationType, y=count,fill=mutation))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Mutation Types")+ylab("Count")
  p <- p+scale_x_discrete(limits = as.character(mutation_order$MutationType))+ggtitle(outputname)
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(#axis.text.x=element_blank(),
    axis.text.x=element_text(size=5,angle=90,hjust=0.9,colour = "black"),
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
  
  o <- ggplot(data=indels_basis_melt, aes(x=indelsubtype, y=freq,fill=indeltype))+ geom_bar(stat="identity",position="dodge", width=.8)+xlab("Indel Types")+ylab("Count")
  o <- o+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  o <- o+scale_fill_manual(values=indel_mypalette_fill)
  o <- o+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5, hjust=0.9,colour = "black"),
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




cos_similarity <- function(v1,v2){
  v1v2 <- sum(v1*v2)
  v1_length <- sqrt(sum(v1*v1))
  v2_length <- sqrt(sum(v2*v2))
  return(v1v2/v1_length/v2_length)
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


is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))}

bootstrapGenomesfun <- function(genomes){
  
  return(apply(genomes, 2, function(x) rmultinom(1, sum(x), x)))
}

bootstrapGenomesfun2 <- function(genomes,n){
  
  return(apply(genomes, 2, function(x) rmultinom(1, n, x)))
}


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

# user-define channel number
RemoveBackground_vector_single_channel <- function(background_profile, sig_profile,sampling_number, start_num,boundary,channel_number){
  
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
  exposure <- rep(0,channel_number)
  origArrayIndex <- 1
  for(i in 1:channel_number){
    if(! i %in% mutationTypesToRemoveSet){
      exposure[i] <- diff_all_save[origArrayIndex]
      origArrayIndex=origArrayIndex+1
    }
  }
  
  # Add Weak mutations for background
  background_exposure<- rep(0,channel_number)
  origArrayIndex <- 1
  for(i in 1:channel_number){
    if(! i %in% mutationTypesToRemoveSet){
      background_exposure[i] <- centroid_background[origArrayIndex]
      origArrayIndex=origArrayIndex+1
    }
  }
  
  
  return(data.frame("KO_exposure"=exposure,"background_exposure"=background_exposure))
  
  
}

Wrap_KOSig_indel <- function(MutCatalogue,bg_column,ko_column,sampling_number, start_num,boundary,channel_number,outputname){
  
  KOSig <- RemoveBackground_vector_single_channel(MutCatalogue[,bg_column], MutCatalogue[,ko_column],sampling_number, start_num,boundary,channel_number)
  KOSig$indelsubtype <- MutCatalogue[,"indelsubtype"]
  write.table(KOSig,paste0(outputname,".txt"),sep = "\t",col.names = T, row.names = F, quote = F)
  plotCountbasis_indel_15types_6(KOSig,1,3,5,paste0(outputname,"_count"))
  plotbasis_indel_15types_6(KOSig,1,3,5,paste0(outputname,"_percentage"))
  
}

#########################################
#
# Topography analysis
#
#########################################
regufea_length <- read.table("/nfs/cancer_archive04/xz3/a_1242/00_common/RegulatoryFeatures_length.txt", sep = "\t", header = T, as.is = T)
replitime_length <- read.table("/nfs/cancer_archive04/xz3/a_1242/00_common/MCF7_RepliTime_length.txt", sep = "\t", header = T, as.is = T)
transtrand_length <- read.table("/nfs/cgp_signatures/Analysis/xz3/00_Resource/GenomeFeature/Strandbias_Polynuc_distribution/Transcriptional/polynucleic_transcription.txt", sep = "\t", header = T, as.is = T)
replitrand_length <- read.table("/nfs/cgp_signatures/Analysis/xz3/00_Resource/GenomeFeature/Strandbias_Polynuc_distribution/Replicative/polynucleic_replication.txt", sep = "\t", header = T, as.is = T)
NucleoPos_length <- read.table("/nfs/cancer_archive04/xz3/a_1242/00_common/NucleosomePos_length41.txt", sep = "\t", header = T, as.is = T)

# bed file is 0-based, half-closed-half-open 
Tab2Bed <- function(muts, outputname){ # bed file is 0-based, half-closed-half-open 
  muts_bed <- muts[,c("Chrom","Pos","Pos","VariantID")]
  muts_bed$Chrom <- paste0("chr",muts_bed$Chrom)
  muts_bed[muts_bed$Chrom=="chr23","Chrom"] <- "chrX"
  muts_bed[muts_bed$Chrom=="chr24","Chrom"] <- "chrY"
  muts_bed[,2] <- muts_bed[,2]-1
  write.table(muts_bed,outputname,sep="\t",col.names = F, row.names = F, quote = F)
}
AddStrandInfo_indel <- function(mutfile1, mutfile2,muts_context,outputname){
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
  
  CTsubs$Type <- NULL
  CTsubs[CTsubs$ref.length==1,"Type"] <- "Ins"
  CTsubs[CTsubs$alt.length==1,"Type"] <- "Del"
  CTsubs[CTsubs$alt.length!=1 & CTsubs$ref.length!=1,"Type"] <- "Complex"
  
  CTsubs$change <- NULL
  CTsubs[CTsubs$Type=="Complex","change"] <- substr( as.character(CTsubs[CTsubs$Type=='Complex',"Ref"]),1,1e5)
  CTsubs[CTsubs$Type=="Ins","change"] <- substr( as.character(CTsubs[CTsubs$Type=='Ins',"Alt"]),2,1e5)
  CTsubs[CTsubs$Type=="Del","change"] <- substr( as.character(CTsubs[CTsubs$Type=='Del',"Ref"]),2,1e5)
  
  CTsubs$change2 <- CTsubs$change
  CTsubs[CTsubs$change %in% c("G","A"),]$change2 <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$change %in% c("G","A"),]$change)))
 
  CTsubs[(CTsubs$change2!=CTsubs$change & CTsubs$Strand_original=="leading_uts"),]$Strand <- "lagging_ts"
  CTsubs[(CTsubs$change2!=CTsubs$change & CTsubs$Strand_original=="lagging_ts"),]$Strand <- "leading_uts"
  
  
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
AddStrandInfo_intersect_indel <- function(mutlist, mutBedfile,featureBedfile_strand1,featureBedfile_strand2,intersectResultfile,outputfilename){
  Tab2Bed(mutlist,mutBedfile)
  intersectBed_command_1 <- paste0("intersectBed -a ",mutBedfile," -b ",featureBedfile_strand1," -wo > ", paste0(intersectResultfile,"_1.txt"))
  intersectBed_command_2 <- paste0("intersectBed -a ",mutBedfile," -b ",featureBedfile_strand2," -wo > ", paste0(intersectResultfile,"_2.txt"))
  
  # nohup intersectBed -a denovo_muts_control_mutagen_bed.txt -b ./00_data/TranscribStrand.uts.txt -wo > subs_control_mutagen_uts.txt &
  # nohup intersectBed -a denovo_muts_control_mutagen_bed.txt -b ./00_data/TranscribStrand.ts.txt -wo > subs_control_mutagen_ts.txt &
  try(system(intersectBed_command_1))
  try(system(intersectBed_command_2))
  
  AddStrandInfo_indel(paste0(intersectResultfile,"_1.txt"),paste0(intersectResultfile,"_2.txt"),mutlist,outputfilename)
}

Rintersect <- function(mutlist, mutBedfile,featureBedfile,intersectResultfile){
  Tab2Bed(mutlist,mutBedfile)
  intersectBed_command <- paste0("intersectBed -a ",mutBedfile," -b ",featureBedfile," -wo > ", paste0(intersectResultfile,".txt"))
  
  try(system(intersectBed_command))
}

# Includding simulation on replication timing distribtution of indel Sig
ReplicationTimingIndel_Sample <- function(muts_replitime,SampleCol,ReptimeCol, DoSimulation="TRUE", outputname){
  
  replitime_observed <- data.frame(table(muts_replitime[,SampleCol], muts_replitime[,ReptimeCol]))
  names(replitime_observed) <- c(SampleCol,"ReplicatingTime","observed_Freq")
  replitime_observed <- replitime_observed[replitime_observed[,ReptimeCol] != "others",]
  # simulate the expected distribution 
  if(DoSimulation){
    bootstrap_num <- 100
    muts_replitime[muts_replitime$Chrom=="23","Chrom"]="X"
    muts_replitime[muts_replitime$Chrom=="24","Chrom"]="Y"
    
    muts_replitime$seq_name <- paste0(muts_replitime$change.pyr,"_",muts_replitime$repcount)
    trinuc_catalogue <- dcast(data.frame(table(muts_replitime$seq_name,muts_replitime$Ko_gene)),Var1~Var2)
    names(trinuc_catalogue)[1] <- "seq_name"
    trinuc_replitime <- read.table("./polynucleic_replitime.txt", sep = "\t", header = T, as.is = T)
    
    # According to distribution of trinuc on replication timing regions (trinuc_replitime), bootstrap this distribution
    # bootstrapping method to evaluate the difference between expected distribution and observed one
    # bootstrapping is used to construct a population of expected distributions
    
    
    # loop for sample
    for(i in 4:dim(trinuc_catalogue)[2]){
      current_sample <- trinuc_catalogue[,c(1,i)]
      replitime_expected <- data.frame("ReplicatingTime"=c("1","2","3","4","5","6","7","8","9","10"))
      
      # loop for bootstrap numbers
      for(j in 1:bootstrap_num){
        bootstrap_replitime_all <- NULL
        
        # loop for trinucleotides
        for(k in 1:dim(current_sample)[1]){
          current_context <- as.character(current_sample[k,1])
          trinuc_replitime_current <- trinuc_replitime[trinuc_replitime$seq_name==current_context,]
          bootstrap_replitime=sample(trinuc_replitime_current$ReplicationTime, current_sample[k,2],replace = T, prob = trinuc_replitime_current$sum/sum(trinuc_replitime_current$sum))
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
StrandBias_indel <- function(denovo_muts_feature, SampleCol,feature,featurelength,h,w,colnum, outputname){
  
  #featurelength$all <-rowSums(featurelength[,2:5])
  
  # 1: uts/leading
  # -1: ts/lagging
  #featurelength_strand <- data.frame("Strand"=c("1","-1"),"C"=c(sum(featurelength[,"C"]), sum(featurelength[,"G"])),"T"=c(sum(featurelength[,"T"]), sum(featurelength[,"A"])))
  featurelength_strand <- dcast(featurelength[,c("seq_name","Strand_pyrimidine","sum")],seq_name~Strand_pyrimidine, value.var="sum")
  names(featurelength_strand) <- c("Ref","ts_lagging_wg","uts_leading_wg")
  denovo_muts_feature <- denovo_muts_feature[denovo_muts_feature[,feature] != "others",]
  denovo_muts_feature$Mutation <- substr(denovo_muts_feature$Subtype,3,11)
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
  gtc_dcast$Ref <- paste0(substr(gtc_dcast$Mutation,3,3),"_",substr(gtc_dcast$Mutation,9,9)) 
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
StrandBias_OddsRatio_indel <- function(denovo_muts_feature, SampleCol,feature,featurelength,h,w,colnum, outputname){
  
  #featurelength$all <-rowSums(featurelength[,2:5])
  
  # 1: uts/leading
  # -1: ts/lagging
  #featurelength_strand <- data.frame("Strand"=c("1","-1"),"C"=c(sum(featurelength[,"C"]), sum(featurelength[,"G"])),"T"=c(sum(featurelength[,"T"]), sum(featurelength[,"A"])))
  featurelength_strand <- dcast(featurelength[,c("seq_name","Strand_pyrimidine","sum")],seq_name~Strand_pyrimidine, value.var="sum")
  names(featurelength_strand) <- c("Ref","ts_lagging_wg","uts_leading_wg")
  denovo_muts_feature <- denovo_muts_feature[denovo_muts_feature[,feature] != "others",]
  denovo_muts_feature$Mutation <- substr(denovo_muts_feature$Subtype,3,11)
  denovo_muts_dis <- data.frame(table(denovo_muts_feature[,SampleCol],denovo_muts_feature[,feature], denovo_muts_feature[,"Mutation"]))
  names(denovo_muts_dis) <- c("Sample",feature,"Mutation", "Freq")
  
  gtc <- denovo_muts_dis
  gtc$targetfeature <- gtc[,feature]
  
  gtc_dcast <- dcast(gtc,Sample+Mutation~Strand, value.var="Freq")
  gtc_dcast$Ref <- paste0(substr(gtc_dcast$Mutation,3,3),"_",substr(gtc_dcast$Mutation,9,9)) 
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
Feature_Sample_indel <- function(muts_feature,SampleCol,FeatureCol,trinuc_feature, FeatureName,outputname){
  
  feature_observed <- data.frame(table(muts_feature[,SampleCol], muts_feature[,FeatureCol]))
  names(feature_observed) <- c(SampleCol,"Feature","observed_Freq")
  # simulate the expected distribution 
  
  bootstrap_num <- 100
  muts_feature[muts_feature$Chrom=="23","Chrom"]="X"
  muts_feature[muts_feature$Chrom=="24","Chrom"]="Y"
  
  muts_feature$seq_name <- paste0(muts_feature$change.pyr,"_",muts_feature$repcount)
  trinuc_catalogue <- dcast(data.frame(table(muts_feature$seq_name,muts_feature$Ko_gene)),Var1~Var2)
  names(trinuc_catalogue)[1] <- "seq_name"
  
  #trinuc_catalogue <- Gen32Catalogue(muts_feature,SampleCol)
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
        trinuc_feature_current <- trinuc_feature[trinuc_feature$seq_name==current_context,]
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
# Signature comparison
#
#########################################
Generate_CossimVector_SingleSample_RepIndel <- function(SingleSample_indels, Sig){
  mut_catalogue <-  gen_indelmuttype_MMRD(SingleSample_indels,"Sample","indeltype_short")
  #mut_catalogue[,3] <- mut_catalogue[,3]/sum(mut_catalogue[,3])
  total_subs_sig <- merge(mut_catalogue,Sig)
  
  
  cossim_del <- abs(cos_similarity(total_subs_sig[total_subs_sig$type=="Del",3],total_subs_sig[total_subs_sig$type=="Del",4]))
  cossim_ins <- abs(cos_similarity(total_subs_sig[total_subs_sig$type=="Ins",3],total_subs_sig[total_subs_sig$type=="Ins",4]))
  cossim <- abs(cos_similarity(total_subs_sig[,3],total_subs_sig[,4]))
  
  return(c(cossim_del, cossim_ins,cossim))
}


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

# compare signature within a data set
Cossimi_CompareSig3 <- function(target_sig,h,w,text_size,x_text,outputname){
  
  cossimil <- as.matrix(proxy::simil(as.matrix(t(target_sig)), diag=FALSE, upper=FALSE, method="cosine",auto_convert_data_frames = FALSE))
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
  
  g1 <-g1 +theme(axis.text.x=element_text(angle=90,hjust=0.9,size=x_text,colour = "black"),
                 axis.text.y=element_text(size=x_text,colour = "black"),
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
