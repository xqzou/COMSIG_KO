library(ggplot2)
library(plyr) 
library(reshape2)
library(gridExtra)
library(scales)
library(VennDiagram) # use for up to 3 sets
library(Rtsne)
#library(factoextra)
library("VariantAnnotation")
library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)



##########################
#
#  Indel classification
#
##########################
prepare.indel.df.tab <- function(indel.data) {
  
  
  if (nrow(indel.data)>0) {
    
    
    indel.data$ref.length <- nchar(indel.data$Ref)
    indel.data$alt.length <- nchar(indel.data$Alt)
    
    indel.data <- indel.data[(indel.data$ref.length+indel.data$alt.length)>2,]
    
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

gen_indelmuttype_15 <- function(muts_list){
  indel_template <- read.table("../00_common/indel_template15.txt",sep = "\t",header = T, as.is = T)
  indel_catalogue <- data.frame(table(muts_list$Sample,muts_list$Subtype))
  names(indel_catalogue) <- c("subclone","indelsubtype","freq")
  indel_catalogue <- dcast(indel_catalogue,indelsubtype~subclone,value.var="freq")
  indel_catalogue <- merge(indel_template,indel_catalogue,by="indelsubtype",all.x=T)
  indel_catalogue[is.na(indel_catalogue)] <- 0
  return(indel_catalogue)
}

indel_classifier <- function(indels){
  
  indel.data <- indels
  indel.data[indel.data$Chrom=="23","Chrom"]="X"
  indel.data[indel.data$Chrom=="24","Chrom"]="Y"
  
  # convert formats, and find context of the indels
  indel.df <- prepare.indel.df.tab(indel.data)
  indel.df.max100 <- indel.df[indel.df$indel.length<=100 & indel.df$indel.length>0,]
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

plot_percentage_type_4_mmrd <- function(muts_basis,colnum,h,w,outputname){
  muts_basis[,-which(colnames(muts_basis)=="type_4")] <- muts_basis[,-which(colnames(muts_basis)=="type_4")]/colSums(muts_basis[,-which(colnames(muts_basis)=="type_4")])[col(muts_basis[,-which(colnames(muts_basis)=="type_4")])]
  indel_template_type_4 <- data.frame("type_4"=c("[+C]NonRep","[+C]ShortRep_leq4","[+C]LongRep_g4",
                                                 
                                                 "A|[+T]Rep=0|A","A|[+T]Rep=0|C","A|[+T]Rep=0|G","C|[+T]Rep=0|A","C|[+T]Rep=0|C","C|[+T]Rep=0|G","G|[+T]Rep=0|A","G|[+T]Rep=0|C","G|[+T]Rep=0|G",
                                                 "A|[+T]Rep=1|A","A|[+T]Rep=1|C","A|[+T]Rep=1|G","C|[+T]Rep=1|A","C|[+T]Rep=1|C","C|[+T]Rep=1|G","G|[+T]Rep=1|A","G|[+T]Rep=1|C","G|[+T]Rep=1|G",
                                                 "A|[+T]Rep=2|A","A|[+T]Rep=2|C","A|[+T]Rep=2|G","C|[+T]Rep=2|A","C|[+T]Rep=2|C","C|[+T]Rep=2|G","G|[+T]Rep=2|A","G|[+T]Rep=2|C","G|[+T]Rep=2|G",
                                                 "A|[+T]Rep=3|A","A|[+T]Rep=3|C","A|[+T]Rep=3|G","C|[+T]Rep=3|A","C|[+T]Rep=3|C","C|[+T]Rep=3|G","G|[+T]Rep=3|A","G|[+T]Rep=3|C","G|[+T]Rep=3|G",
                                                 "A|[+T]Rep=4|A","A|[+T]Rep=4|C","A|[+T]Rep=4|G","C|[+T]Rep=4|A","C|[+T]Rep=4|C","C|[+T]Rep=4|G","G|[+T]Rep=4|A","G|[+T]Rep=4|C","G|[+T]Rep=4|G",
                                                 "A|[+T]Rep=5|A","A|[+T]Rep=5|C","A|[+T]Rep=5|G","C|[+T]Rep=5|A","C|[+T]Rep=5|C","C|[+T]Rep=5|G","G|[+T]Rep=5|A","G|[+T]Rep=5|C","G|[+T]Rep=5|G",
                                                 "A|[+T]Rep=6|A","A|[+T]Rep=6|C","A|[+T]Rep=6|G","C|[+T]Rep=6|A","C|[+T]Rep=6|C","C|[+T]Rep=6|G","G|[+T]Rep=6|A","G|[+T]Rep=6|C","G|[+T]Rep=6|G",
                                                 "A|[+T]Rep=7|A","A|[+T]Rep=7|C","A|[+T]Rep=7|G","C|[+T]Rep=7|A","C|[+T]Rep=7|C","C|[+T]Rep=7|G","G|[+T]Rep=7|A","G|[+T]Rep=7|C","G|[+T]Rep=7|G",
                                                 "A|[+T]Rep=8|A","A|[+T]Rep=8|C","A|[+T]Rep=8|G","C|[+T]Rep=8|A","C|[+T]Rep=8|C","C|[+T]Rep=8|G","G|[+T]Rep=8|A","G|[+T]Rep=8|C","G|[+T]Rep=8|G",
                                                 "A|[+T]Rep=9|A","A|[+T]Rep=9|C","A|[+T]Rep=9|G","C|[+T]Rep=9|A","C|[+T]Rep=9|C","C|[+T]Rep=9|G","G|[+T]Rep=9|A","G|[+T]Rep=9|C","G|[+T]Rep=9|G",
                                                 
                                                 "Ins_NonRep","Ins_nMer_ShortRep_leq4","Ins_nMer_LongRep_g4",
                                                 
                                                 "[-C]NonRep","[-C]ShortRep_leq4","[-C]LongRep_g4",
                                                 
                                                 "A|[-T]Rep=1|A","A|[-T]Rep=1|C","A|[-T]Rep=1|G","C|[-T]Rep=1|A","C|[-T]Rep=1|C","C|[-T]Rep=1|G","G|[-T]Rep=1|A","G|[-T]Rep=1|C","G|[-T]Rep=1|G",
                                                 "A|[-T]Rep=2|A","A|[-T]Rep=2|C","A|[-T]Rep=2|G","C|[-T]Rep=2|A","C|[-T]Rep=2|C","C|[-T]Rep=2|G","G|[-T]Rep=2|A","G|[-T]Rep=2|C","G|[-T]Rep=2|G",
                                                 "A|[-T]Rep=3|A","A|[-T]Rep=3|C","A|[-T]Rep=3|G","C|[-T]Rep=3|A","C|[-T]Rep=3|C","C|[-T]Rep=3|G","G|[-T]Rep=3|A","G|[-T]Rep=3|C","G|[-T]Rep=3|G",
                                                 "A|[-T]Rep=4|A","A|[-T]Rep=4|C","A|[-T]Rep=4|G","C|[-T]Rep=4|A","C|[-T]Rep=4|C","C|[-T]Rep=4|G","G|[-T]Rep=4|A","G|[-T]Rep=4|C","G|[-T]Rep=4|G",
                                                 "A|[-T]Rep=5|A","A|[-T]Rep=5|C","A|[-T]Rep=5|G","C|[-T]Rep=5|A","C|[-T]Rep=5|C","C|[-T]Rep=5|G","G|[-T]Rep=5|A","G|[-T]Rep=5|C","G|[-T]Rep=5|G",
                                                 "A|[-T]Rep=6|A","A|[-T]Rep=6|C","A|[-T]Rep=6|G","C|[-T]Rep=6|A","C|[-T]Rep=6|C","C|[-T]Rep=6|G","G|[-T]Rep=6|A","G|[-T]Rep=6|C","G|[-T]Rep=6|G",
                                                 "A|[-T]Rep=7|A","A|[-T]Rep=7|C","A|[-T]Rep=7|G","C|[-T]Rep=7|A","C|[-T]Rep=7|C","C|[-T]Rep=7|G","G|[-T]Rep=7|A","G|[-T]Rep=7|C","G|[-T]Rep=7|G",
                                                 "A|[-T]Rep=8|A","A|[-T]Rep=8|C","A|[-T]Rep=8|G","C|[-T]Rep=8|A","C|[-T]Rep=8|C","C|[-T]Rep=8|G","G|[-T]Rep=8|A","G|[-T]Rep=8|C","G|[-T]Rep=8|G",
                                                 "A|[-T]Rep=9|A","A|[-T]Rep=9|C","A|[-T]Rep=9|G","C|[-T]Rep=9|A","C|[-T]Rep=9|C","C|[-T]Rep=9|G","G|[-T]Rep=9|A","G|[-T]Rep=9|C","G|[-T]Rep=9|G",
                                                 
                                                 "Del_NonRep","Del_nMer_ShortRep_leq4","Del_nMer_LongRep_g4","Del_Spaced_short_leq5","Del_Spaced_long_g5",
                                                 "Complex"),
                                      
                                      "type_2"=c("[+C]","[+C]","[+C]",
                                                 
                                                 "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",
                                                 "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",
                                                 "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",
                                                 "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",
                                                 "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",
                                                 "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",
                                                 "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",
                                                 "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",
                                                 "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",
                                                 "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",
                                                 
                                                 "Ins_NonRep",
                                                 "Ins_nMer","Ins_nMer",
                                                 
                                                 "[-C]","[-C]","[-C]",
                                                 
                                                 "[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]",
                                                 "[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]",
                                                 "[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]",
                                                 "[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]",
                                                 "[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]",
                                                 "[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]",
                                                 "[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]",
                                                 "[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]",
                                                 "[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]",
                                                 
                                                 "Del_NonRep",
                                                 "Del_nMer","Del_nMer",
                                                 "Del_Spaced","Del_Spaced",
                                                 "Complex")
  )
  
  muts_basis_melt <- melt(muts_basis,"type_4")
  
  muts_basis_melt <- merge(indel_template_type_4, muts_basis_melt,by="type_4",all.x=T)
  muts_basis_melt[is.na(muts_basis_melt)] <- 0
  names(muts_basis_melt) <- c("type_4","type_2","Sample","freq")
  #indel_mypalette_fill <- c("lightcyan", "skyblue","royalblue","lightgoldenrod","goldenrod1","darkorange","deeppink","hotpink","pink","green")
  #indel_mypalette_fill <- c("pink", "deeppink","skyblue","royalblue","goldenrod1","darkorange","green","purple")
  muts_basis_melt$Sample <- as.character(muts_basis_melt$Sample)
  # Fill: [+C]: lightcyan; [+C]Rep<5: skyblue; [+C]Rep>=5: royalblue; [+T]: lightgoldenrod; [+T]Rep<5: goldenrod1; [+T]Rep>=5: darkorgane; 
  #        [+]Mh: deeppink; [+>1]Rep: hotpink; [+>1]Other: pink
  #        Complex:	green
  
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_mypalette_fill <- c("skyblue","orange", "blue","tomato","greenyellow","pink",
                            "grey","purple","darkred","black")
  
  
  indel_positions <- indel_template_type_4$type_4
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis_melt, aes(x=type_4, y=freq,fill=type_2))+ geom_bar(stat="identity",position="dodge", width=.5)+xlab("Indel Types")+ylab("Count")
  p <- p+scale_y_continuous(labels=percent)
  p <- p+scale_x_discrete(limits = indel_positions)+ggtitle(outputname)
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


Generate_CossimVector_SingleSample_RepIndel2 <- function(SingleSample_indels, Sig){
  mut_catalogue <-  gen_indelmuttype_MMRD2(SingleSample_indels,"Sample","indeltype_medium")
  total_subs_sig <- merge(mut_catalogue,Sig)
  
  
  cossim_del <- abs(cos_similarity(total_subs_sig[total_subs_sig$type %in% c("[-C]","[-T]","[->1]Rep"),3],total_subs_sig[total_subs_sig$type %in% c("[-C]","[-T]","[->1]Rep"),4]))
  cossim_ins <- abs(cos_similarity(total_subs_sig[total_subs_sig$type%in% c("[+C]","[+T]","[+>1]Rep"),3],total_subs_sig[total_subs_sig$type%in% c("[+C]","[+T]","[+>1]Rep"),4]))
  
  return(c(cossim_del, cossim_ins))
}
