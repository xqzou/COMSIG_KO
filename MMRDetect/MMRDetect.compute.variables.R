# Calculate variables for MMRDetect


MMRDetect.compute.variables <- function(subs, indels, tissue_type,MMR_subsig96,MMR_sig_indel, tissue_subsig96){
  
  
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
  
  
  # Compute similarity of tissue-specific signatures with MMR KO sigs
  selected_tissueSig <- tissue_subsig96[,c(1,which(sub("_[^_]+$","",names(tissue_subsig96))==tissue_type))]
  
  MMR1Sig <- c("Breast_A", "Colorectal_F", "Liver_E","Stomach_H", "Uterus_C", "Uterus_J")
  MMR2Sig <- c("Biliary_E", "Breast_D", "Colorectal_E", "Ovary_F", "Pancreas_H", "Stomach_B", "Uterus_D", "Kidney_D")
  
  tissue_MMR1Sig <- names(selected_tissueSig)[which(names(selected_tissueSig) %in%MMR1Sig)]
  tissue_MMR2Sig <- names(selected_tissueSig)[which(names(selected_tissueSig) %in%MMR2Sig)]
  
  
  if(length(tissue_MMR1Sig)>0 & length(tissue_MMR2Sig)>0){
    
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }else if(length(tissue_MMR1Sig)>0 & length(tissue_MMR2Sig)==0){
    # add PMS2
    selected_tissueSig <- merge(selected_tissueSig,MMR_subsig96[,c("MutationType","PMS2")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR2")
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }else if(length(tissue_MMR1Sig)==0 & length(tissue_MMR2Sig)>0){
    # add MLH1
    selected_tissueSig <- merge(selected_tissueSig,MMR_subsig96[,c("MutationType","MLH1")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR1")
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
    
  }else if(length(tissue_MMR1Sig)==0 & length(tissue_MMR2Sig)==0){
    
    # add MLH1 and PMS2 
    selected_tissueSig <- merge(selected_tissueSig,MMR_subsig96[,c("MutationType","MLH1","PMS2")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR2")
    names(selected_tissueSig)[dim(selected_tissueSig)[2]-1] <- c("MMR1")
    
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }
  
  a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
  write.table(a$E_median_filtered,paste0("exposure_",tissue_type,".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
  
  
  # Using tissue-specific signatures to fit, but remove the ones that have high similarity with MMR KO signatures.
  #Sig_MMR <- selected_tissueSig[,which(names(selected_tissueSig) %in% ko_pc_cossim[ko_pc_cossim$similarity>=0.85,"CosmicSig"])]
  
  # Choose the samples that have exposure in MMR signatures
  #Exposure <- read.table("exposure_breast741_250breastsamples.txt", sep="\t",header = T, as.is=T)
  Exposure <-a$E_median_filtered
  sample_exposure <- as.data.frame(Exposure)
  MMRsig_sample <- sample_exposure[which(rownames(sample_exposure) %in% c(tissue_MMR1Sig, tissue_MMR2Sig,"MMR2","MMR1")),]
  #MMRsig_sample <- sample_exposure[which(rownames(sample_exposure) %in% ko_pc_cossim_MMR$CosmicSig),]
  # MMRsig_sample <- MMRsig_sample[,which(colSums(MMRsig_sample)>0)]
  
  MMRsig_sample$AndreaSig <- rownames(MMRsig_sample)
  #MMRsig_sample <- merge(MMRsig_sample,ko_pc_cossim_MMR_dcast[,c("CosmicSig","KO_sig")],by="CosmicSig")
  MMRsig_sample_melt <- melt(MMRsig_sample,c("AndreaSig"))
  names(MMRsig_sample_melt) <- c("AndreaSig","Sample","exposure")
  MMRsig_sample_melt_dcast <- dcast(MMRsig_sample_melt,Sample~AndreaSig,value.var="exposure")
  
  if(dim(MMRsig_sample_melt_dcast)[2]==2){
    MMRsig_sample_melt_dcast$MMR_sum <- MMRsig_sample_melt_dcast[,-1]
  }else{
    MMRsig_sample_melt_dcast$MMR_sum <- rowSums(MMRsig_sample_melt_dcast[,-1])
    
  }
  # MMRsig_sample_melt_dcast_realMMR <- MMRsig_sample_melt_dcast[MMRsig_sample_melt_dcast$MMR_sum>=190,]
  #write.table(MMRsig_sample_melt_dcast,"MMRexposure_breast741_175breastsamples.txt",sep = "\t",row.names = F, col.names = T, quote = F)
  
  
  # Step 3. Measure the cosine similarity between RepDel, RepIns and RepIndel
  # Using indel profile to discriminate MLH1/MSH2/MSH6 and PMS2 signatures 
  cossim_allsample <- NULL
  cossim_MMR <- NULL
  cossim_sample <- NULL
  Sample_MMR <- MMRsig_sample_melt_dcast
  for(i in 2:5){
    for(j in 1:dim(Sample_MMR)[1]){
      current_sample_indel <- indel_classied[indel_classied$Sample==as.character(Sample_MMR[j,"Sample"]),]
      if(dim(current_sample_indel)[1]>1){
        cossim <- Generate_CossimVector_SingleSample_RepIndel(current_sample_indel,MMR_sig_indel[,c(1,i)])    
        cossim_allsample <- rbind(cossim_allsample,cossim)
        cossim_MMR <- c(cossim_MMR,names(MMR_sig_indel)[i])
        cossim_sample <- c(cossim_sample,as.character(Sample_MMR[j,"Sample"]))
        
      }
      
    }
    
  }
  cossim_allsample <- data.frame(cossim_allsample)
  names(cossim_allsample) <- c("Del_rep","Ins_rep","Indel_rep")
  cossim_allsample$Sample <- cossim_sample
  cossim_allsample$MMRgene <- cossim_MMR
  
  Del_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Del_rep_mean=mean(Del_rep),Del_rep_sd=sd(Del_rep))
  Ins_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Ins_rep_mean=mean(Ins_rep),Ins_rep_sd=sd(Ins_rep))
  
  Indel_rep_MMR1 <- ddply(cossim_allsample[cossim_allsample$MMRgene!="PMS2",],c("Sample"),summarise,N=length(Sample),Indel_rep_MMR1_mean=mean(Indel_rep),Indel_rep_MMR1_sd=sd(Indel_rep))
  Indel_rep_MMR2 <- cossim_allsample[cossim_allsample$MMRgene=="PMS2",c("Sample","Indel_rep")]
  
  MMRsig_2 <- merge(Del_rep_mean[,-2],Ins_rep_mean[,-2],by="Sample")
  MMRsig_2 <- merge(MMRsig_2,Indel_rep_MMR1[,-2],by="Sample")
  MMRsig_2 <- merge(MMRsig_2,Indel_rep_MMR2,by="Sample")
  MMRsig_2 <- merge(MMRsig_2,indel_classied_rep_summary,by="Sample")
  names(MMRsig_2)[8]="Indel_rep_MMR2"
  
  MMRsig_2 <- merge(MMRsig_2,MMRsig_sample_melt_dcast,by="Sample")
  write.table(MMRsig_2,paste0("sample_feature_summary_",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)
  
  
  return(MMRsig_2) 
}
