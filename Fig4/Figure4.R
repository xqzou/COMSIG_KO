source("../00_common/indel_functions.R")
source("../00_common/sub_functions.R")

##################
#   Figure 4c
#  CMMRD patients subs
##################
# CMMRD patients
muts_patient_devono_all <- read.table("muts_patient_devono_all.txt",sep = "\t",header = T, as.is = T)
muts_patient_devono_all2 <- muts_patient_devono_all
muts_patient_devono_all2$Sample <- muts_patient_devono_all2$Genotype
muts_patient_devono_all2[muts_patient_devono_all2$patient=="MSH77","Sample"] <- "CMMRD_77"
muts_patient_devono_all2[muts_patient_devono_all2$patient=="MSH89","Sample"] <- "CMMRD_89"
muts_patient_devono_all2[muts_patient_devono_all2$patient=="MSH94","Sample"] <- "CMMRD_94"

ko_sigset <- GenCatalogue(muts_patient_devono_all2[muts_patient_devono_all2$celltype=="ips_subclone",],"Sample")
write.table(ko_sigset[,-2],"genotype_catalogue_denovo.txt",sep = "\t",col.names = T, row.names = F, quote = F)
subclones_catalogue <- read.table("genotype_catalogue_denovo.txt",sep = "\t", header = T, as.is = T)
plotCountbasis(subclones_catalogue,1,6,6,"genotype_catalogue_denovo.pdf")

# Remove culture background
mutation_catalogue <- read.table("genotype_catalogue_denovo.txt",sep = "\t", header = T, as.is = T)

Wrap_KOSig(mutation_catalogue,"Control","CMMRD_77",100,150,2,"CMMRD_77_sig")
Wrap_KOSig(mutation_catalogue,"Control","CMMRD_89",100,150,2,"CMMRD_89_sig")
Wrap_KOSig(mutation_catalogue,"Control","CMMRD_94",100,150,2,"CMMRD_94_sig")
Wrap_KOSig(mutation_catalogue,"Control","PMS2",100,150,2,"PMS2_sig")

sig_all <- NULL
sigfilelist <- dir("./sigfiles")
for(i in 1:length(sigfilelist)[1]){
  sig_file <- read.table(paste0("./sigfiles/",sigfilelist[i]), sep = "\t", header = T, as.is = T)
  sig_all <- cbind(sig_all,sig_file[,"KO_exposure"])
  
}
sig_all <- sig_all/colSums(sig_all)[col(sig_all)]
sig_all <- as.data.frame(sig_all)
names(sig_all) <- sub("_sig.txt","",sigfilelist)
sig_all$MutationType <-  sig_file$MutationType

write.table(sig_all,"CMMRD_sig_all.txt",sep = "\t",col.names = T, row.names = F, quote = F)

##################
#   Figure 4d
#  CMMRD patients indels
##################
indel.classified_details <- read.table("./muts_patient_devono_all.txt", sep = "\t", header = T, as.is = T)
indel_templateMMR <- read.table("../00_common/indel_templateMMRD.txt",sep = "\t",header = T, as.is = T)
indel_template2 <- indel_templateMMR
names(indel_template2)[1] <- c("Subtype")
indel.classified_details <- merge(indel.classified_details,indel_template2,by="Subtype")
indel.classified_details <- indel.classified_details[indel.classified_details$patient%in%c("MSH3","MSH77","MSH89","MSH94"),]

ko_sigset <- gen_indelmuttype_MMRD(indel.classified_details[indel.classified_details$celltype=="ips_subclone",],"patient","indeltype_short")
write.table(ko_sigset[,-2],"denovo_Ko_gene_catalogue_indel.txt",sep = "\t",col.names = T, row.names = F, quote = F)

ko_sigset <- read.table("denovo_Ko_gene_catalogue_indel.txt",sep = "\t", header = T, as.is = T)
plotCountbasis_indel_45types_6(ko_sigset,1,10,6,"denovo_Ko_gene_profile_indel")

##################
#   Figure 4e
# RefSigMMR Tissue signature + CMMRD patients
##################
# hcluster
# ko signatures
sig_all <- NULL
sigfilelist <- dir("../Fig2/sigfiles")
for(i in 1:length(sigfilelist)[1]){
  sig_file <- read.table(paste0("../Fig2/sigfiles/",sigfilelist[i]), sep = "\t", header = T, as.is = T)
  sig_all <- cbind(sig_all,sig_file[,"KO_exposure"])
  
}
sig_all <- sig_all/colSums(sig_all)[col(sig_all)]
sig_all <- as.data.frame(sig_all)
names(sig_all) <- sub("\\_.*","",sigfilelist)
sig_all$MutationType <-  sig_file$MutationType

control_sig <- read.table("../Fig2/background_profile.txt",sep = "\t", header = T, as.is = T)
control_sig$Control <- control_sig$mean/sum(control_sig$mean)
sig_total <- merge(sig_all,control_sig[,c("MutationType","Control")],by="MutationType")
sig_total$Mutation <- substr(sig_total$MutationType,3,5)
sig_total <- sig_total[order(sig_total$Mutation),]



# patient signature
sig_cmmrd <- NULL
cmmrd_files <- dir("./sigfiles")
for(i in 1:length(cmmrd_files)[1]){
  sig_file <- read.table(paste0("./sigfiles/",cmmrd_files[i]), sep = "\t", header = T, as.is = T)
  sig_cmmrd <- cbind(sig_cmmrd,sig_file[,"KO_exposure"])
  
}
sig_cmmrd <- sig_cmmrd/colSums(sig_cmmrd)[col(sig_cmmrd)]
sig_cmmrd <- as.data.frame(sig_cmmrd)
names(sig_cmmrd) <- sub("_sig.txt","",cmmrd_files)
sig_cmmrd$MutationType <- sig_file$MutationType
names(sig_cmmrd)[4] <- "CMMRD_3"


PancanSig <- read.table("../00_common/Pancan_signatures_subs_final.txt",sep = "\t",header = T, as.is = T)
names(PancanSig)[1] <- "MutationType"
PancanSig_5 <- PancanSig[,c("MutationType","Biliary_E","Breast_A","Colorectal_E","Colorectal_F",
                            "Kidney_D","Liver_E","Pancreas_G","Stomach_A","Stomach_B","Stomach_H","Uterus_C", "Uterus_J",
                            "Bone_SoftTissue_B","Breast_D","CNS_E","Liver_H","Lymphoid_D","Ovary_F",
                            "Pancreas_H","Skin_H","Uterus_D","Uterus_E")]
kosig_pancansig_all <- merge(PancanSig_5,sig_total[,c("MutationType","MSH2","MSH6","MLH1","PMS2","Mutation")],by="MutationType")
kosig_pancansig_all <- merge(kosig_pancansig_all,sig_cmmrd,by="MutationType")

kosig_pancansig_all <- kosig_pancansig_all[order(kosig_pancansig_all$Mutation),]
kosig_pancansig_all <- kosig_pancansig_all[,-28]
kosig_pancansig_all[,-1] <- kosig_pancansig_all[,-1]/colSums(kosig_pancansig_all[,-1])[col(kosig_pancansig_all[,-1])]
kosig_pancansig_all_2 <- t(kosig_pancansig_all[,-1])

colnames(kosig_pancansig_all_2) <- kosig_pancansig_all$MutationType
kosig_pancansig_all_2 <- data.matrix(kosig_pancansig_all_2)
kosig_pancansig_all_2 <- round(kosig_pancansig_all_2,2)
write.table(kosig_pancansig_all_2,"sigMMR_group_with_cmmrd_patients.txt",sep = "\t",col.names = T, row.names = T, quote = F)


# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white", "blue"))(n = 50)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(0,  # for red
               0.05,           # for yellow
               0.2,max(kosig_pancansig_all_2))             # for green

pdf(file="sigMMR_group_with_cmmrd_patients.pdf", h=8, w=15, onefile=TRUE)
gplots::heatmap.2(kosig_pancansig_all_2,
                  hclustfun=function(x) hclust(x,method = 'complete'),
                  reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean), # Reorder dendrogram by branch means rather than sums
                  main = "Heatmap of substitution signature", # heat map title
                  notecol="black",      # change font color of cell labels to black
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",         # turns off trace lines inside the heat map
                  margins =c(12,9),     # widens margins around plot
                  keysize=1,
                  col=my_palette,       # use on color palette defined earlier
                  key = TRUE, 
                  #      breaks=col_breaks,    # enable color transition at specified limits
                  dendrogram="row",     # only draw a row dendrogram
                  Colv="NA")            # turn off column clustering
dev.off()


##################
#   Figure 4E,F
#   MMRDetect and MSISeq results
##################
source("./Compute.input.variables.R")
source("./find.mono.R")
source("./MSIseq.classify.R")
source("./MSIseq.train.R")
library("RWeka")
library("partykit")
library("Biostrings")
library("C50")
library(ggplot2)
library(scales)

selected_sample <- read.table("../00_common/Pancan_muts_summary_2610_new.txt", sep = "\t", header = T,as.is = T)

########################
#
#   MSIseq EXON
#
########################

pancan_muts_exon <- read.table("../00_common/exon_muts_maf_pancan.txt", sep = "\t", header = T, as.is = T) 
#pancan_muts_exon_2610 <- pancan_muts_exon[pancan_muts_exon$Tumor_Sample_Barcode %in% selected_sample$Sample,]

# Train model 1316
pancan_sample_train <- read.table("../00_common/pancan_sample_classification_train.txt", sep = "\t", header = T, as.is = T)

pancan_muts_train <- pancan_muts_exon[pancan_muts_exon$Tumor_Sample_Barcode %in% pancan_sample_train$Tumor_Sample_Barcode,]
load("Hg19repeats.rda")

pancan.mutationNum.train_all <- NULL
for(i in 1:7){
  print(i)
  if(i*200<=dim(pancan_sample_train)[1]){
    sample_index <- seq(1+(i-1)*200, i*200)
    pancan.mutationNum.train <- Compute.input.variables(pancan_muts_train[pancan_muts_train$Tumor_Sample_Barcode %in% pancan_sample_train[sample_index,"Tumor_Sample_Barcode"],],repeats=Hg19repeats)
    
  }else{
    sample_index <- seq(1+(i-1)*200, dim(pancan_sample_train)[1])
    pancan.mutationNum.train <- Compute.input.variables(pancan_muts_train[pancan_muts_train$Tumor_Sample_Barcode %in% pancan_sample_train[sample_index,"Tumor_Sample_Barcode"],],repeats=Hg19repeats)
    
  }
  pancan.mutationNum.train_all <- rbind(pancan.mutationNum.train_all,pancan.mutationNum.train)
}

#1316
write.table(pancan.mutationNum.train_all,"pancan.mutationNum.train_WES.txt", sep = "\t", col.names = T, row.names = T, quote = F)

pancan.mutationNum.train_all <- read.table("pancan.mutationNum.train_WES.txt",sep="\t" ,header = T, as.is = T)
NGSclassifier_pancan<-MSIseq.train(mutationNum = pancan.mutationNum.train_all,
                                   classification=pancan_sample_train)


pdf(file=paste0("pancan_train_model",".pdf"), onefile=TRUE,width = 8,height = 8)
g <- plot(NGSclassifier_pancan)
print(g)
dev.off()
#write_to_dot(NGSclassifier_breast)


# Test model: 1294 sample
# MSS: 1278
# MSI: 16
pancan_sample_test <- selected_sample[!selected_sample$Sample %in% pancan_sample_train$Tumor_Sample_Barcode,]
pancan_muts_test <- pancan_muts_exon[!pancan_muts_exon$Tumor_Sample_Barcode %in% pancan_sample_train$Tumor_Sample_Barcode,]

pancan.mutationNum.test_all <- NULL
for(i in 1:7){
  print(i)
  if(i*200<=dim(pancan_sample_test)[1]){
    sample_index <- seq(1+(i-1)*200, i*200)
    pancan.mutationNum.test <- Compute.input.variables(pancan_muts_test[pancan_muts_test$Tumor_Sample_Barcode %in% pancan_sample_test[sample_index,"Sample"],],repeats=Hg19repeats)
    
  }else{
    sample_index <- seq(1+(i-1)*200, dim(pancan_sample_train)[1])
    pancan.mutationNum.test <- Compute.input.variables(pancan_muts_test[pancan_muts_test$Tumor_Sample_Barcode %in% pancan_sample_test[sample_index,"Sample"],],repeats=Hg19repeats)
    
  }
  pancan.mutationNum.test_all <- rbind(pancan.mutationNum.test_all,pancan.mutationNum.test)
}

#2194
write.table(pancan.mutationNum.test_all,"pancan.mutationNum.test_WES.txt", sep = "\t", col.names = T, row.names = T, quote = F)

pancan.mutationNum.test_all <- read.table("pancan.mutationNum.test_WES.txt",sep="\t" ,header = T, as.is = T)
Pancan_test_result <- MSIseq.classify(mutationNum = pancan.mutationNum.test_all, classifier =NGSclassifier_pancan)

names(Pancan_test_result)[2] <- "Predicted" 
Pancan_test_result$Real <- "Non-MSI-H"
Pancan_test_result[Pancan_test_result$Tumor_Sample_Barcode %in% MSI_samples,"Real"] <- "MSI-H"

write.table(Pancan_test_result,"MSIseq_Pancan_test_result_WES_retrain.txt", sep = "\t", col.names = T, row.names = F,quote = F)


load("sysdata.rda")
Pancan_test_result_default <- MSIseq.classify(mutationNum = pancan.mutationNum.test_all)
names(Pancan_test_result_default)[2] <- "Predicted" 
Pancan_test_result_default$Real <- "Non-MSI-H"
Pancan_test_result_default[Pancan_test_result_default$Tumor_Sample_Barcode %in% MSI_samples,"Real"] <- "MSI-H"

write.table(Pancan_test_result_default,"MSIseq_Pancan_test_result_WES_default.txt", sep = "\t", col.names = T, row.names = F,quote = F)



########################
#
#   MSIseq WGS
#
########################
pancan_muts_maf_2610 <- read.table("pancan_maf_muts.txt", sep = "\t", header = T, as.is = T) # This dataset is too big to  put in github, it will be avaible on Mendeley
pancan_sample_train <- read.table("pancan_sample_classification_train.txt", sep = "\t", header = T, as.is = T)

load("Hg19repeats.rda")

# Train model 1316

pancan_muts_train <- pancan_muts_maf_2610[pancan_muts_maf_2610$Tumor_Sample_Barcode %in% pancan_sample_train$Tumor_Sample_Barcode,]

pancan.mutationNum.train_all <- NULL
for(i in 1:14){
  print(i)
  if(i*100<=dim(pancan_sample_train)[1]){
    sample_index <- seq(1+(i-1)*100, i*100)
    pancan.mutationNum.train <- Compute.input.variables(pancan_muts_train[pancan_muts_train$Tumor_Sample_Barcode %in% pancan_sample_train[sample_index,"Tumor_Sample_Barcode"],],repeats=Hg19repeats,uniform.seq.len=2861)
    
  }else{
    sample_index <- seq(1+(i-1)*100, dim(pancan_sample_train)[1])
    pancan.mutationNum.train <- Compute.input.variables(pancan_muts_train[pancan_muts_train$Tumor_Sample_Barcode %in% pancan_sample_train[sample_index,"Tumor_Sample_Barcode"],],repeats=Hg19repeats,uniform.seq.len=2861)
    
  }
  pancan.mutationNum.train_all <- rbind(pancan.mutationNum.train_all,pancan.mutationNum.train)
}

#1316
write.table(pancan.mutationNum.train_all,"pancan.mutationNum.train_all_WGS_new.txt", sep = "\t", col.names = T, row.names = T, quote = F)

pancan.mutationNum.train_all <- read.table("./pancan.mutationNum.train_all_WGS_new.txt",sep="\t" ,header = T, as.is = T)
NGSclassifier_pancan<-MSIseq.train(mutationNum = pancan.mutationNum.train_all,
                                   classification=pancan_sample_train)


pdf(file=paste0("pancan_train_model_MSIseq_WGS",".pdf"), onefile=TRUE,width = 8,height = 8)
g <- plot(NGSclassifier_pancan)
print(g)
dev.off()
#write_to_dot(NGSclassifier_breast)


# Test model: 1294 sample
# MSS: 1278
# MSI: 16
pancan_sample_test <- selected_sample[!selected_sample$Sample %in% pancan_sample_train$Tumor_Sample_Barcode,]
pancan_muts_test <- pancan_muts_maf_2610[!pancan_muts_maf_2610$Tumor_Sample_Barcode %in% pancan_sample_train$Tumor_Sample_Barcode,]

pancan.mutationNum.test_all <- NULL
for(i in 1:13){
  print(i)
  
  if(i*100<=dim(pancan_sample_test)[1]){
    sample_index <- seq(1+(i-1)*100, i*100)
    pancan.mutationNum.test <- Compute.input.variables(pancan_muts_test[pancan_muts_test$Tumor_Sample_Barcode %in% pancan_sample_test[sample_index,"Sample"],],repeats=Hg19repeats,uniform.seq.len=2861)
    
  }else{
    sample_index <- seq(1+(i-1)*100, dim(pancan_sample_train)[1])
    pancan.mutationNum.test <- Compute.input.variables(pancan_muts_test[pancan_muts_test$Tumor_Sample_Barcode %in% pancan_sample_test[sample_index,"Sample"],],repeats=Hg19repeats,uniform.seq.len=2861)
    
  }
  pancan.mutationNum.test_all <- rbind(pancan.mutationNum.test_all,pancan.mutationNum.test)
}

#2194
write.table(pancan.mutationNum.test_all,"pancan.mutationNum.test_all_WGS_new.txt", sep = "\t", col.names = T, row.names = T, quote = F)

pancan.mutationNum.test_all <- read.table("pancan.mutationNum.test_all_WGS_new.txt",sep="\t" ,header = T, as.is = T)
Pancan_test_result <- MSIseq.classify(mutationNum = pancan.mutationNum.test_all, classifier =NGSclassifier_pancan)

names(Pancan_test_result)[2] <- "Predicted" 
Pancan_test_result$Real <- "Non-MSI-H"
Pancan_test_result[Pancan_test_result$Tumor_Sample_Barcode %in% MSI_samples,"Real"] <- "MSI-H"

write.table(Pancan_test_result,"MSIseq_Pancan_test_result_WGS_retrain.txt", sep = "\t", col.names = T, row.names = F,quote = F)


load("sysdata.rda")
Pancan_test_result_default <- MSIseq.classify(mutationNum = pancan.mutationNum.test_all)
names(Pancan_test_result_default)[2] <- "Predicted" 
Pancan_test_result_default$Real <- "Non-MSI-H"
Pancan_test_result_default[Pancan_test_result_default$Tumor_Sample_Barcode %in% MSI_samples,"Real"] <- "MSI-H"

write.table(Pancan_test_result_default,"MSIseq_Pancan_test_result_default_WGS.txt", sep = "\t", col.names = T, row.names = F,quote = F)







########################
#
#   Generate variables for MMRDetect
#
########################

dir.create("./MMRDetect_variables")
setwd("./MMRDetect_variables/")
# Used Andrea's signature to fit,but pickup real tissue-specific MMR signatures and short indel signatures
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

PancanSig <- read.table("../../00_common/Pancan_signatures_subs_final.txt",sep = "\t",header = T, as.is = T)
names(PancanSig)[1] <- "MutationType"
#tissue_subsig96 <- PancanSig

MMR_indelsig <- read.table("../../00_common/MMRD_indel_sig.txt", sep = "\t", header = T, as.is = T)
MMR_subsig <- read.table("../../00_common/MMR_sig_96.txt",sep = "\t", header = T, as.is = T)

Pancan_muts_summary_2610 <- read.table("../../00_common/Pancan_muts_summary_2610_info_summary.txt",sep = "\t", header = T, as.is = T)


organ_list <- data.frame(table(Pancan_muts_summary_2610$Organ_system))
for(i in 4:21){
  print(i)
  organ <- as.character(organ_list[i,1])
  tissue_subs <- all_subs[all_subs$Sample %in% Pancan_muts_summary_2610[Pancan_muts_summary_2610$Organ_system==organ,"Sample"],]
  tissue_indels <- all_indels[all_indels$Sample %in% Pancan_muts_summary_2610[Pancan_muts_summary_2610$Organ_system==organ,"Sample"],]
  
  muts_variables <- MMRDetect.compute.variables(tissue_subs, tissue_indels, organ, MMR_subsig,MMR_indelsig,PancanSig)
  
}


########################
#
#   MMRDetect
#
########################

MMRDetect_var <- dir("./Pancan_MMRDetect_variables")

MMRDetect_var_all <- NULL
for(i in 1:length(MMRDetect_var)){
  MMRDetect_var_i <- read.table(paste0("./Pancan_MMRDetect_variables/",MMRDetect_var[i]), sep = "\t", header = T, as.is = T)
  MMRDetect_var_i <- MMRDetect_var_i[,c("Sample","Del_rep_mean","Ins_rep_mean","Indel_rep_MMR1_mean","Indel_rep_MMR2","RepIndel_num","MMR_sum")]
  MMRDetect_var_all <- rbind(MMRDetect_var_all,MMRDetect_var_i)
}
write.table(MMRDetect_var_all,"MMRDetect_var_all.txt", sep = "\t", col.names = T, row.names = F,quote = F) # 2448

# Train
MMRDetect_var_all <- read.table("MMRDetect_var_all.txt", sep = "\t", header = T, as.is = T)
pancan_sample_train <- read.table("pancan_sample_classification_train.txt", sep = "\t", header = T, as.is = T)
names(pancan_sample_train) <- c("Sample","MSI_status")
MMRDetect_var_train <- merge(pancan_sample_train, MMRDetect_var_all, by="Sample", all.x=T)
MMRDetect_var_train[is.na(MMRDetect_var_train)] <- 0 # 1316
MMRDetect_var_train <- MMRDetect_var_train[,-2]
MMRDetect.train <- function(mutationVariable, classification, cancerType = NULL) {
  
  
  ## match the data with classification
  trainset = mutationVariable
  trainset = merge(trainset, classification[,c("Sample","MSI_status")], by="Sample")
  
  if(nrow(trainset)<50){
    warning('training set size < 50')  
  }
  
  ## build model with trainset
  trainset$MSI_status<-as.factor(trainset$MSI_status)
  tree.model = J48(MSI_status ~ . , data=trainset[,c("Del_rep_mean","Ins_rep_mean","MMR_sum","MSI_status")])
  
  tree.eval = evaluate_Weka_classifier(tree.model, numFolds = 5, seed = 1)
  cat('5 fold cross validation result: ',unname(tree.eval$detail[1]), '\n', sep='')
  tree.model
  
}

MMRDetectclassifier_pancan<-MMRDetect.train(MMRDetect_var_train, pancan_sample_train)
pdf(file=paste0("pancan_train_model_C48",".pdf"), onefile=TRUE,width = 8,height = 8)
g <- plot(MMRDetectclassifier_pancan)
print(g)
dev.off()



# Test
selected_sample <- read.table("Pancan_muts_summary_2610_new.txt", sep = "\t", header = T,as.is = T) # 2614
pancan_sample_test <- selected_sample[!selected_sample$Sample %in% pancan_sample_train$Sample,] # 1297
MMRDetect_var_test <- merge(pancan_sample_test, MMRDetect_var_all, by="Sample", all.x=T)
MMRDetect_var_test[is.na(MMRDetect_var_test)] <- 0 # 1294
row.names(MMRDetect_var_test) <- MMRDetect_var_test$Sample
MMRDetect.classify <- function(mutationVariable, classifier = MMRDclassifier) {
  classifyset = mutationVariable[,c("Del_rep_mean","Ins_rep_mean","MMR_sum")]
  
  
  classify.result = predict(classifier, newdata=classifyset, type="class")
  result = as.data.frame(cbind(rownames(classifyset), 
                               as.character(classify.result)))
  colnames(result) = c("Sample", "MSI_status")
  result
} 

MMRDetect_Pancan_test_result <- MMRDetect.classify(mutationVariable = MMRDetect_var_test, classifier =MMRDetectclassifier_pancan)
names(MMRDetect_Pancan_test_result)[2] <- "Predicted" 
MMRDetect_Pancan_test_result$Real <- "Non-MSI-H"
MMRDetect_Pancan_test_result[MMRDetect_Pancan_test_result$Sample %in% MSI_samples,"Real"] <- "MSI-H"

write.table(MMRDetect_Pancan_test_result,"MMRDetect_Pancan_test_result.txt", sep = "\t", col.names = T, row.names = F,quote = F)



MMRDetect_var_all_classify <- MMRDetect_var_all
MMRDetect_var_all_classify$MSI_status <- "Non-MSI-H"
MMRDetect_var_all_classify[MMRDetect_var_all_classify$Sample %in% MSI_samples,"MSI_status"] <- "MSI-H"
MMRDetect_var_all_classify$Group <- "Test"
MMRDetect_var_all_classify[MMRDetect_var_all_classify$Sample %in% pancan_sample_train$Sample,"Group"] <- "Train"

write.table(MMRDetect_var_all_classify,"MMRDetect_var_all_classify.txt", sep = "\t", col.names = T, row.names = F,quote = F)



