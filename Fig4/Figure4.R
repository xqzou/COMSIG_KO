source("../00_common/indel_functions.R")
source("../00_common/sub_functions.R")

##################
#   Figure 4B
# RefSigMMR Tissue signature + CMMRD patients
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

