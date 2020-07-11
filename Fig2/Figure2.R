source("../00_common/indel_functions.R")
source("../00_common/sub_functions.R")

##################
#   Figure 2A
##################
subs <- read.table("../00_common/total_subs_43genes.txt", sep = "\t", header = T, as.is = T)
mutation_catalogue <- GenCatalogue(subs,"Sample_ko")

control_catalogue <- read.table("background_profile.txt",sep = "\t", header = T, as.is = T)
mutation_catalogue <- merge(mutation_catalogue,control_catalogue[,c("MutationType","mean")],by="MutationType")
plotPercentagebasis(mutation_catalogue[,c("MutationType","mean","mean")],1,6,9,paste0("background",".pdf"))
Wrap_KOSig(mutation_catalogue,"mean",c("RNF168_116_s1","RNF168_116_s2","RNF168_14_s1","RNF168_14_s2"),100,150,2,"RNF168_116_14_sig")
Wrap_KOSig(mutation_catalogue,"mean",c("EXO1_71_s2","EXO1_71_s3","EXO1_71_s4"),100,150,2,"EXO1_71_sig")
Wrap_KOSig(mutation_catalogue,"mean",c("OGG1_106_s1","OGG1_106_s2","OGG1_25_s1","OGG1_25_s2"),100,150,2,"OGG1_106_25_sig")
Wrap_KOSig(mutation_catalogue,"mean",c("UNG_6_s3","UNG_6_s4","UNG_6_s5","UNG_6_s6"),100,150,2,"UNG_0.6_sig")

######################################################
#   Figure 2B
#   Cosine similarity between Andrea's reference sig
######################################################


sig_all <- NULL
sigfilelist <- dir("./sigfiles")
for(i in 1:length(sigfilelist)[1]){
  sig_file <- read.table(paste0("./sigfiles/",sigfilelist[i]), sep = "\t", header = T, as.is = T)
  sig_all <- cbind(sig_all,sig_file[,"KO_exposure"])
  
}
sig_all <- sig_all/colSums(sig_all)[col(sig_all)]
sig_all <- as.data.frame(sig_all)
names(sig_all) <- sub("\\_.*","",sigfilelist)
sig_all$MutationType <-  sig_file$MutationType

control_sig <- read.table("./background_profile.txt",sep = "\t", header = T, as.is = T)
control_sig$Control <- control_sig$mean/sum(control_sig$mean)
sig_total <- merge(sig_all,control_sig[,c("MutationType","Control")],by="MutationType")


cosmic_sig <- read.table("../00_common/Andrea_RefSig.txt",sep = "\t",header = T, as.is = T)
names(cosmic_sig)[1] <- "MutationType"
kosig_cosmicsig <- merge(cosmic_sig,sig_total,by="MutationType")

cosmicsig <- kosig_cosmicsig[,2:42]
kosig <- kosig_cosmicsig[,c("Control","OGG1","UNG","RNF168","EXO1")]

cossimil <- as.matrix(proxy::simil(as.matrix(t(kosig)), as.matrix(t(cosmicsig)),diag=FALSE, upper=FALSE, method="cosine",auto_convert_data_frames = FALSE))
#cossimil[lower.tri(cossimil,diag=TRUE)]=NA
cossimil <- as.data.frame(as.table(cossimil))
cossimil=na.omit(cossimil)
names(cossimil) <- c("sample1","sample2","simil")
cossimil$simil <- round(cossimil$simil,2)

sig_order <- c("Control","OGG1","UNG","RNF168","EXO1")
filename=paste0("cossim_KOsig_Refsig",".pdf")
pdf(file=filename, onefile=TRUE,width = 10,height =3)
g1 <-ggplot(cossimil, aes(x=sample2, y=sample1)) 
g1 <-g1 +scale_fill_gradient2(high="dodgerblue", low="white",limits=c(0, 1))
g1 <-g1 +scale_y_discrete(limits = sig_order)+ geom_tile(aes(fill=simil),colour="black")#+geom_text(aes(label=paste(simil)),size=2)
g1 <-g1 +scale_x_discrete(limits = c("Ref.Sig.1","Ref.Sig.2","Ref.Sig.3","Ref.Sig.4", "Ref.Sig.5","Ref.Sig.7","Ref.Sig.8","Ref.Sig.9","Ref.Sig.10",
                                     "Ref.Sig.11","Ref.Sig.13","Ref.Sig.16","Ref.Sig.17","Ref.Sig.18","Ref.Sig.19","Ref.Sig.22","Ref.Sig.24","Ref.Sig.30","Ref.Sig.33",
                                     "Ref.Sig.36","Ref.Sig.38","Ref.Sig.51","Ref.Sig.52","Ref.Sig.PLATINUM","Ref.Sig.MMR1","Ref.Sig.MMR2"))

g1 <-g1 +theme(axis.text.x=element_text(angle=90,size=10,colour = "black",hjust=0.9),
               axis.text.y=element_text(size=10,colour = "black"),
               plot.title = element_text(size=10),
               # panel.grid.minor.x=element_blank(),
               # panel.grid.major.x=element_blank(),
               # panel.grid.major.y = element_blank(),
               # panel.grid.minor.y = element_blank(),
               # panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
g1
dev.off()



######################################################
#   Figure 2C
######################################################
# OGG1
# odds ratio of mutations in 16 trinucleotide for C>T

ung_sig <- read.table(paste0("./sigfiles/","OGG1_106_25_sig.txt"), sep = "\t", header = T, as.is = T)
ung_sig$Mutation <- substr(ung_sig$MutationType,3,5)
ung_sig_CT <- ung_sig[ung_sig$Mutation=="C>A",]
ung_sig_CT[ung_sig_CT$KO_exposure==0,"KO_exposure"] <- 1

ung_sig_CT$trinuc_pyrimidine <- paste0(substr(ung_sig_CT$MutationType, 1,1), substr(ung_sig_CT$MutationType, 3,3), substr(ung_sig_CT$MutationType, 7,7))
tri_ct <- read.table("../00_common/trinucleotides_CT_wg.txt", sep = "\t", header = T, as.is = T)
#tri_ct_C <- tri_ct[tri_ct$Ref=="C",]

ung_sig_CT <- merge(ung_sig_CT,tri_ct, by="trinuc_pyrimidine")
ung_sig_CT$totalmutation <- 0
ung_sig_CT[ung_sig_CT$Ref=="C",]$totalmutation <- sum(ung_sig_CT[ung_sig_CT$Ref=="C",]$KO_exposure)

ung_sig_CT$totalCT <- 0
ung_sig_CT[ung_sig_CT$Ref=="C",]$totalCT <- sum(ung_sig_CT[ung_sig_CT$Ref=="C",]$total_freq)


ung_sig_CT[,"oddsratio"] <- 1
ung_sig_CT[,"lowerconfint"] <- 1
ung_sig_CT[,"higherconfint"] <- 1


for(i in 1:dim(ung_sig_CT)[1]){
  M <- matrix(c(ung_sig_CT[i,"KO_exposure"], ung_sig_CT[i,"total_freq"], (ung_sig_CT[i,"totalmutation"]-ung_sig_CT[i,"KO_exposure"]), (ung_sig_CT[i,"totalCT"]-ung_sig_CT[i,"total_freq"])), ncol = 2)
  b <- data.frame(t(calcOddsRatio(M)))
  names(b) <- c("oddsratio","lowerconfint","higherconfint")
  ung_sig_CT[i,"oddsratio"] <- b$oddsratio
  ung_sig_CT[i,"lowerconfint"] <- b$lowerconfint
  ung_sig_CT[i,"higherconfint"] <- b$higherconfint
  
}

ung_sig_CT <- ung_sig_CT[order(ung_sig_CT$MutationType),]
filename=paste0("OGG1","_oddsratio.pdf")
pdf(file=filename, onefile=TRUE,height=2.5,width = 4.5, useDingbats=FALSE) 
d1 <- ggplot(ung_sig_CT,aes(x=oddsratio,y=MutationType))+geom_vline(aes(xintercept = 1), size = 1, linetype = "dashed") + xlim(0,16)
d1 <- d1 + geom_errorbarh(aes(xmax = higherconfint, xmin = lowerconfint), size = 1, height = .2, color = "gray50")
d1 <- d1 + geom_point(size = 2,colour="blue") +ylab("") +xlab("Odds Ratio")
d1 <- d1 + scale_y_discrete(breaks = ung_sig_CT$MutationType, labels = ung_sig_CT$trinuc_pyrimidine,limits = ung_sig_CT$MutationType)+ coord_flip()
d1 <- d1+theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5,colour = "black"),
        axis.text.y=element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA))
# d1 <- d1+facet_wrap(~indel,ncol=5,scales="free")
print(d1)
dev.off()

write.table(ung_sig_CT,paste0("OGG1","_oddsratio.txt"),sep = "\t", col.names = T, row.names = F, quote = F)



# MUTYH
# odds ratio of mutations in 16 trinucleotide for C>T


cosmic_sig <- read.table("../00_common/sigProfiler_SBS_signatures_2019_05_22.csv",sep = ",", header = T, as.is = T)
cosmic_sig <- cosmic_sig[order(cosmic_sig$Type),]
cosmic_sig$MutationType <- paste0(substr(cosmic_sig$SubType,1,1),"[",cosmic_sig$Type,"]",substr(cosmic_sig$SubType,3,3))
cosmic_sig <- cosmic_sig[,c("SBS36","MutationType")]
ung_sig <- cosmic_sig
ung_sig$Mutation <- substr(ung_sig$MutationType,3,5)
ung_sig_CT <- ung_sig[ung_sig$Mutation=="C>A",]
names(ung_sig_CT) <- c("KO_exposure","MutationType","Mutation")
ung_sig_CT$KO_exposure <- 1000*ung_sig_CT$KO_exposure

ung_sig_CT$trinuc_pyrimidine <- paste0(substr(ung_sig_CT$MutationType, 1,1), substr(ung_sig_CT$MutationType, 3,3), substr(ung_sig_CT$MutationType, 7,7))
tri_ct <- read.table("./00_common/trinucleotides_CT_wg.txt", sep = "\t", header = T, as.is = T)
#tri_ct_C <- tri_ct[tri_ct$Ref=="C",]

ung_sig_CT <- merge(ung_sig_CT,tri_ct, by="trinuc_pyrimidine")
ung_sig_CT$totalmutation <- 0
ung_sig_CT[ung_sig_CT$Ref=="C",]$totalmutation <- sum(ung_sig_CT[ung_sig_CT$Ref=="C",]$KO_exposure)

ung_sig_CT$totalCT <- 0
ung_sig_CT[ung_sig_CT$Ref=="C",]$totalCT <- sum(ung_sig_CT[ung_sig_CT$Ref=="C",]$total_freq)


ung_sig_CT[,"oddsratio"] <- 1
ung_sig_CT[,"lowerconfint"] <- 1
ung_sig_CT[,"higherconfint"] <- 1


for(i in 1:dim(ung_sig_CT)[1]){
  M <- matrix(c(ung_sig_CT[i,"KO_exposure"], ung_sig_CT[i,"total_freq"], (ung_sig_CT[i,"totalmutation"]-ung_sig_CT[i,"KO_exposure"]), (ung_sig_CT[i,"totalCT"]-ung_sig_CT[i,"total_freq"])), ncol = 2)
  b <- data.frame(t(calcOddsRatio(M)))
  names(b) <- c("oddsratio","lowerconfint","higherconfint")
  ung_sig_CT[i,"oddsratio"] <- b$oddsratio
  ung_sig_CT[i,"lowerconfint"] <- b$lowerconfint
  ung_sig_CT[i,"higherconfint"] <- b$higherconfint
  
}

ung_sig_CT <- ung_sig_CT[order(ung_sig_CT$MutationType),]
filename=paste0("MUTYH","_oddsratio.pdf")
pdf(file=filename, onefile=TRUE,height=2.5,width = 4.5, useDingbats=FALSE) 
d1 <- ggplot(ung_sig_CT,aes(x=oddsratio,y=MutationType))+geom_vline(aes(xintercept = 1), size = 1, linetype = "dashed")+ xlim(0,16)
d1 <- d1 + geom_errorbarh(aes(xmax = higherconfint, xmin = lowerconfint), size = 1, height = .2, color = "gray50")
d1 <- d1 + geom_point(size = 2,colour="blue") +ylab("") +xlab("Odds Ratio")
d1 <- d1 + scale_y_discrete(breaks = ung_sig_CT$MutationType, labels = ung_sig_CT$trinuc_pyrimidine,limits = ung_sig_CT$MutationType)+ coord_flip()
d1 <- d1+theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5,colour = "black"),
        axis.text.y=element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA))
# d1 <- d1+facet_wrap(~indel,ncol=5,scales="free")
print(d1)
dev.off()

write.table(ung_sig_CT,paste0("MUTYH","_oddsratio.txt"),sep = "\t", col.names = T, row.names = F, quote = F)


######################################################
#   Figure 2D
######################################################

# UNG
# odds ratio of mutations in 16 trinucleotide for C>T
ung_sig <- read.table(paste0("./sigfiles/","UNG_0.6_sig.txt"), sep = "\t", header = T, as.is = T)
ung_sig$Mutation <- substr(ung_sig$MutationType,3,5)
ung_sig_CT <- ung_sig[ung_sig$Mutation=="C>T",]
ung_sig_CT[ung_sig_CT$KO_exposure==0,"KO_exposure"] <- 1

ung_sig_CT$trinuc_pyrimidine <- paste0(substr(ung_sig_CT$MutationType, 1,1), substr(ung_sig_CT$MutationType, 3,3), substr(ung_sig_CT$MutationType, 7,7))
tri_ct <- read.table("../00_common/trinucleotides_CT_wg.txt", sep = "\t", header = T, as.is = T)
#tri_ct_C <- tri_ct[tri_ct$Ref=="C",]

ung_sig_CT <- merge(ung_sig_CT,tri_ct, by="trinuc_pyrimidine")
ung_sig_CT$totalmutation <- 0
ung_sig_CT[ung_sig_CT$Ref=="C",]$totalmutation <- sum(ung_sig_CT[ung_sig_CT$Ref=="C",]$KO_exposure)

ung_sig_CT$totalCT <- 0
ung_sig_CT[ung_sig_CT$Ref=="C",]$totalCT <- sum(ung_sig_CT[ung_sig_CT$Ref=="C",]$total_freq)


ung_sig_CT[,"oddsratio"] <- 1
ung_sig_CT[,"lowerconfint"] <- 1
ung_sig_CT[,"higherconfint"] <- 1


for(i in 1:dim(ung_sig_CT)[1]){
  M <- matrix(c(ung_sig_CT[i,"KO_exposure"], ung_sig_CT[i,"total_freq"], (ung_sig_CT[i,"totalmutation"]-ung_sig_CT[i,"KO_exposure"]), (ung_sig_CT[i,"totalCT"]-ung_sig_CT[i,"total_freq"])), ncol = 2)
  b <- data.frame(t(calcOddsRatio(M)))
  names(b) <- c("oddsratio","lowerconfint","higherconfint")
  ung_sig_CT[i,"oddsratio"] <- b$oddsratio
  ung_sig_CT[i,"lowerconfint"] <- b$lowerconfint
  ung_sig_CT[i,"higherconfint"] <- b$higherconfint
  
}

ung_sig_CT <- ung_sig_CT[order(ung_sig_CT$Mutation),]
filename=paste0("UNG","_oddsratio2.pdf")
pdf(file=filename, onefile=TRUE,height=2.5,width = 4.5, useDingbats=FALSE) 
d1 <- ggplot(ung_sig_CT,aes(x=oddsratio,y=MutationType))+geom_vline(aes(xintercept = 1), size = 1, linetype = "dashed") + xlim(0,6)
d1 <- d1 + geom_errorbarh(aes(xmax = higherconfint, xmin = lowerconfint), size = 1, height = .2, color = "gray50")
d1 <- d1 + geom_point(size = 2,colour="red") +ylab("") +xlab("Odds Ratio")
d1 <- d1 + scale_y_discrete(breaks = ung_sig_CT$MutationType, labels = ung_sig_CT$trinuc_pyrimidine,limits = ung_sig_CT$MutationType)+ coord_flip()
d1 <- d1+theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5,colour = "black"),
        axis.text.y=element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA))
# d1 <- d1+facet_wrap(~indel,ncol=5,scales="free")
print(d1)
dev.off()

write.table(ung_sig_CT,paste0("UNG","_oddsratio.txt"),sep = "\t", col.names = T, row.names = F, quote = F)

# NTHL1 
# odds ratio of mutations in 16 trinucleotide for C>T
cosmic_sig <- read.table("../00_common/sigProfiler_SBS_signatures_2019_05_22.csv",sep = ",", header = T, as.is = T)
cosmic_sig <- cosmic_sig[order(cosmic_sig$Type),]
cosmic_sig$MutationType <- paste0(substr(cosmic_sig$SubType,1,1),"[",cosmic_sig$Type,"]",substr(cosmic_sig$SubType,3,3))
cosmic_sig <- cosmic_sig[,c("SBS30","MutationType")]
ung_sig <- cosmic_sig
ung_sig$Mutation <- substr(ung_sig$MutationType,3,5)
ung_sig_CT <- ung_sig[ung_sig$Mutation=="C>T",]
names(ung_sig_CT) <- c("KO_exposure","MutationType","Mutation")
ung_sig_CT$KO_exposure <- 4000*ung_sig_CT$KO_exposure

ung_sig_CT$trinuc_pyrimidine <- paste0(substr(ung_sig_CT$MutationType, 1,1), substr(ung_sig_CT$MutationType, 3,3), substr(ung_sig_CT$MutationType, 7,7))
tri_ct <- read.table("../00_common/trinucleotides_CT_wg.txt", sep = "\t", header = T, as.is = T)
#tri_ct_C <- tri_ct[tri_ct$Ref=="C",]

ung_sig_CT <- merge(ung_sig_CT,tri_ct, by="trinuc_pyrimidine")
ung_sig_CT$totalmutation <- 0
ung_sig_CT[ung_sig_CT$Ref=="C",]$totalmutation <- sum(ung_sig_CT[ung_sig_CT$Ref=="C",]$KO_exposure)

ung_sig_CT$totalCT <- 0
ung_sig_CT[ung_sig_CT$Ref=="C",]$totalCT <- sum(ung_sig_CT[ung_sig_CT$Ref=="C",]$total_freq)


ung_sig_CT[,"oddsratio"] <- 1
ung_sig_CT[,"lowerconfint"] <- 1
ung_sig_CT[,"higherconfint"] <- 1


for(i in 1:dim(ung_sig_CT)[1]){
  M <- matrix(c(ung_sig_CT[i,"KO_exposure"], ung_sig_CT[i,"total_freq"], (ung_sig_CT[i,"totalmutation"]-ung_sig_CT[i,"KO_exposure"]), (ung_sig_CT[i,"totalCT"]-ung_sig_CT[i,"total_freq"])), ncol = 2)
  b <- data.frame(t(calcOddsRatio(M)))
  names(b) <- c("oddsratio","lowerconfint","higherconfint")
  ung_sig_CT[i,"oddsratio"] <- b$oddsratio
  ung_sig_CT[i,"lowerconfint"] <- b$lowerconfint
  ung_sig_CT[i,"higherconfint"] <- b$higherconfint
  
}

ung_sig_CT <- ung_sig_CT[order(ung_sig_CT$Mutation),]
filename=paste0("NTHL1","_oddsratio.pdf")
pdf(file=filename, onefile=TRUE,height=2.5,width = 4.5, useDingbats=FALSE) 
d1 <- ggplot(ung_sig_CT,aes(x=oddsratio,y=MutationType))+geom_vline(aes(xintercept = 1), size = 1, linetype = "dashed") + xlim(0,6)
d1 <- d1 + geom_errorbarh(aes(xmax = higherconfint, xmin = lowerconfint), size = 1, height = .2, color = "gray50")
d1 <- d1 + geom_point(size = 2,colour="red") +ylab("") +xlab("Odds Ratio")
d1 <- d1 + scale_y_discrete(breaks = ung_sig_CT$MutationType, labels = ung_sig_CT$trinuc_pyrimidine,limits = ung_sig_CT$MutationType)+ coord_flip()
d1 <- d1+theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5,colour = "black"),
        axis.text.y=element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA))
# d1 <- d1+facet_wrap(~indel,ncol=5,scales="free")
print(d1)
dev.off()

write.table(ung_sig_CT,paste0("NTHL1","_oddsratio.txt"),sep = "\t", col.names = T, row.names = F, quote = F)






######################################################
#   Figure 2E 2F
######################################################
transtrand_length <- read.table("../00_common/TranscriptionalStrand_length.txt", sep = "\t", header = T, as.is = T)
Trans_32length <- read.table("../00_common/TransStrand_32length.txt",sep = "\t", header = T, as.is = T)


total_muts <- read.table("../00_common/total_subs_43genes.txt", sep = "\t", header = T, as.is = T)
total_muts[total_muts$Chrom==23,"Chrom"]="X"
total_muts[total_muts$Chrom==24,"Chrom"]="Y"
total_muts$VariantID <- paste0(total_muts$Sample,"_",total_muts$Chrom,"_",total_muts$Pos)
total_muts <- total_muts[total_muts$Ko_gene %in%c("RNF168","EXO1"),]

#total_muts$pre_context <- as.character(getSeq(Hsapiens, paste0('chr',total_muts$Chrom), total_muts$Pos-1, total_muts$Pos-1))
#total_muts$rear_context <- as.character(getSeq(Hsapiens, paste0('chr',total_muts$Chrom), total_muts$Pos+1, total_muts$Pos+1))

# Transcriptional 
Subs_TranStrand <- AddStrandInfo_intersect(total_muts,"denovo_muts_bed.txt","../00_common/TranscribStrand.uts.txt","../00_common/TranscribStrand.ts.txt","subs_TranscriptionalStrand","total_muts_transtrand.txt")
Subs_TranStrand <- read.table("total_muts_transtrand.txt", sep = "\t", header = T, as.is = T)
StrandBias_6muttype(Subs_TranStrand, "Sample_ko","Strand",transtrand_length,12,12,6,"TranSB")
StrandBias_6muttype(Subs_TranStrand, "Ko_gene","Strand",transtrand_length,5,10,5,"TranSB_Gene")

StrandBias_96muttype(Subs_TranStrand, "Sample_ko","Strand",Trans_32length,5,12,3,"TranSB96")
StrandBias_96muttype(Subs_TranStrand, "Ko_gene","Strand",Trans_32length,2,10,2,"TranSB_Gene96")

# odds ratio
StrandBias_oddsratio(Subs_TranStrand, "Ko_gene","Strand",transtrand_length,5,10,5,"TranSB_Gene")


