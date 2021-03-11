source("../00_common/indel_functions.R")
source("../00_common/sub_functions.R")

##################
#   Figure 3A
##################
subs <- read.table("../00_common/total_subs_43genes.txt", sep = "\t", header = T, as.is = T)
mutation_catalogue <- GenCatalogue(subs,"Sample_ko")

Wrap_KOSig(mutation_catalogue,"mean",c("PMS2_170_s1","PMS2_170_s2","PMS2_171_s1","PMS2_171_s2"),100,150,2,"PMS2_170_171_sig")
Wrap_KOSig(mutation_catalogue,"mean",c("MLH1_172_s1","MLH1_172_s2","MLH1_173_s1","MLH1_173_s2"),100,150,2,"MLH1_172_173_sig")
Wrap_KOSig(mutation_catalogue,"mean",c("MSH6_3_s4","MSH6_3_s6","MSH6_3_s5","MSH6_3_s8","MSH6_4_s2","MSH6_4_s3","MSH6_4_s4","MSH6_4_s7"),100,150,2,"MSH6_0.3_0.4_sig")
Wrap_KOSig(mutation_catalogue,"mean",c("MSH2_120_s1","MSH2_120_s2","MSH2_120_s3"),100,150,2,"MSH2_120_sig")
Wrap_KOSig(mutation_catalogue,"mean",c("PMS1_123_s1","PMS1_123_s2","PMS1_130_s1","PMS1_130_s2"),100,150,2,"PMS1_123_130_sig")

##################
#   Figure 3B
##################
indels <- read.table("../00_common/total_indels_43genes.txt", sep = "\t", header = T, as.is = T)

# 15 channels
indel.classified <- indel_classifier15(indels)
indel.classified$Sample <- indel.classified$Genotype
mutation_catalogue <- gen_indelmuttype_15(indel.classified)
write.table(mutation_catalogue,"indel_catalogue_15channels",sep = "\t", col.names = T, row.names = F, quote = F)

# 45 channels
indel.classified <- indel_classifier(indels)
indel.classified$Sample <- indel.classified$Genotype
mutation_catalogue <- gen_indelmuttype_MMRD(indel.classified,"Sample","indeltype_short")
write.table(mutation_catalogue,"indel_catalogue_45channels",sep = "\t", col.names = T, row.names = F, quote = F)

#ko_indels_type_4_catalouge <- read.table("Ko_gene_type_4_mmrd_channeled.txt", sep = "\t", header = T, as.is = T)
#plot_percentage_type_4_mmrd(ko_indels_type_4_catalouge, 1,26, 13, "Ko_gene_type_4_mmrd_percentage") # for figure 3
Wrap_KOSig_indel(mutation_catalogue,"mean",c("EXO1_71_s2","EXO1_71_s3","EXO1_71_s4"),100,10,2,15,"EXO1_71_indelsig")
Wrap_KOSig_indel(mutation_catalogue,"mean",c("MSH6_3_s4","MSH6_3_s6","MSH6_3_s5","MSH6_3_s8","MSH6_4_s2","MSH6_4_s3","MSH6_4_s4","MSH6_4_s7"),100,10,2,15,"MSH6_3_4_indelsig")
Wrap_KOSig_indel(mutation_catalogue,"mean",c("MSH2_120_s1","MSH2_120_s2","MSH2_120_s3"),100,10,2,15,"MSH2_120_indelsig")
Wrap_KOSig_indel(mutation_catalogue,"mean",c("PMS1_123_s1","PMS1_123_s2","PMS1_130_s1","PMS1_130_s2"),100,8,2,15,"PMS1_123_130_indelsig")
Wrap_KOSig_indel(mutation_catalogue,"mean",c("PMS2_170_s1","PMS2_170_s2","PMS2_171_s1","PMS2_171_s2"),100,10,2,15,"PMS2_170_171_indelsig")
Wrap_KOSig_indel(mutation_catalogue,"mean",c("MLH1_172_s1","MLH1_172_s2","MLH1_173_s1","MLH1_173_s2"),100,10,2,15,"MLH1_172_173_indelsig")


##################
#   Figure 3C
##################

replitrand_length <- read.table("../00_common/ReplicativeStrand_length.txt", sep = "\t", header = T, as.is = T)
Repli_32length <- read.table("../00_common/RepliStrand_32length.txt",sep = "\t", header = T, as.is = T)


total_muts <- read.table("../00_common/total_subs_43genes.txt", sep = "\t", header = T, as.is = T)
total_muts[total_muts$Chrom==23,"Chrom"]="X"
total_muts[total_muts$Chrom==24,"Chrom"]="Y"
total_muts$VariantID <- paste0(total_muts$Sample,"_",total_muts$Chrom,"_",total_muts$Pos)
total_muts <- total_muts[total_muts$Ko_gene %in%c("MSH2","MSH6","MLH1","PMS2"),]

Subs_ReplicStrand <- AddStrandInfo_intersect(total_muts,"denovo_muts_bed.txt","/nfs/cgp_signatures/Analysis/xz3/00_Resource/GenomeFeature/00_common/MCF7_RepliStrand.lagging","/nfs/cgp_signatures/Analysis/xz3/00_Resource/GenomeFeature/00_common/MCF7_RepliStrand.leading","subs_ReplicativeStrand","total_muts_replictrand.txt")
Subs_ReplicStrand <- read.table("total_muts_replictrand.txt", sep = "\t", header = T, as.is = T)

StrandBias_6muttype(Subs_ReplicStrand, "Sample_ko","Strand",replitrand_length,12,12,6,"ReplicSB")
StrandBias_6muttype(Subs_ReplicStrand, "Ko_gene","Strand",replitrand_length,5,10,5,"ReplicSB_Gene")

StrandBias_96muttype(Subs_ReplicStrand, "Sample_ko","Strand",Repli_32length,12,12,4,"ReplicSB96")
StrandBias_96muttype(Subs_ReplicStrand, "Ko_gene","Strand",Repli_32length,12,12,2,"ReplicSB96_Gene")

muts_Strand <- read.table("ReplicSB_Gene_chisq_adjust.txt", sep = "\t", header = T, as.is = T)
muts_Strand[,"oddsratio"] <- 1
muts_Strand[,"lowerconfint"] <- 1
muts_Strand[,"higherconfint"] <- 1


for(i in 1:dim(muts_Strand)[1]){
  M <- matrix(c(muts_Strand[i,"lagging_ts"], muts_Strand[i,"leading_uts"], (muts_Strand[i,"ts_lagging_wg"]-muts_Strand[i,"lagging_ts"]), (muts_Strand[i,"uts_leading_wg"]-muts_Strand[i,"leading_uts"])), ncol = 2)
  b <- data.frame(t(calcOddsRatio(M)))
  names(b) <- c("oddsratio","lowerconfint","higherconfint")
  muts_Strand[i,"oddsratio"] <- b$oddsratio
  muts_Strand[i,"lowerconfint"] <- b$lowerconfint
  muts_Strand[i,"higherconfint"] <- b$higherconfint
  
}

filename=paste0("MMR4_genes","_oddsratio.pdf")
pdf(file=filename, onefile=TRUE,height=2,width = 10, useDingbats=FALSE) 
d1 <- ggplot(muts_Strand,aes(x=oddsratio,y=Sample, color = "red"))+geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") #+ xlim(0.5,2)
d1 <- d1 + geom_errorbarh(aes(xmax = higherconfint, xmin = lowerconfint), size = .5, height = .2, color = "gray50")
d1 <- d1 + geom_point(size = 1.5) +ylab("") +xlab("Odds Ratio")
#  d1 <- d1 + coord_trans(x = scales:::exp_trans(10)) + scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1),
#                   limits = log10(c(0.09,2.5)))
d1 <- d1+theme_bw() +
  theme(axis.text.x=element_text(colour = "black"),
        axis.text.y=element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA))
# d1 <- d1+facet_wrap(~indel,ncol=5,scales="free")
d1 <- d1+facet_grid(.~Mutation)
print(d1)
dev.off()

write.table(muts_Strand,paste0("MMR4_genes","_oddsratio.txt"),sep = "\t", col.names = T, row.names = F, quote = F)




##################
#   Figure 3E
# C>A mutation in polyC/polyG tracks
##################
# Annotate mutations if it is in Polynucleotide tracks
total_muts <- read.table("../00_common/total_subs_43genes.txt", sep = "\t", header = T, as.is = T)
total_muts[total_muts$Chrom==23,"Chrom"]="X"
total_muts[total_muts$Chrom==24,"Chrom"]="Y"
total_muts <- total_muts[total_muts$Ko_gene %in%c("MSH2","MSH6","MLH1","PMS2"),]

total_muts_context <-  Find_SeqContext(total_muts,10)
Find_RepeatSeq(total_muts_context)

total_muts_repeat <- read.table("./CTsubs_repeat.txt",sep = "\t",header = T, as.is = T)
total_muts_repeat[total_muts_repeat$rear_pos==Inf, "rear_repeat"] <- total_muts_repeat[total_muts_repeat$rear_pos==Inf, "Ref_py"]
total_muts_repeat[total_muts_repeat$rear_pos==Inf, "rear_pos"] <- 11

total_muts_repeat[total_muts_repeat$pre_pos==-Inf, "pre_repeat"] <- total_muts_repeat[total_muts_repeat$pre_pos==-Inf, "Ref_py"]
total_muts_repeat[total_muts_repeat$pre_pos==-Inf, "pre_pos"] <- 0

total_muts_repeat$rep_len <- total_muts_repeat$rear_pos-total_muts_repeat$pre_pos+10
total_muts_repeat$mut_pos_to5 <- 10-total_muts_repeat$pre_pos+1
total_muts_repeat$mut_pos_to3 <- total_muts_repeat$rear_pos

muts_rep_summary <- data.frame(table(total_muts_repeat$Ko_gene,total_muts_repeat$rep_len,total_muts_repeat$Ref_py,total_muts_repeat$Alt_py))
names(muts_rep_summary) <- c("Ko_gene","Rep_len","Ref","Alt","Freq")
muts_rep_summary$Mutation <- paste0(muts_rep_summary$Ref,">",muts_rep_summary$Alt)
muts_rep_summary$MutationType <- paste0(muts_rep_summary$Mutation,"_",muts_rep_summary$Rep_len)
muts_rep_summary <- muts_rep_summary[! muts_rep_summary$Mutation %in% c("C>C","T>T"),]
# Only C>A and rep>=2 and MLH1, MSH2, MSH6
total_muts_repeat_CA <- total_muts_repeat[total_muts_repeat$Mutation=="C>A"& total_muts_repeat$rep_len>1 & total_muts_repeat$Ko_gene %in% c("MLH1", "MSH2","MSH6"),]

# rep_length
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
b <- data.frame(table(total_muts_repeat_CA$Ko_gene,total_muts_repeat_CA$rear_repeat,total_muts_repeat_CA$mut_pos_to3,total_muts_repeat_CA$rep_len))
names(b) <- c("Ko_gene","rear_base","mut_pos_to3","rep_len","Freq")
b$seq_context <- paste0(b$rep_len,"C","_",b$rear_base)
pdf(file="MMR_CA_replength.pdf", onefile=TRUE,width=10,height=8)
p <- ggplot(data=b, aes(x=rep_len, y=Freq,fill=mut_pos_to3))+ geom_bar(stat="identity", width=0.8)+xlab("Repeat length")+ylab("Count")
#p <- p+scale_x_discrete(limits = as.character(muts_rep_summary$MutationType))
p <- p+scale_fill_manual(values=cbPalette)
p <- p+theme(#axis.text.x=element_blank(),
  axis.text.x=element_text(size=10,colour = "black"),
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
p <- p+facet_grid(Ko_gene~rear_base,scales = "free")

print(p)
dev.off()

write.table(b,paste0("MMR4_genes","MMR_CA_replength.txt"),sep = "\t", col.names = T, row.names = F, quote = F)


# rep_length percentage
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
b <- data.frame(table(total_muts_repeat_CA$Ko_gene,total_muts_repeat_CA$rear_repeat,total_muts_repeat_CA$mut_pos_to3,total_muts_repeat_CA$rep_len))
names(b) <- c("Ko_gene","rear_base","mut_pos_to3","rep_len","Freq")
b$seq_context <- paste0(b$rep_len,"C","_",b$rear_base)
pdf(file="MMR_CA_replength_percentage.pdf", onefile=TRUE,width=10,height=8)
p <- ggplot(data=b[b$rep_len%in%c(2,3,4,5),], aes(x=rep_len, y=Freq,fill=mut_pos_to3))+ geom_bar(stat="identity", width=0.8,position="fill")+xlab("Repeat length")+ylab("Count")
#p <- p+scale_x_discrete(limits = as.character(muts_rep_summary$MutationType))
p <- p+scale_fill_manual(values=cbPalette)+ scale_y_continuous(labels=percent)
p <- p+theme(#axis.text.x=element_blank(),
  axis.text.x=element_text(size=10,colour = "black"),
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
p <- p+facet_grid(Ko_gene~rear_base,scales = "free")

print(p)
dev.off()





##################
#   Figure 3F
##################
subs <- read.table("../00_common/total_subs_43genes.txt", sep = "\t", header = T, as.is = T)
mmr_muts <- subs[subs$Ko_gene %in%c("MSH2","MSH6","MLH1","PMS2"),]

mut_catalogue <- Gen1536Catalogue(mmr_muts,"Ko_gene")
plot1536Countbasis(mut_catalogue[,-2],1,10,40,"denovo_muts_1536catalogue_KOgene.pdf")
write.table(mut_catalogue,"mut_catalogue1536.txt", sep = "\t", col.names = T, row.names = F, quote = F)

mut_catalogue <- Gen1536Catalogue(mmr_muts,"Sample_ko")
plot1536Countbasis(mut_catalogue[,-2],1,40,40,"denovo_muts_1536catalogue_KOSample.pdf")

##################
#   Figure 3G
##################
subs <- read.table("../00_common/total_subs_43genes.txt", sep = "\t", header = T, as.is = T)
mmr_muts <- subs[subs$Ko_gene %in%c("MSH2","MSH6","MLH1","PMS2"),]

mmr_muts_seqcontext <- Find_SeqContext(mmr_muts,10)

TtoA <- mmr_muts_seqcontext[mmr_muts_seqcontext$Mutation=="T>A",]

TtoA$pre_context_py_1 <- substr(TtoA$pre_context_py,10,10)
TtoA$rear_context_py_1 <- substr(TtoA$rear_context_py,1,1)
TtoA$pre_repeat <- 1
TtoA$rear_repeat <- 1


for(i in 1:dim(TtoA)[1]){
  print(i)
  TtoA[i,"pre_repeat"] <- which(strsplit(stri_reverse(TtoA[i,"pre_context_py"]), "")[[1]]!=TtoA[i,"pre_context_py_1"])[1]-1
  TtoA[i,"rear_repeat"] <- which(strsplit(TtoA[i,"rear_context_py"], "")[[1]]!=TtoA[i,"rear_context_py_1"])[1]-1
  
}
TtoA[is.na(TtoA)] <- 10
write.table(TtoA,"TtoA_MMR_seqcontext.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# 12 doubling time and clonal
TtoA <- TtoA[!TtoA$Sample %in% c("MSK0.170_s2","MSK0.171_s2","MSK0.120_s3","MSK0.3_s5","MSK0.3_s6","MSK0.3_s8","MSK0.4_s3","MSK0.4_s4","MSK0.4_s7"),]

# AT
TtoA_preA <- TtoA[TtoA$pre_context_py_1=="A" ,]
TtoA_preA$flag <- ""
TtoA_preA[TtoA_preA$rear_context_py_1=="T","flag"] <- "AATT"
TtoA_preA[TtoA_preA$rear_context_py_1!="T","flag"] <- "AT"
a <- data.frame(table(TtoA_preA$Sample,TtoA_preA$flag))
names(a) <- c("Sample","flag","sub_num")

sample_ko <- data.frame(table(TtoA_preA$Sample,TtoA_preA$Ko_gene))
sample_ko <- sample_ko[sample_ko$Freq>0,]
names(sample_ko) <- c("Sample","Ko_gene")
a <- merge(a,sample_ko[,c(1,2)],by="Sample")

#/nfs/cgp_signatures/Analysis/xz3/00_Resource/GenomeFeature/Polynucleotides/AT_palidrom
AT_TA_wg <- data.frame("flag"=c("AT","TA","AATT","TTAA"),"Freq"=c(199604821,168805669,21456762,19045923))
a <- merge(a,AT_TA_wg,by="flag")
a$odds <- a$sub_num*1000000000/a$Freq
a_odds_ddply <- ddply(a,c("Ko_gene","flag"),summarise,NChild=length(odds),mean=mean(odds),se=sd(odds)/sqrt(NChild))

# Number
a_numb_ddply <- ddply(a,c("Ko_gene","flag"),summarise,NChild=length(sub_num),mean=mean(sub_num),se=sd(sub_num)/sqrt(NChild))

# TA
TtoA_rearA <- TtoA[TtoA$rear_context_py_1=="A",]
TtoA_rearA$flag <- ""
TtoA_rearA[TtoA_rearA$pre_context_py_1=="T","flag"] <- "TTAA"
TtoA_rearA[TtoA_rearA$pre_context_py_1!="T","flag"] <- "TA"
b <- data.frame(table(TtoA_rearA$Sample,TtoA_rearA$flag))
names(b) <- c("Sample","flag","sub_num")

sample_ko <- data.frame(table(TtoA_rearA$Sample,TtoA_rearA$Ko_gene))
sample_ko <- sample_ko[sample_ko$Freq>0,]
names(sample_ko) <- c("Sample","Ko_gene")
b <- merge(b,sample_ko[,c(1,2)],by="Sample")

#/nfs/cgp_signatures/Analysis/xz3/00_Resource/GenomeFeature/Polynucleotides/AT_palidrom
AT_TA_wg <- data.frame("flag"=c("AT","TA","AATT","TTAA"),"Freq"=c(199604821,168805669,21456762,19045923))
b <- merge(b,AT_TA_wg,by="flag")
b$odds <- b$sub_num*1000000000/b$Freq
b_odds_ddply <- ddply(b,c("Ko_gene","flag"),summarise,NChild=length(odds),mean=mean(odds),se=sd(odds)/sqrt(NChild))

# Number
b_numb_ddply <- ddply(b,c("Ko_gene","flag"),summarise,NChild=length(sub_num),mean=mean(sub_num),se=sd(sub_num)/sqrt(NChild))

c <- rbind(a_odds_ddply,b_odds_ddply)
c_all <- rbind(a,b)
write.table(c_all,"4MMR_TtoA_AnTn_count.txt",sep = "\t",col.names = T, row.names = F)

pdf(file=paste0("4MMR_TtoA_AnTn_odds",".pdf"), onefile=TRUE,height=3,width=5, useDingbats=FALSE)
p <- ggplot(data=c, aes(x=flag, y=mean,fill=Ko_gene))+ geom_bar(stat="identity",position="dodge",width=.8,alpha=0.5)+xlab("Mutations")+ylab("Number")
p <- p+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",position=position_dodge(.8),width=.1)
p <- p+scale_x_discrete(limits = c("AATT","AT","TTAA","TA"))
p <- p+scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
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

print(p)
dev.off()


pdf(file=paste0("4MMR_TtoA_AnTn_odds2",".pdf"), onefile=TRUE,height=3,width=5, useDingbats=FALSE)
p <- ggplot(c, aes(x=flag, y=mean, fill=Ko_gene)) +
  geom_bar(position=position_dodge(), stat="identity",width=0.8,alpha=0.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),color="black",position=position_dodge(.8),width=.1)+ geom_point(data=c_all,aes(x=flag, y=odds),position=position_jitterdodge(0.1))
p <- p+scale_x_discrete(limits = c("AATT","AT","TTAA","TA"))
p <- p+scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))
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

print(p)
dev.off()

