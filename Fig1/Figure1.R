source("../00_common/indel_functions.R")
source("../00_common/sub_functions.R")

dir.create("./Fig1")
setwd("./Fig1/")

##################
#   Figure 1C
##################
subs <- read.table("../00_common/total_subs_43genes.txt", sep = "\t", header = T, as.is = T)

mutation_catalogue <- GenCatalogue(subs,"Genotype")
ko_info <- data.frame(table(subs$Sample_ko))
names(ko_info) <-c("Sample_ko","sub_num")
ko_info$ko <- sub("_s.*","", ko_info$Sample_ko)
ko_info2 <- ddply(ko_info[,c("ko","sub_num")],c("ko"),summarise,NChild=length(sub_num),sum=sum(sub_num))


# read in control profile
control_catalogue <- read.table("background_profile.txt",sep = "\t", header = T, as.is = T)

# calculate cosine similarity of ko profiles and control profile
simi_num <- Cossimi_CompareSig2_noplot(mutation_catalogue[,-2],control_catalogue[,c("MutationType","mean")],"simi")
names(simi_num) <- c("control","ko","simi_mean")
simi_num <- merge(simi_num,ko_info2,by="ko")
names(simi_num) <- c("ko","control","simi_mean","n","count")

# simulate control profile variation at different mutation number using bootstraping method
control_simi <- Bootstrap_profile_similarity(control_catalogue$mean,20,10000,50,"control_simi")

# plot the result
pdf(file="simi_control.pdf", onefile=TRUE,height=4,width=5, useDingbats=FALSE)
q <- ggplot(data=simi_num, aes(x=count, y=simi_mean)) + geom_point(size=3)  
q <- q+geom_point(data=control_simi,aes(x=count, y=simi_mean),colour="green",size=0.5)
q <- q+geom_errorbar(data=control_simi,aes(ymin=simi_mean-3*simi_sd,ymax=simi_mean+3*simi_sd),color="green",position=position_dodge(.9),width=.1,alpha=0.8)+scale_y_continuous()
#q <- q+geom_text(aes(label=KO_sample),color="black", size=2)
q <- q+theme(axis.text.x=element_text(size=15,colour = "black"),
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

q
dev.off()

##################
#   Figure 1D
##################
indels <- read.table("../00_common/total_indels_43genes.txt", sep = "\t", header = T, as.is = T)

indel.classified <- indel_classifier15(indels)
indel.classified$Sample <- indel.classified$Genotype
mutation_catalogue <- gen_indelmuttype_15(indel.classified)


mutation_catalogue <- read.table("denovo_KO_catalogue.txt",sep = "\t", header = T, as.is = T)
control_catalogue <- read.table("background_aggragated_indel.txt",sep = "\t", header = T, as.is = T)
pilot_catalogue <- read.table("../101_indel_pilot/indel_catalogue_subclones_aggregated.txt",sep = "\t", header = T, as.is = T)
mutation_catalogue <- merge(mutation_catalogue,pilot_catalogue,by="indelsubtype")
mutation_catalogue <- merge(mutation_catalogue[,-3],control_catalogue[,c("indelsubtype","aggragate")],by="indelsubtype")
simi_all <- NULL
for(i in 3:(dim(mutation_catalogue)[2]-1)){
  a <- cos_similarity(mutation_catalogue[,i],mutation_catalogue$aggragate)
  simi_all <- c(simi_all,a)
}
simi_num <- data.frame("ko"=colnames(mutation_catalogue[,3:(dim(mutation_catalogue)[2]-1)]),"simi"=simi_all,"num"=colSums(mutation_catalogue[,3:(dim(mutation_catalogue)[2]-1)]))


ko_info <- data.frame(table(indel.classified$Sample_ko))
names(ko_info) <-c("Sample_ko","indel_num")
ko_info$ko <- sub("_s.*","", ko_info$Sample_ko)
ko_info2 <- ddply(ko_info[,c("ko","indel_num")],c("ko"),summarise,NChild=length(indel_num),sum=sum(indel_num))


simi_num <- merge(simi_num,ko_info2,by="ko")
write.table(simi_num[order(simi_num$num,decreasing = T),],"simi_num_indel.txt",sep = "\t",col.names = T, row.names = F, quote = F)

control_simi <- Bootstrap_profile_similarity(mutation_catalogue$aggragate,20,7000,50,"control_simi_7000")
simi_num <- read.table("simi_num_indel.txt",sep = "\t", header = T, as.is = T)
names(simi_num) <- c("Sample.Name","simi_mean","count")
pdf(file="Similarity_number_withcontrol_indel.pdf", onefile=TRUE,height=4,width=4, useDingbats=FALSE)
q <- ggplot(data=simi_num, aes(x=count, y=simi_mean)) + geom_point()  #+ ylim(0.4, 1) + geom_smooth(method=lm)
q <- q+geom_point(data=control_simi,aes(x=count, y=simi_mean),colour="lightblue",size=0.5)
q <- q+geom_errorbar(data=control_simi,aes(ymin=simi_mean-3*simi_sd,ymax=simi_mean+3*simi_sd),color="lightblue",position=position_dodge(.9),width=.1,alpha=0.8)+scale_y_continuous()
#q <- q+geom_text(aes(label=KO_sample),color="black", size=2)
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

q
dev.off()

control_simi <- Bootstrap_profile_similarity(mutation_catalogue$aggragate,20,100,2,"control_simi_100")
simi_num <- read.table("simi_num_indel.txt",sep = "\t", header = T, as.is = T)
names(simi_num) <- c("Sample.Name","simi_mean","count")
pdf(file="Similarity_number_withcontrol_100_indel.pdf", onefile=TRUE,height=4,width=4, useDingbats=FALSE)
q <- ggplot(data=simi_num[simi_num$count<=100,], aes(x=count, y=simi_mean)) + geom_point()  #+ ylim(0.4, 1) + geom_smooth(method=lm)
q <- q+geom_point(data=control_simi,aes(x=count, y=simi_mean),colour="lightblue",size=1)
q <- q+geom_errorbar(data=control_simi,aes(ymin=simi_mean-3*simi_sd,ymax=simi_mean+3*simi_sd),color="lightblue",position=position_dodge(.9),width=.1,alpha=0.8)+scale_y_continuous()
#q <- q+geom_text(aes(label=KO_sample),color="black", size=2)
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

q
dev.off()

##################
#   Figure 1E
##################
sample_info <- read.table("../00_common/allsamples43genes.txt",sep = "\t",header = T, as.is = T)
sample_info$Sample_ko <- paste0(sample_info$Ko_gene,"_",sub("MSK0.","",sample_info$Sample))


sub_indel_burden <- merge(data.frame(table(subs$Sample_ko)), data.frame(table(indels$Sample_ko)), by="Var1",all=T)
sub_indel_burden[is.na(sub_indel_burden)] <- 0
names(sub_indel_burden) <- c("Sample_ko","sub","indel")
write.table(sub_indel_burden,"sub_indel_burden.txt", sep = "\t", col.names = T, row.names = F, quote = F)

#c("MSH2","MLH1","PMS2", "PMS1","EXO1","RNF168","UNG","OGG1")
sample_sel <- sample_info[sample_info$Ko_gene %in% c("ATP2B4", "OGG1", "UNG", "EXO1", "RNF168", "MSH2", "MSH6", "MLH1", "PMS2","PMS1") & sample_info$Doublings==12 & sample_info$Clonality=="Clonal","Sample_ko"]

sub_indel_burden_9sg <- sub_indel_burden[sub_indel_burden$Sample_ko %in%sample_sel,]
sub_indel_burden_9sg$Ko_gene <- sub("\\_.*","",sub_indel_burden_9sg$Sample_ko)
all_good_melt <- melt(sub_indel_burden_9sg,c("Sample_ko","Ko_gene"))
names(all_good_melt) <- c("Sample","Ko_gene","type","freq")
all_good_summary <- ddply(all_good_melt, c("Ko_gene","type"), summarise, NChild=length(freq),mean=mean(freq),sd=sd(freq))

pdf(file=paste0("sig_sub_indel",".pdf"), onefile=TRUE,height=4,width=5.5, useDingbats=FALSE)
p <- ggplot(all_good_summary, aes(x=Ko_gene, y=mean, fill=type)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.8) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                width=.2,                    # Width of the error bars
                position=position_dodge(.8))
p <- p+scale_x_discrete(limits = c("MSH2","MLH1","MSH6","PMS2","EXO1","RNF168","OGG1","PMS1","UNG","ATP2B4"))+scale_fill_manual(values=c("#33FF99", "#0099FF"))
p <- p+theme(axis.text.x=element_text(angle=45,size=10,colour = "black",hjust=0.9,vjust=1),
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



