#tip rates in relation to nodulation and type of nodulator

library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

setwd("~/Dropbox/postdoc_MSU/linear.models.sep21/")
getwd()
DR <- read.csv("tip_DR_sep21.csv", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)


colnames(DR) <- c("species", "DR")
DR <- distinct(DR, species, .keep_all= TRUE)
rownames(DR) <- DR$species


DR$genus <- rownames(DR)

DR$genus <- str_replace(DR$genus, "_.*", "")

nodulation <- read.csv("new.nod.status.sep21.csv")
change_list <- setNames(as.character(nodulation$Genus_Nodulation_Status), as.character(nodulation$Genus))
DR$nodulation <- lapply(DR$genus, function(x) if(any(!is.na(x))) change_list[x] else NA)

change_list <- setNames(as.character(nodulation$Type), as.character(nodulation$Genus))
DR$nodulationtype <- lapply(DR$genus, function(x) if(any(!is.na(x))) change_list[x] else NA)

change_list <- setNames(as.character(nodulation$Subfamily), as.character(nodulation$Genus))
DR$subfamily <- lapply(DR$genus, function(x) if(any(!is.na(x))) change_list[x] else NA)

change_list <- setNames(as.character(nodulation$Subanalysis), as.character(nodulation$Genus))
DR$subanalysis <- lapply(DR$genus, function(x) if(any(!is.na(x))) change_list[x] else NA)


DR.fixed <- na.omit(DR)
DR.fixed <- DR.fixed[!(DR.fixed$nodulation == "Unknown"),]
DR.fixed <- DR.fixed[!(DR.fixed$nodulation == "variable"),]
DR.fixed <- DR.fixed[!(DR.fixed$nodulation == "NULL"),]
DR.fixed <- DR.fixed[!is.null(DR.fixed$nodulation),]
DR.fixed <- DR.fixed[!(DR.fixed$nodulationtype == "Unknown"),]
DR.fixed <- DR.fixed[!(DR.fixed$nodulationtype == "variable"),]
DR.fixed <- DR.fixed[!(DR.fixed$nodulationtype == "NULL"),]
DR.fixed <- DR.fixed[!is.null(DR.fixed$nodulationtype),]
DR.fixed <- DR.fixed[!(DR.fixed$subfamily == "Unknown"),]
DR.fixed <- DR.fixed[!(DR.fixed$subfamily == "variable"),]
DR.fixed <- DR.fixed[!(DR.fixed$subfamily == "NULL"),]
DR.fixed <- DR.fixed[!is.null(DR.fixed$subfamily),]
DR.fixed <- DR.fixed[!(DR.fixed$subanalysis == "Unknown"),]
DR.fixed <- DR.fixed[!(DR.fixed$subanalysis == "variable"),]
DR.fixed <- DR.fixed[!(DR.fixed$subanalysis == "NULL"),]
DR.fixed <- DR.fixed[!is.null(DR.fixed$subanalysis),]


DR.fixed <- as.data.frame(lapply(DR.fixed, unlist))
DR.fixed$nodulation <- as.factor(DR.fixed$nodulation)
DR.fixed$nodulationtype <- as.factor(DR.fixed$nodulationtype)
rownames(DR.fixed) <- DR.fixed$species
DR.fixed$genus <- NULL
names(DR.fixed)[names(DR.fixed) == 'x'] <- 'DR'


aes(width=10, height=10)
DR.type <- ggplot(DR.fixed, aes(x = nodulationtype, y = log(DR), fill = nodulationtype)) + geom_violin(trim = TRUE) + 
  ylim(min(DR.fixed$DR), log(100)) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")

DR.nod <- ggplot(DR.fixed, aes(x = nodulation, y = log(DR), fill = nodulation)) + geom_violin(trim = TRUE) + 
  ylim(min(DR.fixed$DR), log(100)) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")

subfamily<-ggplot(DR.fixed, aes(x = subfamily, y = log(DR), fill = subfamily)) + geom_violin(trim = TRUE) + 
  ylim(min(DR.fixed$DR), log(100)) + scale_fill_brewer(palette="RdYlBu") + theme(legend.position="none")

subanalysis<-ggplot(DR.fixed, aes(x = subanalysis, y = log(DR), fill = subanalysis)) + geom_violin(trim = TRUE) + 
  ylim(min(DR.fixed$DR), log(100)) + scale_fill_brewer(palette="RdYlBu") + theme(legend.position="none")

pdf(file="DR.composite.pdf", width = 24, height = 10)
ggarrange(DR.type, DR.nod, subanalysis, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

pdf(file="comp2.pdf", width = 20, height = 7)
ggarrange(nodstatus, nodtype, subfamily, 
          #labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

tree = read.tree("nitfix.finalgenbank_Feb2021_onlygenusspecies_renamed.noduplicates.tre")
svl<-read.csv("tip_DR_sep21.csv",header=TRUE,row.names=1)
log.svl <-log1p(svl)
log.svl<-setNames(log.svl[,1],rownames(log.svl))
obj<-contMap(tree,log.svl,plot=FALSE)
obj<-setMap(obj,invert=TRUE)
pdf(file="reconstruction.pdf", width = 500, height = 500)
plot(obj,fsize=c(0.1,0.1),outline=FALSE,lwd=c(1,3),leg.txt="log(SVL)")
dev.off()

# Print highest diversification rates with taxa
DR.fixed[DR.fixed$DR > 5 & DR.fixed$nodulation == "No",]
DR.fixed[DR.fixed$DR > 5 & DR.fixed$nodulation == "Yes",]
