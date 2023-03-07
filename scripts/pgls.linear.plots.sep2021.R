library(caper)
library(geiger)
library(stringr)
library(phytools)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stargazer)

setwd("~/Desktop/NitFixSoil-main/") # Change to top-level directory of GitHub repository
getwd()



###########
## Load tree
##########

tree = read.tree("./phylogenetic_trees/nitfix.finalgenbank_Feb2021_onlygenusspecies_renamed.noduplicates.tre")



###########
## Fetch datasets
##########


temperature <- read.csv("./environment_final_species_averages/kew_names_carolcorrected/BIOCLIM_1.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(temperature) <- c("species", "temperature")
temperature <- distinct(temperature, species, .keep_all= TRUE)

bio3 <- read.table("./environment_final_species_averages/kew_names_carolcorrected/BIOCLIM_3.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio3) <- c("species", "bio3")
bio3 <- distinct(bio3, species, .keep_all= TRUE)

bio4 <- read.table("./environment_final_species_averages/kew_names_carolcorrected/BIOCLIM_4.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio4) <- c("species", "bio4")
bio4 <- distinct(bio4, species, .keep_all= TRUE)

precipitation <- read.table("./environment_final_species_averages/kew_names_carolcorrected/BIOCLIM_12.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(precipitation) <- c("species", "precipitation")
precipitation <- distinct(precipitation, species, .keep_all= TRUE)

bio15 <- read.table("./environment_final_species_averages/kew_names_carolcorrected/BIOCLIM_15.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio15) <- c("species", "bio15")
bio15 <- distinct(bio15, species, .keep_all= TRUE)

bio17 <- read.table("./environment_final_species_averages/kew_names_carolcorrected/BIOCLIM_17.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio17) <- c("species", "bio17")
bio17 <- distinct(bio17, species, .keep_all= TRUE)

elevation <- read.table("./environment_final_species_averages/kew_names_carolcorrected/GTOPO30_ELEVATION.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(elevation) <- c("species", "elevation")
elevation <- distinct(elevation, species, .keep_all= TRUE)

nitrogen <- read.table("./environment_final_species_averages/kew_names_carolcorrected/ISRICSOILGRIDS_new_average_nitrogen_reduced.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(nitrogen) <- c("species", "nitrogen")
nitrogen <- distinct(nitrogen, species, .keep_all= TRUE)

carbon <- read.table("./environment_final_species_averages/kew_names_carolcorrected/ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(carbon) <- c("species", "carbon")
carbon <- distinct(carbon, species, .keep_all= TRUE)

ph <- read.table("./environment_final_species_averages/kew_names_carolcorrected/ISRICSOILGRIDS_new_average_phx10percent_reduced.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(ph) <- c("species", "ph")
ph <- distinct(ph, species, .keep_all= TRUE)

coarsefragment <- read.table("./environment_final_species_averages/kew_names_carolcorrected/ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(coarsefragment) <- c("species", "coarsefragment")
coarsefragment <- distinct(coarsefragment, species, .keep_all= TRUE)

sand <- read.table("./environment_final_species_averages/kew_names_carolcorrected/ISRICSOILGRIDS_new_average_sandpercent_reduced.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(sand) <- c("species", "sand")
sand <- distinct(sand, species, .keep_all= TRUE)

needleleaf <- read.table("./environment_final_species_averages/kew_names_carolcorrected/LandCover_1_Needleleaf.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(needleleaf) <- c("species", "needleleaf")
needleleaf <- distinct(needleleaf, species, .keep_all= TRUE)

deciduousbroadleaf <- read.table("./environment_final_species_averages/kew_names_carolcorrected/LandCover_3_Deciduousbroadleaf.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(deciduousbroadleaf) <- c("species", "deciduousbroadleaf")
deciduousbroadleaf <- distinct(deciduousbroadleaf, species, .keep_all= TRUE)

herbaceous <- read.table("./environment_final_species_averages/kew_names_carolcorrected/LandCover_6_Herbaceous.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(herbaceous) <- c("species", "herbaceous")
herbaceous <- distinct(herbaceous, species, .keep_all= TRUE)

aridity <- read.table("./environment_final_species_averages/kew_names_carolcorrected/aridity_index_UNEP.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(aridity) <- c("species", "aridity")
aridity <- distinct(aridity, species, .keep_all= TRUE)

combined = merge(temperature, bio3, by = "species")
combined = merge(combined, bio4, by = "species")
combined = merge(combined, precipitation, by = "species")
combined = merge(combined, bio15, by = "species")
combined = merge(combined, bio17, by = "species")
combined = merge(combined, elevation, by = "species")
combined = merge(combined, nitrogen, by = "species")
combined = merge(combined, carbon, by = "species")
combined = merge(combined, ph, by = "species")
combined = merge(combined, sand, by = "species")
combined = merge(combined, coarsefragment, by = "species")
combined = merge(combined, needleleaf, by = "species")
combined = merge(combined, deciduousbroadleaf, by = "species")
combined = merge(combined, herbaceous, by = "species")
combined = merge(combined, aridity, by = "species")


combined$genus <- combined$species
combined$genus <- str_replace(combined$genus, "_.*", "")

#add nodulation status, type and subfamily
nodulation <- read.csv("./misc/new.nod.status.sep21.csv")
change_list <- setNames(as.character(nodulation$Genus_Nodulation_Status), as.character(nodulation$Genus))
combined$nodulation <- lapply(combined$genus, function(x) if(any(!is.na(x))) change_list[x] else NA)

change_list <- setNames(as.character(nodulation$Type), as.character(nodulation$Genus))
combined$nodulationtype <- lapply(combined$genus, function(x) if(any(!is.na(x))) change_list[x] else NA)

change_list <- setNames(as.character(nodulation$Subfamily), as.character(nodulation$Genus))
combined$subfamily <- lapply(combined$genus, function(x) if(any(!is.na(x))) change_list[x] else NA)


combined.fixed <- na.omit(combined)
combined.fixed <- combined.fixed[!(combined.fixed$nodulation == "Unknown"),]
combined.fixed <- combined.fixed[!(combined.fixed$nodulation == "variable"),]
combined.fixed <- combined.fixed[!(combined.fixed$nodulation == "NULL"),]
combined.fixed <- combined.fixed[!is.null(combined.fixed$nodulation),]
combined.fixed <- combined.fixed[!(combined.fixed$nodulationtype == "Unknown"),]
combined.fixed <- combined.fixed[!(combined.fixed$nodulationtype == "variable"),]
combined.fixed <- combined.fixed[!(combined.fixed$nodulationtype == "NULL"),]
combined.fixed <- combined.fixed[!is.null(combined.fixed$nodulationtype),]


combined.fixed <- as.data.frame(lapply(combined.fixed, unlist))
combined.fixed$nodulation <- as.factor(combined.fixed$nodulation)
combined.fixed$nodulationtype <- as.factor(combined.fixed$nodulationtype)
rownames(combined.fixed) <- combined.fixed$species
combined.fixed$genus <- NULL
names(combined.fixed)[names(combined.fixed) == 'x'] <- 'combined'

row.names(combined.fixed) <- combined.fixed$species
treedata_object = treedata(tree, combined.fixed)

tree.reduced <- treedata_object$phy
data.reduced <- treedata_object$data
data.reduced <- data.frame(matrix(unlist(data.reduced), nrow=nrow(data.reduced)), stringsAsFactors=FALSE)
colnames(data.reduced) = c(colnames(combined.fixed))
data.reduced$nodulation <- as.factor(data.reduced$nodulation)
data.reduced$nodulationtype <- as.factor(data.reduced$nodulationtype)

write.csv(data.reduced, file = "./results/all.layers.combined.sep21.csv")



###### violin plots ######

#change y axis to numeric
data.reduced$aridity <- as.numeric(data.reduced$aridity)
data.reduced$bio3 <- as.numeric(data.reduced$bio3)
data.reduced$bio4 <- as.numeric(data.reduced$bio4)
data.reduced$bio15 <- as.numeric(data.reduced$bio15)
data.reduced$bio17 <- as.numeric(data.reduced$bio17)
data.reduced$carbon <- as.numeric(data.reduced$carbon)
data.reduced$coarsefragment <- as.numeric(data.reduced$coarsefragment)
data.reduced$deciduousbroadleaf <- as.numeric(data.reduced$deciduousbroadleaf)
data.reduced$elevation <- as.numeric(data.reduced$elevation)
data.reduced$herbaceous <- as.numeric(data.reduced$herbaceous)
data.reduced$needleleaf<- as.numeric(data.reduced$needleleaf)
data.reduced$nitrogen <- as.numeric(data.reduced$nitrogen)
data.reduced$ph <- as.numeric(data.reduced$ph)
data.reduced$precipitation <- as.numeric(data.reduced$precipitation)
data.reduced$sand <- as.numeric(data.reduced$sand)
data.reduced$temperature <- as.numeric(data.reduced$temperature)

#plots envVSnodstatus
carb.plot <- ggplot(data.reduced, aes(x = nodulation, y = carbon, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$carbon), (max(data.reduced$carbon))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
ari.plot <- ggplot(data.reduced, aes(x = nodulation, y = aridity, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$aridity), (max(data.reduced$aridity))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
bio3.plot <- ggplot(data.reduced, aes(x = nodulation, y = bio3, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$bio3), (max(data.reduced$bio3))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
bio4.plot <- ggplot(data.reduced, aes(x = nodulation, y = bio4, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$bio4), (max(data.reduced$bio4))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
bio15.plot <- ggplot(data.reduced, aes(x = nodulation, y = bio15, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$bio15), (max(data.reduced$bio15))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
bio17.plot <- ggplot(data.reduced, aes(x = nodulation, y = bio17, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$bio17), (max(data.reduced$bio17))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
coar.plot <- ggplot(data.reduced, aes(x = nodulation, y = coarsefragment, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$coarsefragment), (max(data.reduced$coarsefragment))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
dec.plot <- ggplot(data.reduced, aes(x = nodulation, y = deciduousbroadleaf, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$deciduousbroadleaf), (max(data.reduced$deciduousbroadleaf))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
ele.plot <- ggplot(data.reduced, aes(x = nodulation, y = elevation, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$elevation), (max(data.reduced$elevation))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
herb.plot <- ggplot(data.reduced, aes(x = nodulation, y = herbaceous, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$herbaceous), (max(data.reduced$herbaceous))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
need.plot <- ggplot(data.reduced, aes(x = nodulation, y = needleleaf, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$needleleaf), (max(data.reduced$needleleaf))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
nit.plot <- ggplot(data.reduced, aes(x = nodulation, y = nitrogen, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$nitrogen), (max(data.reduced$nitrogen))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
ph.plot <- ggplot(data.reduced, aes(x = nodulation, y = ph, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$ph), (max(data.reduced$ph))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
prec.plot <- ggplot(data.reduced, aes(x = nodulation, y = precipitation, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$precipitation), (max(data.reduced$precipitation))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
sand.plot <- ggplot(data.reduced, aes(x = nodulation, y = sand, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$sand), (max(data.reduced$sand))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")
tem.plot <- ggplot(data.reduced, aes(x = nodulation, y = temperature, fill = nodulation)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$temperature), (max(data.reduced$temperature))) + scale_fill_manual(values=c("#F1A340","#998EC3")) + theme(legend.position="none")

pdf(file="envVSnod.pdf", width = 20, height = 20)
ggarrange(ari.plot, bio3.plot, bio4.plot, bio15.plot, bio17.plot, carb.plot, coar.plot, dec.plot, ele.plot, herb.plot, need.plot, nit.plot, ph.plot, prec.plot, sand.plot, tem.plot, 
          ncol = 4, nrow = 4)
dev.off()

#plots envVSnodtype
carbtype.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = carbon, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$carbon), (max(data.reduced$carbon))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
aritype.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = aridity, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$aridity), (max(data.reduced$aridity))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
bio3type.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = bio3, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$bio3), (max(data.reduced$bio3))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
bio4type.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = bio4, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$bio4), (max(data.reduced$bio4))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
bio15type.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = bio15, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$bio15), (max(data.reduced$bio15))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
bio17type.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = bio17, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$bio17), (max(data.reduced$bio17))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
coartype.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = coarsefragment, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$coarsefragment), (max(data.reduced$coarsefragment))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
dectype.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = deciduousbroadleaf, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$deciduousbroadleaf), (max(data.reduced$deciduousbroadleaf))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
eletype.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = elevation, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$elevation), (max(data.reduced$elevation))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
herbtype.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = herbaceous, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$herbaceous), (max(data.reduced$herbaceous))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
needtype.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = needleleaf, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$needleleaf), (max(data.reduced$needleleaf))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
nittype.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = nitrogen, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$nitrogen), (max(data.reduced$nitrogen))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
phtype.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = ph, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$ph), (max(data.reduced$ph))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
prectype.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = precipitation, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$precipitation), (max(data.reduced$precipitation))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
sandtype.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = sand, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$sand), (max(data.reduced$sand))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")
temtype.plot <- ggplot(data.reduced, aes(x = nodulationtype, y = temperature, fill = nodulationtype)) + geom_violin(trim = TRUE) + ylim(min(data.reduced$temperature), (max(data.reduced$temperature))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none")


pdf(file="./results/envVSnodtype.pdf", width = 20, height = 20)
ggarrange(aritype.plot, bio3type.plot, bio4type.plot, bio15type.plot, bio17type.plot, carbtype.plot, coartype.plot, dectype.plot, eletype.plot, herbtype.plot, needtype.plot, nittype.plot, phtype.plot, prectype.plot, sandtype.plot, temtype.plot, ncol = 4, nrow = 4)
dev.off()

#### linear model ####

aridity.lin <- (lm(formula = aridity ~ nodulation, data = data.reduced))
bio15.lin <- (lm(formula = bio15 ~ nodulation, data = data.reduced))
bio17.lin <- (lm(formula = bio17 ~ nodulation, data = data.reduced))
bio3.lin <- (lm(formula = bio3 ~ nodulation, data = data.reduced))
bio4.lin <- (lm(formula = bio4 ~ nodulation, data = data.reduced))
carbon.lin <- (lm(formula = carbon ~ nodulation, data = data.reduced))
coarsefragment.lin <- (lm(formula = coarsefragment ~ nodulation, data = data.reduced))
deciduousbroadleaf.lin <- (lm(formula = deciduousbroadleaf ~ nodulation, data = data.reduced))
elevation.lin <- (lm(formula = elevation ~ nodulation, data = data.reduced))
herbaceous.lin <- (lm(formula = herbaceous ~ nodulation, data = data.reduced))
needleleaf.lin <- (lm(formula = needleleaf ~ nodulation, data = data.reduced))
nitrogen.lin <- (lm(formula = nitrogen ~ nodulation, data = data.reduced))
ph.lin <- (lm(formula = ph ~ nodulation, data = data.reduced))
precipitation.lin <- (lm(formula = precipitation ~ nodulation, data = data.reduced))
sand.lin <- (lm(formula = sand ~ nodulation, data = data.reduced))
temperature.lin <- (lm(formula = temperature ~ nodulation, data = data.reduced))

#summarize models
stargazer(aridity.lin, bio15.lin, bio17.lin, bio3.lin, bio4.lin, carbon.lin, coarsefragment.lin, deciduousbroadleaf.lin, elevation.lin, herbaceous.lin, needleleaf.lin, nitrogen.lin, ph.lin, precipitation.lin, sand.lin, temperature.lin, type ="text", title = "Linear models", align = TRUE, out = "linear.txt")

#calculate AICs
AIClin <- AIC(aridity.lin, bio15.lin, bio17.lin, bio3.lin, bio4.lin, carbon.lin, coarsefragment.lin, deciduousbroadleaf.lin, elevation.lin, herbaceous.lin, needleleaf.lin, nitrogen.lin, ph.lin, precipitation.lin, sand.lin, temperature.lin)
extractAIC(aridity.lin)

# logistic models
aridity.log <- (glm(formula = nodulation ~ aridity, data = data.reduced, family = binomial))
bio15.log <- (glm(formula = nodulation ~ bio15, data = data.reduced, family = binomial))
bio17.log <- (glm(formula = nodulation ~ bio17, data = data.reduced, family = binomial))
bio3.log <- (glm(formula = nodulation ~ bio3, data = data.reduced, family = binomial))
bio4.log <- (glm(formula = nodulation ~ bio4, data = data.reduced, family = binomial))
carbon.log <- (glm(formula = nodulation ~ carbon, data = data.reduced, family = binomial))
coarsefragment.log <- (glm(formula = nodulation ~ coarsefragment, data = data.reduced, family = binomial))
deciduousbroadleaf.log <- (glm(formula = nodulation ~ deciduousbroadleaf, data = data.reduced, family = binomial))
elevation.log <- (glm(formula = nodulation ~ elevation, data = data.reduced, family = binomial))
herbaceous.log <- (glm(formula = nodulation ~ herbaceous, data = data.reduced, family = binomial))
needleleaf.log <- (glm(formula = nodulation ~ needleleaf, data = data.reduced, family = binomial))
nitrogen.log <- (glm(formula = nodulation ~ nitrogen, data = data.reduced, family = binomial))
ph.log <- (glm(formula = nodulation ~ ph, data = data.reduced, family = binomial))
precipitation.log <- (glm(formula = nodulation ~ precipitation, data = data.reduced, family = binomial))
sand.log <- (glm(formula = nodulation ~ sand, data = data.reduced, family = binomial))
temperature.log <- (glm(formula = nodulation ~ temperature, data = data.reduced, family = binomial))
#summarize models
stargazer(aridity.log, bio15.log, bio17.log, bio3.log, bio4.log, carbon.log, coarsefragment.log, deciduousbroadleaf.log, elevation.log, herbaceous.log, needleleaf.log, nitrogen.log, ph.log, precipitation.log, sand.log, temperature.log, type ="text", title = "Logistic models", align = TRUE, out = "logistic.txt")
#calculate AICs
AIClog <- AIC(aridity.log, bio15.log, bio17.log, bio3.log, bio4.log, carbon.log, coarsefragment.log, deciduousbroadleaf.log, elevation.log, herbaceous.log, needleleaf.log, nitrogen.log, ph.log, precipitation.log, sand.log, temperature.log)

#AIC table
rows <- c("aridity" ,"bio15" ,"bio17" ,"bio3" ,"bio4" ,"carbon" ,"coarsefragment" ,"deciduousbroadleaf" ,"elevation" ,"herbaceous" ,"needleleaf" ,"nitrogen" ,"ph" ,"precipitation" ,"sand" ,"temperature")
linear <- c(AIClin$AIC)
logistic <- c(AIClog$AIC)
aic.table <- data.frame(rows, linear, logistic)
aic.table <- data.frame(rows, linear)
write.csv(aic.table, file = "nitfix.aic.table.csv")

#some plot tests

plot(data.reduced$aridity, data.reduced$nodulation)
abline((lm(formula = aridity ~ nodulation, data = data.reduced)))

range(data.reduced$aridity)
xari <- seq(0, 6, 0.01)
yari <- predict(aridity.log, list(aridity=xari, type = "response" ))
plot(data.reduced$aridity, data.reduced$nodulation, pch = 16, xlab = "aridity", ylab = "nodulation")
lines(xari, yari)

range(data.reduced$ph)
xph <- seq(21, 82, 0.01)
yph <- predict(ph.log, list(ph=xph, type = "response" ))
plot(data.reduced$ph, data.reduced$nodulation, pch = 16, xlab = "ph", ylab = "nodulation")
lines(xph, yph)

range(data.reduced$nitrogen)
xnit <- seq(0, 1386, 0.01)
ynit <- predict(nitrogen.log, list(nitrogen=xnit, type = "response" ))
plot(data.reduced$nitrogen, data.reduced$nodulation, pch = 16, xlab = "nitrogen", ylab = "nodulation")
lines(xnit, ynit)

###########
## Build PGLS model
##########


# Make caper data object, pre-calculate phylogenetic covariance matrix
#data_object <- comparative.data(tree.reduced, data.reduced, species, vcv=FALSE) # To check structure of data is right
data_object <- comparative.data(tree.reduced, data.reduced, species, vcv=TRUE)
# Save this object -- takes a long time to calculate
save(data_object, file = "./results/combined_pgls_object.Rdata")
# load("~/Desktop/dating_quick/treePL_ultrametric_smoothing100/caper_object_all.robject")

# Test whether nodulators differ in several environmental factors
#temperature
model1 <- pgls(formula = as.numeric(temperature) ~ nodulation, data = data_object)
save(model1, file = "./results/temperature_pgls_model.Rdata")
summary(model1)

#precipitation
model2 <- pgls(formula = as.numeric(precipitation) ~ nodulation, data = data_object)
save(model2, file = "./results/precipitation_pgls_model.Rdata")
summary(model2)

#elevation
model3 <- pgls(formula = as.numeric(elevation) ~ nodulation, data = data_object)
save(model3, file = "./results/elevation_pgls_model.Rdata")
summary(model3)

#nitrogen
model4 <- pgls(formula = as.numeric(nitrogen) ~ nodulation, data = data_object)
save(model4, file = "./results/nitrogen_pgls_model.Rdata")
summary(model4)

#carbon
model5 <- pgls(formula = as.numeric(carbon) ~ nodulation, data = data_object)
save(model5, file = "./results/carbon_pgls_model.Rdata")
summary(model5)

#ph
model6 <- pgls(formula = as.numeric(ph) ~ nodulation, data = data_object)
save(model6, file = "./results/ph_pgls_model.Rdata")
summary(model6)

#aridity
model7 <- pgls(formula = as.numeric(aridity) ~ nodulation, data = data_object)
save(model7, file = "aridity_pgls_model.Rdata")
summary(model7)

#coarsefragment
model8 <- pgls(formula = as.numeric(coarsefragment) ~ nodulation, data = data_object)
save(model8, file = "./results/coarsefragment_pgls_model.Rdata")
summary(model8)

#sand
model9 <- pgls(formula = as.numeric(sand) ~ nodulation, data = data_object)
save(model9, file = "./results/sand_pgls_model.Rdata")
summary(model9)

#needleleaf
model10 <- pgls(formula = as.numeric(needleleaf) ~ nodulation, data = data_object)
save(model10, file = "./results/needleleaf_pgls_model.Rdata")
summary(model10)

#deciduousbroadleaf
model11 <- pgls(formula = as.numeric(deciduousbroadleaf) ~ nodulation, data = data_object)
save(model11, file = "./results/deciduousbroadleaf_pgls_model.Rdata")
summary(model11)

#herbaceous
model12 <- pgls(formula = as.numeric(herbaceous) ~ nodulation, data = data_object)
save(model12, file = "./results/herbaceous_pgls_model.Rdata")
summary(model12)

#bio3
model13 <- pgls(formula = as.numeric(bio3) ~ nodulation, data = data_object)
save(model13, file = "./results/bio3_pgls_model.Rdata")
summary(model13)

#bio4
model14 <- pgls(formula = as.numeric(bio4) ~ nodulation, data = data_object)
save(model14, file = "./results/bio4_pgls_model.Rdata")
summary(model14)

#bio15
model15 <- pgls(formula = as.numeric(bio15) ~ nodulation, data = data_object)
save(model15, file = "./results/bio15_pgls_model.Rdata")
summary(model15)

#bio17
model16 <- pgls(formula = as.numeric(bio17) ~ nodulation, data = data_object)
save(model16, file = "./results/bio17_pgls_model.Rdata")
summary(model16)



