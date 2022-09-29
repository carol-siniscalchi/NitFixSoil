library(caper)
library(geiger)
library(stringr)
library(phytools)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stargazer)


#######################################################
## Load tree
######################################################

tree = read.tree("nitfix.finalgenbank_Feb2021_onlygenusspecies_renamed.noduplicates.tre")
tree <- force.ultrametric(tree, method = "extend")

#######################################################
## Fetch datasets
######################################################


temperature <- read.csv("./../environment_final_species_averages/kew_names_carolcorrected/BIOCLIM_1.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(temperature) <- c("species", "temperature")
temperature <- distinct(temperature, species, .keep_all= TRUE)

bio3 <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/BIOCLIM_3.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio3) <- c("species", "bio3")
bio3 <- distinct(bio3, species, .keep_all= TRUE)

bio4 <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/BIOCLIM_4.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio4) <- c("species", "bio4")
bio4 <- distinct(bio4, species, .keep_all= TRUE)

precipitation <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/BIOCLIM_12.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(precipitation) <- c("species", "precipitation")
precipitation <- distinct(precipitation, species, .keep_all= TRUE)

bio15 <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/BIOCLIM_15.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio15) <- c("species", "bio15")
bio15 <- distinct(bio15, species, .keep_all= TRUE)

bio17 <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/BIOCLIM_17.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio17) <- c("species", "bio17")
bio17 <- distinct(bio17, species, .keep_all= TRUE)

elevation <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/GTOPO30_ELEVATION.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(elevation) <- c("species", "elevation")
elevation <- distinct(elevation, species, .keep_all= TRUE)

nitrogen <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/ISRICSOILGRIDS_new_average_nitrogen_reduced.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(nitrogen) <- c("species", "nitrogen")
nitrogen <- distinct(nitrogen, species, .keep_all= TRUE)

carbon <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(carbon) <- c("species", "carbon")
carbon <- distinct(carbon, species, .keep_all= TRUE)

ph <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/ISRICSOILGRIDS_new_average_phx10percent_reduced.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(ph) <- c("species", "ph")
ph <- distinct(ph, species, .keep_all= TRUE)

coarsefragment <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(coarsefragment) <- c("species", "coarsefragment")
coarsefragment <- distinct(coarsefragment, species, .keep_all= TRUE)

sand <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/ISRICSOILGRIDS_new_average_sandpercent_reduced.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(sand) <- c("species", "sand")
sand <- distinct(sand, species, .keep_all= TRUE)

needleleaf <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/LandCover_1_Needleleaf.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(needleleaf) <- c("species", "needleleaf")
needleleaf <- distinct(needleleaf, species, .keep_all= TRUE)

deciduousbroadleaf <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/LandCover_3_Deciduousbroadleaf.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(deciduousbroadleaf) <- c("species", "deciduousbroadleaf")
deciduousbroadleaf <- distinct(deciduousbroadleaf, species, .keep_all= TRUE)

herbaceous <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/LandCover_6_Herbaceous.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(herbaceous) <- c("species", "herbaceous")
herbaceous <- distinct(herbaceous, species, .keep_all= TRUE)

aridity <- read.table("./../environment_final_species_averages/kew_names_carolcorrected/aridity_index_UNEP.average.final.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
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
nodulation <- read.csv("../misc/new.nod.status.sep21.csv")
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


####################################################################################
# Measure phylogenetic signal        
####################################################################################

aridity <- as.numeric(combined$aridity)
names(aridity) <- combined$species
phylosig(tree, aridity, method = "lambda", test = TRUE)

temperature <- as.numeric(combined$temperature)
names(temperature) <- combined$species
phylosig(tree, temperature, method = "lambda", test = TRUE)

precipitation <- as.numeric(combined$precipitation)
names(precipitation) <- combined$species
phylosig(tree, precipitation, method = "lambda", test = TRUE)

nitrogen <- as.numeric(combined$nitrogen)
names(nitrogen) <- combined$species
phylosig(tree, nitrogen, method = "lambda", test = TRUE)

carbon <- as.numeric(combined$carbon)
names(carbon) <- combined$species
phylosig(tree, carbon, method = "lambda", test = TRUE)

tempseasonality <- as.numeric(combined$bio4
names(tempseasonality) <- combined$species
phylosig(tree, tempseasonality, method = "lambda", test = TRUE)

#precipseasonality <- combined$bio15
#names(precipseasonality) <- combined$species
#phylosig(tree, precipseasonality, method = "lambda", test = TRUE)

bio17 <- combined$bio17
names(bio17) <- combined$species
phylosig(tree, bio17, method = "lambda", test = TRUE)


####################################################################################
# Fit evolutionary models of nodulating and non-nodulating taxa separately         
####################################################################################

# Prepare segregated datasets

data.reduced.nod <- data.reduced[which(data.reduced$nodulation=='Yes'),]
row.names(data.reduced.nod) <- data.reduced.nod$species
treedata_object = treedata(tree, data.reduced.nod)
tree.nod <- treedata_object$phy

data.reduced.nonnod <- data.reduced[which(data.reduced$nodulation=='No'),]
row.names(data.reduced.nonnod) <- data.reduced.nonnod$species
treedata_object = treedata(tree, data.reduced.nonnod)
tree.nonnod <- treedata_object$phy




# Variables

nitrogen.nod <- as.numeric(data.reduced.nod$nitrogen)
names(nitrogen.nod) <- data.reduced.nod$species
nodulation_nitrogen_BM <- fitContinuous(tree.nod, nitrogen.nod, model = "BM")
nodulation_nitrogen_OU <- fitContinuous(tree.nod, nitrogen.nod, model = "OU") # If OU is taking a long time and huge amounts of RAM, force.ultrametric was not applied
nodulation_nitrogen_EB <- fitContinuous(tree.nod, nitrogen.nod, model = "EB")

nitrogen.nonnod <- as.numeric(data.reduced.nonnod$nitrogen)
names(nitrogen.nonnod) <- data.reduced.nonnod$species
nonnodulation_nitrogen_BM <- fitContinuous(tree.nonnod, nitrogen.nonnod, model = "BM")
nonnodulation_nitrogen_OU <- fitContinuous(tree.nonnod, nitrogen.nonnod, model = "OU") # If OU is taking a long time and huge amounts of RAM, force.ultrametric was not applied
nonnodulation_nitrogen_EB <- fitContinuous(tree.nonnod, nitrogen.nonnod, model = "EB")

nodulation_nitrogen_OU_constrained <- fitContinuous(tree.nod, nitrogen.nod, model = "OU", bounds=list(sigsq=c(272774.351383,272774.351383)))
library(lmtest)
lrtest(nodulation_nitrogen_OU, nodulation_nitrogen_OU_constrained)

## Results of above
## Suggests minimal difference between non- and nodulators in nitrogen niche

#nodulation_nitrogen_BM
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 205353.581662
#	z0 = 193.324738
#
# model summary:
#	log-likelihood = -59477.933063
#	AIC = 118959.866127
#	AICc = 118959.867607
#	free parameters = 2
#
#nonnodulation_nitrogen_BM
#GEIGER-fitted comparative model of continuous dataontinuous data
# fitted ‘BM’ model parameters:
#	sigsq = 237003.429622
#	z0 = 199.672115
#
# model summary:
#	log-likelihood = -48163.290507
#	AIC = 96330.581013
#	AICc = 96330.582930
#	free parameters = 2

#nodulation_nitrogen_OU
#GEIGER-fitted comparative model of continuous datadel of continuous data
# fitted ‘OU’ model parameters:s:
#	alpha = 2.718282
#	sigsq = 227254.419268
#	z0 = 138.587935
#
# model summary:
#	log-likelihood = -53339.66097540818
#	AIC = 106685.321951
#	AICc = 106685.324911
#	free parameters = 3
#	
#nonnodulation_nitrogen_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 272774.351383
#	z0 = 187.233221
#
# model summary:
#	log-likelihood = -42178.288151
#	AIC = 84362.576303
#	AICc = 84362.580137
#	free parameters = 3
#
#nodulation_nitrogen_EB
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 205312.680559
#	z0 = 193.324803
#
# model summary:
#	log-likelihood = -59477.935274
#	AIC = 118961.870547
#	AICc = 118961.873507
#	free parameters = 3
#
#nonnodulation_nitrogen_EB
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 237027.445202
#	z0 = 199.672159
#
# model summary:
#	log-likelihood = -48163.294360
#	AIC = 96332.588721
#	AICc = 96332.592554
#	free parameters = 3

#lrtest(nodulation_nitrogen_OU, nodulation_nitrogen_OU_constrained)
#Likelihood ratio test
#
#Model 1: nodulation_nitrogen_OU
#Model 2: nodulation_nitrogen_OU_constrained
#  #Df LogLik Df  Chisq Pr(>Chisq)    
#1   3 -53340                         
#2   3 -53404  0 128.56  < 2.2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


aridity.nod <- as.numeric(data.reduced.nod$aridity)
names(aridity.nod) <- data.reduced.nod$species
nodulation_aridity_BM <- fitContinuous(tree.nod, aridity.nod, model = "BM")
nodulation_aridity_OU <- fitContinuous(tree.nod, aridity.nod, model = "OU")
nodulation_aridity_EB <- fitContinuous(tree.nod, aridity.nod, model = "EB")

aridity.nonnod <- as.numeric(data.reduced.nonnod$aridity)
names(aridity.nonnod) <- data.reduced.nonnod$species
nonnodulation_aridity_BM <- fitContinuous(tree.nonnod, aridity.nonnod, model = "BM")
nonnodulation_aridity_OU <- fitContinuous(tree.nonnod, aridity.nonnod, model = "OU")
nonnodulation_aridity_EB <- fitContinuous(tree.nonnod, aridity.nonnod, model = "EB")

nodulation_aridity_OU_constrained <- fitContinuous(tree.nod, aridity.nod, model = "OU", bounds=list(sigsq=c(7.718612,7.718612)))
library(lmtest)
lrtest(nodulation_aridity_OU, nodulation_aridity_OU_constrained)

## Results
## Suggests aridity is a stable niche for nodulators but not nonnodulators
#
#nodulation_aridity_BM
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 4.307529
#	z0 = 1.084237
#
# model summary:
#	log-likelihood = -15786.197983
#	AIC = 31576.395966
#	AICc = 31576.397446
#	free parameters = 2
#
#nonnodulation_aridity_BM
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 6.514966
#	z0 = 1.125218
#
# model summary:
#	log-likelihood = -15271.878500
#	AIC = 30547.756999
#	AICc = 30547.758916
#	free parameters = 2

#nodulation_aridity_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 4.969321
#	z0 = 0.695261
#
# model summary:
#	log-likelihood = -9816.582468
#	AIC = 19639.164937
#	AICc = 19639.167897
#	free parameters = 3
#	
#nonnodulation_aridity_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 7.718612
#	z0 = 1.042896
#
# model summary:
#	log-likelihood = -9377.586680
#	AIC = 18761.173360
#	AICc = 18761.177194
#	free parameters = 3
#
#
#nodulation_aridity_EB
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 4.308504
#	z0 = 1.084238
#
# model summary:
#	log-likelihood = -15786.199683
#	AIC = 31578.399367
#	AICc = 31578.402327
#	free parameters = 3
#	
#nonnodulation_aridity_EB
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 6.515586
#	z0 = 1.125217
#
# model summary:
#	log-likelihood = -15271.882240
#	AIC = 30549.764481
#	AICc = 30549.768315
#	free parameters = 3

#lrtest(nodulation_aridity_OU, nodulation_aridity_OU_constrained)
#Likelihood ratio test
#
#Model 1: nodulation_aridity_OU
#Model 2: nodulation_aridity_OU_constrained
#  #Df   LogLik Df Chisq Pr(>Chisq)    
#1   3  -9816.6                        
#2   3 -16359.5  0 13086  < 2.2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


temp.nod <- as.numeric(data.reduced.nod$temp)
names(temp.nod) <- data.reduced.nod$species
nodulation_temp_BM <- fitContinuous(tree.nod, temp.nod, model = "BM")
nodulation_temp_OU <- fitContinuous(tree.nod, temp.nod, model = "OU")
nodulation_temp_EB <- fitContinuous(tree.nod, temp.nod, model = "EB")

temp.nonnod <- as.numeric(data.reduced.nonnod$temp)
names(temp.nonnod) <- data.reduced.nonnod$species
nonnodulation_temp_BM <- fitContinuous(tree.nonnod, temp.nonnod, model = "BM")
nonnodulation_temp_OU <- fitContinuous(tree.nonnod, temp.nonnod, model = "OU")
nonnodulation_temp_EB <- fitContinuous(tree.nonnod, temp.nonnod, model = "EB")

nodulation_temp_OU_constrained <- fitContinuous(tree.nod, temp.nod, model = "OU", bounds=list(sigsq=c(51845.489330,51845.489330)))
library(lmtest)
lrtest(nodulation_temp_OU, nodulation_temp_OU_constrained)


### Results
#
#nodulation_temp_BM
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 51924.418088
#	z0 = 157.268437
#
# model summary:
#	log-likelihood = -53901.160051
#	AIC = 107806.320101
#	AICc = 107806.321581
#	free parameters = 2
#	
#nonnodulation_temp_BM 
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 33969.898751
#	z0 = 176.770350
#
# model summary:
#	log-likelihood = -42079.067878
#	AIC = 84162.135755
#	AICc = 84162.137672
#	free parameters = 2
#
#nodulation_temp_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 65823.442222
#	z0 = 167.130984
#
# model summary:
#	log-likelihood = -48313.895123
#	AIC = 96633.790247
#	AICc = 96633.793207
#	free parameters = 3
#	
#nonnodulation_temp_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 51845.489330
#	z0 = 163.376907
#
# model summary:
#	log-likelihood = -36977.988443
#	AIC = 73961.976886
#	AICc = 73961.980720
#	free parameters = 3
#
#nodulation_temp_EB
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 51933.554212
#	z0 = 157.268557
#
# model summary:
#	log-likelihood = -53901.162365
#	AIC = 107808.324729
#	AICc = 107808.327689
#	free parameters = 3
#	
#nonnodulation_temp_EB 
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 33970.209856
#	z0 = 176.770189
#
# model summary:
#	log-likelihood = -42079.071768
#	AIC = 84164.143536
#	AICc = 84164.147370
#	free parameters = 3

#lrtest(nodulation_temp_OU, nodulation_temp_OU_constrained)
#Likelihood ratio test
#
#Model 1: nodulation_temp_OU
#Model 2: nodulation_temp_OU_constrained
#  #Df LogLik Df Chisq Pr(>Chisq)    
#1   3 -48314                        
#2   3 -53901  0 11175  < 2.2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

precipitation.nod <- as.numeric(data.reduced.nod$precipitation)
names(precipitation.nod) <- data.reduced.nod$species
nodulation_precipitation_BM <- fitContinuous(tree.nod, precipitation.nod, model = "BM")
nodulation_precipitation_OU <- fitContinuous(tree.nod, precipitation.nod, model = "OU")
nodulation_precipitation_EB <- fitContinuous(tree.nod, precipitation.nod, model = "EB")

precipitation.nonnod <- as.numeric(data.reduced.nonnod$precipitation)
names(precipitation.nonnod) <- data.reduced.nonnod$species
nonnodulation_precipitation_BM <- fitContinuous(tree.nonnod, precipitation.nonnod, model = "BM")
nonnodulation_precipitation_OU <- fitContinuous(tree.nonnod, precipitation.nonnod, model = "OU")
nonnodulation_precipitation_EB <- fitContinuous(tree.nonnod, precipitation.nonnod, model = "EB")

nodulation_precipitation_OU_constrained <- fitContinuous(tree.nod, precipitation.nod, model = "OU", bounds=list(sigsq=c(8254270.462273,8254270.462273)))
library(lmtest)
lrtest(nodulation_precipitation_OU, nodulation_precipitation_OU_constrained)


### Results
#
#nodulation_precipitation_BM
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 4249612.957849
#	z0 = 1384.618324
#
# model summary:
#	log-likelihood = -71767.004624
#	AIC = 143538.009248
#	AICc = 143538.010727
#	free parameters = 2
#
#nonnodulation_precipitation_BM
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 6005019.798384
#	z0 = 1479.155346
#
# model summary:
#	log-likelihood = -58286.780244
#	AIC = 116577.560488
#	AICc = 116577.562404
#	free parameters = 2
#
#nodulation_precipitation_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 5592190.527506
#	z0 = 966.097752
#
# model summary:
#	log-likelihood = -66331.256619
#	AIC = 132668.513238
#	AICc = 132668.516198
#	free parameters = 3
#
#nonnodulation_precipitation_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 8254270.462273
#	z0 = 1350.970372
#
# model summary:
#	log-likelihood = -52857.909404
#	AIC = 105721.818808
#	AICc = 105721.822642
#	free parameters = 3
#
#nodulation_precipitation_EB
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 4250187.027691
#	z0 = 1384.620863
#
# model summary:
#	log-likelihood = -71767.004757
#	AIC = 143540.009514
#	AICc = 143540.012474
#	free parameters = 3
#
#nonnodulation_precipitation_EB
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 6005978.756292
#	z0 = 1479.153793
#
# model summary:
#	log-likelihood = -58286.784020
#	AIC = 116579.568040
#	AICc = 116579.571874
#	free parameters = 3	

#lrtest(nodulation_precipitation_OU, nodulation_precipitation_OU_constrained)
#Likelihood ratio test
#
#Model 1: nodulation_precipitation_OU
#Model 2: nodulation_precipitation_OU_constrained
#  #Df LogLik Df Chisq Pr(>Chisq)    
#1   3 -66331                        
#2   3 -72417  0 12172  < 2.2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


tempseasonality.nod <- as.numeric(data.reduced.nod$bio4)
names(tempseasonality.nod) <- data.reduced.nod$species
nodulation_tempseasonality_BM <- fitContinuous(tree.nod, tempseasonality.nod, model = "BM")
nodulation_tempseasonality_OU <- fitContinuous(tree.nod, tempseasonality.nod, model = "OU")
nodulation_tempseasonality_EB <- fitContinuous(tree.nod, tempseasonality.nod, model = "EB")

tempseasonality.nonnod <- as.numeric(data.reduced.nonnod$bio4)
names(tempseasonality.nonnod) <- data.reduced.nonnod$species
nonnodulation_tempseasonality_BM <- fitContinuous(tree.nonnod, tempseasonality.nonnod, model = "BM")
nonnodulation_tempseasonality_OU <- fitContinuous(tree.nonnod, tempseasonality.nonnod, model = "OU")
nonnodulation_tempseasonality_EB <- fitContinuous(tree.nonnod, tempseasonality.nonnod, model = "EB")

nodulation_tempseasonality_OU_constrained <- fitContinuous(tree.nod, tempseasonality.nod, model = "OU", bounds=list(sigsq=c(95547446.675058,95547446.675058)))
library(lmtest)
lrtest(nodulation_tempseasonality_OU, nodulation_tempseasonality_OU_constrained)


### Results
### Evolutionary rates vastly higher in nodulators (unstable) and ancestrally much more seasonal overall
#
#nodulation_tempseasonality_BM
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 155989423.932920
#	z0 = 3774.021932
#
# model summary:
#	log-likelihood = -86380.611316
#	AIC = 172765.222631
#	AICc = 172765.224111
#	free parameters = 2
#
#nonnodulation_tempseasonality_BM 
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 61750448.076385
#	z0 = 2852.522623
#
# model summary:
#	log-likelihood = -65585.923391
#	AIC = 131175.846782
#	AICc = 131175.848698
#	free parameters = 2
#
#nodulation_tempseasonality_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 176117049.543342
#	z0 = 3253.349361
#
# model summary:
#	log-likelihood = -80323.555331
#	AIC = 160653.110662
#	AICc = 160653.113622
#	free parameters = 3
#	
#nonnodulation_tempseasonality_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 95547446.675058
#	z0 = 3322.948201
#
# model summary:
#	log-likelihood = -60527.841335
#	AIC = 121061.682670
#	AICc = 121061.686504
#	free parameters = 3
#
#nodulation_tempseasonality_EB
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 155967596.046498
#	z0 = 3774.017683
#
# model summary:
#	log-likelihood = -86380.612691
#	AIC = 172767.225382
#	AICc = 172767.228342
#	free parameters = 3
#
#nonnodulation_tempseasonality_EB 
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 61760824.416610
#	z0 = 2852.525238
#
# model summary:
#	log-likelihood = -65585.927313
#	AIC = 131177.854626
#	AICc = 131177.858460
#	free parameters = 3

#lrtest(nodulation_tempseasonality_OU, nodulation_tempseasonality_OU_constrained)
#Likelihood ratio test
#
#Model 1: nodulation_tempseasonality_OU
#Model 2: nodulation_tempseasonality_OU_constrained
#  #Df LogLik Df Chisq Pr(>Chisq)    
#1   3 -80324                        
#2   3 -86638  0 12629  < 2.2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
	


carbon.nod <- as.numeric(data.reduced.nod$carbon)
names(carbon.nod) <- data.reduced.nod$species
nodulation_carbon_BM <- fitContinuous(tree.nod, carbon.nod, model = "BM")
nodulation_carbon_OU <- fitContinuous(tree.nod, carbon.nod, model = "OU")
nodulation_carbon_EB <- fitContinuous(tree.nod, carbon.nod, model = "EB")

carbon.nonnod <- as.numeric(data.reduced.nonnod$carbon)
names(carbon.nonnod) <- data.reduced.nonnod$species
nonnodulation_carbon_BM <- fitContinuous(tree.nonnod, carbon.nonnod, model = "BM")
nonnodulation_carbon_OU <- fitContinuous(tree.nonnod, carbon.nonnod, model = "OU")
nonnodulation_carbon_EB <- fitContinuous(tree.nonnod, carbon.nonnod, model = "EB")

nodulation_carbon_OU_constrained <- fitContinuous(tree.nod, carbon.nod, model = "OU", bounds=list(sigsq=c(14935.411586,14935.411586)))
library(lmtest)
lrtest(nodulation_carbon_OU, nodulation_carbon_OU_constrained)


### Results
#
#
#nodulation_carbon_BM
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 9666.907878
#	z0 = 51.655148
#
# model summary:
#	log-likelihood = -47082.696787
#	AIC = 94169.393575
#	AICc = 94169.395055
#	free parameters = 2
#
#nonnodulation_carbon_BM
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 12564.233451
#	z0 = 48.878039
#
# model summary:
#	log-likelihood = -38963.916184
#	AIC = 77931.832369
#	AICc = 77931.834286
#	free parameters = 2
#
#nodulation_carbon_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 10879.783404
#	z0 = 30.683454
#
# model summary:
#	log-likelihood = -41012.812333
#	AIC = 82031.624667
#	AICc = 82031.627627
#	free parameters = 3
#	
#nonnodulation_carbon_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 14935.411586
#	z0 = 44.668001
#
# model summary:
#	log-likelihood = -33080.110706
#	AIC = 66166.221411
#	AICc = 66166.225245
#	free parameters = 3
#	
#nodulation_carbon_EB
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 9667.618337
#	z0 = 51.655218
#
# model summary:
#	log-likelihood = -47082.698881
#	AIC = 94171.397762
#	AICc = 94171.400722
#	free parameters = 3
#
#nonnodulation_carbon_EB
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 12566.025225
#	z0 = 48.878082
#
# model summary:
#	log-likelihood = -38963.920088
#	AIC = 77933.840176
#	AICc = 77933.844010
#	free parameters = 3

#lrtest(nodulation_carbon_OU, nodulation_carbon_OU_constrained)
#Likelihood ratio test
#
#Model 1: nodulation_carbon_OU
#Model 2: nodulation_carbon_OU_constrained
#  #Df LogLik Df Chisq Pr(>Chisq)    
#1   3 -41013                        
#2   3 -47416  0 12807  < 2.2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


bio17.nod <- as.numeric(data.reduced.nod$bio17)
names(bio17.nod) <- data.reduced.nod$species
nodulation_bio17_BM <- fitContinuous(tree.nod, bio17.nod, model = "BM")
nodulation_bio17_OU <- fitContinuous(tree.nod, bio17.nod, model = "OU")
nodulation_bio17_EB <- fitContinuous(tree.nod, bio17.nod, model = "EB")

bio17.nonnod <- as.numeric(data.reduced.nonnod$bio17)
names(bio17.nonnod) <- data.reduced.nonnod$species
nonnodulation_bio17_BM <- fitContinuous(tree.nonnod, bio17.nonnod, model = "BM")
nonnodulation_bio17_OU <- fitContinuous(tree.nonnod, bio17.nonnod, model = "OU")
nonnodulation_bio17_EB <- fitContinuous(tree.nonnod, bio17.nonnod, model = "EB")

nodulation_bio17_OU_constrained <- fitContinuous(tree.nod, bio17.nod, model = "OU", bounds=list(sigsq=c(275668.573155,275668.573155)))
library(lmtest)
lrtest(nodulation_bio17_OU, nodulation_bio17_OU_constrained)


### Results
#
#nodulation_bio17_BM
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 146225.290742
#	z0 = 157.773021
#
# model summary:
#	log-likelihood = -58100.577175
#	AIC = 116205.154350
#	AICc = 116205.155830
#	free parameters = 2
#
#nonnodulation_bio17_BM
#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 215567.383006
#	z0 = 147.474331
#
# model summary:
#	log-likelihood = -47866.373570
#	AIC = 95736.747140
#	AICc = 95736.749056
#	free parameters = 2
#
#nodulation_bio17_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 173226.605508
#	z0 = 78.973867
#
# model summary:
#	log-likelihood = -52238.577847
#	AIC = 104483.155694
#	AICc = 104483.158655
#	free parameters = 3
#
#nonnodulation_bio17_OU
#GEIGER-fitted comparative model of continuous data
# fitted ‘OU’ model parameters:
#	alpha = 2.718282
#	sigsq = 275668.573155
#	z0 = 135.462281
#
# model summary:
#	log-likelihood = -42211.344886
#	AIC = 84428.689773
#	AICc = 84428.693607
#	free parameters = 3
#	
#nodulation_bio17_EB
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 146240.813461
#	z0 = 157.773319
#
# model summary:
#	log-likelihood = -58100.579116
#	AIC = 116207.158233
#	AICc = 116207.161193
#	free parameters = 3
#	
#nonnodulation_bio17_EB
#GEIGER-fitted comparative model of continuous data
# fitted ‘EB’ model parameters:
#	a = -0.000001
#	sigsq = 215584.285742
#	z0 = 147.474141
#
# model summary:
#	log-likelihood = -47866.377386
#	AIC = 95738.754772
#	AICc = 95738.758606
#	free parameters = 3

#lrtest(nodulation_bio17_OU, nodulation_bio17_OU_constrained)
#Likelihood ratio test
#
#Model 1: nodulation_bio17_OU
#Model 2: nodulation_bio17_OU_constrained
#  #Df LogLik Df Chisq Pr(>Chisq)    
#1   3 -52239                        
#2   3 -58768  0 13058  < 2.2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#############################################################################
## Niche breadth
############################################################################

table <- read.csv("../environment_final_niche_breadth/niche_breadth_kewnames.csv")
table[table<=0] <- NA
breadth_table <- merge(combined.fixed, table, by = "species")

# Model fit
linear.model <-lm(nichebreadth ~ as.numeric(nitrogen), breadth_table)
exp.model <-lm(nichebreadth ~ !is.infinite(exp(as.numeric(nitrogen))), data = breadth_table)

# Model choice
AIC(linear.model)
AIC(exp.model)

summary(linear.model)

# Model fit
linear.model <-lm(nichebreadth ~ as.numeric(aridity), breadth_table)
exp.model <-lm(nichebreadth ~ !is.infinite(exp(as.numeric(aridity))), data = breadth_table)

# Model choice
AIC(linear.model)
AIC(exp.model)

summary(linear.model)