# Based on FiSSE example file

library(ape)
library(phangorn)
library(diversitree)
library(geiger)
library(phytools)


setwd("~/Desktop/NitFixSoil-main/") # Change to top-level directory of GitHub repository
getwd()


source("traitDependent_functions.R") # Get from https://github.com/macroevolution/fisse

tree = read.tree("./phylogenetic_trees/nitfix.finalgenbank_Feb2021_onlygenusspecies_renamed.noduplicates.tre")


# Check for ultrametricity and if binary; FiSSE has its own function for ultrametricity but this should be faster and fine with slight precision errors 
if(is.binary(tree)) {
	print("TRUE")
	} else {
	tree <- multi2di(tree)
	}

if(is.ultrametric(tree)) {
	print("TRUE")
	} else {
	tree <- force.ultrametric(tree, method="extend")
	}


data = as.data.frame(tree$tip.label)
names(data) = c("species")
data$genus <- data$species
library(stringr)
data$genus <- str_replace(data$genus, "_.*", "")


nodulation <- read.csv("./misc/new.nod.status.sep21.csv")
change_list <- setNames(as.character(nodulation$Genus_Nodulation_Status), as.character(nodulation$Genus))

data$nodulation <- lapply(data$genus, function(x) if(any(!is.na(x))) change_list[x] else NA)
data$genus <- NULL


data <- na.omit(data)
data <- data[!(data$nodulation == "Unknown"),]
data <- data[!(data$nodulation == "variable"),]
data <- data[!(data$nodulation == "NULL"),]
data <- data[!is.null(data$nodulation),]


data <- as.data.frame(lapply(data, unlist))

lookup <- c("No" = 0, "Yes" = 1)
data$nodulation_binary <- lookup[data$nodulation]


traits <- as.numeric(data[,3])
names(traits) <- as.character(data[,1])

treedata_object <- treedata(tree, traits)
tree <- treedata_object$phy


traits <- traits[tree$tip.label]



res <- FISSE.binary(tree, traits)

pval_1tailed <- min(res$pval, 1-res$pval)
pval_1tailed

save(pval_1tailed, file = "./results/pval_1tailed.Rdata")

# pval_2tailed <- min(res$pval, 1-res$pval)*2
