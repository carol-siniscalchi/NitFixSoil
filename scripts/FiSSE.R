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
# 0.05094905
# in Jun 15:  0.00999001
# in Sep 2021: 0.04695305

save(pval_1tailed, file = "./results/pval_1tailed.Rdata")

# pval_2tailed <- min(res$pval, 1-res$pval)*2


######
# HiSSE

library(hisse)

percentsampled <- length(tree$tip.label)/38564 # The denominator The Plant List species total


# Back to dataframe
traits.frame <- as.data.frame(traits)
traits.frame$names <- names(traits)
# Reorder column
traits.frame <- traits.frame[c("names", "traits")]

# Set up
turnover.anc = c(1,2,3,4) # Full HiSSE
turnover.anc.null = c(1,1,2,2) # HiSSE with turnover rates the same net turnover (CID-2)
# eps.anc = c(0,0,0,0) # Yule pure-birth
eps.anc = c(1,1,2,2) # Different extinction parameters
# Transition matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10)) # Prevent dual transitions between hidden and observed
trans.rates.nodual.equal16 = ParEqual(trans.rates.nodual, c(1,6))
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8)) # Equal rates model
trans.rates.bisse <- TransMatMaker(hidden.states=FALSE)

full_hisse = hisse(tree, traits.frame, f = c(percentsampled, percentsampled), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal, output.type="net.div")
bisse <- hisse(tree, traits.frame, hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type="net.div")
cid_2 <- hisse(tree, traits.frame, hidden.states=FALSE, turnover.anc=turnover.anc.null, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type="net.div")
cid_4 <- hisse.null4(tree, traits.frame, f = c(percentsampled, percentsampled), eps.anc=rep(1,8))

df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) = c("model", "aicc")
df[1,] = c("full_hisse", full_hisse$AICc)
df[2,] = c("bisse", bisse$AICc)
df[3,] = c("cid_2", cid_2$AICc)
df[4,] = c("cid_4", cid_4$AICc)

library(qpcR) # For quick akaike weights
df$weights <- akaike.weights(as.numeric(df$aicc))$weights


## LRT
null.logL <- cid_2$loglik
alternative.logL <- bisse$loglik

# Check degrees of freedom
pchisq(-2*(null.logL - alternative.logL) ,df=1, lower.tail=FALSE)
# 0
