# NitFixSoil
Data from "Testing the evolutionary drivers of nitrogen-fixing symbioses in challenging soil environments" (https://doi.org/10.1101/2022.09.27.509719).

## File Structure:

### environment_final_niche_breadth
This folder contains niche breadth measures. Taxa are named according to WCVP nomenclature.

* niche_breadth.csv: final niche breadth calculated per species
* environment_final_species_ranges: contains environmental files with the range value of each variable per species

### environment_final_species_averages: 
Files containing the average values of environmental variables per taxa. These files were used to generate the plots shown in Figure 1 and other analysis in the paper. Notables file include:

* [all.layers.combined.matching.tree.csv](https://github.com/carol-siniscalchi/NitFixSoil/blob/main/environment_final_species_averages/all.layers.combined.matching.tree.csv): final product of environmental variables matched to the phylogenetic sampling, containing average values of all environmental variables per taxa, with assigned nodulation status and type and subfamily assignment for legume species. 
* [all_nod_status_subfam_final.csv](https://github.com/carol-siniscalchi/NitFixSoil/blob/main/environment_final_species_averages/all_nod_status_subfam_final.csv): contains nodulation status and type per genus, plus subfamily assignment for legume species.

### misc
* [dropped.tips.after.environdata.csv](https://github.com/carol-siniscalchi/NitFixSoil/blob/main/misc/dropped.tips.after.environdata.csv): file containing the names of tips that were sampled in the phylogeny but dropped from analysis because they lacked environmental data
* [new.nod.status.sep21.csv](https://github.com/carol-siniscalchi/NitFixSoil/blob/main/misc/new.nod.status.sep21.csv): genus-level nodulation assignment, with subfamily and subanalysis assignment for legumes. Subfamily assignment follows [LPGW 2017](https://doi.org/10.12705/661.3)

### phylogenetic_trees:
Contains several versions of the phylogenetic tree used in the study.
* backbone_100taxa.tre: phylogenetic backbone used to graft the subtrees
* nitfix.finalgenbank_Feb2021_onlygenusspecies_renamed.noduplicates.tre and [nitfix.finalgenbank_Feb2021_onlygenusspecies_renamed.noduplicates.tre.pdf](https://github.com/carol-siniscalchi/NitFixSoil/blob/main/phylogenetic_trees/nitfix.finalgenbank_Feb2021_onlygenusspecies_renamed.noduplicates.tre.pdf): final dated tree used for analyses
* treepl_nitfix_genbank.txt: configuration file for TreePL
* nitfix_total_beforeTreePL_kewnames.tre: final tree before callibration