# 2022.06.12. Anna Lovrics
# filter drugs obtained from the CLUE database

library(ggplot2) # for nice plots
library(gridExtra) # to merge plots in a grid
library(tibble) # simple data frames
library(dplyr) # grammar of data manipulation
library(here) # to handle file paths

# folders
data_phase <-  file.path("Data", "CLUE")
data_calc <- file.path("Data","CLUE_calculated")
if (!file.exists(here(data_calc))) dir.create(here(data_calc))
figfolder <- "Figures"

# utility functions
source(here('Rscripts', 'utils.R'))

#########################################
#           read in data files               #
#########################################
# read in signature info
siginfo_file <- here(data_phase, "siginfo_beta.txt")
siginfo_phase2 <- read.table(siginfo_file, sep="\t", header=T, quote = "") # 1202656*35
print(dim(siginfo_phase2))

# read in structure file and for each drug, select the smiles
# with most information
pertinfo_file <- file.path(data_phase, "compoundinfo_beta.txt")
pertinfo_phase2 <- read.table(pertinfo_file, sep="\t", header=T, quote = "",
                              comment.char="")
# for each of these compounds, decide which BRD-ID to use for 
# structure determination
pertinfo_phase2_filtered <- retain_best_BRDIDs(pertinfo_phase2,
                                            common_name="cmap_name",
                                            important_columns=c("pert_id", "cmap_name"),
                                            dupl_folder=getwd(),
                                            additional_ordering=NULL,
                                            writeDuplicates=FALSE)

# add a new column "drug_name" which combines cmap_name and compound_aliases
# and does not contain capitals or special characters
pertinfo_phase2_filtered$drug_name_orig <- ifelse(
  nchar(pertinfo_phase2_filtered$compound_aliases) > 0 &
    pertinfo_phase2_filtered$compound_aliases != "\"\"",
  pertinfo_phase2_filtered$compound_aliases, pertinfo_phase2_filtered$cmap_name)
pertinfo_phase2_filtered$drug_name_small <- sapply(pertinfo_phase2_filtered$drug_name, drug2lower)
# note - in some cases, compounds with different cmap_name have the same compound_aliases
# name and hence will be duplicates for drug_name
#.########################################### #
# add structure information to signature information
# IMPORTANT - merge by "cmap_name" only (the pert_id may become different, because
# of different batches of the same compound)
siginfo_merged <- left_join(siginfo_phase2,
                        pertinfo_phase2_filtered, 
                        by=c("cmap_name"))
# rename appropriate columns
colnames(siginfo_merged)[
  colnames(siginfo_merged) == "pert_id.x"] <- "pert_id"
colnames(siginfo_merged)[
  colnames(siginfo_merged) == "pert_id.y"] <- "pert_id.structure"

# print number of instances and drugs
print("number of instances:")
print(nrow(siginfo_merged))
print("number of drugs:")
print(length(unique(siginfo_merged$drug_name_orig)))
# [1] "number of instances:"
# [1] 1202656
# [1] "number of drugs:"
# [1] 33552


# save for later use
merged_file <- here(data_calc, "merged_instances.rds")
saveRDS(siginfo_merged, file=merged_file)
###################################

# finally keep signatures with good quality and high enough tas score & measured in a core cell line
tasthr_str <- '0.2'
tasthr <- as.numeric(tasthr_str)
tasthr_name <- gsub(".", "", tasthr_str, fixed=T)
core_cells <- c('A375', 'A549', 'HA1E', 'HCC515', 'HT29', 'HEPG2',
                'MCF7', 'PC3', 'VCAP')
# keep a signature only if 
# 1) tas >= 0.2
# 2) high quality
# 3) measured in a core cell line
siginfo_reduced <- siginfo_merged[which(siginfo_merged$tas>=tasthr &
                                      siginfo_merged$is_hiq>=1 &
                                      siginfo_merged$cell_iname %in% core_cells),]

# remove compounds with unknown structures
ui <- which(siginfo_reduced$canonical_smiles 
              %in% c("restricted", "", "\"\""))
siginfo_reduced <- siginfo_reduced[-ui,]

# print number of retained instances and drugs
print("number of instances:")
print(nrow(siginfo_reduced))
print("number of drugs:")
print(length(unique(siginfo_reduced$drug_name_orig)))
# check for canononical smiles as well - there is a bit of discrepancy -> check
print("number of canonical smiles")
print(length(unique(siginfo_reduced$canonical_smiles)))

# NEW IN 2023 January - keep drugs with at least 5 instances
# count the number of occurrences
num_occurences <- unlist(sort(table(siginfo_reduced$drug_name_orig)))
keep_drugs <- names(num_occurences)[which(num_occurences >=5 )]
siginfo_reduced <- siginfo_reduced[siginfo_reduced$drug_name_orig %in% keep_drugs,]


# print number of retained instances and drugs
print("number of instances:")
print(nrow(siginfo_reduced))
print("number of drugs:")
print(length(unique(siginfo_reduced$drug_name_orig)))
# check for canononical smiles as well - there is a bit of discrepancy -> check
print("number of canonical smiles")
print(length(unique(siginfo_reduced$canonical_smiles)))

############
# find signatures where drug_name is the same but canonical smiles differ
tmp_can <- siginfo_reduced
tmp_can <- tmp_can[-which(duplicated(tmp_can$canonical_smiles)),]
ddif <- which(duplicated(tmp_can$drug_name_orig))
name_comp <- c()
for (ddi in ddif){
  dname <- tmp_can$drug_name_orig[ddi]
  if (!dname %in% name_comp) name_comp <- c(name_comp, dname)
  ddii <- which(tmp_can$drug_name_orig == dname)
  # print(dname)
  # print(tmp_can[ddii, c("cmap_name", "canonical_smiles", "pert_id.structure",
  #                       "compound_aliases", "drug_name_orig")])
  # print("############")
}
# print the compounds with stereoisomers
print(name_comp)
print("percentage of stereoisomers:")
print(100*length(name_comp)/length(unique(siginfo_reduced$drug_name_orig)))

# which smiles string is a duplicate??? - none
si <- which(!duplicated(siginfo_reduced$pert_id.structure))
tmp <- siginfo_reduced[si,]
ssi <- which(duplicated(tmp$canonical_smiles))
ssi2 <- which(tmp$canonical_smiles %in% tmp$canonical_smiles[ssi])
tmp_final <- tmp[ssi2, c("cmap_name", "canonical_smiles", "pert_id.structure",
             "compound_aliases")]
print(tmp_final[order(tmp_final$canonical_smiles),])

#############################################################################
# for the final dataset, choose from the BRD-IDs again 
# for each duplicate, find best BRD-ID

# find duplicates to decide what to do:
di <- which(duplicated(siginfo_reduced$cmap_name))
tt <- siginfo_reduced[-di,]
ttdupl <- which(duplicated(tt$drug_name_small))
ttdrug <- unique(tt$drug_name_small[ttdupl])
siginfo_reduced_tmp <- siginfo_reduced
for (drug in ttdrug){
  tt_tmp <- tt[tt$drug_name_small == drug,]
  tt_red <- retain_best_BRDIDs(tt_tmp,
                               common_name="drug_name_orig",
                               important_columns=c("pert_id", "drug_name_orig"),
                               dupl_folder=getwd(),
                               additional_ordering=NULL,
                               writeDuplicates=FALSE)
  winner <- which(siginfo_reduced_tmp$pert_id == tt_red$pert_id)
  loser <- which(siginfo_reduced_tmp$pert_id %in% setdiff(tt_tmp$pert_id, tt_red$pert_id))
  # replace the canonical smiles but keep the cmap_name-s & pert_id-s
  new_smile <- siginfo_reduced_tmp$canonical_smiles[winner[1]]
  siginfo_reduced_tmp$canonical_smiles[loser] <- new_smile
}
siginfo_reduced <- siginfo_reduced_tmp
#############################################################################
# write smiles to file for structure based comparisons
if (!file.exists(data_calc)) dir.create(data_calc)
# drug sturctures
smiles_file <- file.path(data_calc, "reduced_canonical_smiles.smiles")
smiles_drugs2 <- siginfo_reduced[,
                    c('canonical_smiles', 'pert_id.structure', 'drug_name_orig')]
#di <- which(duplicated(smiles_drugs2$pert_id.structure))
di <- which(duplicated(smiles_drugs2$drug_name_orig))
smiles_drugs2 <- smiles_drugs2[-di,]
write.table(smiles_drugs2, file=smiles_file, row.names = F, quote=F,
            col.names = F, sep="\t")

# also write to file for oeinput
smiles_oeinput_file <- file.path(data_calc,
                              "reduced_canonical_smiles_oeinput.smi")
write.table(siginfo_reduced[-di,c('canonical_smiles','drug_name_small')],
            file=smiles_oeinput_file ,
            sep=" ", row.names = F, col.names = F, quote=F)

# also save the reduced
instances_file <- file.path(data_calc, "reduced_instances.rds")
saveRDS(siginfo_reduced, file=instances_file)

################################################
# table of cell types
list_ctypes <- list()
for (ctype in core_cells){
  siginfo_tmp <- siginfo_reduced[siginfo_reduced$cell_iname == ctype,]
  list_ctypes[[ctype]] <- c("num_instances" = nrow(siginfo_tmp),
                            "num_drugs" = length(unique(siginfo_tmp$drug_name_orig)))
}
table_ctypes <- as.data.frame(list_ctypes)
print("number of instances and drugs measured in each core cell line")
print(table_ctypes)
# A375 A549 HA1E HCC515 HT29 HEPG2 MCF7  PC3 VCAP
# num_instances 5574 5431 4535   2543 4099  2349 8474 6697 1448
# num_drugs     1444 1185 1320    976 1215   844 1494 1516  623

# ###############################################################################
# Define anti-malarial drugs based on WHO and tse2019past
WHO <- c("Amodiaquine", "Artemether", "Lumefantrine", "Artesunate", 
                   "Mefloquine", "Pyronaridine", "Chloroquine", "Artenimol",
                   "Piperaquine", "Doxycycline", "Primaquine", "Quinine", 
                   "Sulfadoxine", "Pyrimethamine", "Proguanil")
WHO <- drug2lower(WHO)
pastAM <- c("Artemisinin", "Artemotil", "Atovaquone",  "Halofantrine", 
           "Hydroxychloroquine", "Mepacrine" , "Quinidine",  "Tafenoquine",
            "cinchonidine", "cinchonine")
pastAM <- drug2lower(pastAM)
antimal_drugs <- sort(union(WHO, pastAM))
# check the number of antimalarial drugs present in the database
antimal_drugs <- antimal_drugs[which(antimal_drugs %in% 
                                 siginfo_reduced$drug_name_small)]

# write AM drugs smiles seperately to files
am_smiles_file <- file.path(data_calc, "AMdrugs_reduced.smiles")
ami <- match(antimal_drugs, drug2lower(smiles_drugs2$drug_name_orig))
am_smiles_drugs2 <- smiles_drugs2[ami,]
write.table(am_smiles_drugs2, file=am_smiles_file, row.names = F, quote=F,
            col.names = F, sep="\t")
