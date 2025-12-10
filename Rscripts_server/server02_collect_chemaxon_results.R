# 2022.06.18. Anna Lovrics
# collect results from rocs, put into rds file
# 
library(ggplot2)
library(gridExtra) # to merge plots in a grid
library(tibble)
library(dplyr)
library(readr) # to read files into tibble
library(cinf) # to read sd files
# library(TopKLists) # for rank aggregations
# library(RCSM) # for GSEA1Scores

server <- T
if (server){
  shome <- '/bigdisk/users/alovrics'
} else {
  shome <- '/media/anna/bigdisk/Anna/Enzimhome_local' 
}


# data structure
EpiGene <- file.path(shome,'Projects/EpiGene')
Malaria <- file.path(shome,'Projects/Malaria')
# subfolders
data_phase <-  file.path(Malaria,'Data/CLUE')
data_calc <- file.path(Malaria,'Data/CLUE_calculated')
figfolder <- file.path(Malaria,'Figures/Final')
#######

# functions in external files
rscriptfolder <-  file.path(Malaria,'Rscripts')
source(file.path(rscriptfolder, 'utils.R'))

#######################################
# read in drug instances
instances_file <- file.path(data_calc, "reduced_instances.rds")
siginfo_reduced <- readRDS(file=instances_file)
drugs <- unique(siginfo_reduced$drug_name_small)
instances_overlap_file <- file.path(data_calc, "reduced_overlap_instances.rds")
instances_overlap <- readRDS(file=instances_overlap_file)
drugs_overlap <- unique(instances_overlap$drug_name_small)
#######################################

# read in results
ql <-  list.files(path =file.path(data_calc, "chemaxon_input"),
                  pattern = glob2rx("query*.txt"))
qfile1 <- ql[1]
resdf <- read_delim(file=file.path(data_calc, "chemaxon_input", qfile1),
                     delim=";")
colnames(resdf)[1] <- "target"
for (qfile in ql[2:length(ql)]){
  resdf_ext <-  read_delim(file=file.path(data_calc, "chemaxon_input", qfile),
                           delim=";")
  colnames(resdf_ext)[1] <- "target"
  resdf <- dplyr::full_join(resdf, resdf_ext)
}

# compare with CFP_simil_merged.rds

# now select the drugs that have been included in the final selection round
colnames(resdf) <- drug2lower(colnames(resdf))
resdf$target <- drug2lower(resdf$target)
length(setdiff(drugs, resdf$target))

final_file <- file.path(data_calc, "chemaxon_simil.rds")
saveRDS(resdf, file=final_file)
#nai <- which(is.na(resdf), arr.ind = TRUE)
