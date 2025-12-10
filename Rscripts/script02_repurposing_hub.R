# 2022.06.30. Anna Lovrics
# filter Repurposing Hub compounds
# modified: 2024.01.17.

#############################################
library(tibble)
library(dplyr)
library(tidyr)
library(here) # to handle file paths
#############################################

# folders
data_phase <-  file.path("Data", "CLUE")
data_repurpose <-  file.path("Data", "RepurposingHub")
data_calc <- file.path("Data","CLUE_calculated")
figfolder <- "Figures"

# utility functions
source(here('Rscripts', 'utils.R'))
############################
#          part 1         #
############################

# 1) 
# read MoA downloaded from repurposing hub
moa_file <- here(data_repurpose, 'repurposing_drugs_20200324.txt')
moa <- read.table(moa_file, header=T, sep='\t', quote="", comment.char = "!")
moa <- as_tibble(moa)
moa <- moa[, c('pert_iname', 'moa', 'target')]


# keep compounds that have info on both MOA & targets
moa_keep <- which(moa$moa != "" & moa$target != "")
moa_both <- moa[moa_keep,]

# create a longer tibble, where there is only one moa or target in each cell
moa_long_list <- list()
for (ri in 1:nrow(moa_both)){
  drug <- moa_both$pert_iname[ri]
  # test if duplicates exist
  di <- which(moa_both$pert_iname == drug)
  if (length(di) != 1){
    print (paste("oops", drug))
    next
  }
  moas <- strsplit(moa_both$moa[ri], "|", fixed=T)[[1]]
  mnum <- length(moas)
  targets <- strsplit(moa_both$target[ri], "|", fixed=T)[[1]]
  tnum <- length(targets)
  moa_long_list[[drug]] <- as_tibble(cbind('pert_iname'=rep(drug, times=mnum*tnum),
                   'moa'=rep(moas, each=tnum),
                   'target'=rep(targets, times=mnum)))
}
moa_long <- do.call(rbind, moa_long_list)
moa_long$both <- paste(moa_long$moa, moa_long$target, sep="_")

# keep compounds for which at least one other compound
# with the same MoA and target exists
keep <- c()
for (drug in unique(moa_long$pert_iname)){
  di <- which(moa_long$pert_iname == drug)
  drug_both <- moa_long$both[di]
  other_both <- unique(moa_long$both[-di])
  if (length(intersect(drug_both, other_both))>0){
    keep <- c(keep, drug)
  }
}

moa_both <- moa_long[moa_long$pert_iname %in% keep,]
moa_both$drug_name_small <- sapply(moa_both$pert_iname, drug2lower)

############################
#          part 2          #
############################
# create an overlapping dataframe

# read in drug instances
instances_file <- here(data_calc, "reduced_instances.rds")
siginfo_reduced <- readRDS(file=instances_file)

drug_overlap <- intersect(moa_both$drug_name_small, 
                      siginfo_reduced$drug_name_small)
instances_overlap <- siginfo_reduced[siginfo_reduced$drug_name_small %in% drug_overlap,]
instances_overlap_file <- file.path(data_calc, "reduced_overlap_instances.rds")
saveRDS(instances_overlap, file=instances_overlap_file)
moa_overlap <- moa_both[moa_both$drug_name_small %in% drug_overlap, ]
moa_overlap_file <- file.path(data_calc, "reduced_overlap_moa.rds")
saveRDS(moa_overlap, file=moa_overlap_file)

# save moa_overlap with structures
moa_overlap_file <- file.path(data_calc,
                              "moa_overlap_structures_oeinput.smi")
# add smiles, but use the pre-calculated siginfo_reduced
if (!file.exists(moa_overlap_file)){
  dupl_structures <- which(duplicated(siginfo_reduced$drug_name_orig))
  siginfo_structures <- siginfo_reduced[-dupl_structures, c('drug_name_orig', 
                                        'drug_name_small', 'canonical_smiles')]
  
  moa_overlap_structures <- inner_join(moa_overlap, 
                                      siginfo_structures, by="drug_name_small")
  # remove duplicate rows, also the moa, target info are not needed
  moa_overlap_structures <- moa_overlap_structures[,c('drug_name_orig', 
                                           'drug_name_small', 'canonical_smiles')]
  moa_dupl <- which(duplicated(moa_overlap_structures))
  moa_overlap_structures <- moa_overlap_structures[-moa_dupl,]
  
  # write to file as smi (rocs input)
  write.table(moa_overlap_structures[,c('canonical_smiles','drug_name_small')],
              file=moa_overlap_file,
              sep=" ", row.names = F, col.names = F, quote=F)
}

#############################################

# create a (long) similarity tibble
# similar if shares 'both', dissimilar if does not
# moa_simil <- tibble('drug1'=character(), 'drug2'=character(), 
#                      'simil'=numeric())
# save the obtained list
moa_file <- file.path(data_calc, "moa_simil.rds")
if (file.exists(moa_file)){
  moa_simil <- readRDS(file=moa_file)
} else {
  numrows <- length(unique(moa_overlap$drug_name_small)) * 
    (length(unique(moa_overlap$drug_name_small))-1)/2
  moa_simil_list <- list() 
  list_index <- 0
  for (d1 in 1:(length(unique(moa_overlap$drug_name_small))-1)){
    for (d2 in (d1+1):length(unique(moa_overlap$drug_name_small))){
      list_index <- list_index + 1
      if (list_index %% 10000 == 0){
        done <- paste(as.character(list_index), '/', as.character(numrows))
        print(done)
      } 
      drug1 <- unique(moa_overlap$drug_name_small)[d1]
      drug2 <- unique(moa_overlap$drug_name_small)[d2]
      drug1_tbl <- moa_overlap[moa_overlap$drug_name_small==drug1,]
      drug2_tbl <- moa_overlap[moa_overlap$drug_name_small==drug2,]
      # drugs are "moa_similar" if at least one pair of moa & target are the same
      simil <- length(intersect(drug1_tbl$both, drug2_tbl$both)) > 0
      moa_simil_list[[list_index]] <- tibble('drug1'=drug1, 'drug2'=drug2,
                                             'moa_simil'=as.numeric(simil))
    }
  }
  moa_simil <- do.call(rbind, moa_simil_list)
  moa_simil <- as_tibble(moa_simil) %>%
    rowwise() %>%
    mutate(both_drugs = paste0(sort(c(drug1, drug2)), collapse = '_'))
  saveRDS(moa_simil, file=moa_file)
}

###########################################
# now check the numbers for each core cell line
# table of cell types
core_cells <- c('A375', 'A549', 'HA1E', 'HCC515', 'HT29', 'HEPG2',
                'MCF7', 'PC3', 'VCAP')
list_ctypes <- list()
for (ctype in core_cells){
  instances_tmp <- instances_overlap[instances_overlap$cell_iname == ctype,]
  list_ctypes[[ctype]] <- c("num_instances" = nrow(instances_tmp),
                            "num_drugs" = length(unique(instances_tmp$drug_name_orig)))
}
table_ctypes <- as.data.frame(list_ctypes)
print("number of instances and drugs measured in each core cell line and also in connectivity hub")
print(table_ctypes)

# A375 A549 HA1E HCC515 HT29 HEPG2 MCF7  PC3 VCAP
# num_instances 3682 3774 2978   1507 2642  1326 5834 4493  816
# num_drugs      665  535  616    409  534   315  662  694  250
