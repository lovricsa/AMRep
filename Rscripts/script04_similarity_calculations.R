# 2022.06.18. Anna Lovrics
# collect results from rocs, put into rds file
# 
# 2022.06.27. - now work with phase 2 data only

library(tibble) # for data matrices (better than dataframe)
library(tidyr)
library(ggplot2)
library(gridExtra) # to merge plots in a grid
library(dplyr)
library(cmapR) # to read the cmap gctx files
library("easyPubMed") # for automated PubMed search
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

opt <- list()
opt$genetype <- "BING"

###################################
# read in antimalarial drugs
am_smiles_file <- file.path(data_calc, "AMdrugs_reduced.smiles")
am_smiles <- read.table(am_smiles_file, header=F, sep="\t")
am_drugs <- am_smiles[,3]

# read in overlap with moa drugs (this is needed for the ROC analysis)
instances_overlap_file <- file.path(data_calc, "reduced_overlap_instances.rds")
instances_overlap <- readRDS(file=instances_overlap_file)
overlap_drugs <- unique(instances_overlap$drug_name_small)

# read in aggregated drug instances
rank_files <- list.files(path=data_calc, pattern="reduced_aggregated_drug_ranks.rds")

# define chunk size
chunk_size <- 50
num_genes <- 150
for (rfile in rank_files){
  ctype <- gsub("_reduced_aggregated_drug_ranks.rds", "", rfile)
  # check whether calculation was already done
  simil_file <- file.path(data_calc, 
                paste(ctype, "_reduced_aggregated_drug_simil.rds", sep=""))
  if (file.exists(simil_file)) next
  # otherwise start similarity calculations
  rank_drug_df <- readRDS(file.path(data_calc, rfile))
  # select indices of overlap drugs
  overlap_cols <- which(colnames(rank_drug_df) %in% overlap_drugs)
  # select indices of antimalarial drugs
  am_cols <- which(colnames(rank_drug_df) %in% am_drugs) # note that halofantrine is missing
  # union of columns
  select_cols <- union(overlap_cols, am_cols) # calculate similarites to both am_drugs and drugs with known MOA

  # calculate similarity between select_cols and rest of the matrix - do in parts
  chunk_col_index <- 1
  simil_res <- list()
  while ((chunk_col_index-1)*chunk_size < length(select_cols)){ # loop for final columns
    firstcol <- (chunk_col_index-1)*chunk_size + 1
    lastcol <- min(chunk_col_index*chunk_size, length(select_cols))
    scols <- select_cols[firstcol:lastcol]
    print("start of a new set of column indices similarity calculation")
    print(c(chunk_col_index, firstcol, lastcol))
    simil_res[[chunk_col_index]] <- list()
    chunk_row_index <- 1
    while ((chunk_row_index-1)*chunk_size < ncol(rank_drug_df)){ # loop for final rows
      firstrow <- (chunk_row_index-1)*chunk_size + 1
      lastrow <- min(chunk_row_index*chunk_size, ncol(rank_drug_df))
      print("start of a new set of row indices similarity calculation")
      print(c(chunk_row_index, firstrow, lastrow))
      simil_res[[chunk_col_index]][[chunk_row_index]] <- simES(
                                         rank_drug_df[,firstrow:lastrow], 
                                         rank_drug_df[,scols],  num_genes = num_genes)
      simil_res[[chunk_col_index]][[chunk_row_index]] <- as.data.frame(
                                simil_res[[chunk_col_index]][[chunk_row_index]])
      print("a row similarity calculation finished")
      print("--------")
      chunk_row_index <- chunk_row_index + 1
    }
    print("a column similarity calculation finished")
    print("###########")
    # put together the row chunks
    simil_res[[chunk_col_index]] <- bind_rows(simil_res[[chunk_col_index]])
    # print dimension of the chunk
    print(dim(simil_res[[chunk_col_index]]))
    chunk_col_index <- chunk_col_index + 1
  }
  # put together the column chunks
  simil_res <- bind_cols(simil_res)
  # save the obtained similarity matrix
  saveRDS(simil_res, file=simil_file)
}

simil_files <- list.files(path=data_calc, pattern="reduced_aggregated_drug_simil.rds")
sfile <- simil_files[3]
simil_res <- readRDS(file.path(data_calc, sfile))

# save similarity among AM drugs
g <- ggcorrplot::ggcorrplot(simil_res[am_drugs, am_drugs], lab=T, hc.order=T,
                            legend.title = "Zhang-similarity")
figname <- file.path(figfolder, 
                     paste("Zhang_similarity_AM_drugs_", "all",  ".pdf", sep=""))
ggsave(figname, height=7, width=7)

# use TS compounds
# read in TS drugs
TS_file <- file.path(data_phase, "touchstone.csv")
TS <- read.csv(TS_file, header=TRUE, sep=';')
TS$name_lower <- sapply(TS$Name, drug2lower)
TS$id_lower <- sapply(TS$ID, drug2lower)
# find duplicates
name_dupl <- TS$name_lower[which(duplicated(TS$name_lower))]
dupl_ind <- which(TS$name_lower %in% name_dupl)
# TS[dupl_ind,]

TS_orig <- TS
TS <- TS[-which(duplicated(TS$name_lower)),]
for (dname in name_dupl){
  brd_vector <- TS_orig$ID[which(TS_orig$name_lower == dname)]
  new_brd <- paste(brd_vector, collapse="|")
  TS$ID[which(TS$name_lower == dname)] <- new_brd
}

# calculate heatmap for each TS class
heatmapfolder <- file.path(figfolder, "Heatmaps")
if (!file.exists(heatmapfolder)) dir.create(heatmapfolder)
class_simil <- c()
ii <- 0
for (cmap.class in unique(TS$CMap.Class)){
  ii <- ii+1
  if (cmap.class=="") next()
  sel_drugs <- TS$name_lower[TS$CMap.Class==cmap.class]
  sel_drugs <- sel_drugs[which(sel_drugs %in% colnames(simil_res))]
  if (length(sel_drugs)<=1){
    # print(ii)
    # print(cmap.class)
    # print("not in colnames")
    next
  }
  g <- ggcorrplot::ggcorrplot(simil_res[sel_drugs, sel_drugs], lab=T, hc.order=T,
                              legend.title = "Zhang-similarity")
  cmap.class.str <- paste0(cmap.class, sep="_")
  cmap.class.str <- gsub("/", "_or_", cmap.class.str)
  figname <- file.path(heatmapfolder, 
          paste("Zhang_similarity_",cmap.class.str, "_drugs_", "all",  ".pdf", sep=""))
  ggsave(figname, height=7, width=7)
  class_simil <- c(class_simil, 
          convert_dist2vector(as.matrix(simil_res[sel_drugs, sel_drugs])))
}
#mean(simil_res)

# do we get similar results based on the repurposing hub?
moa_file <- file.path(data_calc, "moa_simil.rds")
moa_simil <- readRDS(moa_file)
# moa_simil_wide <- pivot_wider(moa_simil[,1:3], names_from=drug2, values_from=moa_simil)
# select drug pairs, where moa_simil is one
rowdrugs <- moa_simil$drug1[which(moa_simil$moa_simil==1)]
rowinds <- match(rowdrugs, rownames(simil_res))
coldrugs <- moa_simil$drug2[which(moa_simil$moa_simil==1)]
colinds <- match(coldrugs, colnames(simil_res))
#mean(as.matrix(simil_res)[rowinds, colinds])

aggr_all1 <- convert_dist2vector(as.matrix(simil_res)[colnames(simil_res), colnames(simil_res)])
aggr_all2 <- as.vector(as.matrix(simil_res[setdiff(rownames(simil_res), colnames(simil_res)), ]))
aggr_all <- as.data.frame(c(aggr_all1, aggr_all2))
colnames(aggr_all) <- "Zhang"
aggr_moa <- as.data.frame(as.vector(as.matrix(simil_res)[cbind(rowinds, colinds)]))
colnames(aggr_moa) <- "Zhang"
aggr_pcl <- as.data.frame(class_simil)
colnames(aggr_pcl) <- "Zhang"
#sim_thr <- quantile(aggr_all$Zhang, probs=c(0.05,0.95))
sim_thr <- quantile(aggr_all$Zhang, probs=c(0.99))
g <- ggplot(aggr_all, aes(Zhang))+
  geom_density(aes(color='all similarity values'))+
  geom_density(data=aggr_moa, aes(Zhang, color='shared MoA'))+
  geom_density(data=aggr_pcl, aes(Zhang, color='shared PCL class'))+
  xlab("Zhang similarity")+
  scale_color_manual(name='',
                     breaks=c('all similarity values', 'shared MoA', 'shared PCL class'),
                     values=c('all similarity values'='black', 'shared MoA'='darkgreen', 'shared PCL class'='red'))+
  geom_vline(xintercept = sim_thr, linetype="dashed")

figname <- file.path(figfolder, 
                     paste("Density_all_vs_moa_simil",  ".pdf", sep=""))
ggsave(figname, height=5, width=5)


# for each AM drug, select top hits
am_hits_all <- list()
am_hits <- list()
for (am_drug in am_drugs){
  am_hits_all[[am_drug]] <- rownames(simil_res)[which(simil_res[,am_drug, drop=T]>=sim_thr)]
  am_hits_all[[am_drug]] <- am_hits_all[[am_drug]][!(am_hits_all[[am_drug]] %in% am_drugs)]
  am_hits[[am_drug]] <- am_hits_all[[am_drug]][am_hits_all[[am_drug]] %in% TS$name_lower]
  am_hits[[am_drug]] <- TS$Name[match(am_hits[[am_drug]], TS$name_lower)]
}
am_hits_vec <- unique(Reduce(c, am_hits_all))
am_hits_TS <- unique(Reduce(c, am_hits))
# save hits
hitfile_all <- file.path(data_calc, "am_hits_all.rds")
saveRDS(am_hits_vec, file=hitfile_all)
hitfile <- file.path(data_calc, "am_hits.rds")
saveRDS(am_hits, file=hitfile)
# output all hits including Zhang similarities & whether it was TS
hits_df <- simil_res[am_hits_vec, am_drugs]
hits_df <- as_tibble(round(hits_df, digits=2), rownames="hit_drug")
hits_df$TS <- "no"
hits_df$TS[drug2lower(hits_df$hit_drug) %in% TS$name_lower] <- "yes"
# add moa information as well
instances_hits_df <- instances_orig[which(instances_orig$drug_name_small 
                                       %in% drug2lower(am_hits_vec)),
                                 c("cmap_name", "canonical_smiles", "drug_name_small", "moa")]
instances_hits_df <- instances_hits_df[!duplicated(instances_hits_df),]
hits_df <- inner_join(instances_hits_df, hits_df,  by=c("drug_name_small"="hit_drug"))
hits_df <- hits_df[order(hits_df$cmap_name), c("canonical_smiles", "cmap_name", "moa", "TS", am_drugs)]
hit_suppfile <- file.path(data_calc, "am_hits_supplement.csv")
write.csv(hits_df, file=hit_suppfile, row.names=F)
# also save to open with mvieww
hit_smilesfile <- gsub(".csv", ".smiles", hit_suppfile)
write.table(hits_df, file=hit_smilesfile, row.names=F, col.names = F,
             quote=F, sep="\t")

################################
# now finally create a PubMed search
drugs_found <- c()
for (drug in am_hits_TS){
  my_query <- paste(drug, "AND malaria", sep=" ")
  my_entrez_id <- get_pubmed_ids(my_query)
  if (as.numeric(my_entrez_id$Count) > 0) drugs_found <- c(drugs_found, drug)
  print(drug)
}

# # join them into their PCL-s 
# TS_found <- TS[TS$Name %in% drugs_found,]
# for (pcl.class in unique(TS_found$CMap.Class)){
#   print(pcl.class)
#   print(TS_found[TS_found$CMap.Class==pcl.class,])
#   invisible(readline(prompt="Press [enter] to continue"))
# }

abstracts_dir <- file.path(data_calc, "pubmed_abstracts_final")
if (!file.exists(abstracts_dir)) dir.create(abstracts_dir)
# write abstract to text file
for (drug in drugs_found){
  abstracts_file <- file.path(abstracts_dir,
                              paste("pubmed_abstracts_", drug, ".txt", sep=""))
  if (file.exists(abstracts_file)) unlink(abstracts_file)
  sink(abstracts_file, append=TRUE, split=TRUE)
  my_query <- paste(drug, "AND malaria", sep=" ")
  my_entrez_id <- get_pubmed_ids(my_query)
  my_abstracts_txt <- fetch_pubmed_data(my_entrez_id, format = "abstract")
  my_abstracts_txt[my_abstracts_txt == ""] <- "\n"
  cat("\n", file=abstracts_file, append=T)
  cat(drug, file=abstracts_file, append=T)
  cat("\n===========================\n", file=abstracts_file, append=T)
  #cat(my_abstracts_txt, file=abstracts_file, append=T)
  cat(paste(my_abstracts_txt, collapse = ""), file=abstracts_file, append=T)
  cat("\n===========================", file=abstracts_file, append=T)
  cat("\n===========================", file=abstracts_file, append=T)
  closeAllConnections()  
} 

# save structures using the wfpaper_final10_save_molecule2Dstructure.R 
# this means this script needs to be run in two seperate halves!
# start a new script instead!
instances_file <- file.path(data_calc, "reduced_instances.rds")
instances_orig  <- readRDS(file=instances_file)
am_hits_TS <- setdiff(am_hits_TS, am_drugs)
instances_hits <- instances_orig[which(instances_orig$drug_name_small 
                                       %in% drug2lower(am_hits_TS)),
                                 c("cmap_name", "drug_name_small", "moa")]
instances_hits <- instances_hits[!duplicated(instances_hits),]

# add MoA based on Drug Repurposing Hub - check whether this is the same
moa_file <- file.path(Malaria, 'Data/RepurposingHub',
                      'repurposing_drugs_20200324.txt')
moa <- read.table(moa_file, header=T, sep='\t', quote="", comment.char = "!")
moa <- as_tibble(moa)
moa$name_lower <- sapply(moa$pert_iname, drug2lower)

# join with previous results
instances_hits <- left_join(instances_hits, moa,
                                  join_by(drug_name_small==name_lower))
hit_file <- file.path(data_calc, 'AMdrugs_hits_final.txt')
write.table(instances_hits , file = hit_file, sep=";",
            row.names=F, col.names=T, quote=T)

####################################################################

# create a latex table using these information
# question - which moa to use???
outfile <- gsub(".txt", ".tex", hit_file, fixed=T)
instances_hits$structure <- paste('\\includegraphics[width=0.2\\textwidth]{Structures/',
                      instances_hits$drug_name_small, '.png}', sep="")
if (file.exists(outfile)){
  print("latex table already printed")
} else {
  latex_table <- instances_hits[, c("structure",
                                            "cmap_name",
                                            "moa.x",
                                            "moa.y",
                                            "target",
                                            "disease_area",
                                            "clinical_phase")]
  print(xtable(latex_table, type = "latex"), file = outfile,
        include.rownames = F,sanitize.text.function=function(x){x})
}

###############################################################################
# this could only be done later
published <- c("bithionol", "cyclosporin-a", "gossypol", "ivermectin", "lonidamine",
               "methotrexate", "nelfinavir", "sorafenib", "thapsigargin",
               "thiazolopyrimidine", "triamterene", "tunicamycin")
am_publ <- list()
for (dname in names(am_hits)){
  am_publ[[dname]] <- am_hits[[dname]][am_hits[[dname]] %in% published]
}
###############################################################################
# save these results as well
am_hits_publ <- Reduce(c, am_publ)
instances_publ <- instances_orig[which(instances_orig$drug_name_small 
                                       %in% drug2lower(am_hits_publ)),
                                 c("cmap_name", "drug_name_small", "moa")]
instances_publ <- instances_publ[!duplicated(instances_publ),]

# add MoA based on Drug Repurposing Hub - check whether this is the same
moa_file <- file.path(Malaria, 'Data/RepurposingHub',
                      'repurposing_drugs_20200324.txt')
moa <- read.table(moa_file, header=T, sep='\t', quote="", comment.char = "!")
moa <- as_tibble(moa)
moa$name_lower <- sapply(moa$pert_iname, drug2lower)

# join with previous results
instances_publ <- left_join(instances_publ, moa,
                            join_by(drug_name_small==name_lower))
publ_file <- file.path(data_calc, 'AMdrugs_hits_final_published.txt')
write.table(instances_publ , file = publ_file, sep=";",
            row.names=F, col.names=T, quote=T)

####################################################################

# create a latex table using these information
# question - which moa to use???
outfile <- gsub(".txt", ".tex", publ_file, fixed=T)
instances_publ$structure <- paste('\\includegraphics[width=0.2\\textwidth]{Structures/',
                                  instances_publ$drug_name_small, '.png}', sep="")
if (file.exists(outfile)){
  print("latex table already printed")
} else {
  latex_table <- instances_publ[, c("structure",
                                    "cmap_name",
                                    "moa.x",
                                    "moa.y",
                                    "target",
                                    "disease_area",
                                    "clinical_phase")]
  print(xtable(latex_table, type = "latex"), file = outfile,
        include.rownames = F,sanitize.text.function=function(x){x})
}