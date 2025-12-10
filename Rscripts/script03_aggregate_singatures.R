# 2022.06.18. Anna Lovrics - finalized 2024.02.26.
# rank aggregation
# the input file is 33.1Gb, you may choose to work with the output files
# of this script: <ctype>_reduced_aggregated_drug_ranks.rds,
# where <ctype> is a core cell type or "all"

library(tibble)
library(dplyr)
library(readr) # to read files into tibble
library(TopKLists) # for rank aggregations
library(cmapR) # to read the cmap gctx files
library(here) # to handle file paths
#############################################

# folders
data_phase <-  file.path("Data", "CLUE")
data_repurpose <-  file.path("Data", "RepurposingHub")
data_calc <- file.path("Data","CLUE_calculated")
figfolder <- "Figures"
resfolder <- "Results"

# utility functions
source(here('Rscripts', 'utils.R'))
############################

opt <- list()
opt$genetype <- "BING"

# read in metadata
datafile <- here(data_phase,'level5_beta_trt_cp_n720216x12328.gctx')
col_meta <- read_gctx_meta(datafile, dim="column")
row_meta <- read_gctx_meta(datafile, dim="row")
geneinfo_file <- here(data_phase, "geneinfo_beta.txt") 
geneinfo <- read.table(geneinfo_file, sep="\t", header=T, quote = "") # 12328*7
# note - some ensembl id is missing which(geneinfo$ensembl_id=="\"\"")

#find the indices for the genes used
try(if (!opt$genetype %in% (c("landmark", "BING"))){
  stop('genetype should be "landmark" or "BING", please choose from these')
})

if (opt$genetype == "landmark"){
  landmark_indices <- which(geneinfo$feature_space == "landmark")
  landmark_genes <- geneinfo$gene_id[landmark_indices]
  signature_rows <- match(landmark_genes, row_meta[,1])
} 
if (opt$genetype == "BING"){
  bing_indices <- which(geneinfo$feature_space == "best inferred")
  bing_genes <- geneinfo$gene_id[bing_indices]
  signature_rows <- match(bing_genes, row_meta[,1])
}

# find the column indices of the specific signatures
instances_file <- file.path(data_calc, "reduced_instances.rds")
#siginfo_reduced <- readRDS(file=instances_file)
instances <- list()
instances[["all"]] <- readRDS(file=instances_file)

# repeat calculation for each core cell type
core_cells <- c('A375', 'A549', 'HA1E', 'HCC515', 'HT29', 'HEPG2',
                'MCF7', 'PC3', 'VCAP')
for (ctype in core_cells){
  instances[[ctype]] <- instances[["all"]][instances[["all"]]$cell_iname==ctype,]
  num_drugs <- table(instances[[ctype]]$drug_name_small)
  keep_drugs <- names(which(num_drugs > 0))
  instances[[ctype]] <- instances[[ctype]][instances[[ctype]]$drug_name_small %in% keep_drugs,]
}

# # check
# for (ctype in names(instances)){
#   print(ctype)
#   print(length(unique(instances[[ctype]]$drug_name_small)))
# }

for (ctype in names(instances)){ # loop for type of core cell or "all"
  # save file if does not exist yet
  merged_ranks_file <- here(data_calc, 
              paste(ctype, "_reduced_aggregated_drug_ranks.rds", sep=""))
  # create log file
  logfile <- file.path(data_calc, paste(ctype, "_reduced.log", sep=""))
  if (!file.exists(merged_ranks_file)){
    drug_index <- 0
    res_drug <- list()
    for (drug in unique(instances[[ctype]]$drug_name_small)){ # loop for drugs
      drug_index <- drug_index + 1
      print(c(drug_index, length(unique(instances[[ctype]]$drug_name_small))))
      print(drug)
      di <- which(instances[[ctype]]$drug_name_small == drug)
      sig_vector <- unique(instances[[ctype]]$sig_id[di])
      di_colmeta <- which(col_meta[,1] %in% sig_vector)
      # remove missing values
      if (length(which(is.na(di_colmeta)))>0){
        di_colmeta <- di_colmeta[-which(is.na(di_colmeta))]
      }
      # if no result is available, move on to the next drug
      if (length(di_colmeta)==0){
        #print(paste(drug, " signatrue is missing", paste=""))
        cat(paste(drug, " signature is missing\n", sep=""), 
            file = logfile, append = TRUE)
        next
      }
      drug_matrix <- parse_gctx(datafile, cid = di_colmeta, rid = signature_rows)
      drug_annot <- col_meta[di_colmeta, 1]
      cat(paste(drug, ': \n', sep=""), file = logfile, append = TRUE)
      cat(paste(nrow(drug_matrix@mat), " ", ncol(drug_matrix@mat),
                '\n', sep=""), file = logfile, append = TRUE)
      
      # we want top overexpressed gene on the top (rank 1)
      rankmx <- apply(-drug_matrix@mat, 2, rank)
      rank_geo_mean <- apply(rankmx, 1, function(x) exp(mean(log(x))))
      rank_gene_id <- names(rank_geo_mean)[order(rank_geo_mean)]
      res_drug[[drug]] <- rank_gene_id
    } # loop for drugs
    # do not write result into file yet
    # saveRDS(res_drug, file=merged_ranks_file)
   
    # create a dataframe instead
    rank_drug <- list()
    gene_id <- rownames(rankmx)
    for (drug in names(res_drug)){
      rank_drug[[drug]] <- match(gene_id, res_drug[[drug]])
    }
    rank_drug_df <- do.call(cbind, rank_drug)
    # replace gene id with gene symbol
    gene_symbol <- geneinfo$gene_symbol[match(gene_id, geneinfo$gene_id)]
    rownames(rank_drug_df) <- gene_symbol
    
    # write result into file
    saveRDS(rank_drug_df, file=merged_ranks_file)
    
    # also write into a text file - TO DO
  }
} # loop for type of cell


# ctype <- "all"
# merged_ranks_file <- here(data_calc, 
#                      paste(ctype, "_reduced_aggregated_drug_ranks.rds", sep=""))
# rank_drug_df <- readRDS(merged_ranks_file)

# write out top and bottom genes from the antimalarials
am_smiles_file <- file.path(data_calc, "AMdrugs_reduced.smiles")
am_smiles <- read.table(am_smiles_file, header=F, sep="\t")
am_drugs <- am_smiles[,3]

# for each am_drug, select top and bottom 150 genes
num_genes <- 150
am_drug_top <- list()
am_drug_bottom <- list()
for (am_drug in am_drugs){
  gene_vector <- rank_drug_df[,  which(colnames(rank_drug_df)==am_drug)]
  gene_vector <- sort(gene_vector, decreasing = F)
  am_drug_top[[am_drug]] <- names(gene_vector[1:num_genes])
  am_drug_bottom[[am_drug]] <- names(rev(gene_vector)[1:num_genes])
}

am_drug_both <- list("upregulated"=as.data.frame(am_drug_top),
                     "downregulated"=as.data.frame(am_drug_bottom))
am_both_file <- here(resfolder, "AMdrugs_top_bottom_genes.xlsx")
openxlsx::write.xlsx(am_drug_both, file=am_both_file)

