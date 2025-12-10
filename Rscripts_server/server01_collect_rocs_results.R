# 2022.06.18. Anna Lovrics
# collect results from rocs, put into rds file
# this only needs to be run once
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
#######################################

# read in results
# all .rpt files are needed
fl <- list.files(path =file.path(data_calc, "rocs_result"),
                    pattern = "_1.rpt$")
queries <- sapply(fl, function(x) gsub("_1.rpt", "", x, fixed=T))
query_res <- list() # list of one column matrices
query_num <- 0
query_problem <- c()
for (query in queries){
  query_num <- query_num + 1
  print(paste("query_num: ", as.character(query_num), sep=""))
  # read in all report files corresponding to this query
  ql <-  list.files(path =file.path(data_calc, "rocs_result"),
                    pattern = glob2rx(paste(query, "*.rpt", sep="")))
  qfile1 <- ql[1]
  resdf <- read.csv(file=file.path(data_calc, "rocs_result", qfile1), 
                    sep="\t", stringsAsFactors=FALSE, na.strings='NULL')
  # last column is an artefact
  resdf$X <- NULL
  if (nrow(resdf) == 0 ){
    print( paste('problem with ', query, sep=""))
    query_problem <- c(query_problem, query)
    next
  }
  if (length(ql) > 1){
    for (qfile in ql[2:length(ql)]){
      resdf_ext <- read.csv(file=file.path(data_calc, "rocs_result", qfile), 
                        sep="\t", stringsAsFactors=FALSE, na.strings='NULL')
      resdf_ext$X <- NULL
      resdf <- rbind(resdf, resdf_ext)
    }
  } # end of reading in rocs result
  # now find best value for the final similarity matrix
  # obtain name of the query (once)
  query_name <- obtain_target(resdf$ShapeQuery[1])
  # obtain name of the targets
  resdf$target <- sapply(resdf$Name, obtain_target)
  target_res <- list() # list of size one elements (the similarity score)
  for (target in unique(resdf$target)){
    ti <- which(resdf$target == target)
    target_res[[target]] <- max(resdf$TanimotoCombo[ti])
  } # end of loop finding best target
  # convert results matrix of one column
  query_res_tmp <- as.data.frame(t(as_tibble_row(target_res)))
  query_res_tbl <- rownames_to_column(query_res_tmp, 
                                      var = "target") %>% as_tibble()
  colnames(query_res_tbl)[2] <- query_name
  query_res[[query_name]] <- query_res_tbl
} # end of cycling through each query
# now put together results in a tibble / dataframe and save for further use
final_res <- purrr::reduce(query_res, dplyr::full_join, by = "target")

colnames(final_res) <- drug2lower(colnames(final_res))
final_res$target <- drug2lower(final_res$target)
length(setdiff(drugs, final_res$target))

# now select the drugs that have been included in the final selection round


final_file <- file.path(data_calc, "rocs_simil.rds")
saveRDS(final_res, file=final_file)
# check for missing values
# nai <- which(is.na(final_res), arr.ind = TRUE)
