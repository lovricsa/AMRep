# outdated in Bioconductor

##This function computes enrichment score for both input 'geneList' 
##and its permutations for one gene set.

##This function computes enrichment scores for GSEA, running score and 
##position of hits for a gene set.
gseaScores <- function(geneList, geneSet, exponent=1, mode="score") {
  # ##check arguments
  # paraCheck("genelist", geneList)
  # paraCheck("exponent", exponent)
  # paraCheck("gs", geneSet)
  # paraCheck("gseaScore.mode", mode)
  ##The geneSet should be a subset of the gene universe, i.e. we keep 
  ##only those element of the gene set that appear in the geneList		
  geneSet<-intersect(names(geneList), geneSet)
  ##Compute the size of the gene set and of the genelist	
  nh <- length(geneSet)
  N <- length(geneList)
  ##Initialize the ES, runningES and the Phit and Pmiss by position 
  ##(the actual values of Phit and Pmiss are actually cumulative sums 
  ##of these 'by position' values)	
  ES <- 0
  Phit <- rep(0, N)
  Pmiss <- rep(0, N)
  runningES <- rep(0, N)
  ##Stop if the geneSet is larger than the gene universe	
  if(nh > N) {
    stop("Gene Set is larger than Gene List")
  } else {
    ##Compute the positions of the hits in the geneList (0 if there 
    ##is no match, 1 if there is a match)	
    hits <- rep(FALSE, N)
    hits[which(!is.na(match(names(geneList), geneSet)))] <- TRUE
    ##If sum(hits)=0 then there is no match between geneList and 
    ##geneSet, and all scores stay at 0.		
    if(sum(hits)!=0) {
      ##Fill the Phit by position		
      Phit[which(hits)]<-abs(geneList[which(hits)])^exponent
      NR=sum(Phit)
      ##Fill the Pmiss by positions			
      Pmiss[which(!hits)]<-1/(N-nh)
      ##Do the cumulative sums	and compute the runningES		
      Phit=cumsum(Phit/NR)
      Pmiss=cumsum(Pmiss)
      runningES<-Phit-Pmiss
      ##Compute the maximal (positive) and minimal (or maximal 
      ##negative) values of the ES, and choose which one is kept			
      ESmax<-max(runningES)
      ESmin<-min(runningES)
      ES<-ifelse(abs(ESmin)>abs(ESmax), ESmin, ESmax)
    }
  }
  ##Return the relevant information according to mode		
  if(mode=="score")
    return(ES)	
  if(mode=="graph")	
    return(list("enrichmentScore"=ES, "runningScore"=runningES, 
                "positions"=as.integer(hits)))
}


# my own function to convert distance matrices
convert_dist2vector <- function(distancemx){
  distancemx[upper.tri(distancemx, diag = TRUE)] <- NA
  distvector <- as.vector(distancemx)
  distvector <- distvector[!is.na(distvector)]
  return(distvector)
}

# quick check whether a BRD-ID corresponds to a fully known structure

# convert dose to M (we know that dose can be uM only)
singleconvert <- function(dose, unit){
  if (unit == 'uM'){
    Mdose <- 10^(-6)*dose
  } else {
    print('check unit manually')
    Mdose <- NULL
  }
  return(Mdose)
}

convert2M <- function(inlist){
  dose <- inlist[['dose']]
  unit <- inlist[['unit']]
  output <- list()
  for (ii in 1:length(dose)){
    output[[ii]] <- singleconvert(dose[ii], unit[ii])
  }
  return(unlist(output))
}

# distES <- function(drug_matrix, num_genes){
#   es.dist <- matrix(data=0, nrow=dim(drug_matrix)[2], ncol=dim(drug_matrix)[2])
#   if (is.null(colnames(drug_matrix))) stop("the drug matrix needs column names")
#   # redifine list of drugs
#   for (test_drug in colnames(drug_matrix)){
#     ci <- which(colnames(drug_matrix) == test_drug) # column index
#     # select top and bottom genes - remember, it is rank ordered!
#     top_thr <- sort(drug_matrix[,test_drug], decreasing = F)[num_genes]
#     bottom_thr <- sort(drug_matrix[,test_drug], decreasing = T)[num_genes]
#     top_ind <- which(drug_matrix[,test_drug] <= top_thr)
#     bottom_ind <- which(drug_matrix[,test_drug] >= bottom_thr)
#     gs_up <- rownames(drug_matrix)[top_ind]
#     gs_down <- rownames(drug_matrix)[bottom_ind]
#     for (other_drug in colnames(drug_matrix)){
#       ri <- which(colnames(drug_matrix) == other_drug)
#       gl <- drug_matrix[,other_drug][order(drug_matrix[,other_drug], decreasing=F)]
#       names(gl) <- rownames(drug_matrix)[order(drug_matrix[,other_drug], decreasing=F)]
#       # calculate ES scores for the top (es_up) and bottom (es_down) genesets
#       es_up <- gseaScores(geneList=gl, geneSet=gs_up)
#       es_down <- gseaScores(geneList=gl, geneSet=gs_down)
#       # this would give same as GSEAweight1Score
#       # tes <- ifelse(es_up * es_down > 0, 0, es_up - es_down)
#       tes <- ifelse(es_up * es_down > 0, 0, (es_up - es_down)/2)
#       es.dist[ri, ci] <- 1-tes
#     }
#   } # end of loop for calculating ES-based distance 
#   es.dist <- (es.dist + t(es.dist))/2
#   colnames(es.dist) <- colnames(drug_matrix)
#   rownames(es.dist) <- colnames(drug_matrix)
#   return(es.dist)
# }


distES <- function(dmx1, dmx2=NULL, num_genes=150){
  if (is.null(dmx2)) dmx2 <- dmx1
  if (is.null(colnames(dmx1))) stop("the drug matrix needs column names")
  if (is.null(colnames(dmx2))) stop("the drug matrix needs column names")
  # redifine list of drugs
  # calculations needed twice reversing the matrices
  calc_matrices <- list('orig' = list('dmx1'=dmx1, 'dmx2'=dmx2),
               'reversed'= list('dmx1'=dmx2, 'dmx2'=dmx1))
  res_list <- list()
  for (morder in c('orig', 'reversed')){
    drug_matrix1 <- calc_matrices[[morder]][['dmx1']]
    drug_matrix2 <- calc_matrices[[morder]][['dmx2']]
    es.dist <- matrix(data=0, 
                      nrow=dim(drug_matrix1)[2], 
                      ncol=dim(drug_matrix2)[2])
    for (test_drug in colnames(drug_matrix1)){
      ri <- which(colnames(drug_matrix1) == test_drug) # row index
      # select top and bottom genes - remember, it is rank ordered!
      top_thr <- sort(drug_matrix1[,test_drug], decreasing = F)[num_genes]
      bottom_thr <- sort(drug_matrix1[,test_drug], decreasing = T)[num_genes]
      top_ind <- which(drug_matrix1[,test_drug] <= top_thr)
      bottom_ind <- which(drug_matrix1[,test_drug] >= bottom_thr)
      gs_up <- rownames(drug_matrix1)[top_ind]
      gs_down <- rownames(drug_matrix1)[bottom_ind]
      for (other_drug in colnames(drug_matrix2)){
        ci <- which(colnames(drug_matrix2) == other_drug)
        gl <- drug_matrix2[,other_drug][
          order(drug_matrix2[,other_drug], decreasing=F)]
        names(gl) <- rownames(drug_matrix2)[
          order(drug_matrix2[,other_drug], decreasing=F)]
        # calculate ES scores for the top (es_up) and bottom (es_down) genesets
        es_up <- gseaScores(geneList=gl, geneSet=gs_up)
        es_down <- gseaScores(geneList=gl, geneSet=gs_down)
        # this would give same as GSEAweight1Score
        # tes <- ifelse(es_up * es_down > 0, 0, es_up - es_down)
        tes <- ifelse(es_up * es_down > 0, 0, (es_up - es_down)/2)
        es.dist[ri, ci] <- 1-tes
        # # TEMPORARY
        # es.dist[ri, ci] <- tes
      }
    } # end of loop for calculating ES-based distance 
    res_list[[morder]] <- es.dist
  }
  es.dist <- (res_list[['orig']] + t(res_list[['reversed']]))/2
  rownames(es.dist) <- colnames(dmx1)
  colnames(es.dist) <- colnames(dmx2)
  return(es.dist)
}


# use  GSEAweight1Score instead
distES2 <- function(dmx1, dmx2=NULL, num_genes=150){
  if (is.null(dmx2)) dmx2 <- dmx1
  if (is.null(colnames(dmx1))) stop("the drug matrix needs column names")
  if (is.null(colnames(dmx2))) stop("the drug matrix needs column names")
  # redifine list of drugs
  # calculations needed twice reversing the matrices
  calc_matrices <- list('orig' = list('dmx1'=dmx1, 'dmx2'=dmx2),
                        'reversed'= list('dmx1'=dmx2, 'dmx2'=dmx1))
  res_list <- list()
  for (morder in c('orig', 'reversed')){
    drug_matrix1 <- calc_matrices[[morder]][['dmx1']]
    drug_matrix2 <- calc_matrices[[morder]][['dmx2']]
    es.dist <- matrix(data=0, 
                      nrow=dim(drug_matrix1)[2], 
                      ncol=dim(drug_matrix2)[2])
    for (test_drug in colnames(drug_matrix1)){
      ri <- which(colnames(drug_matrix1) == test_drug) # row index
      # select top and bottom genes - remember, it is rank ordered!
      top_thr <- sort(drug_matrix1[,test_drug], decreasing = F)[num_genes]
      bottom_thr <- sort(drug_matrix1[,test_drug], decreasing = T)[num_genes]
      top_ind <- which(drug_matrix1[,test_drug] <= top_thr)
      bottom_ind <- which(drug_matrix1[,test_drug] >= bottom_thr)
      gs_up <- rownames(drug_matrix1)[top_ind]
      gs_down <- rownames(drug_matrix1)[bottom_ind]
      # not use GSEAweight to calculate a row from the result matrix
      # uses up and down sets in the reverse order
      # cres <- GSEAweight1Score(drug_matrix2, gs_down, gs_up, permuteNum = 1,
      #                          pAdjMethod = "BH", mcCore = 1)
      cres <- ZhangScore(drug_matrix2, gs_down, gs_up, permuteNum = 1,
                               pAdjMethod = "BH", mcCore = 1)
      es.dist[ri,] <- 1 - cres[,1] / 2
      # # TEMPORARY
      # es.dist[ri, ] <- cres[,1]
    } # end of loop for calculating ES-based distance 
    res_list[[morder]] <- es.dist
  }
  es.dist <- (res_list[['orig']] + t(res_list[['reversed']]))/2
  rownames(es.dist) <- colnames(dmx1)
  colnames(es.dist) <- colnames(dmx2)
  return(es.dist)
}

simES <- function(dmx1, dmx2=NULL, num_genes=50,
                  simfunction=RCSM::ZhangScore){
  # input matrices: rank matrices!
  if (is.null(dmx2)) dmx2 <- dmx1
  if (is.null(colnames(dmx1))) stop("the drug matrix needs column names")
  if (is.null(colnames(dmx2))) stop("the drug matrix needs column names")
  # redifine list of drugs
  # calculations needed twice reversing the matrices
  calc_matrices <- list('orig' = list('dmx1'=dmx1, 'dmx2'=dmx2),
                        'reversed'= list('dmx1'=dmx2, 'dmx2'=dmx1))
  res_list <- list()
  for (morder in c('orig', 'reversed')){
    drug_matrix1 <- calc_matrices[[morder]][['dmx1']]
    drug_matrix2 <- calc_matrices[[morder]][['dmx2']]
    es.sim <- matrix(data=0, 
                      nrow=dim(drug_matrix1)[2], 
                      ncol=dim(drug_matrix2)[2])
    for (test_drug in colnames(drug_matrix1)){
      ri <- which(colnames(drug_matrix1) == test_drug) # row index
      # select top and bottom genes - remember, it is rank ordered!
      top_thr <- sort(drug_matrix1[,test_drug], decreasing = F)[num_genes]
      bottom_thr <- sort(drug_matrix1[,test_drug], decreasing = T)[num_genes]
      top_ind <- which(drug_matrix1[,test_drug] <= top_thr)
      bottom_ind <- which(drug_matrix1[,test_drug] >= bottom_thr)
      gs_up <- rownames(drug_matrix1)[top_ind]
      gs_down <- rownames(drug_matrix1)[bottom_ind]
      # not use GSEAweight to calculate a row from the result matrix
      # it uses a different method to define rank matrices
      # -> create mock log fold changes
      # uses up and down sets in the reverse order
      cres <- simfunction(rankmatrixTomockLFC(drug_matrix2), 
                          gs_up, gs_down, permuteNum = 1,
                         pAdjMethod = "BH", mcCore = 1)
      # first column is the score
      es.sim[ri,] <- cres[,1]
    } # end of loop for calculating ES-based similarity
    res_list[[morder]] <- es.sim
  }
  es.sim <- (res_list[['orig']] + t(res_list[['reversed']]))/2
  rownames(es.sim) <- colnames(dmx1)
  colnames(es.sim) <- colnames(dmx2)
  return(es.sim)
}


# temp
# retain_best_BRDIDs(pertinfo_phase2,
#                    common_name="cmap_name",
#                    important_columns=c("pert_id", "cmap_name"),
#                    dupl_folder=getwd(),
#                    additional_ordering=NULL,
#                    writeDuplicates=FALSE)
# end of temp

retain_best_BRDIDs <- function(allinfo,
                        common_name="pert_iname",
                        important_columns=c("pert_id", "pert_iname"),
                        dupl_folder=getwd(),
                        additional_ordering="pubchem_cid",
                        writeDuplicates=FALSE){
  # instead of the previous reduction method, select structures in the order
  # K (fully known), A (ambiguous), and U (unknown). M (mixture) should not
  # appear for the same compound, but if it does, flag it!
  allinfo$str1 <- substr(allinfo$pert_id, 5,5) # str1: K, A, U or M
  # decide which columns to use to detect duplicates
  if (important_columns[1] == "all") important_columns <- colnames(allinfo)
  allinfo$all <- apply(allinfo[, important_columns], 1,
                        paste , collapse = "_" )
  # delete "total duplicates"
  tdi <- which(duplicated(allinfo[, 'all']))
  if (length(tdi)>0)  allinfo <- allinfo[-tdi,]
  # now choose from the available common_name-s
  di <- duplicated(allinfo[, common_name])
  tmp <- tibble(str1 = c('K', 'A', 'U', 'M'))
  omit <- c()
  for (cname in unique(allinfo[di, common_name])){
    ci <- which(allinfo[, common_name] == cname)
    # this will keep the order given in tmp ('K', 'A', 'U', 'M')
    tmp_tbl <- dplyr::left_join(tmp, tibble::as_tibble(allinfo[ci,]), 
                                  by = c("str1" = "str1")) 
    # delete NA columns
    ni <- which(is.na(tmp_tbl[, common_name]))
    if (length(ni)>0) tmp_tbl <- tmp_tbl[-ni,]
    # if there is an additional ordering column, order accordingly
    if (!is.null(additional_ordering)){
      # put all rows not "-666" to the front
      ri <- which(tmp_tbl[,additional_ordering] != "-666")
      new_order <- c(ri, setdiff(c(1:nrow(tmp_tbl)), ri))
      tmp_tbl <- tmp_tbl[new_order,]
    }
    # also flag all but the first to omit
    ci_omit <- which(allinfo[, "all"] %in% unlist(tmp_tbl[2:dim(tmp_tbl)[1], "all"]))
    omit <- c(omit, ci_omit)
    # save structures for later 2D comparison calculations
    # reduce to canonical_smiles
    si <- which(duplicated(tmp_tbl$canonical_smiles))
    if (length(si) > 0) tmp_tbl <- tmp_tbl[-si,]
    if (dim(tmp_tbl)[1]>1 & writeDuplicates){
      dupl_file <- file.path(dupl_folder, paste(cname, '.smiles', sep=""))
      dupl_smiles <- tmp_tbl[,c('canonical_smiles', 'pert_id', common_name)]
      write.table(dupl_smiles, file=dupl_file, row.names = F, quote=F, 
                  col.names = F, sep="\t")
    }
  }
  # return the final dataframe
  if (length(omit) == 0){
    return(allinfo)
  } else {
    return(allinfo[-omit,])
  }
}

## util functions for rocs results
obtain_target <- function(extended_target){
  tgl <- strsplit(extended_target, split="_")[[1]]
  output <- paste(tgl[1:(length(tgl)-2)], collapse="_")
  return(output)
}

drug2lower <- function(drug){
  # converts the name lowercase and also removes special characters
  lower <- tolower(gsub("[^[:alnum:] ]", "", drug))
  lower <- gsub(" ", "", lower, fixed=T)
  return(lower)
}

# to convert 'Anna-type' rank matrix to 'LFC-type matrix
rankmatrixTomockLFC <- function(rankMatrix) {
  refLFC <- matrix(data=NA, nrow=nrow(rankMatrix), ncol=ncol(rankMatrix),
                   dimnames = dimnames(rankMatrix))
  for (i in 1:ncol(rankMatrix)) {
    refLFC[,i] <- rankvectorTomockLFC(rankMatrix[,i])
  }
  return(refLFC)
}

# rankvectorTomockLFC <- function(rankVector, topperc = 0.33,
#                                   bottomperc = 0.33){
#   mockLFC <- rep(NA, times=length(rankVector))
#   for (ii in 1:length(rankVector)){
#     if (rankVector[ii] <= topperc * length(rankVector)){
#       mockLFC[ii] <- length(rankVector)/2 - rankVector[ii] + 1
#     } else {
#       if (rankVector[ii] > (1-bottomperc)* length(rankVector)){
#         mockLFC[ii] <- -length(rankVector)/2 + rankVector[ii] - 1
#       } else {
#         mockLFC[ii] <- 0
#       }
#     }
#   }
#   return(mockLFC)
# }


rankvectorTomockLFC <- function(rankVector, topperc = 0.33,
                                bottomperc = 0.33){
  mockLFC <- rep(NA, times=length(rankVector))
  for (ii in 1:length(rankVector)){
    if (rankVector[ii] <= topperc * (length(rankVector)+1) |
          rankVector[ii] >=  (1-bottomperc) * (length(rankVector)+1) ){
      mockLFC[ii] <- (length(rankVector)+1)/2 - rankVector[ii]
    } else {
      mockLFC[ii] <- 0
    }
  }
  return(mockLFC)
}

distmatrix2vector <- function(distmatrix){
  keep <- lower.tri(distmatrix, diag = FALSE)
  return(unlist(distmatrix[keep]))
}

#https://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
readROCSoutput <- function(fpath) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one 
      # R expression in the "try" part then you'll have to 
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      # message("This is the 'try' part")
      header <- scan(fpath, what=character(),
                     sep="\t", nlines=1, quiet = TRUE)
      # ryCatch returns the value associated to executing expr unless there's an error or a warning
      tmp <- read.table(file.path(opt$infolder, frfile), sep="\t", 
                        row.names=NULL, comment.char="",
                        header=F, skip=1)
      indata <- as_tibble(tmp[,c(1,2,4)])
      colnames(indata) <- header[c(1,2,4)]
      return(indata)
      # You don't need to state the return value via `return()` as code 
      # in the "try" part is not wrapped inside a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message("No lines in file")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NULL)
    },
    warning=function(cond) {
      message("reading the file caused a warning")
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      # message(paste("Processed URL:", url))
      # message("Some other message at the end")
    }
  )    
  return(out)
}

# specifically for ROCS runs
#########################################
obtain_run_num <- function(filename, fsplit="_", fpattern="run"){
  file_prefix <- strsplit(filename, split=fsplit, fixed=T)[[1]][1] 
  run_num <- gsub(fpattern, "", file_prefix, fixed=T)
  return(run_num)
}
remove_wart_num <- function(inname){
  inname_short <- strsplit(inname, split="_", fixed=T)[[1]][1] 
  return(inname_short)
}

#########################################
select_intersection <- function(intible, inmolecules, colname){
  outtible <- intible[which(intible[,colname,drop=T] %in% inmolecules),]
  return(outtible)
}

# #########################################
# # from https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
# library(gtable)
# library(cowplot)
# 
# shift_legend <- function(p){
#   
#   # check if p is a valid object
#   if(!"gtable" %in% class(p)){
#     if("ggplot" %in% class(p)){
#       gp <- ggplotGrob(p) # convert to grob
#     } else {
#       message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
#       return(p)
#     }
#   } else {
#     gp <- p
#   }
#   
#   # check for unfilled facet panels
#   facet.panels <- grep("^panel", gp[["layout"]][["name"]])
#   empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
#   empty.facet.panels <- facet.panels[empty.facet.panels]
#   if(length(empty.facet.panels) == 0){
#     message("There are no unfilled facet panels to shift legend into. Returning original plot.")
#     return(p)
#   }
#   
#   # establish extent of unfilled facet panels (including any axis cells in between)
#   empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
#   empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
#                              max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
#   names(empty.facet.panels) <- c("t", "l", "b", "r")
#   
#   # extract legend & copy over to location of unfilled facet panels
#   guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
#   if(length(guide.grob) == 0){
#     message("There is no legend present. Returning original plot.")
#     return(p)
#   }
#   gp <- gtable_add_grob(x = gp,
#                         grobs = gp[["grobs"]][[guide.grob]],
#                         t = empty.facet.panels[["t"]],
#                         l = empty.facet.panels[["l"]],
#                         b = empty.facet.panels[["b"]],
#                         r = empty.facet.panels[["r"]],
#                         name = "new-guide-box")
#   
#   # squash the original guide box's row / column (whichever applicable)
#   # & empty its cell
#   guide.grob <- gp[["layout"]][guide.grob, ]
#   if(guide.grob[["l"]] == guide.grob[["r"]]){
#     gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
#   }
#   if(guide.grob[["t"]] == guide.grob[["b"]]){
#     gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
#   }
#   gp <- gtable_remove_grobs(gp, "guide-box")
#   
#   return(gp)
# }

obtain_heat_tibble <- function(simil_tibble,
                              am_drugs,
                              hit_drugs,
                              simil_name = "CFP"){
  # Ensure both directions are covered (drug1 vs drug2)
  cfp_trans <-  mutate(simil_tibble, drug_tmp = drug1, drug1 = drug2,  
                       drug2 = drug_tmp) %>% 
    select(-drug_tmp)
  cfp_long <- bind_rows(simil_tibble, cfp_trans)
  
  # Filter only rows where the "row drug" is in am_drugs 
  # and the "column drug" is in hit_drugs
  cfp_filtered <- cfp_long %>%
    filter(drug1 %in% am_drugs,
           drug2 %in% hit_drugs)
  
  
  # Pivot to wide format
  cfp_wide <- cfp_filtered %>%
    pivot_wider(
      names_from = drug2,
      values_from = !!sym(simil_name),
      names_prefix = ""
    ) %>%
    rename(am_drug = drug1)
  
  # report missing drugs
  missing_hit <- hit_drugs[(!hit_drugs %in% colnames(cfp_wide))]
  missing_am <- am_drugs[(!am_drugs %in% cfp_wide$am_drug)]
  
  return(list("heat_tibble"=cfp_wide,
              "missing_hit"=missing_hit,
               "missing_am"=missing_am))
}

