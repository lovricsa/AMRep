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
library(fgsea)
library(msigdbr) # for the gene sets
#############################################

# folders
data_phase <-  file.path("Data", "CLUE")
data_repurpose <-  file.path("Data", "RepurposingHub")
data_calc <- file.path("Data","CLUE_calculated")
figfolder <- "Figures"
resfolder <- "Results"
folderFA <- file.path(resfolder, "Enrichment")
if (!file.exists(folderFA)) dir.create(folderFA)

# utility functions
source(here('Rscripts', 'utils.R'))
############################

# read in genes universe
geneinfo_file <- here(data_phase, "geneinfo_beta.txt") 
geneinfo <- read.table(geneinfo_file, sep="\t", header=T, quote = "") # 12328*7
# note - some ensembl id is missing which(geneinfo$ensembl_id=="\"\"")

#find the indices for the genes used
bing_indices <- which(geneinfo$feature_space == "best inferred")
bing_genes_ID <- geneinfo$ensembl_id[bing_indices]

am_both_file <- here(resfolder, "AMdrugs_top_bottom_genes.xlsx")
am_drug_both <- list()
am_drug_both[["upregulated"]] <- openxlsx::read.xlsx(am_both_file, sheet="upregulated")
am_drug_both[["downregulated"]] <- openxlsx::read.xlsx(am_both_file, sheet="downregulated")

#GSEA analysis for the contrasts with selected gene sets
# selected gene sets (GO:BP, CP:KEGG, hallmark)
category_names <- c("C5", "C2", "H")
subcategory_names <- c("GO:BP","CP:KEGG", "")
prefix <- c("gobp_", "kegg_", "hallmark_")
# repeat for different coefficients
FA_lists <- list()

am_drug <- "artesunate"
for (gsi in c(1:length(category_names))){
  broad_name <- subcategory_names[gsi]
  folder_name <- ifelse(nchar(broad_name) > 0, broad_name, category_names[gsi])
  # create a folder
  folderenrres <- here(folderFA, folder_name)
  if (!file.exists(folderenrres)) dir.create(folderenrres)
  # use the actual genesets
  GO_db <- msigdbr(species= "Homo sapiens",
                   category=category_names[gsi], subcategory = subcategory_names[gsi])
  GO_db_list <- unique(GO_db[, c("ensembl_gene","gs_name")])
  GO_db_list <- split(x=GO_db_list$ensembl_gene, f=GO_db_list$gs_name)
  print(paste("GO_db_list ", folder_name, " created", sep = ""))
  
  # replace gene symbol with ensembl ID
  upreg_ID <- geneinfo$ensembl_id[match(
    am_drug_both[["upregulated"]][[am_drug]], geneinfo$gene_symbol
             )] 
  upreg_rank <- c(1:length(upreg_ID))
  names(upreg_rank) <- upreg_ID
  
  # remove genes not in universe 
  for (gsname in names(GO_db_list)){
    # keepgene <- which(GO_db_list[[gsname]] %in% names(fite$Amean))
    keepgene <- which(GO_db_list[[gsname]] %in% names(upreg_rank))
    GO_db_list[[gsname]] <- GO_db_list[[gsname]][keepgene]
  }
  # now the gene enrichment analysis
  # genes_universe <- names(fite$Amean)
  genes_universe <- bing_genes_symbol
  set.seed(99)
  system.time( # keep track of elapsed time
    fgsea_res <- fgseaMultilevel(pathways = GO_db_list,
                                 stats = upreg_rank, 
                                 minSize = 0, 
                                 maxSize = 500, 
                                 eps = 0, 
                                 nPermSimple = 50000)
  )
  fgsea_res <- arrange(fgsea_res, -NES)
  # remove na rows
  omit <- which(is.na(fgsea_res$pval))
  if (length(omit) > 0) fgsea_res <- fgsea_res[-omit,]
  # write out leading edge genes' symbol as well
  fgsea_res$leadingSymbols <- fgsea_res$leadingEdge
  for (ri in 1:nrow(fgsea_res)){
    gsymbols <- geneinfo$gene_symbol[match(
      fgsea_res$leadingEdge[[ri]], geneinfo$ensembl_id)]
    fgsea_res$leadingEdge[[ri]] <- gsymbols
  }
    
    # add number of leading edge genes
    fgsea_res$count = lengths(fgsea_res$leadingEdge)
    
    # add gene ratio & background ratio 
    Nval <- nrow(fite$coefficients)
    nval <- length(which(toptt_orig[[mcf]]$adj.P.Val<=0.05 & 
                           toptt_orig[[mcf]]$logFC>=1))
    fgsea_res$GeneRatio <- fgsea_res$count/nval
    fgsea_res$BgRatio <- fgsea_res$size/Nval
    
    # write to file (TO DO: for each antimalarial drug, write to a seperate tab)
    outfile <- file.path(folderenrres, 
                         paste("enriched_genesets_fgsea_", am_drug, ".xlsx", sep=""))
    openxlsx::write.xlsx(fgsea_res, file = outfile, colNames=T, rowNames=F)
  }




# put into the list
FA_lists[[mcf]][[folder_name]] <- fgsea_res


# Visual representation of gene set enrichment results of "H", GO:BP and KEGG pathways. Included graphics for Figure 4 L, M (dotplots for contrasts TIS vs CTR and REPOP vs TIS in malignant cells for hallmark gene sets), Supplementary Figure S3E (dotplot for contrast TIS vs CTR in HFF for hallmark gene sets), Figure S3F (dotplot for contrast TIS vs CTR in malignant cells for KEGG gene sets), Figure S3G (dotplot for contrast TIS vs CTR in HFF for KEGG gene sets)
for (mcf in mcf_vector){
  geneset_file <- file.path(folderFA, "GO_category_enrichments.rds")
  for (gsi in c(1:length(category_names))){
    y <- FA_lists[[mcf]][[gsi]]
    
    # select activated and supressed genesets
    y$direction <- "activated"
    y$direction[y$NES<=0] <- "suppressed"
    numact <- length(which(y$direction=="activated"))
    numsupp <- length(which(y$direction=="suppressed"))
    
    # pick top and bottom genesets
    act_Nvals <- y$NES[y$direction=="activated"]
    act_thr <- sort(act_Nvals, decreasing = T)[min(10, numact)]
    act_keep <- which(y$direction=="activated" & y$NES >= act_thr)
    sup_Nvals <- y$NES[y$direction=="suppressed"]
    sup_thr <- sort(sup_Nvals, decreasing = F)[min(10, numsupp)]
    sup_keep <- which(y$direction=="suppressed" & y$NES <= sup_thr)
    y <- y[c(act_keep, sup_keep),]
    
    # sort by normalized enrichment
    y <- y[order(y$NES, decreasing = F),]
    # change the names of the gene sets
    y$pathway <- tolower(y$pathway)
    y$pathway <- gsub(prefix[gsi], "", y$pathway)
    y$pathway <- gsub("_", " ", y$pathway)
    y$pathway <- str_to_sentence(y$pathway)
    y$pathway <- gsub("dna", "DNA", y$pathway)
    y$pathway <- gsub("Dna", "DNA", y$pathway)
    #y$pathway <- nice_names$name[match(y$pathway, nice_names$compare)]
    
    y$pathway <- factor(y$pathway, levels=y$pathway)
    
    # plot together
    p <- ggplot(y, 
                aes(x = NES, y = pathway)) + 
      geom_point(aes(size = count, color = padj)) +
      guides(size=guide_legend(title="Count", ncol=1,byrow=TRUE,
                               title.vjust=0.9, order=1))+
      theme_bw(base_size = 12) +
      scale_color_continuous(low="red", high="blue", name = "p.adjust",
                             guide=guide_colorbar(reverse=TRUE))+
      theme(legend.position="right", 
            legend.box = "vertical",
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            #strip.background =element_rect(fill=cfill),
            strip.text = element_text(colour = 'black'))+
      ylab(NULL) + xlab("NES")+ 
      scale_size(range=c(3, 8))+
      #scale_y_discrete(labels = function(x) str_wrap(x, width = 24))+
      geom_vline(xintercept = 0, linetype="dashed")
    figfile <- file.path(figfolder, paste(mcf, '_',
                                          names(FA_lists[[mcf]])[gsi], '_', postfix, '_', 'both',
                                          "_enrich_dotplot.pdf", sep=""))
    ggsave(figfile, plot=p, width=7, height=6)
    
  }
}


