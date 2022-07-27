library('tidyverse')
library('fgsea')

common_proteins <- read.csv("data/common_proteins.csv")
common_proteins <- common_proteins %>% 
  as_tibble() %>% 
  pivot_wider(names_from = X, values_from = counts)
              
head(common_proteins)
# rnk_list <- setNames(counts, common_proteins)
# head(rnk_list)

pathways_fgsea <- fgsea::gmtPathways('data/c2.cp.v7.5.1.symbols.gmt')

fgsea_results <- fgsea(pathways_fgsea, 
                       common_proteins, 
                       minSize = 15, 
                       maxSize=500)

fgsea_results <- fgsea_results %>% as_tibble()

#' Function to run fgsea on DESeq2 results
#'
#' @param labeled_results (tibble): the labeled results from DESeq2
#' @param gmt (str): the path to the GMT file
#' @param min_size: the threshold for minimum size of the gene set
#' @param max_size: the threshold for maximum size of the gene set
#'
#' @return tibble containing the results from running fgsea using descending
#' log2foldchange as a ranking metric
#' @export
#'
#' @examples fgsea_results <- run_gsea(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500)
run_gsea <- function(labeled_results, min_size, max_size) {
  
  labeled_results <- labeled_results %>% separate(X, sep='\\.', into='X', remove=TRUE)
  gene_ids <- labeled_results %>% pull(X)
  
  
  hgnc_symbols <- read_csv("data/biomart_ensmusg_to_hgnc.csv", show_col_types = FALSE)
  
  hgnc_results <- labeled_results %>% left_join(hgnc_symbols, by=c('X' = 'HGNC.symbol'))
  
  rnks <- hgnc_results %>% 
    drop_na(X, counts) %>% 
    distinct(X, counts, .keep_all=TRUE) %>%
    arrange(desc(counts)) %>% 
    dplyr::select(X, counts) %>% 
    deframe()
  
  c2_pathways <- gmtPathways("data/c2.cp.v7.5.1.symbols.gmt")
  
  fgsea_results <- fgsea(c2_pathways, rnks, 
                         minSize = min_size, maxSize = max_size) %>% as_tibble()
  
  return(fgsea_results)
}


########## FUNCTION CALLS ##############

common_proteins <- read.csv("data/common_proteins.csv")


fgsea_results <- run_gsea(common_proteins, 0, 500)
fgsea_results
fgsea_results <- write_csv(fgsea_results, "data/fgsea.csv")

####################
