library('tidyverse')
library('fgsea')
library(biomaRt)


#' Function to run fgsea on results
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
run_common_gsea <- function(labeled_results, min_size, max_size) {
  
  labeled_results <- labeled_results %>% separate(First.Protein.Name, sep='\\.', into='First.Protein.Name', remove=TRUE)
  gene_ids <- labeled_results %>% pull(First.Protein.Name)
  
  
  hgnc_symbols <- read_csv("data/biomart_ensmusg_to_hgnc.csv", show_col_types = FALSE)
  
  hgnc_results <- labeled_results %>% left_join(hgnc_symbols, by=c('First.Protein.Name' = 'HGNC.symbol'))
  
  rnks <- hgnc_results %>% 
    drop_na(First.Protein.Name, fc_betas) %>% 
    distinct(First.Protein.Name, fc_betas, .keep_all=TRUE) %>%
    arrange(desc(fc_betas)) %>% 
    dplyr::select(First.Protein.Name, fc_betas) %>% 
    deframe()
  
  c2_pathways <- gmtPathways("data/c2.cp.v7.5.1.symbols.gmt")
  
  fgsea_results <- fgsea(c2_pathways, rnks, 
                         minSize = min_size, maxSize = max_size) %>% as_tibble()
  
  return(fgsea_results)
}


########## FUNCTION CALLS ##############

common_proteins <- read.csv("data/common_proteins.csv")

mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
  distinct() %>%
  as_tibble() %>%
  na_if("") %>% 
  na.omit()

res <- inner_join(common_proteins, bm, by=c("First.Protein.Name"="hsapiens_homolog_associated_gene_name"))

common_proteins <- res %>% 
  dplyr::select(First.Protein.Name, fc_betas) %>% 
  na.omit()

common_fgsea_results <- run_common_gsea(common_proteins, 0, 500)
common_fgsea_results
common_fgsea_results <- write_csv(common_fgsea_results, "data/common_proteins_fgsea.csv")

####################

fgseaResTidy <- common_fgsea_results %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

### Attempted Plotting ###
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA") + 
  theme_minimal()

# common_fgsea_results %>%
#   mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
#   ggplot() +
#   geom_bar(aes(x=pathway, y=NES, fill = padj < .95), stat='identity') +
#   scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
#   theme_minimal() +
#   ggtitle('fgsea results for protein sets') +
#   ylab('Normalized Enrichment Score (NES)') +
#   xlab('') +
#   coord_flip()
