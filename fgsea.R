######
# Author: Kyra Griffin
# Project: BMSIP Internship 2022
# PI: Dr. Francesca Seta
# Mentor: Joey Orofino
# Internship Director: Adam Labadorf
######


library('tidyverse')
library('fgsea')
library(biomaRt)
library(data.table)


#' Function to run fgsea on results
#'
#' @param labeled_results (tibble): the results
#' @param min_size: the threshold for minimum size of the gene set
#' @param max_size: the threshold for maximum size of the gene set
#'
#' @return tibble containing the results from running fgsea using descending
#' beta values as a ranking metric
#' @export
#'
#' @examples fgsea_results <- run_gsea(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500)
run_gsea <- function(labeled_results, min_size, max_size) {
  
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
                         minSize = min_size, maxSize = max_size) %>% as_tibble() # , nperm=10000
  
  return(fgsea_results)
}


########## FUNCTION CALLS ##############

common_proteins <- read.csv("data/common_proteins.csv")
ack_proteins <- read.csv("data/ack_proteins.csv")
total_proteins <- read.csv("data/total_proteins.csv")

mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
  distinct() %>%
  as_tibble() %>%
  na_if("") %>% 
  na.omit()

common <- inner_join(common_proteins, bm, by=c("First.Protein.Name"="hsapiens_homolog_associated_gene_name"))
ack <- inner_join(ack_proteins, bm, by=c("First.Protein.Name"="hsapiens_homolog_associated_gene_name"))
total <- inner_join(total_proteins, bm, by=c("First.Protein.Name"="hsapiens_homolog_associated_gene_name"))

common_proteins <- common %>% 
  dplyr::select(First.Protein.Name, fc_betas) %>% 
  na.omit() %>% 
  group_by(First.Protein.Name) %>% 
  filter(abs(fc_betas) == max(abs(fc_betas))) %>% 
  distinct()

ack_proteins <- ack %>% 
  dplyr::select(First.Protein.Name, fc_betas) %>% 
  na.omit() %>% 
  group_by(First.Protein.Name) %>% 
  filter(abs(fc_betas) == max(abs(fc_betas))) %>% 
  distinct()

total_proteins <- total %>% 
  dplyr::select(First.Protein.Name, fc_betas) %>% 
  na.omit() %>% 
  group_by(First.Protein.Name) %>% 
  filter(abs(fc_betas) == max(abs(fc_betas))) %>% 
  distinct()


common_fgsea_results <- run_gsea(common_proteins, 0, 500)
common_fgsea_results
common_fgsea_results <- write_csv(common_fgsea_results, "data/common_proteins_fgsea.csv")
fwrite(common_fgsea_results, file="data/common_fgsea_results.csv", sep=",", sep2=c("", " ", ""))

ack_fgsea_results <- run_gsea(ack_proteins, 0, 500)
ack_fgsea_results
ack_fgsea_results <- write_csv(ack_fgsea_results, "data/ack_proteins_fgsea.csv")
fwrite(ack_fgsea_results, file="data/ack_fgsea_results.csv", sep=",", sep2=c("", " ", ""))

total_fgsea_results <- run_gsea(total_proteins, 0, 500)
total_fgsea_results
total_fgsea_results <- write_csv(total_fgsea_results, "data/total_proteins_fgsea.csv")
fwrite(total_fgsea_results, file="data/total_fgsea_results.csv", sep=",", sep2=c("", " ", ""))
###################

common<- common_fgsea_results %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
common <- common %>% 
  dplyr::select(-ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

c <- common %>% arrange(pval)%>% 
  slice(1:25) 

fwrite(c, file="data/common_fgsea_results_25.csv", sep=",", sep2=c("", " ", ""))

ack <- ack_fgsea_results %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
ack %>% 
  dplyr::select(-ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

a <- ack %>% arrange(pval)%>% 
  slice(1:25) 

fwrite(a, file="data/ack_fgsea_results_25.csv", sep=",", sep2=c("", " ", ""))


total <- total_fgsea_results %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
total %>% 
  dplyr::select(-ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

t <- total %>% arrange(pval)%>% 
  slice(1:25) 

fwrite(t, file="data/total_fgsea_results_25.csv", sep=",", sep2=c("", " ", ""))








### Plotting ###

c %>%
  mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot() +
  geom_bar(aes(x=pathway, y=NES), fill="blue", stat='identity') +
  # scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
  theme_minimal() +
  ggtitle('Top 25 Pathways from GSEA for Common Proteins') +
  ylab('Normalized Enrichment Score (NES)') +
  xlab('') +
  coord_flip()

a %>%
  mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot() +
  geom_bar(aes(x=pathway, y=NES), fill = 'blue', stat='identity') +
  # scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
  theme_minimal() +
  ggtitle('Top 25 Pathways from GSEA for Ack Proteins') +
  ylab('Normalized Enrichment Score (NES)') +
  xlab('') +
  coord_flip()

t %>%
  mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot() +
  geom_bar(aes(x=pathway, y=NES), fill = 'blue', stat='identity') +
  # scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
  theme_minimal() +
  ggtitle('Top 25 Pathways from GSEA for Total Proteins') +
  ylab('Normalized Enrichment Score (NES)') +
  xlab('') +
  coord_flip()

