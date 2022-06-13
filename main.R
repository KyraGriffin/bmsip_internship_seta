
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")
# BiocManager::install("fgsea")

library('tidyverse')
library(readxl)
# library('SummarizedExperiment')
# library('DESeq2')
# library('biomaRt')
# library('testthat')
# library('fgsea')

#' Function to read in the proteomics data
#'
#' @param p_path (str): path to the file proteomics_ack.xlsx
#' 
#' @return 
#' @export
#'
read_data <- function(p_path) {

  proteomics_ack <- read_excel(p_path, 
                               sheet = "Summary", skip = 10)
  return(proteomics_ack)
}


select_columns <- function(proteomics_ack) {
  
  proteomics_ack <- proteomics_ack %>% 
    dplyr::select(c("Index","Index in Detail","Sample 2 : Sample 1...3", 
                "Sample 3 : Sample 1...4", "Sample 3 : Sample 2...5", 
                "Max Abundance", "Max % CV", "Gene Name", "Protein Name",
                "Site", "Description", "Accession", "kD", "Peptide",
                "Charge","Calc. m/z", "Count In Details", "Average RT", "Species"))

  
  proteomics_ack <- proteomics_ack %>% 
    dplyr::rename("Sample2vSample1" = "Sample 2 : Sample 1...3", 
                  "Sample3vSample1" = "Sample 3 : Sample 1...4",
                  "Sample3vSample2" = "Sample 3 : Sample 2...5")
  
  return(proteomics_ack)
}

id_protein_groups <- function(proteomics_ack){
  
  prot_list <- proteomics_ack[is.na(proteomics_ack$Sample2vSample1),]
  
  protein_groups <- prot_list %>% dplyr::select("Index in Detail") %>% 
    dplyr::rename("ProteinType" = "Index in Detail")
  
  return(protein_groups)
}


consolidate_data <- function(proteomics_ack, protein_groups) {
  
  proteomics_ack <- proteomics_ack %>% 
    drop_na("Sample2vSample1") %>% 
    dplyr::select(-c("Index", "Index in Detail"))
  
  return(proteomics_ack)
}




