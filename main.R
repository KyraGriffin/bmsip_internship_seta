
library('tidyverse')
library(readxl)
library(rstatix)
library(ggplot2)
library(stringr)
install.packages('GGally')
library('GGally')

#' Function to read in the proteomics data
#'
#' @param p_path (str): path to the file proteomics_ack.xlsx
#' 
#' @return 
#'
read_data <- function(p_path) {

  proteomics_ack <- read_excel(p_path, 
                               sheet = "Summary", skip = 10)
  return(proteomics_ack)
}

#' Function to read in the proteomics data
#'
#' @param p_path (str): path to the file proteomics_ack.xlsx
#' 
#' @return 
#'
select_columns <- function(proteomics_ack) {
  
  # proteomics_ack <- proteomics_ack %>% 
  #   dplyr::select(c("Index","Index in Detail","Sample 2 : Sample 1...3", 
  #               "Sample 3 : Sample 1...4", "Sample 3 : Sample 2...5", 
  #               "Max Abundance", "Max % CV", "Gene Name", "Protein Name",
  #               "Site", "Description", "Accession", "kD", "Peptide",
  #               "Charge","Calc. m/z", "Count In Details", "Average RT", "Species"))

  
  proteomics_ack <- proteomics_ack %>% 
    dplyr::rename("Sample2vSample1" = "Sample 2 : Sample 1...3", 
                  "Sample3vSample1" = "Sample 3 : Sample 1...4",
                  "Sample3vSample2" = "Sample 3 : Sample 2...5")
  
  return(proteomics_ack)
}

#' Function to read in the proteomics data
#'
#' @param p_path (str): path to the file proteomics_ack.xlsx
#' 
#' @return 
#'
id_protein_groups <- function(proteomics_ack){
  
  prot_list <- proteomics_ack[is.na(proteomics_ack$Sample2vSample1),]
  
  protein_groups <- prot_list %>% dplyr::select("Index in Detail") %>% 
    dplyr::rename("ProteinType" = "Index in Detail")
  
  return(protein_groups)
}

#' Function to read in the proteomics data
#'
#' @param p_path (str): path to the file proteomics_ack.xlsx
#' 
#' @return 
#'
consolidate_data <- function(proteomics_ack, protein_groups) {
  
  proteomics_ack <- proteomics_ack %>% 
    drop_na("Sample2vSample1") 
  #%>% 
    # dplyr::select(-c("Index", "Index in Detail"))
  
  return(proteomics_ack)
}

#' Function to read in the proteomics data
#'
#' @param p_path (str): path to the file proteomics_ack.xlsx
#' 
#' @return 
#'
form_t_test_data <- function(proteomics_ack) {
  
  proteomics_ack$Sample3vSample1 <- str_replace(proteomics_ack$Sample3vSample1,"-", "0") 
  proteomics_ack$Sample3vSample1 <- as.double(proteomics_ack$Sample3vSample1)
  
  #replace NA values with zero in rebs column only
  proteomics_ack <- proteomics_ack %>% mutate(Sample3vSample1 = ifelse(is.na(Sample3vSample1), 0, Sample3vSample1))
  
  t_test_data <- proteomics_ack %>% 
    dplyr::select(c("Sample2vSample1", "Sample3vSample1", "Sample3vSample2")) %>% 
    pivot_longer(everything(),
                 names_to = "sampleVsample",
                 values_to = "normFoldChange")
  
  return(t_test_data)
}


#' Function to read in the proteomics data
#'
#' @param p_path (str): path to the file proteomics_ack.xlsx
#' 
#' @return 
#'
summary_t_test_data <- function(t_test_data){
  
  summary <- group_by(t_test_data, sampleVsample) %>%
    summarise(
    count = n(),
    mean = mean(normFoldChange, na.rm = TRUE),
    sd = sd(normFoldChange, na.rm = TRUE)
  )
  
  return(summary)
}


#' Function to read in the proteomics data
#'
#' @param p_path (str): path to the file proteomics_ack.xlsx
#' 
#' @return 
#'
visulaize_t_test_data <- function(t_test_data){
  
  p <- ggplot(t_test_data, 
                 aes(x=sampleVsample, y=normFoldChange, color=normFoldChange)) +
    geom_boxplot() +
    scale_color_manual(values=c("#81CBD5", "#993342", "#D2B279")) +
    labs(title="Preliminary Visualization",x="Samples", y = "Normalized Fold Change") +
    theme(legend.position="bottom",
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 16))
  
  return(p)
}



#' Function to plot the unadjusted p-values as a histogram
#'
#' @param proteomics_ack (tibble): Tibble with proteomics results
#'
#' @return ggplot: a histogram of the raw max abundance values 
#' from the proteomics results
#' @export
#'
#' @examples abundance_plot <- plot_abundance(proteomics_ack)
plot_abundance <- function(proteomics_ack) {
  # Making ggplot object to return
  a_plot <- proteomics_ack %>% 
    # Making aesthetics based in max abundance column of data
    ggplot(aes(proteomics_ack$`Max Abundance`)) + 
    # coloring plot to look like the example
    geom_histogram(color='black', fill='lightblue') + 
    # setting theme to be simple squares
    theme_minimal() + 
    labs(title="Distribution of Max Abundance",x="x", y = "Abundance") +
    # Setting the location of the plot title
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10))
  
  return(a_plot)
}










#' Perform and plot PCA using processed data.
#' 
#' PCA is performed over genes, and samples should be colored by time point.
#' Both `y` and `x` axis should have percent of explained variance included.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param meta tibble: sample-level meta information (_S_ x 3)
#' @param title string: title for plot
#'
#' @return ggplot: scatter plot showing each sample in the first two PCs. 
#'
#' @examples
#' `plot_pca(data, meta, "Raw Count PCA")`

plot_pca <- function(data, meta, title="") {
  #Performing PCA on a transposed version of the filtered count data
  pca <- prcomp(t(data))
  
  # making a cop of the metadata and assigning new columns names PC1 and PC2
  # while extracting out the PCA results
  metadata <- meta
  PC1 = pca$x[ , 1]
  PC2 = pca$x[ , 2]
  # metadata <- metadata %>% dplyr::mutate(PC1 = PC1)
  # metadata <- metadata %>% dplyr::mutate(PC2 = PC2)
  metadata$PC1 = PC1
  metadata$PC2 = PC2
  
  # Calculating variance using the pca standard deviation
  percent_var <- pca$sdev^2 / sum(pca$sdev^2)
  
  #Constructing / Plotting PCA results
  pca_plot <- ggplot(metadata, aes(x=PC1, y=PC2, col=timepoint)) +
    geom_point() +
    xlab(paste0("PC1 (",round(percent_var[1] * 100),"% variance)")) +
    ylab(paste0("PC2 (",round(percent_var[2] * 100),"% variance)")) +
    ggtitle(title)
  
  return(pca_plot)
}


