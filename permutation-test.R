##############################################################################
##
## Permutation test to compare observed overlap of 
## alternatively spliced genes and ASD risk gene list (SFARI)
## to overlap expected by chance
##
##############################################################################

# Data #######################################################################
#
# SFARI_Genes: https://gene.sfari.org/database/gene-scoring/
# Background_Genes: all genes observed in this study
# Spliced_Genes: all genes with exonic alternative splicing events in PTEN-/- 
#                observed in this study (supplementary tables 6-11)
#                or in NSC-study (Duan et. al 2015)
# RBP_Genes: RBP-database (Cook et al, 2011)
#
# Mm = murine gene nomenclature
# Hs = human gene nomenclature
# 
##############################################################################




library(readr)
library(ggplot2)


## Relevant functions

# simulate random overlap to ASD risk gene list by randomly picking a specified number of genes from whole genome
simulate_overlap <- function(Background_Genes, SFARI_Genes, num_sims, list_size) {
  
  overlap_score <- data.frame(round=1:num_sims, overlap=NA)
  
  for(repetition in overlap_score$round){
    
    Overlap_Genes <- sample(Background_Genes, list_size)
    overlap_score$overlap[repetition] <- length(intersect(SFARI_Genes, Overlap_Genes))
  }
  
  return(overlap_score$overlap)
}


# Calcualte p-value and visualize observed versus simulated overlap by chance
compare_random_to_observed <- function(Background_Genes, SFARI_Genes, Spliced_Genes, num_sims, list_size) {
  
  # The observed overlap of alternatively spliced genes with the ASD dataset
  overlap_observed <- length(intersect(SFARI_Genes, Spliced_Genes))
  
  # Array of overlap by chance
  overlap_random <- data.frame(overlap=simulate_overlap(Background_Genes, SFARI_Genes, num_sims, length(Spliced_Genes)))
  
  Differences <- overlap_random - overlap_observed
  
  # The fraction of obtaining a larger random overlap than the observed overlap
  pvalue <- sum(Differences>0) / num_sims
  
  # plot the simulated overlap by chance and the observed overlap (red line)
  print(plot <- ggplot(data = overlap_random, aes(x = overlap)) + 
      geom_histogram(binwidth = 1) +
      geom_vline(aes(xintercept = overlap_observed), color = "red") + 
      labs(x = "Number of overlapping genes per iteration",
           y = "Count") + theme_minimal())
  
  # Return p-value
  if(pvalue == 0) paste("p < ", 1/num_sims)
  else paste("p = ", pvalue)


}




## load dataframes

# ASD risk genes (SFARI database)
Autism_Genes <- read_delim(file.path("data", "Autism_Genes.txt"), delim = "\n", col_names = F)

# Alternatively spliced in PTEN-KO (mapped against human or mouse or human genome)
Spliced_Hs <- read_delim(file.path("data", "Spliced_Genes_NSC_Hs.txt"), delim = "\n", col_names = F)
Spliced_Mm <- read_delim(file.path("data", "Spliced_Genes_Neurons_Mm.txt"), delim = "\n", col_names = F)


# All detected genes in samples (as a control group)
Background_Genes_Hs <- read_delim(file.path("data", "All_Genes_Hs.txt"), delim = "\n", col_names = F)
Background_Genes_Mm <- read_delim(file.path("data", "All_Genes_Mm.txt"), delim = "\n", col_names = F)


# RBPs
RBP_Hs <- read_delim(file.path("data", "RBP_Genes_Hs.txt"), delim = "\n", col_names = F)
RBP_Mm <- read_delim(file.path("data", "RBP_Genes_Mm.txt"), delim = "\n", col_names = F)


## Harmonize case of gene names
Background_Genes_Hs <- tolower(Background_Genes_Hs$X1)
Background_Genes_Mm <- tolower(Background_Genes_Mm$X1)
SFARI_Genes <- tolower(Autism_Genes$X1)
Spliced_Hs <- tolower(Spliced_Hs$X1)
Spliced_Mm <- tolower(Spliced_Mm$X1)
RBP_Hs <- tolower(RBP_Hs$X1)
RBP_Mm <- tolower(RBP_Mm$X1)



## Final calculations
num_sim = 100000

# Overlap of altered splicing in neurons and ASD risk gene dataset
compare_random_to_observed(Background_Genes_Mm, SFARI_Genes, Spliced_Mm, num_sim, length(Spliced_Mm))

# Overlap of altered splicing in NSCs (Duan et al., 2015) and ASD risk gene dataset
compare_random_to_observed(Background_Genes_Hs, SFARI_Genes, Spliced_Hs, num_sim, length(Spliced_Hs))


# Overlap of altered splicing in neurons and RBP dataset
compare_random_to_observed(Background_Genes_Mm, RBP_Mm, Spliced_Mm, num_sim, length(Spliced_Mm))

# Overlap of altered splicing in NSCs (Duan et al., 2015) and RBP dataset
compare_random_to_observed(Background_Genes_Hs, RBP_Hs, Spliced_Hs, num_sim, length(Spliced_Hs))
