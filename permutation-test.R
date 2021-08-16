## Permutation test to compare observed overlap of 
## alternatively spliced genes and ASD risk gene list (SFARI)
## to overlap expected by chance
##
################################################################

# SFARI_Genes: https://gene.sfari.org/database/gene-scoring/
# Background_Genes: all genes observed in this study
# Spliced_Genes: all genes with exonic alternative splicing events in PTEN-/- 
#                observed in this study (supplementary tables 6-11)


# simulate random overlap to ASD risk gene list by randomly picking a specified number of genes from whole genome
simulate_overlap <- function(Background_Genes, SFARI_Genes, num_sims, list_size) {
  
  overlap_score <- data.frame(round=1:num_sims, overlap=NA)
  
  for(repetition in overlap_score$round){
    
    Overlap_Genes <- sample(Background_Genes, list_size)
    overlap_score$overlap[repetition] <- length(intersect(SFARI_Genes, Overlap_Genes))
  }
  
  return(overlap_score$overlap)
}

# The observed overlap of alternatively spliced genes with the ASD dataset
overlap_observed <- length(intersect(SFARI_Genes, Spliced_Genes))

# Array of overlap by chance
overlap_random <- data.frame(overlap=simulate_overlap(Background_Genes, SFARI_Genes, num_sims, length(Spliced_Genes)))

Differences <- overlap_random - overlap_observed

# The fraction of obtaining a larger random overlap than the observed overlap
pvalue <- sum(Differences>0)/length(Differences)
