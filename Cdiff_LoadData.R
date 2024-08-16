## Characterizing the spatio-temporal phylogenomic patterns of C. difficile 
## in two densely sampled critical care units 
## Lindsay Keegan
## March 2023


# LOAD LIBRARIES ----------------------------------------------------------

## Install Packages --------------------------------------------------------
{
  library(RColorBrewer)
  library(tidyverse)
  library(ggtext)
  library(ggtree)
  library(ggnewscale)
  library(gdata)
  library(dplyr)
  library(ape)
  library(ggplot2)
  library(stringr)
  library(grid)
  library(gridExtra)
  library(cowplot)
  library(Polychrome)
  library(GGally)
  library(network)
  library(sna)
  library(ggplot2)
  library(igraph)
  library(wesanderson)
  library(ggraph)
  library(tidygraph)
  library(lubridate)
  library(readxl)
  library(Polychrome)
  library(readxl)
  library(RColorBrewer)
  
}

# DATA PRE-PROCESSING -----------------------------------------------------

## Define Colors -----------------------------------------------------------
{
  Glasbey = as.character(glasbey.colors(32)) # Color palette
  Glasbey[32] <-"#ebc38d" # change the colors I don't like 
  Glasbey[1] <-"khaki2"
  Glasbey[5] <-"#3c32a8"
  
  
  
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  # Order the colors by hue
  Ordered.Glasbey.extra <- c(  "#C9CC3F", Glasbey[23], "#baffba", Glasbey[4], Glasbey[22],
                               "#669e05", Glasbey[7], Glasbey[13],"#6d8567", 
                               Glasbey[11], Glasbey[28], Glasbey[31], Glasbey[9], 
                               Glasbey[2],  Glasbey[5], Glasbey[29], Glasbey[21], 
                               "#DECBE4","#C9ACF9", "#7c5496", Glasbey[12], 
                               Glasbey[20], "#808080", "#FED9A6",
                               Glasbey[25], "#f58884", Glasbey[14], Glasbey[6], 
                               Glasbey[30], Glasbey[16],  Glasbey[3],"#E41A1C", Glasbey[27], 
                               Glasbey[17], "#f56a0a", "#A65628", "black")
  

  Ordered.Glasbey <- c( "#C9CC3F", Glasbey[23], "#baffba", Glasbey[4], Glasbey[22],
                        "#669e05", Glasbey[7], Glasbey[13],"#6d8567", 
                        Glasbey[11], Glasbey[28], Glasbey[31], Glasbey[9], 
                        Glasbey[2],  Glasbey[5], Glasbey[29], Glasbey[21], 
                        "#DECBE4","#C9ACF9", "#7c5496", Glasbey[12], 
                        Glasbey[20], "#808080", "#FED9A6",
                        Glasbey[25], "#f58884", Glasbey[14], Glasbey[6], 
                        Glasbey[30], Glasbey[16],  Glasbey[3], Glasbey[27], 
                        Glasbey[17], "#A65628")
  swatch(Ordered.Glasbey)
  
  
  Ordered.Glasbey.plus1 <-  c( "#C9CC3F", Glasbey[23], "#baffba", Glasbey[4], Glasbey[22],
                               "#669e05", Glasbey[7], Glasbey[13],"#6d8567", 
                               Glasbey[11], Glasbey[28], Glasbey[31], Glasbey[9], 
                               Glasbey[2],  Glasbey[5], Glasbey[29], Glasbey[21], 
                               "#DECBE4","#C9ACF9", "#7c5496", Glasbey[12], 
                               Glasbey[20], "#808080", "#FED9A6",
                               Glasbey[25], "#f58884", Glasbey[14], Glasbey[6], 
                               Glasbey[30], Glasbey[16],  Glasbey[3], Glasbey[27], 
                               Glasbey[17],  "#A65628", "black")
}

## 1.2 Read in data ------------------------------------------------------------
{
  gm.metadata= read.csv(file = "gm_metadata_final.csv") # Read in the Metadata
  
  Room.dist.metadata = read.csv(file = "Input_data/Distance_Between_Rooms.csv", 
                                header = TRUE) # Distance between rooms
  
    ## MLST Data
  mlst <- read.csv("C_diffMLST.csv", header = TRUE)
  gm.metadata <- merge(gm.metadata, mlst, by = "Short.WGS", all.x = TRUE)
  
  snp.dsts <- read.csv("Distances/clade1_distances.csv")
  
  ## Phylogenetic tree data 
  tree.clade1 = read.tree("Genome/clade1-13Apr.full.final_bootstrapped_tree.tre")
  tree.clade2 = read.tree("Genome/clade2-13Apr.full.final_bootstrapped_tree.tre")
  tree.clade4 = read.tree("Genome/clade4_C4.core.full.final_bootstrapped_tree.tre")
  tree.all.boot = read.tree("Genome/clade4.core.full.final_bootstrapped_tree.tre")
  tree.clade5 = tree.all.boot
}



## Filter out everything but C. diff
temp.cdiff <- gm.metadata %>%  
  filter(Bacterial.Species == "C.difficile") 

pat.list <- unique(temp.cdiff$PATID)
pat.list = pat.list[!is.na(pat.list)]

cdiff.metadata2 <- gm.metadata %>%  
  filter(PATID %in% pat.list | Collection.Sample.ID == "A-SN-014-1")

cdiff.metadata<-gm.metadata %>%  
  filter(Bacterial.Species == "C.difficile")  
