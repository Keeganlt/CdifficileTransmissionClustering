## Characterizing the spatio-temporal phylogenomic patterns of C. difficile 
## in two densely sampled critical care units 
## Lindsay Keegan
## March 2023


# LOAD LIBRARIES ----------------------------------------------------------

## Install Packages --------------------------------------------------------
{

  library(ape)
  library(cowplot)
  library(dplyr)
  library(GGally)
  library(gdata)
  library(ggnewscale)
  library(ggplot2)
  library(ggraph)
  library(ggtext)
  library(ggtree)
  library(grid)
  library(gridExtra)
  library(igraph)
  library(lubridate)
  library(network)
  library(png)
  library(Polychrome)
  library(RColorBrewer)
  library(readxl)
  library(sna)
  library(stringr)
  library(tidygraph)
  library(tidyverse)
  library(wesanderson)
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
  filter(Bacterial.Species == "C.difficile") %>% 
  group_by(Short.ID) %>% 
  mutate(Patient.Day = as.numeric(as.Date(Collection.date) - as.Date(Admit) +1))


## Make the data set we need
cdiff.descriptive <- cdiff.metadata2 %>% 
  select(Collection.date, Collection.Sample.ID, Room, Short.Room, Sample.ID, 
         PATID, Short.ID, sample.env, Admit, Disch, Cdiff.Toxin, C.diff, 
         Hospital, Room.Start, Room.End, Study.Day, Room.Start.Day, 
         Room.End.Day, contamination, sample.loc, ToxinA, MLST.ST, Short.WGS)


cdiff.descriptive$Admit <- as.Date(cdiff.descriptive$Admit)
cdiff.descriptive$Disch <- as.Date(cdiff.descriptive$Disch)
cdiff.descriptive$Room.Start <- format(parse_date_time(cdiff.descriptive$Room.Start, "YmdHMS"), "%Y-%m-%d")
cdiff.descriptive$Room.End <- format(parse_date_time(cdiff.descriptive$Room.End, "YmdHMS"), "%Y-%m-%d")


### Figure out who is on contact precautions

# Function to calculate contact.prec based on conditions
cdiff.descriptive$contact.prec <- ifelse(cdiff.descriptive$sample.loc %in% c('TO', 'BR', 'CC', 'TB', 'DH'), 1,
                                         ifelse(cdiff.descriptive$sample.loc %in% c('Z1', 'Z2', 'Z3'), 0,
                                                ifelse(cdiff.descriptive$sample.loc == 'SN', 0, NA)))



