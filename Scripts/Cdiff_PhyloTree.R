source("Scripts/Cdiff_GlobalFunctions.R")
source("Scripts/Cdiff_LoadData.R")

## Process data for phylogenetic tree --------------------------------------

### Rename the tip labels in the tree from WGS Sample ID to Collection Sample ID
rename.tips <- function(tree){
  for(i in 1:dim(cdiff.metadata)[1]){
    tree$tip.label[which((tree$tip.label) == 
                           (cdiff.metadata$WGS.Sample.ID[i]))] <- 
      cdiff.metadata$Collection.Sample.ID[i]
  }
  return(tree)
}

tree.clade1 <- rename.tips(tree = tree.clade1)
tree.clade2 <- rename.tips(tree = tree.clade2)
tree.clade4 <- rename.tips(tree = tree.clade4)
tree.all.boot <- rename.tips(tree = tree.all.boot)
tree.clade5 <- rename.tips(tree = tree.clade5)


### Simplify the  metadata
#### Make a smaller metadata data frame with only the columns we need
cdiff.metadata.simple <- cdiff.metadata %>% 
  as_tibble %>% 
  select(WGS.Sample.ID, Cdiff.Toxin, PATID, Sample.ID, Collection.Sample.ID, 
         Collection.date, Collection.time, Room, Short.ID, Short.Room, 
         sample.env, sample.loc, Unit, Hospital, Study.Day, Room.Start.Day, 
         Room.End.Day, contamination, Patient.Day)

cdiff.metadata.simple.X <- cdiff.metadata.simple
cdiff.metadata.simple.X$WGS.Sample.ID <- paste("X",
                                               cdiff.metadata.simple.X$WGS.Sample.ID, sep="")


### Long data set with duplicates
### Remove upper half of the triangular matrix and flip from wide to long

create.snp.lower.matrix.all <- function(snp.distance){
  dat <- snp.distance %>%
    as_tibble %>% 
    pivot_longer(cols = -"snp.dists.0.8.2",
                 names_to = 'sample2',
                 values_to = 'dist') %>% 
    left_join(cdiff.metadata.simple %>% rename_with(.fn = ~paste0(.,'1')), 
              by = c('snp.dists.0.8.2'='WGS.Sample.ID1')) %>% 
    left_join(cdiff.metadata.simple.X %>% rename_with(.fn = ~paste0(.,'2')), 
              by = c('sample2'='WGS.Sample.ID2'))
  return(dat)
}

## Distance between clades
lower.all.dups <- create.snp.lower.matrix.all(snp.dsts)

### Create a single data set from all 3 clades
dist.lower.meta.tall.all.dup <- lower.all.dups

### Drop distance 0
#### Pivot from wide to long
create.snp.lower.matrix.drop.self <- function(snp.distance){
  dat <- snp.distance %>%
    as_tibble %>% 
    pivot_longer(cols = -snp.dists.0.8.2,
                 names_to = 'sample2',
                 values_to = 'dist') %>% 
    filter(snp.dists.0.8.2!=sample2) %>% 
    left_join(cdiff.metadata.simple %>% rename_with(.fn = ~paste0(.,'1')), 
              by = c('snp.dists.0.8.2'='WGS.Sample.ID1')) %>% 
    left_join(cdiff.metadata.simple.X %>% rename_with(.fn = ~paste0(.,'2')), 
              by = c('sample2'='WGS.Sample.ID2'))
  return(dat)
}


## Isolate distances 
lower.all <- create.snp.lower.matrix.drop.self(snp.dsts)
dist.lower.meta.tall.all <- lower.all


### filter to dist 0:  no longer all pairwise comparisons

dist.lower.meta.tall.filtered.specific <- dist.lower.meta.tall.all %>% 
  filter(dist==0) %>%
  filter(PATID1==PATID2) %>% 
  filter(Cdiff.Toxin1 == Cdiff.Toxin2) %>% 
  filter(sample.loc1==sample.loc2)

specific.unique.samples <- 
  unique(c(dist.lower.meta.tall.filtered.specific$Collection.Sample.ID1))

specific.unique.samples.df <- tibble(specific.unique.samples,group = NA_integer_)
for (i in 1:length(specific.unique.samples)){
  if(i==1){
    j=1 #count for group
    specific.unique.samples.df$group[i] <- j
  }
  else{
    #find which rows ith person belongs to in tall filtered distance matrix
    matched_rows_index <- 
      which(dist.lower.meta.tall.filtered.specific$Collection.Sample.ID1 == 
              specific.unique.samples.df$specific.unique.samples[i])
    #of those rows, find corresponding partner
    matched_comparison_ids <- dist.lower.meta.tall.filtered.specific$
      Collection.Sample.ID2[matched_rows_index]
    #find distinct groups assigned to partners -- if no assignment, then length == 0
    matched_comparison_groups <- specific.unique.samples.df %>% 
      filter(specific.unique.samples %in% matched_comparison_ids) %>% 
      filter(!is.na(group)) %>% .$group %>% unique
    if(length(matched_comparison_groups)==0L){
      #assign new group
      j=j+1
      specific.unique.samples.df$group[i] <- j
    } else if(length(matched_comparison_groups)>1L){
      #throw error because there is more than one group
      stop("uhh ohh spaghetti-os")
    } else if(length(matched_comparison_groups)==1){
      #assign as matched_comparison_groups
      specific.unique.samples.df$group[i] <- matched_comparison_groups
    }
  }
}

specific.unique.samples.df <- specific.unique.samples.df %>% 
  group_by(group) %>% 
  mutate(keep_me = row_number()==1L) %>% 
  ungroup

## List of which to drop 
specific.unique.samples.drop <- specific.unique.samples.df %>% 
  filter(!keep_me)

### Repeat for dropping any environmental match 
dist.lower.meta.tall.filtered.general <- dist.lower.meta.tall.all %>% 
  filter(dist==0) %>%
  filter(PATID1==PATID2) %>% 
  filter(Cdiff.Toxin1 == Cdiff.Toxin2) %>% 
  filter(sample.env1==sample.env2)

general.unique.samples <- 
  unique(c(dist.lower.meta.tall.filtered.general$Collection.Sample.ID1))
general.unique.samples.df <- tibble(general.unique.samples,group = NA_integer_)
for (i in 1:length(general.unique.samples)){
  if(i==1){
    j=1 #count for group
    general.unique.samples.df$group[i] <- j
  }
  else{
    #find which rows ith person belongs to in tall filtered distance matrix
    matched_rows_index <- 
      which(dist.lower.meta.tall.filtered.general$Collection.Sample.ID1 ==
              general.unique.samples.df$general.unique.samples[i])
    #of those rows, find corresponding partner
    matched_comparison_ids <- dist.lower.meta.tall.filtered.general$
      Collection.Sample.ID2[matched_rows_index]
    #find distinct groups assigned to partners -- if no assignment, then length == 0
    matched_comparison_groups <- general.unique.samples.df %>% 
      filter(general.unique.samples %in% matched_comparison_ids) %>% 
      filter(!is.na(group)) %>% .$group %>% unique
    if(length(matched_comparison_groups)==0L){
      #assign new group
      j=j+1
      general.unique.samples.df$group[i] <- j
    } else if(length(matched_comparison_groups)>1L){
      #throw error because there is more than one group
      stop("uhh ohh spaghetti-os")
    } else if(length(matched_comparison_groups)==1){
      #assign as matched_comparison_groups
      general.unique.samples.df$group[i] <- matched_comparison_groups
    }
  }
}

general.unique.samples.df <- general.unique.samples.df %>% 
  group_by(group) %>% 
  mutate(keep_me = row_number()==1L) %>% 
  ungroup

## Tips to drop if you want to drop all room environment not just exact room environment
general.unique.samples.drop <- general.unique.samples.df %>% 
  filter(!keep_me)



phylo.tree <- tree.clade5

## Define colors
global.color.shape <- tree.color.shape(phylo.tree, cdiff.metadata) 


### 3.1.1 Figure: all isolate plot ----------
p.all.iso <- ggtree(all.iso.tree) +
  geom_tippoint(aes(color = Room.Loc, shape = factor(Sample.Loc)), stroke = 1, size = 1.5) +
  scale_color_manual(values = Ordered.Glasbey.plus1, 
                     name = "Unique patient ID") +
  scale_shape_manual(breaks = unique(all.iso.dat$Sample.Loc), 
                     values = c("HCW Hand" = 18, "Patient" = 19, "Room Env." =15, 
                                "Shared Env." = 17, "Reference" = 82), 
                     name = "Sampling location")  +
  # theme_cowplot()  +
  geom_nodelab(hjust = 1.5, vjust = 0) +
  theme(axis.text.y = element_text(colour="white"), 
        axis.ticks.y = element_line(colour="white"), 
        axis.line.y=element_line(colour="white")) +
  geom_treescale(offset=1) +
  geom_cladelab(node=179, label="Clade 2", angle=270,
                hjust=0.6, offset.text=1,  align=TRUE, vjust =0) +
  geom_cladelab(node=182, label="Clade 1", angle=270,
                hjust='center', offset.text=1,  align=TRUE, vjust =0) +
  geom_cladelab(node= 335, label="Clade 4", angle=270,
                hjust=0.35,  align=TRUE, vjust =0) 

# xlab("SNPs") 

p.all.iso


tree.clade5.tox.dat <- all.iso.dat %>% 
  select(label, tox.pos) %>%
  filter(!is.na(tox.pos))
tree.clade5.tox.dat$tox.pos <- ifelse(tree.clade5.tox.dat$tox.pos == "-", "Negative", "Positive")
tree.clade5.tox.dat <-as.data.frame(tree.clade5.tox.dat)
colnames(tree.clade5.tox.dat) <-c("label", "Toxin")
row.names(tree.clade5.tox.dat) <- tree.clade5.tox.dat$label
tree.clade5.tox.dat$label <- as.character(tree.clade5.tox.dat$label)
col <- "Toxin"

df <- tree.clade5.tox.dat[tree.clade5$tip.label,][col]

p2 <- p.all.iso + new_scale_fill()

p.all.iso.tox <- gheatmap(p2, df, 
                          colnames=FALSE, 
                          colnames_position = "top", width=0.01, offset = 5) + 
  guides(color=guide_legend(nrow=5,byrow=FALSE)) + 
  # new_scale_fill()  + 
  scale_fill_manual(values=c("#B2C1DC","#F8766D"), 
                    name= expression(paste(italic("C. difficile"), " Toxin")),
                    guide = guide_legend(title.position = "top", order = 0)) + 
  ggtree::vexpand(0.05, 1) 

p1 <- p.all.iso.tox + 
  theme(legend.position = "none")

p.legend<- cowplot::get_plot_component(p.all.iso.tox, "guide-box", return_all = TRUE)[[1]] 


p.tree <- ggdraw(p1) +
  draw_plot(p.legend, 0.15, 0.6, .5, .5) # +
# draw_plot_label(
#   c("A", "B"),
#   c(0, 0.45),
#   c(1, 0.95),
#   size = 12
# )