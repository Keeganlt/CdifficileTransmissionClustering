## Define Global Functions  ----------------------------------------------

### Not in function -----------

`%!in%` <- Negate(`%in%`)


### Color extractor function --------------

### Extracts colors from ordered.glassby that correspond to the order they are 
### in the phylo tree
legend.color.picker <- function(id, global.sorted.id){
  legend.colors <- NULL
  for(i in 1:length(id)){
    legend.colors[i] <- Ordered.Glasbey.extra[which(id[i] == global.sorted.id)]
  }
  return(legend.colors)
}



### Tip order, shape, and color extractor function --------------

#### Creates a data frame of colors and shapes for room location or granular
#### room location 

tree.color.shape <- function(phylo.tree, metadata){
  p<- ggtree(phylo.tree)
  tree.tib <- as_tibble(phylo.tree)
  tree.tib <- tree.tib %>%
    mutate("sample.site" = label) %>% 
    mutate_at("sample.site", str_extract, "PAT|Z1|Z2|Z3|HH|TO|BR|CC|TB|DH|SN|Reference|Node") %>%
    mutate("Room.Loc" = NA) %>%
    mutate("Sample.Loc" = NA) %>%
    mutate("granular.sample.site" = NA) %>%
    mutate("tox.pos" = NA) 
  
  for(i in 1:dim(tree.tib)[1]){
    if(str_detect(tree.tib$label[i], "Node") == TRUE){tree.tib$label[i]}
  }
  
  for(i in 1:dim(tree.tib)[1]){
    print(tree.tib$label[i])
    if(!is.na(tree.tib$sample.site[i])){
      if(tree.tib$sample.site[i] == "Reference"){
        tree.tib$Room.Loc[i]<- "Reference"
      }else if(tree.tib$sample.site[i] == "Node"){tree.tib$Room.Loc[i] <- "Node"
      }else if(tree.tib$sample.site[i] %in% 
               c("PAT","Z1","Z2","Z3","HH","TO","BR","CC","TB","DH","SN")){tree.tib$Room.Loc[i]<- 
                 metadata$Short.ID[which(metadata$Collection.Sample.ID == tree.tib$label[i])]}
    }
  }
  
  for(i in 1:dim(tree.tib)[1]){
    if(!is.na(tree.tib$sample.site[i])){
      if(tree.tib$sample.site[i] == "Reference"){
        tree.tib$Sample.Loc[i]<- "Reference"
      }else if(tree.tib$sample.site[i] == "Node"){tree.tib$Sample.Loc[i] <- "Node"
      }else if(is.na(tree.tib$Sample.Loc[i])){
        tree.tib$Sample.Loc[i]<- metadata$sample.env[
          which(metadata$Collection.Sample.ID == tree.tib$label[i])]} 
    } 
  }
  
  for(i in 1:dim(tree.tib)[1]){
    if(!is.na(tree.tib$label[i])){tree.tib$granular.sample.site[i] <- tree.tib$sample.site[i]}
  }
  
  toxin.positive.list <-cdiff.metadata$Collection.Sample.ID[which(cdiff.metadata$ToxinA == "Positive")]
  
  for(i in 1:dim(tree.tib)[1]){
    if(!is.na(tree.tib$sample.site[i])){
      if(tree.tib$label[i] %in% toxin.positive.list == TRUE){
        tree.tib$tox.pos[i] <- "+"
      }else{tree.tib$tox.pos[i] <- "-"}
    }
  }
  
  tree.tib <- tree.tib %>% mutate(plot.shape = ifelse(is.na(Sample.Loc), 
                                                      as.integer(NA),paste(Sample.Loc, tox.pos, sep = "")))
  
  
  
  pat.only <- str_extract(tree.tib$granular.sample.site, "PAT")
  
  for(i in 1:length(pat.only)){
    if(grepl("PAT", pat.only[i])){
      tree.tib$granular.sample.site[i] <-word(tree.tib$label[i],sep = "-", start = 3)
      tree.tib$granular.sample.site[i] <-str_extract(tree.tib$granular.sample.site[i], 
                                                     "[:alpha:]")
    }
  }
  return(tree.tib)
}