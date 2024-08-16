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

