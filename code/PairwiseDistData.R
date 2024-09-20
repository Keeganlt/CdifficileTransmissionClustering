source("code/Functions.R")
source("code/LoadData.R")


# WITHIN AND BETWEEN SAMPLE DISTANCES --------------------------------------

## Create the data ------
pairwise.dist.dat <- dist.lower.meta.tall.all.dup %>% 
  filter(snp.dists.0.8.2 != "Reference", sample2 != "Reference")%>% 
  mutate(pairwise.dist.name = NA)

for(i in 1:dim(pairwise.dist.dat)[1]){
  # print(i)
  if(pairwise.dist.dat$Short.ID1[i] == pairwise.dist.dat$Short.ID2[i]){
    if(pairwise.dist.dat$sample.env1[i] == "Patient" & pairwise.dist.dat$sample.env2[i] == "Patient"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "sppp"
    } else if(pairwise.dist.dat$sample.env1[i] == "Patient" & pairwise.dist.dat$sample.env2[i] == "Room Env."){
      pairwise.dist.dat$pairwise.dist.name[i] <- "sppe"
    } else if(pairwise.dist.dat$sample.env1[i] == "Room Env." & pairwise.dist.dat$sample.env2[i] == "Patient"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "sppe"
    } else if(pairwise.dist.dat$sample.env1[i] == "Patient" & pairwise.dist.dat$sample.env2[i] == "HCW Hand"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "spph"
    } else if(pairwise.dist.dat$sample.env1[i] == "HCW Hand" & pairwise.dist.dat$sample.env2[i] == "Patient"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "spph"
    } else if(pairwise.dist.dat$sample.env1[i] == "Room Env." & pairwise.dist.dat$sample.env2[i] == "Room Env."){
      pairwise.dist.dat$pairwise.dist.name[i] <- "spee"
    } else if(pairwise.dist.dat$sample.env1[i] == "Room Env." & pairwise.dist.dat$sample.env2[i] == "HCW Hand"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "speh"
    } else if(pairwise.dist.dat$sample.env1[i] == "HCW Hand" & pairwise.dist.dat$sample.env2[i] == "Room Env."){
      pairwise.dist.dat$pairwise.dist.name[i] <- "speh"
    } else if(pairwise.dist.dat$sample.env1[i] == "HCW Hand" & pairwise.dist.dat$sample.env2[i] == "HCW Hand"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "sphh"
    } else if(pairwise.dist.dat$sample.env1[i] == "Shared Env." & pairwise.dist.dat$sample.env2[i] == "Shared Env."){
      pairwise.dist.dat$pairwise.dist.name[i] <- "spss"
    }
  } else{
    if(pairwise.dist.dat$sample.env1[i] == "Patient" & pairwise.dist.dat$sample.env2[i] == "Patient"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dppp"
    } else if(pairwise.dist.dat$sample.env1[i] == "Patient" & pairwise.dist.dat$sample.env2[i] == "Room Env."){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dppe"
    } else if(pairwise.dist.dat$sample.env1[i] == "Room Env." & pairwise.dist.dat$sample.env2[i] == "Patient"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dppe"
    } else if(pairwise.dist.dat$sample.env1[i] == "Patient" & pairwise.dist.dat$sample.env2[i] == "HCW Hand"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dpph"
    }  else if(pairwise.dist.dat$sample.env1[i] == "HCW Hand" & pairwise.dist.dat$sample.env2[i] == "Patient"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dpph"
    } else if(pairwise.dist.dat$sample.env1[i] == "Patient" & pairwise.dist.dat$sample.env2[i] == "Shared Env."){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dpps"
    } else if(pairwise.dist.dat$sample.env1[i] == "Shared Env." & pairwise.dist.dat$sample.env2[i] == "Patient"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dpps"
    } else if(pairwise.dist.dat$sample.env1[i] == "Room Env." & pairwise.dist.dat$sample.env2[i] == "Room Env."){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dpee"
    } else if(pairwise.dist.dat$sample.env1[i] == "Room Env." & pairwise.dist.dat$sample.env2[i] == "HCW Hand"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dpeh"
    } else if(pairwise.dist.dat$sample.env1[i] == "HCW Hand" & pairwise.dist.dat$sample.env2[i] == "Room Env."){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dpeh"
    } else if(pairwise.dist.dat$sample.env1[i] == "Room Env." & pairwise.dist.dat$sample.env2[i] == "Shared Env."){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dpes"
    } else if(pairwise.dist.dat$sample.env1[i] == "Shared Env." & pairwise.dist.dat$sample.env2[i] == "Room Env."){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dpes"
    } else if(pairwise.dist.dat$sample.env1[i] == "HCW Hand" & pairwise.dist.dat$sample.env2[i] == "HCW Hand"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dphh"
    } else if(pairwise.dist.dat$sample.env1[i] == "HCW Hand" & pairwise.dist.dat$sample.env2[i] == "Shared Env."){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dphs"
    } else if(pairwise.dist.dat$sample.env1[i] == "Shared Env." & pairwise.dist.dat$sample.env2[i] == "HCW Hand"){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dphs"
    } else if(pairwise.dist.dat$sample.env1[i] == "Shared Env." & pairwise.dist.dat$sample.env2[i] == "Shared Env."){
      pairwise.dist.dat$pairwise.dist.name[i] <- "dpss"
    }
  }
}


### How many same pt-pt pairs don't cluster

same.pt.pt <- pairwise.dist.dat %>% filter(pairwise.dist.name == "sppp") 


(sum(same.pt.pt$dist <=2) /length(same.pt.pt$dist))*100

### Figure 2: pairwise distance ----------------------------------

pairwise.dist.dat.one <- pairwise.dist.dat %>% 
  mutate(sample.combination.sorted = if_else(Collection.Sample.ID1<Collection.Sample.ID2,
                                             paste0(Collection.Sample.ID1,":",Collection.Sample.ID2),
                                             paste0(Collection.Sample.ID2,":",Collection.Sample.ID1))) %>% 
  group_by(sample.combination.sorted) %>% 
  slice(1) %>% 
  ungroup
