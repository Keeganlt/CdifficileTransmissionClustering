
source("code/GlobalFunctions.R")
source("code/LoadData.R")

# CLUSTERING (MAIN TEXT) ----------------------------------------------------

## SNPs cluster -------------------------------------------------------

Room.dist.metadata.long1 <-Room.dist.metadata %>% 
  as_tibble %>%
  pivot_longer (cols = -Room,
                names_to = 'room2',
                values_to = 'dist') %>% 
  drop_na()

colnames(Room.dist.metadata.long1) <- c("room1", "room2", "spatial.dist")

Room.dist.metadata.long2 <- Room.dist.metadata.long1[c("room2", "room1", "spatial.dist")]
colnames(Room.dist.metadata.long2) <- c("room1", "room2", "spatial.dist")
Room.dist.metadata.long<- rbind(Room.dist.metadata.long1, Room.dist.metadata.long2)


## Make the 0, 1, and 2 SNP datasets
merged.room.dist.tall.nodup = merge(dist.lower.meta.tall.all.dup, 
                                    Room.dist.metadata.long1, 
                                    by.x=c("Short.Room1", "Short.Room2"), 
                                    by.y=c("room1", "room2"),
                                    all = TRUE)


dat <- merged.room.dist.tall.nodup %>% 
  filter(dist<=2) %>%
  filter(PATID1 != PATID2) %>%
  select(Collection.Sample.ID1, Collection.Sample.ID2, Short.ID1, dist) 

save(dat, file = "cluster_dat.RData")

colnames(dat)<-c("from", "to", "short.ID", "dist")

edge.dat <- dat %>% 
  select(from, to, dist)

node.dat <- dat %>% 
  select(from, short.ID) %>% 
  distinct(from, short.ID) %>% 
  mutate(sample.loc = from) %>%  
  mutate_at("sample.loc", str_extract, "PAT|Z1|Z2|Z3|HH|TO|BR|CC|TB|DH|SN") %>%
  mutate(sample.loc = case_when(
    sample.loc=="PAT"~"Patient",
    sample.loc=="Z1"~"Room Env.",
    sample.loc=="Z2"~"Room Env.",
    sample.loc=="Z3"~"Room Env.",
    sample.loc=="HH"~"HCW Hand",
    sample.loc=="TO"~"Room Env.",
    sample.loc=="BR"~"Room Env.",
    sample.loc=="CC"~"Room Env.",
    sample.loc=="TB"~"Room Env.",
    sample.loc=="DH"~"Room Env.",
    sample.loc=="SN"~"Shared Env."))

## Extract colors
legend.colors <- legend.color.picker(sort(unique(node.dat$short.ID)), 
                                     sorted.id)
## Convert to graph network format
graph.dat <- tbl_graph(nodes = node.dat, edges = edge.dat, directed = FALSE)

graph_tbl <- graph.dat %>% 
  as_tbl_graph() %>% 
  activate(nodes) 


all.cluster.shortIDs <- sort(unique(node.dat$short.ID))

### Network clustering 0-2 SNPs plot ----------------------------

layout <- create_layout(graph_tbl, layout = 'igraph', algorithm = 'nicely',  
                        maxiter = 100)
all.net.plot <-  ggraph(layout) +
  geom_edge_link(aes(color = factor(dist)))+
  geom_node_point(aes(color = short.ID, shape = as.factor(sample.loc)), 
                  size = 4, show.legend = T)  +
  scale_edge_colour_manual(values = c("#000000FF",  "#525252", "#BDBDBD"),
                           name = "Number of SNPs", guide = "none") +
  scale_color_manual(breaks = sort(factor(unique(node.dat$short.ID))), 
                     values = legend.colors, 
                     name = "Unique \npatient-stay ID") +
  scale_shape_manual(values = c(19, 15, 18, 17), 
                     breaks = c("Patient", "Room Env.",
                                "HCW Hand", "Shared Env."), 
                     name = "Sampling \nEnvironment") + 
  guides(shape = guide_legend(nrow = 4, byrow =TRUE)) +
  coord_equal() +
  theme_void() + 
  theme(legend.position = "bottom", legend.box="horizontal", legend.margin=margin()) 


all.net.plot


counter <- 1
dat$cluster <- NA

##### Loop to cluster the data

counter <- 1
dat$cluster <- NA
while(sum(is.na(dat$cluster))>0){
  length.0 <- 0
  ind<-which(is.na(dat$cluster))[1]
  cur.clust <- c(dat$from[ind], dat$to[ind])
  length.1 <- length(cur.clust)
  while(length.1 > length.0){
    length.0 <- length(cur.clust)
    cur.clust <- unique(c(dat$from[which(dat$to %in% cur.clust)],
                          dat$to[which(dat$from %in% cur.clust)], cur.clust))
    dat$cluster[which(dat$to %in% cur.clust)] <-counter
    dat$cluster[which(dat$from %in% cur.clust)] <-counter
    length.1 <-length(cur.clust)
  }
  counter <- counter +1
}
