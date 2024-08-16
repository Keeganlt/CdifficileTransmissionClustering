
source("Scripts/Cdiff_GlobalFunctions.R")
source("Scripts/Cdiff_LoadData.R")


# CLUSTERING (Supplement TEXT) ----------------------------------------------------

## 8.2 0-7 SNPs cluster -------------------------------------------------------

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


## Make the 0-7 SNP datasets
merged.room.dist.tall.nodup = merge(dist.lower.meta.tall.all.dup, 
                                    Room.dist.metadata.long1, 
                                    by.x=c("Short.Room1", "Short.Room2"), 
                                    by.y=c("room1", "room2"),
                                    all = TRUE)


dat <- merged.room.dist.tall.nodup %>% 
  filter(dist<=7) %>%
  filter(PATID1 != PATID2) %>%
  select(Collection.Sample.ID1, Collection.Sample.ID2, Short.ID1, dist) 

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



### Figure Network clustering 0-7 SNPs plot ----------------------------


layout <- create_layout(graph_tbl, layout = 'igraph', algorithm = 'nicely',  
                        maxiter = 100)
all.net.plot.5 <-  ggraph(layout) +
  geom_edge_link(aes(color = factor(dist)))+
  geom_node_point(aes(color = short.ID, shape = as.factor(sample.loc)), 
                  size = 4, show.legend = T)  +
  scale_edge_colour_manual(values = c("#000000", "#252525", "#525252", "#737373", "#969696","#BDBDBD", "#D9D9D9", "#F0F0F0"),
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
  # annotate("text", x=-3, y=10, label= "A.") + 
  # annotate("text", x=4, y=7, label= "B.") + 
  # annotate("text", x=6, y=4, label= "C.") + 
  # annotate("text", x=2.5, y=-4, label= "D.") + 
  # annotate("text", x=-5.2, y=-8, label= "E.") + 
  # annotate("text", x=-5.2, y=2, label= "F.") + 
  # annotate("text", x=-11.4, y=4.1, label= "G.") + 
  theme(legend.position = "bottom", legend.box="horizontal", legend.margin=margin()) 


all.net.plot.5

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


##Function to make network plots
network.plot.fun <- function(cluster.num, dat, lay.out){
  edge.dat <- dat %>% 
    filter(cluster == cluster.num) %>% 
    select(from, to, dist)
  
  node.dat <- dat %>% 
    filter(cluster == cluster.num) %>% 
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
  legend.colors <- legend.color.picker(unique(node.dat$short.ID), 
                                       sorted.id)
  ## Convert to graph network format
  graph.dat <- tbl_graph(nodes = node.dat, edges = edge.dat, directed = FALSE)
  
  graph_tbl <- graph.dat %>% 
    as_tbl_graph() %>% 
    activate(nodes) 
  if(lay.out == "igraph"){
    layout <- create_layout(graph_tbl, layout = "igraph", 
                            algorithm = 'nicely', maxiter = 100)
    
  }else{
    layout <- create_layout(graph_tbl, layout = lay.out)
    
  }
  plot <- ggraph(layout) +
    geom_edge_link(aes(color = factor(dist)))+
    geom_node_point(aes(color = short.ID, shape = as.factor(sample.loc)), 
                    size = 5, show.legend = F)  +
    scale_edge_colour_manual(values = c("#000000", "#252525", "#525252", "#737373", "#969696","#BDBDBD", "#D9D9D9", "#F0F0F0"),
                             #values = c("#000000FF", "#494848", "#909090", "#B4B4B4", "#D4D4D4","#ECECECFF"),
                             name = "Number of SNPs", guide = "none") +
    scale_color_manual(breaks = factor(unique(node.dat$short.ID)), 
                       values = legend.colors, 
                       name = "Unique ID") +
    scale_shape_manual(values = c(19, 15, 18, 17), 
                       breaks = c("Patient", "Room Env.",
                                  "HCW Hand", "Shared Env."), 
                       name = "Sampling Environment") +
    coord_equal() + 
    theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
          legend.position = "none") +
    theme_void()
  return(plot)
}

## Loop over all the plots
num.clusts <- seq(from = 1, to= counter-1, by = 1)



for(i in 1:length(num.clusts)){
  assign(paste("p.net",num.clusts[i], sep ="." ),
         network.plot.fun(num.clusts[i], dat, lay.out = "igraph"))
}

sample.timing.plot.fun <- function(cluster.num, dat, x.label){
  node.dat <- dat %>% 
    filter(cluster == cluster.num) %>% 
    select(from, short.ID) %>% 
    distinct(from, short.ID) %>% 
    mutate(sample.loc = from) 
  
  short.name.1 <- unique(node.dat$from)
  
  ## Make the dataset we need
  cdiff.descriptive.1 <- cdiff.descriptive %>% 
    # filter(str_detect(Collection.Sample.ID, "VA")) %>% 
    filter(Collection.Sample.ID %in% short.name.1, negate = TRUE)
  
  cdiff.descriptive.1
  
  edge.dat <- dat %>% 
    filter(cluster == cluster.num) # %>% 
  # select(from, to, dist)
  
  short.name.2 <- unique(edge.dat$short.ID)
  
  ## Make the dataset we need
  cdiff.descriptive.2 <- cdiff.descriptive %>% 
    # filter(str_detect(Collection.Sample.ID, "VA")) %>% 
    filter(Short.ID %in% short.name.2, negate = TRUE)
  
  if(length(unique(cdiff.descriptive.1$Short.ID)) < 4){
    error.bar.height <- 0.2
  } else{
    error.bar.height <- 0.6
  }
  
  
  legend.colors <- legend.color.picker(unique(cdiff.descriptive.1$Short.ID), 
                                       sorted.id)
  if(x.label == TRUE){
    plot <- ggplot() + 
      geom_point(data = cdiff.descriptive.1 %>% filter(C.diff == 1 &contamination == 1), 
                 aes(x = Study.Day, y= Short.Room, shape = sample.env), 
                 color = "darkgrey", size = 6) +
      geom_point(data = cdiff.descriptive.1 %>% filter(!is.na(Cdiff.Toxin)), 
                 aes(x = Study.Day, y = Short.Room, color = Short.ID, 
                     shape = sample.env), size = 3) + 
      geom_point(data = cdiff.descriptive.2 %>% filter(!is.na(Cdiff.Toxin)), 
                 aes(x = Study.Day, y = Short.Room, color = Short.ID, 
                     shape = sample.env), size = 3, alpha =0.1)  +
      geom_errorbarh(data = cdiff.descriptive.1, aes(xmax = Room.Start.Day, 
                                                     xmin = Room.End.Day, y= Short.Room,
                                                     color = Short.ID), height = error.bar.height) + 
      geom_point(data = cdiff.descriptive.1 %>% filter(ToxinA == "Positive"), 
                 aes(x = Study.Day, y= Short.Room), color = "black", shape = 3, stroke = 1.3, size = 3) + 
      geom_point(data = cdiff.descriptive.2 %>% filter(ToxinA == "Positive"), 
                 aes(x = Study.Day, y= Short.Room), color = "black", shape = 3,alpha =0.1, size = 3) + 
      scale_color_manual(breaks = factor(unique(cdiff.descriptive.1$Short.ID)), 
                         values = legend.colors, 
                         name = "Unique ID") + 
      ylab("Room Number") +
      xlab("Study Day") +
      theme_cowplot() +
      background_grid() +
      scale_shape_manual(values = c(19, 15, 18, 17), 
                         breaks = c("Patient", "Room Env.",
                                    "HCW Hand", "Shared Env."),
                         name = "Sample Location") +
      theme(text = element_text(size=18), 
            axis.text.x = element_text(angle = 15, hjust = 1), 
            legend.position = "none") +
      scale_y_discrete(drop=T, 
                       labels=c("A13" = "1", "A14" = "2","A15" = "3", 
                                "A16" = "4","A17" = "5", "A18" = "6",
                                "A19" = "7", "A20" = "8", "A21" = "9",
                                "A22" = "10", "A23" = "11", "A24" = "12", 
                                "A25" = "13", "A26" = "14", "A27" = "15",
                                "A28" = "16", "A28.5" = "16.5", "A29" = "17",
                                "A30" = "18", "B1" = "1", "B2" = "2", 
                                "B3" = "3", "B4" = "4", "B5" = "5", 
                                "B6" = "6", "B7" = "7", "B8" = "8",
                                "B9" = "9", "B10" = "10"))
    
    return(plot)
  }else{
    plot <- ggplot() + 
      geom_point(data = cdiff.descriptive.1 %>% filter(C.diff == 1 &contamination == 1), 
                 aes(x = Study.Day, y= Short.Room, shape = sample.env), 
                 color = "darkgrey", size = 6) +
      geom_point(data = cdiff.descriptive.1 %>% filter(!is.na(Cdiff.Toxin)), 
                 aes(x = Study.Day, y = Short.Room, color = Short.ID, 
                     shape = sample.env), size = 3) + 
      geom_point(data = cdiff.descriptive.2 %>% filter(!is.na(Cdiff.Toxin)), 
                 aes(x = Study.Day, y = Short.Room, color = Short.ID, 
                     shape = sample.env), size = 3, alpha =0.1) +
      # geom_point(data = cdiff.descriptive.1 %>% filter(!is.na(Cdiff.Toxin)), 
      #            aes(x = Study.Day, y = Short.Room, color = Short.ID, 
      #                shape = sample.env), size = 3, alpha = 0.1) + 
      # geom_point(data = cdiff.descriptive.2 %>% filter(!is.na(Cdiff.Toxin)), 
      #            aes(x = Study.Day, y = Short.Room, color = Short.ID, 
      #                shape = sample.env), size = 3)  +
      geom_errorbarh(data = cdiff.descriptive.2, aes(xmax = Room.Start.Day, 
                                                     xmin = Room.End.Day, y= Short.Room,
                                                     color = Short.ID), height = error.bar.height) + 
      geom_point(data = cdiff.descriptive.1 %>% filter(ToxinA == "Positive"), 
                 aes(x = Study.Day, y= Short.Room), color = "black", shape = 3, stroke = 1.3, size = 3) + 
      geom_point(data = cdiff.descriptive.2 %>% filter(ToxinA == "Positive"), 
                 aes(x = Study.Day, y= Short.Room), color = "black", shape = 3,alpha =0.1, size = 3) + 
      scale_color_manual(breaks = factor(unique(cdiff.descriptive.1$Short.ID)), 
                         values = legend.colors, 
                         name = "Unique ID") + 
      ylab("Room Number") +
      xlab("Study Day") +
      xlab(" ") + 
      theme_cowplot() +
      background_grid() +
      scale_shape_manual(values = c(19, 15, 18, 17), 
                         breaks = c("Patient", "Room Env.",
                                    "HCW Hand", "Shared Env."),
                         name = "Sample Location") +
      theme(text = element_text(size=18), 
            axis.text.x = element_text(angle = 15, hjust = 1), 
            legend.position = "none") +
      scale_y_discrete(drop=T, 
                       labels=c("A13" = "1", "A14" = "2","A15" = "3", 
                                "A16" = "4","A17" = "5", "A18" = "6",
                                "A19" = "7", "A20" = "8", "A21" = "9",
                                "A22" = "10", "A23" = "11", "A24" = "12", 
                                "A25" = "13", "A26" = "14", "A27" = "15",
                                "A28" = "16", "A28.5" = "16.5", "A29" = "17",
                                "A30" = "18", "B1" = "1", "B2" = "2", 
                                "B3" = "3", "B4" = "4", "B5" = "5", 
                                "B6" = "6", "B7" = "7", "B8" = "8",
                                "B9" = "9", "B10" = "10"))
    
    if("A056" %in% cdiff.descriptive.2$Short.ID){
      A056.col <- legend.color.picker("A056", id.list)
      
      plot <- plot +
        annotate("segment", x = 34, 
                 xend =  36, 
                 y = 4, yend = 5, colour = A056.col, 
                 size = 0.8, linetype = "dashed")
    } else if("A034" %in% cdiff.descriptive.2$Short.ID){
      A034.col <- legend.color.picker("A034", id.list)
      plot <- plot +
        annotate("segment", x = 53, ## day 53 is this persons transfer day
                 xend =  53, 
                 y = 1, yend = 2, colour = A034.col, 
                 size = 0.8, linetype = "dashed")
    } else if("A004" %in% cdiff.descriptive.2$Short.ID){
      A004.col <- legend.color.picker("A004", id.list)
      A013.col <- legend.color.picker("A013", id.list)
      plot <- plot + 
        annotate("segment", x = 11, 
                 xend =  11, 
                 y = 4, yend = 2, colour = A004.col, 
                 size = 0.8, linetype = "dashed") +
        annotate("segment", x = 20, 
                 xend =  20, 
                 y = 2, yend = 1, colour = A004.col, 
                 size = 0.8, linetype = "dashed")+ 
        annotate("segment", x = 8, 
                 xend =  8, 
                 y = 1, yend = 5, colour = A013.col, 
                 size = 0.8, linetype = "dashed")
      
    }
    
    return(plot)
  }
  
}

for(i in 1:length(num.clusts)){
  if(num.clusts[i] %in% c(5,6,7)){
    assign(paste("p.net.time",num.clusts[i], sep ="." ),
           sample.timing.plot.fun(num.clusts[i], dat, x.label = TRUE))
  } else{
    assign(paste("p.net.time",num.clusts[i], sep ="." ),
           sample.timing.plot.fun(num.clusts[i], dat, x.label = FALSE))
  }
  
}

### Clusters 
clst1 <- plot_grid(p.net.1, p.net.time.1, ncol = 2) #F
clst2 <- plot_grid(p.net.2, p.net.time.2, ncol = 2) #E
clst3 <- plot_grid(p.net.3, p.net.time.3, ncol = 2) #C
clst4 <- plot_grid(p.net.4, p.net.time.4, ncol = 2) #H
clst5 <- plot_grid(p.net.5, p.net.time.5, ncol = 2) #A
clst6 <- plot_grid(p.net.6, p.net.time.6, ncol = 2) #I
clst7 <- plot_grid(p.net.7, p.net.time.7, ncol = 2) #G
clst8 <- plot_grid(p.net.8, p.net.time.8, ncol = 2) #


space<- ggplot() + theme_void()

figure5.dummy <- cdiff.descriptive %>% filter(Short.ID %in% sort(unique(node.dat$short.ID)))
figure5.dummy$sample.env[3] <-  "Positive"
figure5.dummy$Study.Day[3] <-  0


fig.5.dat <- figure5.dummy %>% filter(C.diff == 1)
legend.colors <- legend.color.picker(sort(unique(fig.5.dat$Short.ID)), 
                                     sorted.id)

f5.legend <- ggplot() + 
  geom_point(data = figure5.dummy %>% filter(C.diff == 1), aes(x = Study.Day, y= Short.Room, 
                                                               color = Short.ID, shape = sample.env), size = 3) + 
  scale_y_discrete(drop=F) + # new line added here
  scale_color_manual(values = legend.colors, 
                     name = "Unique patient ID") + 
  guides(color=guide_legend(nrow=5,byrow=FALSE)) +
  ylab("Room Number") +
  xlab("Study Day") + 
  scale_shape_manual(values = shapes,
                     labels = c("Patient", "Room Env.",
                                "Shared Env.", "HCW Hand"),
                     breaks = c("Patient", "Room Env.",
                                "Shared Env.", "HCW Hand"),
                     name = "Sample Location",
                     guide = guide_legend(title.position = "top", order = 1)) +
  new_scale("shape") +
  geom_point(data = figure5.dummy %>% filter(C.diff == 1), aes(x = Study.Day, y= Short.Room, 
                                                               color = Short.ID, shape = sample.env), size = 3) + 
  scale_shape_manual(values = shapes,
                     labels = c("Positive"),
                     breaks = c("Positive"),
                     name = "ToxinA",
                     guide = guide_legend(title.position = "top", order = 2)) +  
  theme_cowplot() +
  background_grid() + 
  new_scale("color") +
  geom_line(data = figure5.dummy %>% 
              filter(Study.Day <8), aes(x = Study.Day, y= Short.Room, 
                                        color = as.factor(Study.Day))) + # new line added here
  scale_color_manual(values = c("#000000", "#252525", "#525252", "#737373", "#969696","#BDBDBD", "#D9D9D9", "#F0F0F0"),
                     #values = c("#000000FF", "#494848", "#909090", "#B4B4B4", "#D4D4D4","#ECECECFF"), 
                     name = "Genomic Distance (SNPs)") + 
  guides(color=guide_legend(nrow=3,byrow=FALSE)) +
  theme(legend.position = "bottom", 
        legend.direction="vertical")


legend_c <- get_legend(f5.legend)

#legends <- plot_grid(space,  legend_b, space, rel_widths = c(1.6,0.4,1,1), nrow = 1)

individual.net.plot <-plot_grid(clst1, clst2, clst3, clst6, 
                                clst7, clst8, clst4, clst5, 
                                ncol = 3, labels = "AUTO") 

individual.net.plot

legend_c <- plot_grid(space, legend_c, space, nrow= 1, rel_widths = c(0.5,1,0.5))

plot_grid(individual.net.plot, legend_c,  nrow = 2, rel_heights = c(1,0.22))

# legend2 <- ggplot(dat, aes(short.ID, cluster)) + 
#   geom_tile(aes(fill = dist)) + 
#   scale_fill_gradient(low = "#000000FF",
#                       high ="#ECECECFF",
#                       name = "Number of \nSNPs") +
#   theme(legend.position = "bottom")

legend_c <- get_legend(legend2)

