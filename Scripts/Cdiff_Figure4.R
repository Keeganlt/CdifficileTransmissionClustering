
## 0-2 SNPs cluster -------------------------------------------------------
source("Scripts/Cdiff_ClusteringData.R")

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
  
  dat <- dat %>% 
    mutate("sample.loc" = from) %>%
    mutate_at("sample.loc", str_extract, "PAT|Z1|Z2|Z3|HH|TO|BR|CC|TB|DH|SN|SS|SE") %>%
    mutate(sample.env = case_when(
      sample.loc == "PAT" ~ "Patient",
      sample.loc %in% c("Z1", "Z2", "Z3", "TB", "TO", "DH") ~ "Room Env.",
      sample.loc == "HH" ~ "HCP Hand",
      TRUE ~ NA_character_  # If none of the above conditions are met, set to NA
    )) 
  
  from_hospital <- sub("-.*", "", dat$from)
  to_hospital <- sub("-.*", "", dat$to)
  
  # Creating the hospital variable based on the comparison
  dat$hospital <- ifelse(from_hospital == to_hospital, "same", "different")
  
  
  edge.drop <- dat %>%
    group_by(short.ID, sample.env) %>%
    slice(1)  
  
  
  edge.dat <- dat %>% 
    filter(cluster == cluster.num) %>% 
    select(from, to, dist)
  
  names_to_drop <- edge.dat %>%  
    filter(from %!in% edge.drop$from) %>% 
    select(from)
  
  names_to_drop <- unlist(names_to_drop)
  
  
  node.dat <- dat %>% 
    filter(cluster == cluster.num) %>% 
    select(from, short.ID, hospital) %>% 
    distinct(from, short.ID, hospital) %>% 
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
  
  graph_df <- as.data.frame(graph.dat, what = "vertices")
  
  # Identify the nodes you want to drop based on their attributes (e.g., sample.loc)
  vertices_to_drop <- which(graph_df$from %in% names_to_drop)
  
  # Now, delete the vertices
  graph.dat <- delete_vertices(graph.dat, vertices_to_drop)
  
  if(cluster.num == 2){
    subgraph <- delete_vertices(graph.dat, V(graph.dat)[!V(graph.dat) %in% c(1,3, 4,6,7,8)])
    graph.dat <- subgraph
  }
  
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
    geom_edge_link(aes(
      # color = factor(hospital), 
      linetype = factor(dist))) +
    geom_node_point(aes(color = short.ID, shape = as.factor(sample.loc)), 
                    size = 5, show.legend = F) +
    scale_edge_linetype_manual(name = "Distance",
                               values = c("0" = "solid", "1" = "dashed", "2" = "dotted")) + 
    # scale_edge_colour_manual(values =c("#000000FF",  "#525252", "#BDBDBD"),
    #                          name = "Number of SNPs", guide = "none") +
    scale_color_manual(breaks = factor(unique(node.dat$short.ID)), 
                       values = legend.colors, 
                       name = "Unique ID") +
    scale_shape_manual(values = c(19, 15, 18, 17), 
                       breaks = c("Patient", "Room Env.",
                                  "HCW Hand", "Shared Env."), 
                       name = "Sampling Environment") +
    coord_equal() + 
    theme_void() + 
    theme(legend.position = "none",
          panel.background = element_rect(fill = "transparent", color = NA)
    )
  
  return(plot)
}


blank.space<- ggplot() + theme_void() + 
  theme(panel.background = element_rect(fill = "transparent", color = NA))



## UNCOMMENT IF YOU WANT TO CHANGE THE ORIENTATION OF A CLUSTER OTHERWISE LOAD THE RDATA BELOW
# num.clusts <- seq(from = 1, to= counter-1, by = 1)

# num.clusts <- c(4)
# 
# for(i in 1:length(num.clusts)){
#   assign(paste("p.net.short",num.clusts[i], sep ="." ),
#          network.plot.fun(num.clusts[i], dat, lay.out = "igraph"))
# }


load("Data/ClusterA.RData")
load("Data/ClusterB.RData")
load("Data/ClusterC.RData")
load("Data/ClusterD.RData")
load("Data/ClusterE.RData")
load("Data/ClusterF.RData")
load("Data/ClusterG.RData")


p1<- plot_grid(blank.space, p.net.short.1,blank.space, p.net.short.2, blank.space, nrow = 1, rel_widths = c(0.2,1,0.1,1,0.2))

p2<- plot_grid(blank.space, p.net.short.3, p.net.short.5, blank.space, nrow = 1, rel_widths = c(0.1,1,0.5,0.3))

p3<- plot_grid(blank.space, p.net.short.4, blank.space, p.net.short.6, blank.space, nrow = 1, rel_widths = c(0.2,0.2,0.05,0.1,0))

p4<- plot_grid(blank.space, p.net.short.7, blank.space, nrow = 1, rel_widths = c(0.1,0.58,0.6))

plot1 <- plot_grid(p1, p2, rel_widths = c(1, 0.8), nrow =1 )
plot2 <- plot_grid(p3, p4, rel_widths = c(2,1.5), nrow =1 )

full_plot <- plot_grid(blank.space, plot1,blank.space, plot2, blank.space, legend, nrow = 6, rel_heights = c(0.1,0.6,0.05,0.5,0.1,0.2))

full_plot
draw_obj <- ggdraw() +
  draw_plot(full_plot)

draw_obj+
  geom_ellipse(aes(x0 = 0.15, y0 = 0.74, a = 0.09, b = 0.2, angle = 110, m1 = 5), color = "#02401b", fill = NA, size =2) +
  draw_text("ST 15", x = 0.04, y = 0.91, label = "ST 15", size = 20, color = "black")  +
  draw_text("ST 15", x = 0.155, y = 0.89, label = "Cluster A", size = 18, color = "#3f3f3f", angle = 0) + # Add text "Center" at coordinates (5, 3)
  geom_ellipse(aes(x0 = 0.42, y0 = 0.765, a = 0.12, b = 0.14, angle = 10, m2 = 3), color = "#972d14", fill = NA, size =2) +
  draw_text("ST 12", x = 0.295, y = 0.87, label = "ST 12", size = 20, color = "black") +
  draw_text("ST 12", x = 0.435, y = 0.83, label = "Cluster B", size = 18, color = "#3f3f3f", angle = 0) +
  geom_ellipse(aes(x0 = 0.775, y0 = 0.73, a = 0.18, b = 0.22, angle = 0, m1 = 3), color = "#EBDA96", fill = NA, size =2) +
  draw_text("ST 3", x = 0.61, y = 0.91, label = "ST 3", size = 20, color = "black")  +
  draw_text("ST 3", x = 0.64, y = 0.84, label = "Cluster C", size = 18, color = "#3f3f3f", angle = 0) + # Add text "Center" at coordinates (5, 3)
  draw_text("ST 3", x = 0.88, y = 0.84, label = "Cluster D", size = 18, color = "#3f3f3f", angle = 0) +
  geom_ellipse(aes(x0 = 0.32, y0 = 0.355, a = 0.07, b = 0.16, angle = pi/4.5, m2 = 3), color = "#a2a475", fill = NA, size =2)  +
  draw_text("ST 39", x = 0.179, y = 0.45, label = "ST 39", size = 20, color = "black")  +
  draw_text("ST 39", x = 0.27, y = 0.45, label = "Cluster E", size = 18, color = "#3f3f3f", angle = 0) +
  geom_ellipse(aes(x0 = 0.53, y0 = 0.35, a = 0.05, b = 0.17, angle = pi/10, m2 = 3), color = "#c7b19c", fill = NA, size =2) +
  draw_text("ST 110", x = 0.43, y = 0.48, label = "ST 110", size = 20, color = "black")  +
  draw_text("ST 110", x = 0.525, y = 0.43, label = "Cluster F", size = 18, color = "#3f3f3f", angle = 0) + # Add text "Center" at coordinates (5, 3)
  geom_ellipse(aes(x0 = 0.7, y0 = 0.355, a = 0.08, b = 0.18, angle = pi/7, m2 = 3), color = "#7fa6ad", fill = NA, size =2) +
  draw_text("ST 125", x = 0.569, y = 0.49, label = "ST 125", size = 20, color = "black") +
  draw_text("ST 125", x = 0.68, y = 0.38, label = "Cluster G", size = 18, color = "#3f3f3f", angle = 0)



