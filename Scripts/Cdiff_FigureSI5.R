source("Scripts/Cdiff_GlobalFunctions.R")
source("Scripts/Cdiff_LoadData.R")

# WITHIN HOST GENOMIC DISTANCE --------------------------------------------

## Figure out which patients have more than 1 sample 
patient.counts<-dist.lower.meta.tall.all.dup %>% 
  count(PATID1, sort = T) 

## list of unique patients 
pat.id<-unique(dist.lower.meta.tall.all.dup$Short.ID1)
sorted.id <- c(sort(pat.id), "Reference")

## Function to make the dataset needed- within-patient SNP dists
distance.from.start <- function(dataset, id.num){
  temp <- NULL
  temp.start <- NULL
  id.list <- sort(unique(dataset$Short.ID1))
  ## Looking only at one patient, filter out everyone else
  temp<- dataset %>% 
    filter(Short.ID1 %in% id.num, Short.ID2 %in% id.num) 
  if(length(id.num) > 1){
    temp.start <- temp
  }  else{
    ## since we want to compare distances to the start, filter out all the "start"
    ## dates, aka collection.date1 to only be the date of the first sample
    start<- min(temp$Collection.date1)
    temp.start <- temp %>% 
      filter(Collection.date1 == start)
    if(length(unique(temp.start$Collection.Sample.ID1)) >1){
      if(temp.start$sample.loc1[1] == "HH"){
        start <- temp.start$Collection.Sample.ID1[which(temp.start$sample.loc1 != "HH")[1]]
      } else{
        start<- unique(temp.start$Collection.Sample.ID1)[1]
      }
      temp.start <- temp %>% 
        filter(Collection.Sample.ID1 == start)  
    }
  }
  return(temp.start)
}

### 6.0.1 Plot making function  --------------------------------------------
int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}


within.dist.plot <- function(dataset, id.num,  position.jitter = TRUE, 
                             logscale = FALSE, legend.pos = "none", 
                             xaxis =  element_text(size = 12), 
                             yaxis = element_text(size = 12)){
  names.plots <- NULL
  temp.start <- distance.from.start(dataset, id.num)
  legend.colors <- legend.color.picker(id.num, id.list)
  max.snp <- max(temp.start$dist) + 1
  if(position.jitter == TRUE){
    jitter <- position_jitter()
    
    tmp.plot<- ggplot(temp.start) + 
      geom_jitter(aes(x= Study.Day2, y = round(dist), 
                      shape = sample.env2, color = Short.ID1),
                  size = 4, width = 0.01, height = 0.1) + 
      scale_y_continuous(limits = c(0, max.snp), oob = scales::squish,
                         breaks = 
                           function(x) unique(floor(pretty(seq(0, (max(x) + 1) *
                                                                 1.1))))) + 
      xlab("Study Day") + 
      ylab("Genomic Distance (SNPs)") +
      scale_color_manual(breaks = id.num, 
                         values= legend.colors, 
                         name = "Unique ID") + 
      scale_shape_manual(values = c(19, 15, 18, 17), 
                         breaks = c("Patient", "Room Env.",
                                    "HCW Hand", "Shared Env."), 
                         labels = c("Patient", "Room Env.",
                                    "HCP Hand", "Shared Env."), 
                         name = "Sample\nLocation") +
      ggtitle(id.num) +
      theme_half_open() +
      theme(legend.position= legend.pos, #text = element_text(size = 20), 
            #axis.text.x = element_text(angle = 15, hjust = 1),
            axis.title.x = xaxis, axis.title.y = yaxis) +
      xlim(min(temp.start$Study.Day2 -1), max(temp.start$Study.Day2 +1)) +
      scale_x_continuous(breaks = int_breaks)+ 
      background_grid() 
  }else{
    if(logscale == TRUE){
      tmp.plot<- ggplot(temp.start) + 
        geom_point(aes(x= Study.Day2, y = round(dist), 
                       shape = sample.env2, color = Short.ID1),
                   size = 4) + 
        xlab("Study Day") + 
        ylab("Genomic Distance (SNPs)") +
        scale_color_manual(breaks = id.num, 
                           values= legend.colors, 
                           name = "Unique ID") + 
        scale_shape_manual(values = c(19, 15, 18, 17), 
                           breaks = c("Patient", "Room Env.",
                                      "HCW Hand", "Shared Env."), 
                           labels = c("Patient", "Room Env.",
                                      "HCP Hand", "Shared Env."), 
                           name = "Sample\nLocation") +
        ggtitle(id.num) +
        theme_half_open() +
        theme(legend.position= legend.pos, #text = element_text(size = 20), 
              #axis.text.x = element_text(angle = 15, hjust = 1),
              axis.title.x = xaxis, axis.title.y = yaxis) +
        background_grid() + 
        xlim(min(temp.start$Study.Day2 -1), max(temp.start$Study.Day2 +1)) +
        scale_x_continuous(breaks = int_breaks)+ 
        scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
    }else{
      tmp.plot<- ggplot(temp.start) + 
        geom_point(aes(x= Study.Day2, y = round(dist), 
                       shape = sample.env2, color = Short.ID1),
                   size = 4) + 
        xlab("Study Day") + 
        ylab("Genomic Distance (SNPs)") +
        scale_color_manual(breaks = id.num, 
                           values= legend.colors, 
                           name = "Unique ID") + 
        scale_shape_manual(values = c(19, 15, 18, 17), 
                           breaks = c("Patient", "Room Env.",
                                      "HCW Hand", "Shared Env."), 
                           labels = c("Patient", "Room Env.",
                                      "HCP Hand", "Shared Env."), 
                           name = "Sample\nLocation") +
        ggtitle(id.num) +
        theme_half_open() +
        theme(legend.position= legend.pos, #text = element_text(size = 20), 
              #axis.text.x = element_text(angle = 15, hjust = 1),
              axis.title.x = xaxis, axis.title.y = yaxis) +
        background_grid() + 
        xlim(min(temp.start$Study.Day2 -1), max(temp.start$Study.Day2 +1)) +
        scale_x_continuous(breaks = int_breaks)+ 
        scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) 
    }
  }
  # }
  return(tmp.plot)
}

#### Figure: Within host diversity figure ----------


## Make the plots
p.U002 <- within.dist.plot(dist.lower.meta.tall.all.dup, "A002", xaxis = element_blank())
p.U007 <- within.dist.plot(dist.lower.meta.tall.all.dup, "A007", 
                           xaxis = element_blank(), yaxis = element_blank())
p.U010 <- within.dist.plot(dist.lower.meta.tall.all.dup, "A010", 
                           xaxis = element_blank(), yaxis = element_blank())

p.U013 <- within.dist.plot(dist.lower.meta.tall.all.dup, "A013", 
                           xaxis = element_blank(), yaxis = element_blank())
p.U034 <- within.dist.plot(dist.lower.meta.tall.all.dup, "A034", 
                           xaxis = element_blank())
p.U046 <- within.dist.plot(dist.lower.meta.tall.all.dup, "A046", 
                           xaxis = element_blank(), position.jitter = FALSE, 
                           logscale = TRUE, yaxis = element_blank())

p.U048 <- within.dist.plot(dist.lower.meta.tall.all.dup, "A048", 
                           xaxis = element_blank(), yaxis = element_blank())
p.U238 <- within.dist.plot(dist.lower.meta.tall.all.dup, "A238", 
                           xaxis = element_blank(), yaxis = element_blank())
p.U250 <- within.dist.plot(dist.lower.meta.tall.all.dup, "A250", 
                           xaxis = element_blank())
p.U912 <- within.dist.plot(dist.lower.meta.tall.all.dup, "A912", position.jitter = FALSE,
                           yaxis = element_blank(), xaxis = element_blank())

p.V003 <- within.dist.plot(dist.lower.meta.tall.all.dup, "B003",
                           yaxis = element_blank(), xaxis = element_blank())
p.V008 <- within.dist.plot(dist.lower.meta.tall.all.dup, "B008", 
                           xaxis = element_blank(), yaxis = element_blank())
p.V012 <- within.dist.plot(dist.lower.meta.tall.all.dup, "B012")

p.V018 <- within.dist.plot(dist.lower.meta.tall.all.dup, "B018", 
                           yaxis = element_blank())
p.V019 <- within.dist.plot(dist.lower.meta.tall.all.dup, "B019", 
                           yaxis = element_blank())
p.V058 <- within.dist.plot(dist.lower.meta.tall.all.dup, "B058", yaxis = element_blank())
# p.V917 <- within.dist.plot(dist.lower.meta.tall.all.dup, "B917", yaxis = element_blank())


## Make the plot with just the legend to extract
plot.ids <-c("A002", "A007", "A010", "A013", "A034", "A046", "A048", "A238", "A250", "A912",
             "B003", "B008", "B012", "B018", "B019", "B058") #, "B917"
p.legend <- within.dist.plot(dist.lower.meta.tall.all.dup, 
                             plot.ids, legend.pos = "bottom")
legend_b <- get_legend(p.legend)


## Put it all together into one plot 
blank.row <-ggplot() + theme_void()
row.2 <- plot_grid(p.U002, p.U007, p.U010, ncol = 3, label_size = 12)
row.3 <- plot_grid(p.U013,  p.U034, p.U046, ncol = 3, label_size = 12)
row.4 <- plot_grid(p.U048, p.U238,  p.U250, ncol = 3, label_size = 12) ## FIX A222
row.6 <- plot_grid(p.U912, p.V003, p.V008, ncol = 3, label_size = 12)
row.7 <- plot_grid(p.V012, p.V018,  p.V019, ncol = 3, label_size = 12)
row.8 <- plot_grid(p.V058, legend_b,  ncol = 3, label_size = 12) #p.V917,

within.host.plot.A <- plot_grid(blank.row, row.2, row.3, row.4,  
                                rel_heights = c(1,5,5,5), nrow = 4) + 
  draw_text("Hospital A", x = 0.1, y = 0.98, size = 20)

within.host.plot.B <- plot_grid(blank.row, row.6, row.7, row.8,
                                rel_heights = c(1,5,5,5), nrow = 4) + 
  draw_text("Hospital B", x = 0.1, y = 0.98, size = 20)

p.within <- plot_grid(within.host.plot.A, within.host.plot.B)


## Make the plot with just the legend to extract
plot.ids <-c("A002", "A007", "A010", "A013", "A034", "A046", "A048", "A238", "A250", "A912",
             "B003", "B008", "B012", "B018", "B019", "B058") #, "B917"
p.legend <- within.dist.plot(dist.lower.meta.tall.all.dup, 
                             plot.ids, legend.pos = "right")
legend_R <- get_legend(p.legend)


row.1 <- plot_grid(p.U002, p.U007, p.U010, p.U013, ncol = 4, label_size = 12)
row.2 <- plot_grid(p.U034, p.U046, p.U048, p.U238, ncol = 4, label_size = 12)
row.3 <- plot_grid(p.U250, p.U912, p.V003, p.V008, ncol = 4, label_size = 12) ## FIX A222
row.4 <- plot_grid(p.V012, p.V018,  p.V019,p.V058, ncol = 4, label_size = 12)

p.within <- plot_grid(row.1, row.2, row.3, row.4, nrow = 4)

p.within <- plot_grid(p.within, legend_R, rel_widths = c(1,0.15))

p.within