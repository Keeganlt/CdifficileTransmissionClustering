source("code/GlobalFunctions.R")
source("code/LoadData.R")
source("code/PairwiseDistData.R")

### **** MAKE A FUNCTION FOR WITHIN BETWEEN BODY/ENV SITES DISTANCES  --------

pairwise.dist.dat.one <- pairwise.dist.dat.one %>% 
  mutate("granular.sample.site1" = sample.loc1) %>% 
  mutate("granular.sample.site2" = sample.loc2) %>% 
  mutate("pt.status1" = NA) %>% 
  mutate("pt.status2" = NA)

Pt.ID.list <- cdiff.descriptive3 %>% filter(sample.env == "Patient", C.diff == 1) %>% distinct(PATID)

for(i in 1:nrow(pairwise.dist.dat.one)){
  if(pairwise.dist.dat.one$sample.env1[i] == "Patient"){
    pairwise.dist.dat.one$granular.sample.site1[i] <- word(pairwise.dist.dat.one$Sample.ID1[i],sep = "-", start = 3)
    pairwise.dist.dat.one$granular.sample.site1[i] <-str_extract(pairwise.dist.dat.one$granular.sample.site1[i], 
                                                                 "[:alpha:]")    
  }
  if(pairwise.dist.dat.one$sample.env2[i] == "Patient"){
    pairwise.dist.dat.one$granular.sample.site2[i] <- word(pairwise.dist.dat.one$Sample.ID2[i],sep = "-", start = 3)
    pairwise.dist.dat.one$granular.sample.site2[i] <-str_extract(pairwise.dist.dat.one$granular.sample.site2[i], 
                                                                 "[:alpha:]")    
  }
  if(pairwise.dist.dat.one$Short.ID1[i] %in% Pt.ID.list$Short.ID){
    pairwise.dist.dat.one$pt.status1[i] <- "pt.positive"    
  } else {
    pairwise.dist.dat.one$pt.status1[i] <- "pt.not.positive"
  }
  if(pairwise.dist.dat.one$Short.ID2[i] %in% Pt.ID.list$Short.ID){
    pairwise.dist.dat.one$pt.status2[i] <- "pt.positive"    
  } else {
    pairwise.dist.dat.one$pt.status2[i] <- "pt.not.positive"
  }
}




pairwise.dist.dat.one <- tidyr::unite(pairwise.dist.dat.one, granular.pair, 
                                      granular.sample.site1, granular.sample.site2, 
                                      na.rm = TRUE, remove = FALSE)



distances.pt.granular <- pairwise.dist.dat.one %>% 
  filter(Short.ID1 == Short.ID2, 
         sample.env1 %!in% c("HCW Hand","Shared Env."), 
         sample.env2 %!in% c("HCW Hand","Shared Env.")) %>% 
  group_by(granular.pair) %>% 
  summarise(mean = mean(dist), 
            median = median(dist), 
            sd = sd(dist, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n),
            lower = mean - qt(1 - (0.05 / 2), n - 1) * se,
            upper = mean + qt(1 - (0.05 / 2), n - 1) * se, 
            Q1 = quantile(dist, 0.25), 
            Q3 = quantile(dist, 0.75)
  ) %>% 
  na.omit()



distances.by.loc.all <- pairwise.dist.dat.one %>% 
  group_by(pairwise.dist.name) %>% 
  summarise(mean = mean(dist), 
            median = median(dist), 
            sd = sd(dist, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n),
            lower = mean - qt(1 - (0.05 / 2), n - 1) * se,
            upper = mean + qt(1 - (0.05 / 2), n - 1) * se, 
            Q1 = quantile(dist, 0.25), 
            Q3 = quantile(dist, 0.75)) %>% 
  na.omit()


contact_labels <- c("TO", "BR", "CC", "TB", "DH")

contact.dist.dat <- pairwise.dist.dat.one %>%
  filter(Short.ID1 == Short.ID2) %>% 
  group_by(Short.ID1) %>%
  mutate(contact_precautions = if_else(any(sample.loc1[sample.env1 == "Room Env."] 
                                           %in% contact_labels), 1, 0)) %>%
  ungroup()

distances.by.contact.prec <- contact.dist.dat %>% 
  group_by(pairwise.dist.name, contact_precautions) %>% 
  summarise(mean = mean(dist), 
            median = median(dist), 
            Q1 = quantile(dist, 0.25), 
            Q3 = quantile(dist, 0.75)) %>% 
  na.omit()


distances.by.loc.pt.status <- pairwise.dist.dat.one %>% 
  filter(pt.status1 == pt.status2, 
         sample.env1 == "Room Env.", 
         sample.env2 == "Room Env.") %>% 
  group_by(pt.status2) %>% 
  summarise(mean = mean(dist), 
            median = median(dist),
            sd = sd(dist, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n),
            lower = mean - qt(1 - (0.05 / 2), n - 1) * se,
            upper = mean + qt(1 - (0.05 / 2), n - 1) * se) %>% 
  na.omit()




all.granular.plot <- ggplot(distances.pt.granular , aes(x = granular.pair, y = mean), 
                            group_by = granular.pair) +
  # geom_errorbar(aes(ymin = lower, ymax = upper, color = pairwise.dist.name),
  #               width = 0.05, size  = 0.5, position=position_dodge(width=0.4)) +
  geom_point(aes(color = granular.pair), shape = 15, size  = 4, 
             position=position_dodge(width=0.4))  +
  scale_y_log10() 


#

distances.pt.granular <- distances.pt.granular %>% 
  mutate("granular.site1" = NA) %>% 
  mutate("granular.site2" = NA)

for(i in 1:nrow(distances.pt.granular)){
  distances.pt.granular$granular.site1[i] <- word(distances.pt.granular$granular.pair[i],sep = "_", start = 1)
  distances.pt.granular$granular.site2[i] <- word(distances.pt.granular$granular.pair[i],sep = "_", start = 2)
}

distances.by.loc.ppee <- distances.by.loc.all %>% filter(pairwise.dist.name %in% c("sppp", "sppe", "spee", "spph", "speh", "sphh"))

distances.by.loc.ppee <- distances.by.loc.ppee %>% 
  mutate("granular.site1" = NA) %>% 
  mutate("granular.site2" = NA)


distances.by.loc.ppee <- rbind(distances.by.loc.ppee, 
                               distances.by.loc.ppee[which(distances.by.loc.ppee$pairwise.dist.name =="sppe"), ],
                               distances.by.loc.ppee[which(distances.by.loc.ppee$pairwise.dist.name =="spph"), ],
                               distances.by.loc.ppee[which(distances.by.loc.ppee$pairwise.dist.name =="speh"), ])
distances.by.loc.ppee$granular.site1[1] <- "Env."
distances.by.loc.ppee$granular.site2[1] <- "Env."
distances.by.loc.ppee$granular.site1[2] <- "Env."
distances.by.loc.ppee$granular.site2[2] <- "HCP"
distances.by.loc.ppee$granular.site1[3] <- "HCP"
distances.by.loc.ppee$granular.site2[3] <- "HCP"
distances.by.loc.ppee$granular.site1[4] <- "Pat."
distances.by.loc.ppee$granular.site2[4] <- "Env."
distances.by.loc.ppee$granular.site1[5] <- "Pat."
distances.by.loc.ppee$granular.site2[5] <- "HCP"
distances.by.loc.ppee$granular.site1[6] <- "Pat."
distances.by.loc.ppee$granular.site2[6] <- "Pat."
distances.by.loc.ppee$granular.site1[7] <- "Env."
distances.by.loc.ppee$granular.site2[7] <- "Pat."
distances.by.loc.ppee$granular.site1[8] <- "HCP"
distances.by.loc.ppee$granular.site2[8] <- "Pat."
distances.by.loc.ppee$granular.site1[9] <- "HCP"
distances.by.loc.ppee$granular.site2[9] <- "Env."



colnames(distances.by.loc.ppee)[1] <- colnames(distances.pt.granular)[1]
distances.pt.granular <- rbind(distances.pt.granular, distances.by.loc.ppee)

distances.pt.granular$granular.site1 <- factor(distances.pt.granular$granular.site1, 
                                               levels = c("Pat.","P", "G", "A", 
                                                          "Z1", "BR", "TB", "DH", 
                                                          "Z2", "CC", "Z3", "TO", "Env."))

distances.pt.granular$granular.site2 <- factor(distances.pt.granular$granular.site2, 
                                               levels = c("Pat.","P", "G", "A", 
                                                          "Z1", "BR", "TB", "DH", 
                                                          "Z2", "CC", "Z3", "TO", "Env."))



pairwise.dist.dat.one <- pairwise.dist.dat.one %>% 
  filter(Short.ID1 == Short.ID2, 
         sample.env1 %!in% c("Shared Env."), 
         sample.env2 %!in% c("Shared Env.")) %>% 
  group_by(granular.pair) %>% 
  summarise(mean = mean(dist), 
            median = median(dist), 
            sd = sd(dist, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n),
            lower = mean - qt(1 - (0.05 / 2), n - 1) * se,
            upper = mean + qt(1 - (0.05 / 2), n - 1) * se
  ) %>% 
  na.omit()


distances.pt.granular <- distances.pt.granular %>% 
  mutate("granular.site1" = NA) %>% 
  mutate("granular.site2" = NA)

for(i in 1:nrow(distances.pt.granular)){
  distances.pt.granular$granular.site1[i] <- word(distances.pt.granular$granular.pair[i],sep = "_", start = 1)
  distances.pt.granular$granular.site2[i] <- word(distances.pt.granular$granular.pair[i],sep = "_", start = 2)
}

distances.by.loc.ppee <- distances.by.loc.all %>% filter(pairwise.dist.name %in% c("sppp", "sppe", "spee"))

distances.by.loc.ppee <- distances.by.loc.ppee %>% 
  mutate("granular.site1" = NA) %>% 
  mutate("granular.site2" = NA)

distances.by.loc.ppee <- distances.by.loc.all %>% filter(pairwise.dist.name %in% c("sppp", "sppe", "spee", "sphh", "speh", "spph"))

distances.by.loc.ppee <- distances.by.loc.ppee %>% 
  mutate("granular.site1" = NA) %>% 
  mutate("granular.site2" = NA)


distances.by.loc.ppee <- rbind(distances.by.loc.ppee, 
                               distances.by.loc.ppee[which(distances.by.loc.ppee$pairwise.dist.name =="sppe"), ], 
                               distances.by.loc.ppee[which(distances.by.loc.ppee$pairwise.dist.name =="spph"), ], 
                               distances.by.loc.ppee[which(distances.by.loc.ppee$pairwise.dist.name =="speh"), ])
distances.by.loc.ppee$granular.site1[1] <- "Env."
distances.by.loc.ppee$granular.site2[1] <- "Env."
distances.by.loc.ppee$granular.site1[2] <- "Env."
distances.by.loc.ppee$granular.site2[2] <- "HCP"
distances.by.loc.ppee$granular.site1[3] <- "HCP"
distances.by.loc.ppee$granular.site2[3] <- "HCP"
distances.by.loc.ppee$granular.site1[4] <- "Pat."
distances.by.loc.ppee$granular.site2[4] <- "Env."
distances.by.loc.ppee$granular.site1[5] <- "Pat."
distances.by.loc.ppee$granular.site2[5] <- "HCP"
distances.by.loc.ppee$granular.site1[6] <- "Pat."
distances.by.loc.ppee$granular.site2[6] <- "Pat."
distances.by.loc.ppee$granular.site1[7] <- "Env."
distances.by.loc.ppee$granular.site2[7] <- "Pat."
distances.by.loc.ppee$granular.site1[8] <- "HCP"
distances.by.loc.ppee$granular.site2[8] <- "Pat."
distances.by.loc.ppee$granular.site1[9] <- "HCP"
distances.by.loc.ppee$granular.site2[9] <- "Env."




colnames(distances.by.loc.ppee)[1] <- colnames(distances.pt.granular)[1]
distances.pt.granular <- rbind(distances.pt.granular, distances.by.loc.ppee)

distances.pt.granular$granular.site1 <- factor(distances.pt.granular$granular.site1, 
                                               levels = c("Pat.","P", "G", "A", 
                                                          "Z1", "BR", "TB", "DH", 
                                                          "Z2", "CC", "Z3", "TO", "Env.", "HCP", "HH"))

distances.pt.granular$granular.site2 <- factor(distances.pt.granular$granular.site2, 
                                               levels = c("Pat.","P", "G", "A", 
                                                          "Z1", "BR", "TB", "DH", 
                                                          "Z2", "CC", "Z3", "TO", "Env.", "HCP", "HH"))




colours <- c("#ffffd9",
             "#edf8b1",
             "#c7e9b4",
             "#7fcdbb",
             "#41b6c4",
             "#1d91c0",
             "#225ea8",
             "#253494",
             "#081d58")

colour_breaks <- c(0,0.5, 1 ,2,5,10, 110)

colours <- c("#3f9cb3",
             "#79b7c3",
             "#97bc9b",
             "#bec365",
             "#e4ca33",
             "#e7c019",
             "#e29e00",
             "#e58601",
             "#f12300")


p1 <- ggplot(data = distances.pt.granular %>% filter(granular.site2 %in% c("Pat.", "Env.", "HCP"))) + 
  geom_tile(aes(x = granular.site1, y = granular.site2, fill = (mean))) + 
  theme_cowplot() + 
  xlab("Sample Location 1") + 
  ylab("Sample Location 2") +
  scale_fill_gradientn(
    limits  = c(0,110),
    colours = colours[c(1, seq_along(colours), length(colours))],
    values  = c(0, scales::rescale(colour_breaks, from = c(0,100)), 1),
    name = "Mean Distance (SNPs)")  +
  theme(legend.position = 'none') #+ 
# ggtitle("All Patients")

p2<- ggplot(data = distances.pt.granular %>% filter(granular.site1 %in% c("P", "G", "A", "BR", "TB", "DH", "CC", "TO", "HH") &
                                                      granular.site2 %in% c("P", "G", "A", "BR", "TB", "DH", "CC", "TO", "HH"))) + 
  geom_tile(aes(x = granular.site1, y = granular.site2, fill = (mean))) + 
  theme_cowplot() + 
  xlab("Sample Location 1") + 
  ylab("Sample Location 2") +
  scale_fill_gradientn(
    limits  = c(0,110),
    colours = colours[c(1, seq_along(colours), length(colours))],
    values  = c(0, scales::rescale(colour_breaks, from = c(0,100)), 1),
    name = "Mean Distance (SNPs)")  +
  theme(legend.position = 'none') #+ 
# ggtitle("Contact Precautions Patients")


p.3<-ggplot(data = distances.pt.granular %>% filter(granular.site1 %in% c("P", "G", "A", "Z1", "Z2", "Z3", "HH") &
                                                      granular.site2 %in% c("P", "G", "A", "Z1", "Z2", "Z3", "HH"))) + 
  geom_tile(aes(x = granular.site1, y = granular.site2, fill = (mean))) + 
  theme_cowplot() + 
  xlab("Sample Location 1") + 
  ylab("Sample Location 2") +
  scale_fill_gradientn(
    limits  = c(0,110),
    colours = colours[c(1, seq_along(colours), length(colours))],
    values  = c(0, scales::rescale(colour_breaks, from = c(0,100)), 1),
    name = "Mean\nDistance\n(SNPs)")  
#ggtitle("Non-Contact Precautions Patients") 


p3 <- p.3 + 
  theme(legend.position = "none")

legend.p <- get_legend(p.3)

plot_grid(p1, p2, p3, legend.p, nrow = 1, rel_widths = c(1,1,1,0.4), labels = c("A", "B", "C", ""))
