source("Scripts/Cdiff_GlobalFunctions.R")
source("Scripts/Cdiff_LoadData.R")
source("Scripts/Cdiff_PairwiseDistData.R")


# WITHIN AND BETWEEN SAMPLE DISTANCES --------------------------------------

table(table(c(unique(pairwise.dist.dat.one$Collection.Sample.ID1),unique(pairwise.dist.dat.one$Collection.Sample.ID2))))


# Define custom bin edges
bins <- c(-1, 1, 2, 5, 10, 99, 499, 749, 874, 999, 2499, 4999, 9999, Inf)

# Define custom bin labels
bin_labels <- c("0", "1-2", "3-5", "6-10", "11-99", 
                "100-499", "500-749", "750-874","875-999","1,000-2,499","2,500-4,999",
                "5,000-9,999", "10,000+")

pairwise.dist.dat.one$bin <- cut(pairwise.dist.dat.one$dist, breaks = bins, 
                                 labels = bin_labels)



pairwise.dist.dat.one <- pairwise.dist.dat.one %>% 
  mutate(pairwise.dist.ordered = factor(pairwise.dist.name, 
                                        levels = c("spss", "dphs","dpes","dpps", "dphh",
                                                   "dpeh", "dpee", "dpph", "dppe", "dppp", 
                                                   "sphh", "speh", "spee", "spph", "sppe", "sppp")))

p1 <- ggplot(pairwise.dist.dat.one) + 
  geom_bar(aes(x = bin, fill = pairwise.dist.ordered)) +
  scale_y_continuous() + 
  scale_fill_manual(name = "Pairwise Distance", 
                    breaks = c("sppp", "sppe", "spph", "spee", "speh", "sphh",
                               "dppp", "dppe", "dpph", "dpee","dpeh", "dphh",
                               "dpps", "dpes",
                               "dphs", "spss"), 
                    values = c("#b53737", "#c76b6a", "#d08585", 
                               "#d99f9f", "#e3b9b9", "#ecd3d2",
                               "#023858","#045a8d", "#0470b0", 
                               "#3690c0", "#74a9cf", "#a6bddb",
                               "#de7a00", "#ff9b21","#ffc47d", 
                               "#ffe1bd"),
                    labels = c("Patient-Patient pair", 
                               "Patient-Room env. pair",
                               "Patient-HCW hand pair", 
                               "Room env.-Room env. pair",
                               "HCW hand-Room env. pair",
                               "HCW hand-HCW hand pair",
                               "Patient-Patient pair",
                               "Patient-Room env. pair",
                               "Patient-HCW hand pair",
                               "Room env.-Room env. pair",
                               "HCW hand-Room env. pair",
                               "HCW hand-HCW hand pair",
                               "Shared env.-Patient pair", "Shared env.- Room env. pair",
                               "Shared env.-HCW hand pair","Shared env.-Shared env. pair")) +
  theme_cowplot() +
  background_grid() + 
  xlab("Genomic Distance (SNPs)") +
  ylab("Number of Pairs of Isolates") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

p1


# AVERAGE DISTANCES W/IN AND BETWEEN HOSPITALS ----------------------------

distances.by.loc.granular <- pairwise.dist.dat.one %>% 
  filter(Hospital1 == Hospital2) %>% 
  group_by(pairwise.dist.name, Hospital1) %>% 
  summarise(mean = mean(dist), 
            sd = sd(dist, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n),
            lower = mean - qt(1 - (0.05 / 2), n - 1) * se,
            upper = mean + qt(1 - (0.05 / 2), n - 1) * se) %>% 
  na.omit()

distances.by.loc.granular2 <- pairwise.dist.dat.one %>% 
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

distances.by.loc.granular.log <- pairwise.dist.dat.one %>% 
  filter(Hospital1 == Hospital2) %>% 
  group_by(pairwise.dist.name, Hospital1) %>% 
  summarise(mean = mean(dist),
            mean_log = mean(log(dist + 0.00001)), 
            sd_log = sd(log(dist + 0.00001), na.rm = TRUE),
            n = n(),
            se_log = sd_log / sqrt(n),
            lower_log = mean_log - qt(1 - (0.05 / 2), n - 1) * se_log,
            upper_log = mean_log + qt(1 - (0.05 / 2), n - 1) * se_log) %>% 
  na.omit()



distances.by.loc.cross.hosp <- pairwise.dist.dat.one %>% 
  filter(Hospital1 != Hospital2) %>% 
  group_by(pairwise.dist.name) %>% 
  summarise(mean = mean(dist), 
            sd = sd(dist, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n),
            lower = mean - qt(1 - (0.05 / 2), n - 1) * se,
            upper = mean + qt(1 - (0.05 / 2), n - 1) * se) %>% 
  na.omit()


distances.by.loc <- pairwise.dist.dat.one %>% 
  mutate(pairwise.broad = str_sub(pairwise.dist.name, start = -2)) %>% 
  filter(Hospital1 == Hospital2) %>% 
  group_by(pairwise.broad, Hospital1) %>% 
  summarise(mean = mean(dist), 
            sd = sd(dist, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n),
            lower = mean - qt(1 - (0.05 / 2), n - 1) * se,
            upper = mean + qt(1 - (0.05 / 2), n - 1) * se) %>% 
  na.omit() 


colnames(distances.by.loc)[1] <- colnames(distances.by.loc.granular)[1]

all.distances.by.loc <- rbind(distances.by.loc.granular)


U.pat.dist <- all.distances.by.loc %>% filter(pairwise.dist.name %in% c("dppp", "sppp"), 
                                              Hospital1 == "U")
V.pat.dist <- all.distances.by.loc %>% filter(pairwise.dist.name %in% c("dppp", "sppp"), 
                                              Hospital1 == "VA")
distances.by.loc.cross.hosp %>% filter(pairwise.dist.name %in% c("dppp", "sppp"))

wilcox.test(U.pat.dist$mean, distances.by.loc.cross.hosp$mean)
wilcox.test(V.pat.dist$mean, distances.by.loc.cross.hosp$mean)



colors <-c("#b53737", "#c76b6a", "#d08585", 
           "#d99f9f", "#e3b9b9", "#ecd3d2",
           "#023858","#045a8d", "#0470b0", 
           "#3690c0", "#74a9cf", "#a6bddb",
           "#de7a00", "#ff9b21","#ffc47d", 
           "#ffe1bd")

colors <- setNames(colors, c("sppp", "sppe", "spph", "spee", "speh", "sphh",
                             "dppp", "dppe", "dpph", "dpee", "dpeh", "dphh",
                             "dpps", "dpes", "dphs", "spss"))

labels <- c("sppp", "sppe", "spph", "spee", "speh", "sphh",
            "dppp", "dppe", "dpph", "dpee", "dpeh", "dphh",
            "dpps", "dpes", "dphs", "spss")
label.labels <- c("Patient-Patient pair", 
                  "Patient-Room env. pair",
                  "Patient-HCP hand pair", 
                  "Room env.-Room env. pair",
                  "HCP hand-Room env. pair",
                  "HCP hand-HCP hand pair",
                  "Patient-Patient pair",
                  "Patient-Room env. pair",
                  "Patient-HCW hand pair",
                  "Room env.-Room env. pair",
                  "HCP hand-Room env. pair",
                  "HCP hand-HCP hand pair",
                  "Shared env.-Patient pair", "Shared env.- Room env. pair",
                  "Shared env.-HCP hand pair","Shared env.-Shared env. pair")

labels <- setNames(labels, labels)

distances.by.loc.granular.log <- distances.by.loc.granular.log %>% 
  mutate(pairwise.dist.ordered = factor(pairwise.dist.name, 
                                        levels = c("spss", "dphs","dpes","dpps", "dphh",
                                                   "dpeh", "dpee", "dpph", "dppe", "dppp", 
                                                   "sphh", "speh", "spee", "spph", "sppe", "sppp")))


all.distances.plot.log <- ggplot(distances.by.loc.granular.log , aes(x = Hospital1, y = exp(mean_log)), 
                                 group_by = pairwise.dist.ordered) +
  geom_errorbar(aes(ymin = exp(lower_log), ymax = exp(upper_log), color = pairwise.dist.ordered),
                width = 0.05, size  = 0.5, position=position_dodge(width=0.4)) +
  geom_point(aes(color = pairwise.dist.ordered), shape = 15, size  = 4, 
             position=position_dodge(width=0.4)) + 
  scale_y_log10(breaks = c(0.00001, 1,  100, 1000), 
                labels = function(x) sprintf("%g", x)) +
  scale_color_manual(values = colors,
                     breaks = names(colors)[1:6], 
                     labels = c("Patient\u2013Patient pair",
                                "Patient\u2013Room env. pair",
                                "Patient\u2013HCW hand pair",
                                "Room env.\u2013Room env. pair",
                                "HCP hand\u2013Room env. pair",
                                "HCP hand\u2013HCP hand pair"),
                     name = "Within-Patient Stay",
                     guide = guide_legend(title.position = "top", 
                                          order = 1, override.aes = list(linetype = 0))) + 
  new_scale_color() +
  geom_errorbar(aes(ymin = exp(lower_log), ymax = exp(upper_log), color = pairwise.dist.ordered),
                width = 0.05, size  = 0.5, position=position_dodge(width=0.4)) +
  geom_point(aes(color = pairwise.dist.ordered), shape = 15, size  = 4, 
             position=position_dodge(width=0.4)) + 
  scale_color_manual(values = colors, 
                     labels = c("Shared env.\u2013Patient pair", 
                                "Shared env.\u2013Room env. pair",
                                "Shared env.\u2013HCP hand pair",
                                "Shared env.\u2013Shared env. pair"),
                     breaks = names(colors)[13:16], name = "Shared Env.\u2013Patient Stay",
                     guide = guide_legend(title.position = "top", 
                                          order = 0, override.aes = list(linetype = 0)))  + 
  new_scale_color() +
  geom_errorbar(aes(ymin = exp(lower_log), ymax = exp(upper_log), color = pairwise.dist.ordered),
                width = 0.05, size  = 0.5, position=position_dodge(width=0.4)) +
  geom_point(aes(color = pairwise.dist.ordered), shape = 15, size  = 4, 
             position=position_dodge(width=0.4)) + 
  scale_color_manual(values = colors, 
                     labels = c("Patient\u2013Patient pair",
                                "Patient\u2013Room env. pair",
                                "Patient\u2013HCW hand pair",
                                "Room env.\u2013Room env. pair",
                                "HCP hand\u2013Room env. pair",
                                "HCP hand\u2013HCP hand pair"),
                     breaks = names(colors)[7:12], name = "Between-Patient Stays",
                     guide = guide_legend(title.position = "top", 
                                          order = 2, override.aes = list(linetype = 0)))  +
  theme_cowplot() +
  # theme(aspect.ratio = 0.6) + 
  xlab("Hospital\n\n\n") + 
  ylab("Mean Genomic Distance (SNPs)") + 
  labs(colour="Pairwise Sample Location") +
  background_grid() + 
  scale_x_discrete(labels=c("U" = "Hospital A", "VA" = "Hospital B"))



all.distances.plot <- ggplot(all.distances.by.loc , aes(x = Hospital1, y = mean), 
                             group_by = pairwise.dist.name) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = pairwise.dist.name),
                width = 0.05, size  = 0.5, position=position_dodge(width=0.4)) +
  geom_point(aes(color = pairwise.dist.name), shape = 15, size  = 4, 
             position=position_dodge(width=0.4))  +
  scale_color_manual(values = colors,
                     breaks = c("sphh", "speh", "spee", "spph", "sppe","sppp"),
                     #names(colors)[1:6], 
                     labels = c("HCP hand\u2013HCP hand pair",
                                "HCP hand\u2013Room env. pair",
                                "Room env.\u2013Room env. pair",
                                "Patient\u2013HCP hand pair",
                                "Patient\u2013Room env. pair",
                                "Patient\u2013Patient pair"),
                     name = "Within-Patient Stay",
                     guide = guide_legend(title.position = "top", 
                                          order = 1, override.aes = list(linetype = 0))) + 
  new_scale_color() +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = pairwise.dist.name),
                width = 0.05, size  = 0.5, position=position_dodge(width=0.4)) +
  geom_point(aes(color = pairwise.dist.name), shape = 15, size  = 4, 
             position=position_dodge(width=0.4)) + 
  scale_color_manual(values = colors, 
                     labels = c("Shared env.\u2013Shared env. pair",
                                "Shared env.\u2013HCP hand pair",
                                "Shared env.\u2013Room env. pair",
                                "Shared env.\u2013Patient pair"
                     ),
                     breaks = names(colors)[16:13], name = "Shared Env.\u2013Patient Stay",
                     guide = guide_legend(title.position = "top", 
                                          order = 0, override.aes = list(linetype = 0)))  + 
  new_scale_color() +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = pairwise.dist.name),
                width = 0.05, size  = 0.5, position=position_dodge(width=0.4)) +
  geom_point(aes(color = pairwise.dist.name), shape = 15, size  = 4, 
             position=position_dodge(width=0.4)) + 
  scale_color_manual(values = colors, 
                     labels = c("HCP hand\u2013HCP hand pair",
                                "HCP hand\u2013Room env. pair",
                                "Room env.\u2013Room env. pair",
                                "Patient\u2013HCP hand pair",
                                "Patient\u2013Room env. pair",
                                "Patient\u2013Patient pair"
                     ),
                     breaks = names(colors)[12:7], name = "Between-Patient Stays",
                     guide = guide_legend(title.position = "top", 
                                          order = 2, override.aes = list(linetype = 0)))  +
  theme_cowplot() +
  # theme(aspect.ratio = 0.6) + 
  xlab("") + 
  ylab("Log Mean Genomic Distance (SNPs)") + 
  labs(colour="Pairwise Sample Location") +
  background_grid() + 
  scale_x_discrete(labels=c("U" = "Hospital A", "VA" = "Hospital B"))

p2 <- all.distances.plot.log +theme(legend.position = "none")

legend_pa <- cowplot::get_plot_component(all.distances.plot.log, "guide-box", return_all = TRUE)[[1]]


top.p <- plot_grid(p2, p1, nrow = 1, labels = c("A", "B"))
plot_grid(top.p, legend_pa, nrow = 1, rel_widths = c(1,0.4))
