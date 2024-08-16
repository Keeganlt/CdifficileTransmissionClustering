## Everyone
multi.pts <- c("A002", "A007", "A010", "A013", "A034", "A046", "A048", "A238",
               "A250", "A912", "B003", "B008", "B012", "B018", "B019", "B058")

dat <- distance.from.start(dist.lower.meta.tall.all.dup, multi.pts[1])
for(i in 2:length(multi.pts)){
  tmp<- distance.from.start(dist.lower.meta.tall.all.dup, multi.pts[i])
  dat <- rbind(dat, tmp)
}


tmp<-distance.from.start(dist.lower.meta.tall.all.dup, "A046")

max.snp <- max(dat$dist) + 1


legend.colors.dat <- legend.color.picker(c(multi.pts), 
                                      sorted.id)

p1 <- ggplot(dat) + 
  geom_jitter(aes(x=Patient.Day2, y = round(dist) + 0.0000001, 
                  shape = sample.env2, color = Short.ID1),
              size = 4, width = 0.01, height = 0.1) + 
  theme_cowplot() +
  scale_y_log10() 
  scale_color_manual(values = legend.colors.dat, 
                     name = "Unique patient ID") + 
  xlab("Days since Admission") +
  ylab("Genomic Distance (SNPs)") + 
  scale_shape_manual(values = c(19, 15, 18, 17), 
                     breaks = c("Patient", "Room Env.",
                                "HCW Hand", "Shared Env."), 
                     name = "Sample\nLocation")
  
  
## PATIENT 

multi.pts <- c("A002", "A007", "A010", "A013", "A034", "A046", "A048", "A238",
               "A250", "A912", "B003", "B008", "B012", "B018", "B019", "B058")

dat <- distance.from.start(dist.lower.meta.tall.all.dup %>% 
                             filter(sample.env1 == "Patient" & sample.env2 == "Patient"), multi.pts[1])
for(i in 2:length(multi.pts)){
  tmp<- distance.from.start(dist.lower.meta.tall.all.dup %>% 
                              filter(sample.env1 == "Patient" & sample.env2 == "Patient"), multi.pts[i])
  dat <- rbind(dat, tmp)
}


legend.colors.dat <- legend.color.picker(sort(unique(dat$Short.ID1)), 
                                         sorted.id)

max.snp <- max(dat$dist) + 1

p2 <- ggplot(dat) + 
  geom_jitter(aes(x=Patient.Day2, y = round(dist) , 
                   color = Short.ID1), shape = 19,
              size = 3, width = 0.01, height = 0.1) +
  theme_cowplot() +
  scale_y_continuous(limits = c(0, max.snp), oob = scales::squish,
      breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) + 
  scale_color_manual(values = legend.colors.dat, 
                     name = "Unique patient ID") + 
  theme(legend.position = "none") +
  xlab("Days since Admission") +
  ylab("Genomic Distance (SNPs)")

## ENVIRONMENT 

multi.pts <- c("A002", "A007", "A010", "A013", "A034", "A046", "A048", "A238",
               "A250", "A912", "B003", "B008", "B012", "B018", "B019", "B058")

dat <- distance.from.start(dist.lower.meta.tall.all.dup %>% 
                             filter(sample.env1 == "Room Env." & sample.env2 == "Room Env."), multi.pts[1])
for(i in 2:length(multi.pts)){
  tmp<- distance.from.start(dist.lower.meta.tall.all.dup %>% 
                              filter(sample.env1 == "Room Env." & sample.env2 == "Room Env."), multi.pts[i])
  dat <- rbind(dat, tmp)
}


legend.colors.dat <- legend.color.picker(sort(unique(dat$Short.ID1)), 
                                         sorted.id)
max.snp <- max(dat$dist) + 1

p3 <- ggplot(dat) + 
  geom_jitter(aes(x=Patient.Day2, y = round(dist), 
                   color = Short.ID1), shape = 15,
              size = 3, width = 0.01, height = 0.1) + 
  theme_cowplot() +
  scale_y_continuous(limits = c(0, max.snp), oob = scales::squish,
       breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
  scale_color_manual(values = legend.colors.dat, 
                     name = "Unique patient ID") + 
  theme(legend.position = "none") +
  xlab("Days since Admission") +
  ylab("Genomic Distance (SNPs)")


## HCW hands 

multi.pts <- c("A002", "A007", "A010", "A013", "A034", "A046", "A048", "A238",
               "A250", "A912", "B003", "B008", "B012", "B018", "B019", "B058")

dat <- distance.from.start(dist.lower.meta.tall.all.dup %>% 
                             filter(sample.env1 == "HCW Hand", sample.env2 == "HCW Hand"), multi.pts[1])
for(i in 2:length(multi.pts)){
  tmp<- distance.from.start(dist.lower.meta.tall.all.dup %>% 
                              filter(sample.env1 == "HCW Hand", sample.env2 == "HCW Hand"), multi.pts[i])
  dat <- rbind(dat, tmp)
}

max.snp <- max(dat$dist) + 1

legend.colors.dat <- legend.color.picker(sort(unique(dat$Short.ID1)), 
                                         sorted.id)

p4 <- ggplot(dat) + 
  geom_jitter(aes(x=Patient.Day2, y = round(dist), 
                  color = Short.ID1), shape = 18,
              size = 3, width = 0.01, height = 0.1) +
  theme_cowplot()+ 
  scale_y_continuous(limits = c(0, max.snp), oob = scales::squish,
       breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
  scale_color_manual(values = legend.colors.dat, 
                     name = "Unique patient ID") + 
  theme(legend.position= "none") +
  xlab("Days since Admission") +
  ylab("Genomic Distance (SNPs)")


legend.p <- get_legend(p1)

plot_grid(p2, p3, p4, legend.p, nrow = 1, rel_widths = c(1,1,1,0.3))


