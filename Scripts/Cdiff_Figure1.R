source("Scripts/Cdiff_GlobalFunctions.R")
source("Scripts/Cdiff_LoadData.R")


## Panel A: Prevalence ----


gm.metadata %>% filter(Bacterial.Species == "C.difficile") %>% group_by(Study.Day) %>% count() -> Postive.by.day
gm.metadata %>% group_by(Study.Day) %>% count() -> total.samps.day
gm.metadata %>% filter(Bacterial.Species == "C.difficile", sample.env == "Room Env.") %>% group_by(Study.Day) %>% count() -> Postive.Room.by.day
gm.metadata %>% filter(Bacterial.Species == "C.difficile", sample.env == "Patient") %>% group_by(Study.Day) %>% count() -> Postive.pt.by.day
gm.metadata %>% filter(Bacterial.Species == "C.difficile", sample.env == "HCW Hand") %>% group_by(Study.Day) %>% count() -> Postive.HCP.by.day



total.samps.day <- total.samps.day %>%
  left_join(Postive.by.day, by = "Study.Day", suffix = c(".A", ".B")) %>%
  rename(n = n.A,
         "Num.Pos" = n.B)

total.samps.day <- total.samps.day %>%
  left_join(Postive.Room.by.day, by = "Study.Day", suffix = c(".A", ".B")) %>%
  rename(n = n.A,
         "Num.Pos.Room" = n.B)

total.samps.day <- total.samps.day %>%
  left_join(Postive.pt.by.day, by = "Study.Day", suffix = c(".A", ".B")) %>%
  rename(n = n.A,
         "Num.Pos.Pat" = n.B)

total.samps.day <- total.samps.day %>%
  left_join(Postive.HCP.by.day, by = "Study.Day", suffix = c(".A", ".B")) %>%
  rename("Samps" = n.A,
         "Num.Pos.HCP" = n.B)

total.samps.day <- mutate_all(total.samps.day, ~replace(., is.na(.), 0))

missing.days <- data.frame(Study.Day = seq(57, 92, by = 1),
                           Samps = rep(NA, 36),
                           Num.Pos = rep(NA, 36),
                           Num.Pos.Room = rep(NA, 36), 
                           Num.Pos.Pat = rep(NA, 36),
                           Num.Pos.HCP = rep(NA, 36)
)

total.samps.day <- rbind(total.samps.day, missing.days)

long_total_samps_day <- total.samps.day %>%
  mutate(across(starts_with("Num."), ~./Samps)) %>%
  pivot_longer(cols = starts_with("Num."),
               names_to = "Sample_Type",
               values_to = "Proportion")



f1 <- ggplot(long_total_samps_day) + 
  geom_line(aes(x = Study.Day, y = Proportion*100, group = Sample_Type, color = Sample_Type), size =0.8) + 
  scale_color_manual(limits= c("Num.Pos", "Num.Pos.HCP", "Num.Pos.Pat", "Num.Pos.Room"),
                     labels= c("All", "HCP Hands","Patient Body Sites", "Room E"),
                     values = c("#132157FF","#1BB6AFFF","#FFAD0AFF" , "#D72000FF"), 
                     name = "Sample Location") + 
  xlab("Study Day") +
  ylab("Prevalence (Num. Isolates/Num. Samples)") +
  theme_cowplot() +
  theme(legend.position = c(0.65,0.75)) + 
  annotate("text", x = 30, y = 19, label = "Hospital A") + 
  annotate("text", x = 100, y = 19, label = "Hospital B")

## Daily patient count --- 
census.data <- gm.metadata %>% filter(Unit %in% c("CVICU", "VA_SICU"))


counts.census.data <- census.data %>%
  group_by(Study.Day, C.diff) %>%
  summarize(unique_patients = n_distinct(PATID))

census.missing.days <- data.frame(Study.Day = c(seq(57, 92, by = 1), seq(57, 92, by = 1)),
                                  C.diff = c(rep(1, 36), rep(0, 36)),
                                  unique_patients = c(rep(NA, 36),rep(NA, 36))
                                  
)

counts.census.data <- rbind(counts.census.data, census.missing.days)

# Calculate the total unique patients with C.diff == 1 and C.diff == 2 for each day
total_patients <- counts.census.data %>%
  group_by(Study.Day) %>%
  summarize(unique_total_patients = sum(unique_patients[C.diff == 1 | C.diff == 0]))

# Filter the dataframe for C.diff == 1
cdiff_1 <- counts.census.data %>%
  filter(C.diff == 1)

# Calculate the ratio of unique patients with cdiff == 1 to all patients with cdiff == 1 + cdiff == 2 for each day
result <- cdiff_1 %>%
  left_join(total_patients, by = "Study.Day") %>%
  mutate(ratio = unique_patients / unique_total_patients) 

long_result <- result %>%
  select(-ratio) %>%
  pivot_longer(cols = starts_with("Unique"),
               names_to = "PatientType",
               values_to = "Count")

###### Panel B: Toixgenic ---

counts.census.data <- census.data %>%
  group_by(Study.Day, C.diff, ToxinA) %>%
  summarize(unique_patients = n_distinct(PATID))

census.missing.days <- data.frame(Study.Day = c(seq(57, 92, by = 1), seq(57, 92, by = 1)),
                                  C.diff = c(rep(1, 36), rep(0, 36)),
                                  unique_patients = c(rep(NA, 36),rep(NA, 36))
                                  
)

counts.census.data <- rbind(counts.census.data, census.missing.days)

# Calculate the total unique patients with C.diff == 1 and C.diff == 2 for each day
total_patients <- counts.census.data %>%
  group_by(Study.Day) %>%
  summarize(unique_total_patients = sum(unique_patients[C.diff == 1]))

# Filter the dataframe for C.diff == 1
cdiff_1_tox <- counts.census.data %>%
  filter(C.diff == 1, ToxinA == "Positive") %>% 
  rename("unique_toxigenic_patients" = unique_patients) %>% 
  select(-c(ToxinA))

cdiff_1_nontox <- counts.census.data %>%
  filter(C.diff == 1, is.na(ToxinA)) %>%
  mutate(ToxinA = if_else(is.na(ToxinA), "Negative", ToxinA))%>% 
  rename("unique_nontoxigenic_patients" = unique_patients)  %>% 
  select(-ToxinA)

result <- cdiff_1_nontox %>%
  left_join(cdiff_1_tox, by = "Study.Day")


long_result <- result %>%
  pivot_longer(cols = starts_with("Unique"),
               names_to = "PatientType",
               values_to = "Count")

long_result <- long_result  %>%
  mutate(Count = if_else(Study.Day < 57 & is.na(Count), 0, Count))%>%
  mutate(Count = if_else(Study.Day > 92 & is.na(Count), 0, Count))


tox <- ggplot(long_result) + 
  geom_line(aes(x = Study.Day, y = Count, group = PatientType, color = PatientType),  size = 0.8) +
  ylim(0,5) +
  theme_cowplot() +
  ylab("Number of Isolates") +
  xlab("Study Day") + 
  scale_color_manual(values = c("#F0BB7B", "#FC6467"), 
                     name = "Toxin Status", 
                     labels = c("Non-toxigenic", "Toxigenic")) + 
  theme(legend.position = c(0.75, 0.75)) + 
  annotate("text", x = 25, y = 5, label = "Hospital A") + 
  annotate("text", x = 100, y = 5, label = "Hospital B")


## Panel C: MLST -------------------

cdiff.neg <- cdiff.descriptive %>% 
  filter(is.na(ToxinA)) %>% 
  group_by(MLST.ST) %>% 
  count()

cdiff.pos <- cdiff.descriptive %>% 
  filter(ToxinA != "NA") %>% 
  group_by(MLST.ST) %>% 
  count()

mlst.dat <- cdiff.descriptive2 %>% filter(C.diff == 1, !is.na(MLST.ST))

get_study_week <- function(day) {
  return(ceiling(day / 7))
}

mlst.dat$Study.Week <- NULL
for(i in 1:nrow(mlst.dat)){
  mlst.dat$Study.Week[i] <- get_study_week(mlst.dat$Study.Day[i] )
}

f2 <- ggplot(mlst.dat, 
             aes(x = Study.Week, group = MLST.ST, fill = as.factor(MLST.ST))) + 
  geom_bar( alpha=0.8) +
  scale_fill_manual(breaks = c("2", "3", "14", "15", "26", "39", "41", "42", "55", "110", "125"),
                    values = c("#972d14", "#d8b709","#fd6467", 
                               "#02401b", "#f1bb7b","#a2a475",
                               "#d67237", "#446455", "#f3df6c",
                               "#c7b19c", "#7fa6ad")) + 
  guides(fill=guide_legend(title="ST")) +
  xlab("Study Week") + 
  ylab("Number of Isolates") +
  theme_cowplot() + 
  annotate("text", x = 4, y = 70, label = "Hospital A") + 
  annotate("text", x = 16, y = 70, label = "Hospital B") +
  theme(legend.position = c(0.9, 0.6))



# Panel D: MLST by day
mlst.counts <- mlst.dat %>%
  group_by(MLST.ST, Study.Week) %>%
  summarise(count = n())

get_study_week(127)

all.days.mlst <- data.frame(Study.Week = seq(1,19, by = 1), 
                            count = rep(0, 19), 
                            MLST.ST = rep(NA, 19))


all.mlst <- all.days.mlst %>%
  left_join(mlst.counts, by = "Study.Week", suffix = c(".A", ".B")) %>%
  mutate(counts = coalesce(count.B, count.A)) %>%
  mutate(MLST.ST = coalesce(MLST.ST.B, MLST.ST.A )) %>%
  select(Study.Week, MLST.ST, counts)

topf <- plot_grid(f1, tox, f2, labels = c("A", "B", "C"), nrow = 1)


## 
all.cdiff <- rbind(cdiff.descriptiveA, cdiff.descriptiveB)



all.cdiff <- all.cdiff  %>%
  mutate(Short.Room.Ordered = factor(Short.Room.Ordered, 
                                     levels=c("A1", "A2","A3", "A4","A5","A6", "A7", "A8", 
                                              "A9", "A10", "A11", "A12", "A13", "A14", "A15",
                                              "A16", "A16.5", "A17", "A18","A19", "A20", 
                                              "B1", "B2", "B3", "B4", "B5", "B6", "B7",
                                              "B8","B9","B10")))


ST.plot <- 
  ggplot() +
  geom_errorbarh(data = all.cdiff %>% filter(!is.na(MLST.ST)), aes(xmax = (Room.End.Day), 
                                                                   xmin = (Room.Start.Day), y= Short.Room.Ordered, color = factor(MLST.ST))) + 
  annotate("segment", x = 8,
           xend =  8,
           y = 3, yend = 14, colour = "#d8b709",
           size = 0.8, linetype = "dashed") +
  geom_point(data = all.cdiff %>% filter(!is.na(MLST.ST)), 
             aes(x = Study.Day, y = Short.Room.Ordered, group = Short.ID, color = factor(MLST.ST), 
                 shape = sample.env), size = 3) +
  geom_point(data = all.cdiff%>% filter(MLST.ST=="3"), 
             aes(x = Study.Day, y = Short.Room.Ordered, group = Short.ID, color = factor(MLST.ST), 
                 shape = sample.env), size = 3) +
  scale_color_manual(breaks = c("2", "3", "14", "15", "26", "39", "41", "42", "55", "110", "125"),
                     values = c("#972d14", "#d8b709","#fd6467", 
                                "#02401b", "#f1bb7b","#a2a475",
                                "#d67237", "#446455", "#f3df6c",
                                "#c7b19c", "#7fa6ad"),
                     name = "ST") +
  ylab("Room Number") +
  xlab("Study Day") +
  theme_cowplot() +
  background_grid() +
  scale_shape_manual(values = c(19, 15, 18, 17), 
                     breaks = c("Patient", "Room Env.",
                                "HCW Hand", "Shared Env."), 
                     labels = c("Patient", "Room Env.",
                                "HCP Hand", "Shared Env."),
                     name = "Sample\nLocation") +
  theme(legend.position = "none") +
  scale_y_discrete(drop=F, 
                   labels=c("A1" = "A1","A2" = "A2",
                            "A3" = "A3","A4" = "A4",
                            "A5" = "A5","A6" = "A6",
                            "A7" = "A7","A8" = "A8",
                            "A9" = "A9","A10" = "A10",
                            "A11" = "A11","A12" = "A12",
                            "A13" = "A13","A14" = "A14",
                            "A15" = "A15","A16" = "A16", "A16.5" = "A16.5",
                            "A17" = "A17","A18" = "A18",
                            "A19" = "A19","A20" = "A20",
                            "B1" = "B1","B2" = "B2",
                            "B3" = "B3","B4" = "B4",
                            "B5" = "B5","B6" = "B6",
                            "B7" = "B7","B8" = "B8",
                            "B9" = "B9","B10" = "B10"))


ST.legend.plot <- ggplot() +
  geom_point(data = all.cdiff %>% filter(!is.na(MLST.ST)), 
             aes(x = Study.Day, y = Short.Room.Ordered, group = Short.ID, color = factor(MLST.ST), 
                 shape = sample.env), size = 3) +
  scale_color_manual(breaks = c("2", "3", "14", "15", "26", "39", "41", "42", "55", "110", "125"),
                     values = c("#972d14", "#d8b709","#fd6467", 
                                "#02401b", "#f1bb7b","#a2a475",
                                "#d67237", "#446455", "#f3df6c",
                                "#c7b19c", "#7fa6ad"),
                     name = "ST")+
  scale_shape_manual(values = c(19, 15, 18, 17), 
                     breaks = c("Patient", "Room Env.",
                                "HCW Hand", "Shared Env."), 
                     labels = c("Patient", "Room Env.",
                                "HCP Hand", "Shared Env."),
                     name = "Sample Location") +
  theme_cowplot()

st.legend <- get_legend(ST.legend.plot)


## Panel E: 2 ST 1 room -----
### NEEDS STUFF FROM CLUSTERING ANALYSIS 

node.dat1 <- dat %>% 
  filter(cluster == 1) %>% 
  select(from, short.ID) %>% 
  distinct(from, short.ID) %>% 
  mutate(sample.loc = from) %>% 
  filter(short.ID == "A046")

short.name1.1 <- unique(node.dat1$from)

## Make the dataset we need
cdiff.descriptive1.1 <- cdiff.descriptive %>% 
  filter(Collection.Sample.ID %in% short.name1.1, negate = TRUE) %>% 
  mutate(clusternum = "1")

node.dat3 <- dat %>% 
  filter(cluster == 5) %>% 
  select(from, short.ID) %>% 
  distinct(from, short.ID) %>% 
  mutate(sample.loc = from) %>% 
  filter(short.ID == "A046")

short.name3.1 <- unique(node.dat3$from)

## Make the dataset we need
cdiff.descriptive3.1 <- cdiff.descriptive %>% 
  filter(Collection.Sample.ID %in% short.name3.1, negate = TRUE) %>% 
  mutate(clusternum = "3")

cdiff.descriptive.A046 <- rbind(cdiff.descriptive1.1,cdiff.descriptive3.1)


p.A046 <- ggplot() + 
  geom_point(data = cdiff.descriptive.A046 %>% filter(!is.na(Cdiff.Toxin)), 
             aes(x = Study.Day, y = factor(sample.env,levels=c("HCW Hand", "Room Env.", "Patient")), color = clusternum, 
                 shape = sample.env), size = 3) + 
  scale_color_manual(breaks = c("1", "3"), 
                     labels = c("A (15)", "C (3)"),
                     values = c("#02401b","#d8b709"), 
                     name = "Cluster (Sequence Type)") + 
  xlim(18, 42) + 
  ylab("Sample Location") +
  ggtitle("Occupant ID: A046") +
  scale_y_discrete(labels = c("HCP\nHand","Room\nEnv.", "Patient")) +
  xlab("Study Day") +
  theme_cowplot() +
  background_grid() +
  scale_shape_manual(values = c(19, 15, 18, 17), 
                     breaks = c("Patient", "Room Env.",
                                "HCW Hand", "Shared Env."),
                     name = "Sample Location") +
  theme(legend.position = "none") 

cdiff.descriptiveA250 <- cdiff.descriptive %>% filter(Short.ID == "A250", C.diff == "1", !is.na(MLST.ST))

p.A250 <- ggplot() + 
  geom_point(data = cdiff.descriptiveA250, 
             aes(x = Study.Day, 
                 y = factor(sample.env,levels=c("HCW Hand", "Room Env.", "Patient")), 
                 color = as.factor(MLST.ST), 
                 shape = sample.env), size = 3) + 
  scale_color_manual(breaks = c("2", "41"), 
                     labels = c("B (2)", "NA (41)"),
                     values = c("#972d14", "#d67237"), 
                     name = "Cluster (Sequence Type)") + 
  xlim(1, 5) + 
  ylab("Sample Location") +
  ggtitle("Occupant ID: A250") +
  scale_y_discrete(labels = c("HCP\nHand","Room\nEnv.", "Patient")) +
  xlab("Study Day") +
  theme_cowplot() +
  background_grid() +
  scale_shape_manual(values = c(19, 15, 18, 17), 
                     breaks = c("Patient", "Room Env.",
                                "HCW Hand", "Shared Env."),
                     name = "Sample Location") +
  theme(legend.position = "none") 
legend.data <- cdiff.descriptive

cdiff.descriptiveA250$Short.Room.Ordered <- "A22"
cdiff.descriptiveA250$clusternum <- "NA"

f3.1 <- plot_grid(p.A250, p.A046, nrow = 2,  labels = c("E.i", "E.ii"))

## put it all together

bottomf <-plot_grid(ST.plot, f3.1, st.legend, nrow = 1, rel_widths = c(1,0.5,0.15),labels = c("D", "", ""))

plot_grid(topf, bottomf, nrow = 2, rel_heights = c(1, 1.1))
