source("code/GlobalFunctions.R")
source("code/LoadData.R")


## 2.3 Within patient stay -----------------------------------------------------

cdiff.descriptive3 <- cdiff.descriptive2 %>% 
  group_by(Short.ID) %>% 
  mutate(Patient.Day = as.numeric(as.Date(Collection.date) - as.Date(Admit) +1)) %>%
  mutate(Admit.Patient.Day = as.numeric(as.Date(Admit) - as.Date(Admit) +1)) %>%
  mutate(Disch.Patient.Day = as.numeric(as.Date(Disch) - as.Date(Admit) +1)) %>%
  mutate(Cdiff.Toxin = replace_na(Cdiff.Toxin, "none")) %>% 
  filter(Short.ID != "ASN")



### Figure 3

cdiff.descriptive.hcp <- cdiff.descriptive3 %>%
  mutate(sample.env = recode(sample.env, "HCW Hand" = "HCP Hand" ))


y.names<- interaction(cdiff.descriptive.hcp$sample.env, cdiff.descriptive.hcp$Short.ID)

p.acq.A <- ggplot(cdiff.descriptive.hcp %>% filter(Hospital == "A", Short.ID != "A912"), 
                  aes(x = Patient.Day, y = factor(sample.env,levels=c("HCP Hand", "Room Env.", "Patient")))) +
  geom_point(aes(color = as.factor(C.diff), shape = sample.env), size = 3) +
  geom_point(data = filter(cdiff.descriptive.hcp  %>% filter(Hospital == "A", Short.ID != "A912"), C.diff == "1"), aes(shape =sample.env), color = "#fa7367", size = 3) + 
  geom_point(data = filter(cdiff.descriptive.hcp %>% filter(Hospital == "A", Short.ID != "A912"), ToxinA == "Positive"), color = "black", shape = 3, size = 3) + 
  scale_color_manual(values = c("1" = "#fa7367", "0" = "#BDC6DE"), 
                     labels = c("1" = "Positive", "0" = "Negative"),
                     name = "Cdiff") +
  scale_shape_manual(values = c("HCP Hand" = 18, 
                                "Patient" = 19, 
                                "Room Env." = 15),
                     labels = c("HCP Hand" = "HCP Hand", 
                                "Patient" = "Patient", 
                                "Room Env." = "Room Env."),
                     name = "Sample Location") +
  xlab("Days since Admission") + 
  ylab("") +
  scale_y_discrete(labels = c("HCP Hand.A002" = " ",
                              "Patient.A002" = "A002",
                              "Room Env..A002" = "",
                              "HCP Hand.A004" = "",
                              "Patient.A004" = "A004",
                              "Room Env..A004" = "",
                              "HCP Hand.A007" = "",
                              "Patient.A007" = "A007",
                              "Room Env..A007" = "",
                              "HCP Hand.A010" = "",
                              "Patient.A010" = "A010",
                              "Room Env..A010" = "",
                              "HCP Hand.A012" = "",
                              "Patient.A012" = "A012",
                              "Room Env..A012" = "",
                              "HCP Hand.A013" = "",
                              "Patient.A013" = "A013",
                              "Room Env..A013" = "",
                              "HCP Hand.A034" = "",
                              "Patient.A034" = "A034",
                              "Room Env..A034" = "",
                              "HCP Hand.A046" = "",
                              "Patient.A046" = "A046",
                              "Room Env..A046" = "",
                              "HCP Hand.A048" = "",
                              "Patient.A048" = "A048",
                              "Room Env..A048" = "",
                              "HCP Hand.A051" = "",
                              "Patient.A051" = "A051",
                              "Room Env..A051" = "",
                              "HCP Hand.A052" = "",
                              "Patient.A052" = "A052",
                              "Room Env..A052" = "",
                              "HCP Hand.A056" = "",
                              "Patient.A056" = "A056",
                              "Room Env..A056" = "",
                              "HCP Hand.A058" = "",
                              "Patient.A058" = "A058",
                              "Room Env..A058" = "",
                              "HCP Hand.A079" = "",
                              "Patient.A079" = "A079",
                              "Room Env..A079" = "",
                              "HCP Hand.A231" = "",
                              "Patient.A231" = "",
                              "Room Env..A231" = "A231",
                              "HCP Hand.A233" = "",
                              "Patient.A233" = "",
                              "Room Env..A233" = "A233",
                              "HCP Hand.A238" = "",
                              "Patient.A238" = "",
                              "Room Env..A238" = "A238",
                              "HCP Hand.A239" = "",
                              "Patient.A239" = "",
                              "Room Env..A239" = "A239",
                              "HCP Hand.A250" = "",
                              "Patient.A250" = "",
                              "Room Env..A250" = "A250",
                              "HCP Hand.A273" = "",
                              "Patient.A273" = "",
                              "Room Env..A273" ="A273",
                              "HCP Hand.A279" = "",
                              "Patient.A279" = "",
                              "Room Env..A279" = "A279",
                              "HCP Hand.A912" = "",
                              "Patient.A912" = "",
                              "Room Env..A912" = "A912",
                              "HCP Hand.B002" = "",
                              "Patient.B002" = "B002",
                              "Room Env..B002" = "",
                              "HCP Hand.B003" = "",
                              "Patient.B003" = "",
                              "Room Env..B003" = "B003",
                              "HCP Hand.B006" = "",
                              "Patient.B006" = "B006",
                              "Room Env..B006" = "",
                              "HCP Hand.B008" = "",
                              "Patient.B008" = "B008",
                              "Room Env..B008" = "",
                              "HCP Hand.B009" = "",
                              "Patient.B009" = "B009",
                              "Room Env..B009" = "",
                              "HCP Hand.B012" = "",
                              "Patient.B012" = "B012",
                              "Room Env..B012" = "",
                              "HCP Hand.B015" = "",
                              "Patient.B015" = "B015",
                              "Room Env..B015" = "",
                              "HCP Hand.B018" = "",
                              "Patient.B018" = "B018",
                              "Room Env..B018" = "",
                              "HCP Hand.B019" = "",
                              "Patient.B019" = "B019",
                              "Room Env..B019" = "",
                              "HCP Hand.B048" = "",
                              "Patient.B048" = "B048",
                              "Room Env..B048" = "",
                              "HCP Hand.B058" = "",
                              "Patient.B058" = "B058",
                              "Room Env..B058" = "",
                              "HCW Hand.B211" = "",
                              "Patient.B211" = "",
                              "Room Env..B211" = "B211",
                              "HCP Hand.B212" = "",
                              "Patient.B212" =  "",
                              "Room Env..B212" =  "B212")) +
  scale_x_continuous(limits = c(0,60),
                     breaks = c(0,10,20,30,40,50,60))


p.acq.B <- ggplot(cdiff.descriptive.hcp %>% filter(Hospital == "B"),
                  aes(x = Patient.Day, y = factor(sample.env,levels=c("HCP Hand", "Room Env.", "Patient")))) +
  geom_point(aes(color = as.factor(C.diff), shape = sample.env), size = 3) +
  geom_point(data = filter(cdiff.descriptive.hcp  %>% filter(Hospital == "B"), C.diff == "1"), aes(shape =sample.env), color = "#fa7367", size = 3) + 
  geom_point(data = filter(cdiff.descriptive.hcp %>% filter(Hospital == "B"), ToxinA == "Positive"), color = "black", shape = 3, size = 3) + 
  scale_color_manual(values = c("1" = "#fa7367", "0" = "#BDC6DE"), 
                     labels = c("1" = "Positive", "0" = "Negative"),
                     name = "C. difficile") +
  scale_shape_manual(values = c("HCP Hand" = 18, 
                                "Patient" = 19, 
                                "Room Env." = 15),
                     labels = c("HCP Hand" = "HCP Hand", 
                                "Patient" = "Patient", 
                                "Room Env." = "Room Env."),
                     name = "Sample Location") +
  guides(shape = "none") +
  xlab("Days since Admission") + 
  ylab(" ") +
  scale_y_discrete(labels = c("HCP Hand.A002" = " ",
                              "Patient.A002" = "A002",
                              "Room Env..A002" = "",
                              "HCP Hand.A004" = "",
                              "Patient.A004" = "A004",
                              "Room Env..A004" = "",
                              "HCP Hand.A007" = "",
                              "Patient.A007" = "A007",
                              "Room Env..A007" = "",
                              "HCP Hand.A010" = "",
                              "Patient.A010" = "A010",
                              "Room Env..A010" = "",
                              "HCP Hand.A012" = "",
                              "Patient.A012" = "A012",
                              "Room Env..A012" = "",
                              "HCP Hand.A013" = "",
                              "Patient.A013" = "A013",
                              "Room Env..A013" = "",
                              "HCP Hand.A034" = "",
                              "Patient.A034" = "A034",
                              "Room Env..A034" = "",
                              "HCP Hand.A046" = "",
                              "Patient.A046" = "A046",
                              "Room Env..A046" = "",
                              "HCP Hand.A048" = "",
                              "Patient.A048" = "A048",
                              "Room Env..A048" = "",
                              "HCP Hand.A051" = "",
                              "Patient.A051" = "A051",
                              "Room Env..A051" = "",
                              "HCP Hand.A052" = "",
                              "Patient.A052" = "A052",
                              "Room Env..A052" = "",
                              "HCP Hand.A056" = "",
                              "Patient.A056" = "A056",
                              "Room Env..A056" = "",
                              "HCP Hand.A058" = "",
                              "Patient.A058" = "A058",
                              "Room Env..A058" = "",
                              "HCP Hand.A079" = "",
                              "Patient.A079" = "A079",
                              "Room Env..A079" = "",
                              "HCP Hand.A231" = "",
                              "Patient.A231" = "",
                              "Room Env..A231" = "A231",
                              "HCP Hand.A233" = "",
                              "Patient.A233" = "",
                              "Room Env..A233" = "A233",
                              "HCP Hand.A238" = "",
                              "Patient.A238" = "",
                              "Room Env..A238" = "A238",
                              "HCP Hand.A239" = "",
                              "Patient.A239" = "",
                              "Room Env..A239" = "A239",
                              "HCP Hand.A250" = "",
                              "Patient.A250" = "",
                              "Room Env..A250" = "A250",
                              "HCP Hand.A273" = "",
                              "Patient.A273" = "",
                              "Room Env..A273" ="A273",
                              "HCP Hand.A279" = "",
                              "Patient.A279" = "",
                              "Room Env..A279" = "A279",
                              "HCP Hand.A912" = "",
                              "Patient.A912" = "",
                              "Room Env..A912" = "A912",
                              "HCP Hand.B002" = "",
                              "Patient.B002" = "B002",
                              "Room Env..B002" = "",
                              "HCP Hand.B003" = "",
                              "Patient.B003" = "",
                              "Room Env..B003" = "B003",
                              "HCP Hand.B006" = "",
                              "Patient.B006" = "B006",
                              "Room Env..B006" = "",
                              "HCP Hand.B008" = "",
                              "Patient.B008" = "B008",
                              "Room Env..B008" = "",
                              "HCP Hand.B009" = "",
                              "Patient.B009" = "B009",
                              "Room Env..B009" = "",
                              "HCP Hand.B012" = "",
                              "Patient.B012" = "B012",
                              "Room Env..B012" = "",
                              "HCP Hand.B015" = "",
                              "Patient.B015" = "B015",
                              "Room Env..B015" = "",
                              "HCP Hand.B018" = "",
                              "Patient.B018" = "B018",
                              "Room Env..B018" = "",
                              "HCP Hand.B019" = "",
                              "Patient.B019" = "B019",
                              "Room Env..B019" = "",
                              "HCP Hand.B048" = "",
                              "Patient.B048" = "B048",
                              "Room Env..B048" = "",
                              "HCP Hand.B058" = "",
                              "Patient.B058" = "B058",
                              "Room Env..B058" = "",
                              "HCP Hand.B211" = "",
                              "Patient.B211" = "",
                              "Room Env..B211" = "B211",
                              "HCP Hand.B212" = "",
                              "Patient.B212" =  "",
                              "Room Env..B212" =  "B212")) +
  scale_x_continuous(limits = c(0,40),
                     breaks = c(0,10,20,30,40))

admit.p1<- p.acq.A + facet_grid(Short.ID ~ .,
                                switch = "y" # flip the facet labels along the y axis from the right side to the left
) + 
  theme_bw(base_size=13) + # remove the word "values"
  theme(strip.background = element_blank(), # remove the background
        strip.placement = "outside") +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

admit.p2<-p.acq.B + facet_grid(Short.ID ~ .,
                               switch = "y" # flip the facet labels along the y axis from the right side to the left
) + 
  theme_bw(base_size=13) + # remove the word "values"
  theme(strip.background = element_blank(), # remove the background
        strip.placement = "outside") +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

acq.dummy <- cdiff.descriptive3 
acq.dummy$sample.env[3] <-  "Positive"


shapes <- c(18, 19, 15,  17, 3)
shapes <- setNames(shapes, c("HCP Hand","Patient", "Room Env.",
                             "Shared Env.", "Positive"))

p.legend <- ggplot(acq.dummy, aes(x = Patient.Day, y = sample.env)) +
  geom_point(aes(color = as.factor(C.diff), shape = sample.env), size = 3) +
  geom_point(data = filter(acq.dummy  %>% filter(C.diff == "1")), 
             aes(shape =sample.env), color = "#fa7367", size = 3) + 
  geom_point(data = filter(acq.dummy %>% filter(ToxinA == "Positive")), 
             color = "black", shape = 3, size = 3) + 
  theme_cowplot() +
  background_grid() +
  scale_color_manual(values = c("1" = "#fa7367", "0" = "#BDC6DE"), 
                     labels = c("1" = "Isolate", "0" = "Sample"),
                     name = expression(italic("C. difficile"))) +
  scale_shape_manual(values = shapes,
                     labels = c("Toxigenic"),
                     breaks = c("Positive"),
                     name = "Toxin Status",
                     guide = guide_legend(title.position = "top", order = 2))  +
  guides(shape = guide_legend(override.aes = list(color = "black"))) + 
  xlab("Days since Admission") + 
  ylab(" ") +
  xlim(1,60) + 
  theme(legend.position = "bottom", 
        legend.direction="vertical")

acq.legend <- cowplot::get_plot_component(p.legend, "guide-box", return_all = TRUE)[[3]]
space<- ggplot() + theme_void()

legend.spacing <- plot_grid(space, acq.legend, space, rel_widths = c(0.3, 1, 0.1))

left.plot <- plot_grid(space, legend.spacing, admit.p2, rel_heights = c(0.05,0.1,1), nrow = 3)
plot_grid(admit.p1, left.plot, rel_widths = c(1,0.71,0.3), nrow = 1)

