
source("Scripts/Cdiff_GlobalFunctions.R")
source("Scripts/Cdiff_LoadData.R")


###  Supplemental Figure Descriptive Plot ----------

cdiff.descriptive2 <- cdiff.descriptive


cdiff.descriptive <- cdiff.descriptive %>%
  mutate(Short.Room = ifelse(Short.Room == "A28.5", "16.5", Short.Room)) 

cdiff.descriptiveA <-cdiff.descriptive  %>%
  filter(Hospital == "A")  %>%
  arrange(Collection.date) %>%
  mutate(Short.Room.Ordered = factor(Short.Room, 
                                     levels=c("A1", "A2","A3", "A4","A5","A6", "A7", "A8", 
                                              "A9", "A10", "A11", "A12", "A13", "A14", "A15",
                                              "A16", "A16.5", "A17","A18", "A19", "A20"))) 

cdiff.descriptiveB <-cdiff.descriptive  %>%
  filter(Hospital == "B")%>%
  arrange(Collection.date) %>%
  mutate(Short.Room.Ordered = factor(Short.Room, 
                                     levels=c("B1", "B2", "B3", "B4", "B5", "B6", "B7",
                                              "B8","B9","B10"))) 


id.list <- sort(unique(cdiff.descriptive2$Short.ID))

A034.col <- legend.color.picker("A034", id.list)
A056.col <- legend.color.picker("A056", id.list)
A004.col <- legend.color.picker("A004", id.list)
A013.col <- legend.color.picker("A013", id.list)


legend.colorsA <- legend.color.picker(c(sort(unique(cdiff.descriptiveA$Short.ID))), 
                                      sorted.id)
legend.colorsB <- legend.color.picker(c(sort(unique(cdiff.descriptiveB$Short.ID))), 
                                      sorted.id)

f1.A <- ggplot() + 
  geom_point(data = cdiff.descriptiveA %>% filter(C.diff == 1 &contamination == 1), aes(x = Study.Day, y= Short.Room.Ordered,
                                                                                        shape = sample.env), 
             color = "darkgrey", size = 6) +
  geom_point(data = cdiff.descriptiveA %>% filter(C.diff == 1), aes(x = Study.Day, y= Short.Room.Ordered, 
                                                                    color = Short.ID, shape = sample.env), size = 3) + 
  annotate("segment", x = 53, ## day 53 is this persons transfer day
           xend =  53, 
           y = 2, yend = 4, colour = A034.col, 
           size = 0.8, linetype = "dashed") +
  annotate("segment", x = 34, 
           xend =  36, 
           y = 11, yend = 16, colour = A056.col, 
           size = 0.8, linetype = "dashed") +
  annotate("segment", x = 11, 
           xend =  11, 
           y = 12, yend = 5, colour = A004.col, 
           size = 0.8, linetype = "dashed") +
  annotate("segment", x = 11, 
           xend =  11, 
           y = 12, yend = 5, colour = A004.col, 
           size = 0.8, linetype = "dashed") +
  annotate("segment", x = 8, 
           xend =  8, 
           y = 3, yend = 14, colour = A013.col, 
           size = 0.8, linetype = "dashed") +
  # shape = interaction(sample.env, C.diff)), size = 3) +
  geom_errorbarh(data = cdiff.descriptiveA, aes(xmax = (Room.End.Day), 
                                                xmin = (Room.Start.Day), y= Short.Room.Ordered,
                                                color = Short.ID)) + 
  # annotate("rect",xmin = 7, xmax = 24, ## study day 7 and 24 represent start and end of contamination
  #          ymin = 0.5, ymax = 19.5, alpha=0.1,  fill = '#333333FF', 
  #          linewidth=0.1) +
  scale_y_discrete(drop=F, 
                   labels=c("A1" = "1", "A2" = "2","A3" = "3", 
                            "A4" = "4","A5" = "5", "A6" = "6",
                            "A7" = "7", "A8" = "8", "A9" = "9",
                            "A10" = "10", "A11" = "11", "A12" = "12", 
                            "A13" = "13", "A14" = "14", "A15" = "15",
                            "A16" = "16", "A16.5" = "16.5", "A17" = "17",
                            "A18" = "18", "A19" = "19", "A20" = "20")) + # new line added here
  scale_color_manual(values = legend.colorsA, 
                     name = "Unique patient ID") + 
  scale_fill_manual(values = legend.colorsA, 
                    name = "Unique patient ID") + 
  geom_point(data = cdiff.descriptiveA %>% filter(ToxinA == "Positive"), aes(x = Study.Day, y= Short.Room.Ordered), color = "black", shape = 3, stroke = 1.3, size = 3) + 
  #geom_point(data = cdiff.descriptive2 %>% filter(contamination == 1), aes(x = Study.Day, y= Short.Room.Ordered, shape = sample.env), color = "black", shape = 0, stroke = 1.3, size = 3) + 
  ylab("Room Number") +
  xlab("Study Day") + 
  scale_shape_manual(values = c(18, 19, 15,  17),
                     labels = c("HCW Hand","Patient", "Room Env.",
                                "Shared Env."),
                     breaks = c("HCW Hand","Patient", "Room Env.",
                                "Shared Env."),
                     name = "Sample Location") +
  #geom_point(data = cdiff.descriptive2 %>%  filter(contact.prec == 1), aes(x = Study.Day, y= Short.Room.Ordered), shape =8, color = "black") +
  theme_cowplot() +
  background_grid()+ 
  theme(legend.position = "none", 
        aspect.ratio = 0.6)+
  ggtitle("Hospital A") + 
  xlim(0,60)


f1.B <- ggplot() + 
  # geom_point(data = cdiff.descriptiveB %>% filter(contamination == 1), aes(x = Study.Day, y= Short.Room.Ordered,
  #                                                                          shape = sample.env), 
  #            color = "grey", size = 5) +
  geom_point(data = cdiff.descriptiveB %>% filter(C.diff == 1), aes(x = Study.Day, y= Short.Room.Ordered, 
                                                                    color = Short.ID, shape = sample.env), size = 3) + 
  theme(legend.position = "none") +
  # shape = interaction(sample.env, C.diff)), size = 3) +
  geom_errorbarh(data = cdiff.descriptiveB, aes(xmax = (Room.End.Day), 
                                                xmin = (Room.Start.Day), y= Short.Room.Ordered,
                                                color = Short.ID), 
                 height = 0.4) + 
  scale_y_discrete(drop=F, 
                   labels=c("B1" = "1", "B2" = "2", "B3" = "3", 
                            "B4" = "4", "B5" = "5", "B6" = "6",
                            "B7" = "7", "B8" = "8","B9" = "9",
                            "B10" = "10")) + # new line added here
  scale_color_manual(values = legend.colorsB, 
                     name = "Unique patient ID") + 
  scale_fill_manual(values = legend.colorsB, 
                    name = "Unique patient ID") + 
  geom_point(data = cdiff.descriptiveB %>% filter(ToxinA == "Positive"), aes(x = Study.Day, y= Short.Room.Ordered), color = "black", shape = 3, stroke = 1.3, size = 3) + 
  #geom_point(data = cdiff.descriptive2 %>% filter(contamination == 1), aes(x = Study.Day, y= Short.Room.Ordered, shape = sample.env), color = "black", shape = 0, stroke = 1.3, size = 3) + 
  ylab(" ") +
  ylab("Room Number") +
  xlab("Study Day") + 
  scale_shape_manual(values = c(18, 19, 15,  17),
                     labels = c("HCW Hand","Patient", "Room Env.",
                                "Shared Env."),
                     breaks = c("HCW Hand","Patient", "Room Env.",
                                "Shared Env."),
                     name = "Sample Location") +
  #geom_point(data = cdiff.descriptive2 %>%  filter(contact.prec == 1), aes(x = Study.Day, y= Short.Room.Ordered), shape =8, color = "black") +
  theme_cowplot() +
  background_grid() + 
  theme(legend.position = "none", 
        aspect.ratio = 0.6) +
  ggtitle("Hospital B") +
  xlim(90, 130)

figure1.dummy <- cdiff.descriptive

figure1.dummy$sample.env[3] <-  "Positive"
figure1.dummy$Study.Day[3] <-  0

figure1.dummy$sample.env[3] <-  "Positive"

shapes <- c(18, 19, 15,  17, 3)
shapes <- setNames(shapes, c("HCW Hand","Patient", "Room Env.",
                             "Shared Env.", "Positive"))

legend.colors <- legend.color.picker(c(sort(unique(figure1.dummy$Short.ID))), 
                                     sorted.id)

f1.legend <-  ggplot() + 
  geom_point(data = figure1.dummy %>% filter(C.diff == 1), aes(x = Study.Day, y= Short.Room, 
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
  # geom_point(data = figure1.dummy %>% filter(ToxinA == "Positive"), 
  #            aes(x = Study.Day, y= Short.Room.Ordered), color = "black", shape = 3, stroke = 1.3, size = 3) + 
  geom_point(data = figure1.dummy, aes(x = Study.Day, y= Short.Room, 
                                       color = Short.ID, shape = sample.env), size = 3) + 
  scale_shape_manual(values = shapes,
                     labels = c("Positive"),
                     breaks = c("Positive"),
                     name = "Toxin",
                     guide = guide_legend(title.position = "top", order = 2)) +  
  theme_cowplot() +
  background_grid() + 
  theme(legend.position = "bottom", 
        legend.direction="vertical")




f1.legend <- get_legend(f1.legend)

f1 <- plot_grid(f1.A, f1.B, f1.legend, nrow = 1, 
                rel_widths = c(1,1,0.2), labels = c("A", "B", ""))

vert.f1 <- plot_grid(f1.A, f1.B,  nrow = 1,
                     rel_widths = c(1,1), labels = c("A", "B"))

vert.legend.f1 <- plot_grid(space, f1.legend, space, nrow = 1, rel_widths = c(0.5,1,0.5))

plot_grid(vert.legend.f1,vert.f1,  nrow = 2, rel_heights = c(0.3,1))
