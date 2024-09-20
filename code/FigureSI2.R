source("code/GlobalFunctions.R")
source("code/LoadData.R")

# Load the image
img_path <- "~/figures/Sampling_blank.png"
img <- readPNG(img_path)

img_dims <- dim(img)  # Get image dimensions

width <- img_dims[2]
height <- img_dims[1]

asp <- height/width

gm.metadata <- gm.metadata %>%
  mutate(granular.sample.loc = if_else(
    sample.loc == "PAT",
    substr(Sample.ID, nchar(Sample.ID), nchar(Sample.ID)),
    sample.loc
  ))

gm.metadata <- gm.metadata %>%
  mutate(granular.sample.loc = case_when(
    sample.loc == "PAT" & Hospital == "A" ~ substr(Sample.ID, nchar(Sample.ID), nchar(Sample.ID)),
    sample.loc == "PAT" & Hospital == "B" ~ substr(Sample.ID, nchar(Sample.ID) - 1, nchar(Sample.ID)),
    TRUE ~ sample.loc
  ))

# Combine categories and sum the counts
gm.metadata <- gm.metadata %>%
  mutate(granular.sample.loc = case_when(
    granular.sample.loc %in% c("NE", "SE") ~ "NE_SE",
    granular.sample.loc %in% c("NS", "SN") ~ "NS_SN",
    granular.sample.loc %in% c("AX") ~ "A",
    granular.sample.loc %in% c("GR") ~ "G",
    granular.sample.loc %in% c("1P","2P", "3P", 
                               "4P", "5P", "6P",
                               "7P", "8P", "9P",
                               "0P") ~ "P",
    TRUE ~ granular.sample.loc
  ))



# Num Samples

# Count the number of samples by sample location
## THIS IS THE DATA FOR MAKING THE HEAT MAP. 
## FOR DIFFERENT HEAT MAPS, CHANGE THIS
sample_count <- gm.metadata %>% 
  count(granular.sample.loc) %>% 
  na.omit()

# Define bounding box coordinates for each sample location (adjust based on the image layout)
coords <- data.frame(
  granular.sample.loc 
  = c("A",   "BR",  "CC",  "DH",   "G", "HH", "NE_SE", "NS_SN", "P",  "SC", "SS",  "TB",  "TO",  "Z1", "Z2", "Z3"),
  x_min = c(0.019, 0.385,  0.295, 0.425, 0.035, 0.88, 0.63,  0.6,    0.125, 0.65, 0.295, 0.495, 0.49, 0.374, 0.346, 0.478),
  x_max = c(0.052, 0.51,  0.345, 0.55,  0.068, 0.999, 0.79,  0.83,   0.162, 0.78, 0.34,  0.55,  0.553, 0.495, 0.555, 0.558),
  y_min = c(0.56,  0.675, 0.065, 0.055,  0.43,  0.73,  0.71,  0.07,   0.455, 0.45, 0.741, 0.41,  0.84,  0.405,  0.115, 0.75),
  y_max = c(0.66,  0.73,  0.33,  0.11,  0.53,  0.403, 0.91,  0.39,   0.555, 0.62, 0.865, 0.499, 0.94,  0.495, 0.36,  0.84))

# Merge sample_count with coordinates
sample_data <- merge(sample_count, coords, by = "granular.sample.loc")

p1 <- ggplot(sample_data) +
  # Overlay the heatmap with custom dimensions
  geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, fill = n), alpha = 0.9) +
  scale_fill_gradient(low = "#e7e1ef", high = "#980043",
                      name = "Number of\nSamples") +
  # Add the background image
  annotation_custom(rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc")), 
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_cowplot() +
  ylim(0,1) + 
  xlim(0,1) + 
  theme(aspect.ratio = asp,
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) 


## Num isolates


# Sample data

gm.cdiff <- gm.metadata %>% 
  filter(C.diff == 1)

iso_count <- gm.cdiff %>% 
  count(granular.sample.loc) %>% 
  na.omit()

# Merge iso_count with coordinates
sample_data <- merge(iso_count, coords, by = "granular.sample.loc")

p2 <- ggplot(sample_data) +
  # Overlay the heatmap with custom dimensions
  geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, fill = n), alpha = 0.9) +
  scale_fill_gradient(low = "#c7e9b4", high = "#253494",
                      name = "Number of\nIsolates") +
  # Add the background image
  annotation_custom(rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc")), 
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_cowplot() +
  ylim(0,1) + 
  xlim(0,1) + 
  theme(aspect.ratio = asp,
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) 


plot_grid(p1, p2, labels = "AUTO", nrow = 2)
