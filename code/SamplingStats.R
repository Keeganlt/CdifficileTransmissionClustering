source("Scripts/Cdiff_GlobalFunctions.R")
source("Scripts/Cdiff_LoadData.R")



## Patient only

cdiff.descriptive3.pat<-cdiff.descriptive3 %>% filter(sample.env == "Patient")

samples_per_patient <- cdiff.descriptive3 %>%
  group_by(Short.ID) %>%
  summarise(Num_Samples = n())


# Calculate the median, mean, and IQR of the number of samples
summary_samples <- summary(samples_per_patient$Num_Samples)

samples_per_patient_gm <- gm.metadata %>%
  group_by(Short.ID) %>%
  summarise(Num_Samples = n()) 

df_filtered <- samples_per_patient_gm %>%
  filter(as.numeric(gsub("[^0-9]", "", Short.ID)) <= 900)

summary(df_filtered$Num_Samples)

# Calculate the median, mean, and IQR of the number of samples
summary_samples <- summary(samples_per_patient$Num_Samples)


df_filtered2 <- df_filtered %>%
  filter(!grepl("^A", Short.ID))

cdiff.los.data <- cdiff.descriptive3 %>%
  group_by(PATID) %>%
  summarise(
    LOS = as.numeric(max(Disch) - min(Admit))
  )

cdiff.los.data <- cdiff.los.data %>% mutate("group" = "Cdiff")

# Calculate the median, mean, and IQR of the LOS
summary.los <- summary(cdiff.los.data$LOS, na.rm = TRUE)

not.cdiff.los.data <- gm.metadata %>%
  filter(PATID %!in% unique(cdiff.descriptive3$PATID)) %>% 
  group_by(PATID) %>%
  summarise(
    LOS = as.numeric(max(as.Date(Disch)) - min(as.Date(Admit)))
  )
# Calculate the median, mean, and IQR of the LOS
summary.los <- summary(not.cdiff.los.data$LOS, na.rm = TRUE)
summary.los

not.cdiff.los.data <- not.cdiff.los.data %>% mutate("group" = "Other")

LOS.data <- rbind(cdiff.los.data, not.cdiff.los.data)

wilcox.test(cdiff.los.data$LOS, not.cdiff.los.data$LOS)

unique_diff_counts <- cdiff.descriptive3 %>%
  group_by(Short.ID) %>%
  filter(C.diff == 1) %>% 
  summarise(
    Unique_C.diff_Count = n_distinct(ToxinA)
  )
