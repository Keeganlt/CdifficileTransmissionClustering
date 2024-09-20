source("code/GlobalFunctions.R")
source("code/LoadData.R")

## Days from pt positive to room positive ----------------------------------


# Filter only rows with patient samples (sample.env == 1) and positive samples (C.diff == 1)
positive_patient_samples <- cdiff.descriptive3 %>%
  filter(sample.env == "Patient", C.diff == 1)

# Group the data by short.ID (patient) and find the date of their first positive sample
first_positive_dates <- positive_patient_samples %>%
  group_by(Short.ID) %>%
  summarize(first_positive_date = min(as.Date(Collection.date)))

# Join the original data with the first_positive_dates data to calculate days since the first positive sample
cdiff.descriptive4 <- cdiff.descriptive3 %>%
  left_join(first_positive_dates, by = "Short.ID") %>%
  mutate(days_since_first_positive = as.numeric(as.Date(Collection.date) - first_positive_date))

positive_envHH_samples <- cdiff.descriptive4 %>%
  filter(sample.env != "Patient", C.diff == 1)

# Find the first positive environmental sample for each Short.ID
first_positive_env_samples <- positive_envHH_samples %>%
  group_by(Short.ID) %>%
  arrange(as.Date(Collection.date) ) %>%
  slice_head(n = 1) %>% 
  ungroup()


frequency_table <- first_positive_env_samples %>%
  count(days_since_first_positive, wt = !is.na(days_since_first_positive))

# Rename the columns for clarity
colnames(frequency_table) <- c("Days_Since_First_Positive", "Frequency")


# Filter rows for environmental samples
env_samples_neg <- cdiff.descriptive4 %>%
  filter(sample.env != "Patient", C.diff == 0)

env_samples_pos <- cdiff.descriptive4 %>%
  filter(sample.env != "Patient", C.diff == 1)

length(setdiff(env_samples_neg$Short.ID, env_samples_pos$Short.ID))


# Filter rows for patients samples
pat_samples_neg <- cdiff.descriptive4 %>%
  filter(sample.env == "Patient", C.diff == 0)

pat_samples_pos <- cdiff.descriptive4 %>%
  filter(sample.env == "Patient", C.diff == 1)

length(setdiff(pat_samples_neg$Short.ID, pat_samples_pos$Short.ID))

## Number of samples to first positive --------------------------------------------------

patient_samples <- cdiff.descriptive3 %>%
  filter(sample.env == "Patient")

# Group the data by short.ID (patient) and find the date of their first patient sample
first_pt_sample_dates <- patient_samples %>%
  group_by(Short.ID) %>%
  summarize(first_pt_sample_date = min(as.Date(Collection.date)))

# Join the original data with the first_pt_sample_dates data to calculate days since the first positive sample
cdiff.descriptive5 <- cdiff.descriptive3 %>%
  left_join(first_pt_sample_dates, by = "Short.ID") %>%
  mutate(days_since_first_pt_sample = as.numeric(as.Date(Collection.date) - first_pt_sample_date))

# Filter rows for patient samples (sample.env == "patient") and positive samples (C.diff == 1)

# Group data by Short.ID and arrange samples by Collection.date
patient_positive_samples <- cdiff.descriptive5 %>%
  filter(sample.env == "Patient", C.diff == 1) %>%
  group_by(Short.ID) %>%
  arrange(as.Date(Collection.date) ) %>%
  slice_head(n = 1)

frequency_table <- patient_positive_samples %>%
  ungroup() %>% 
  count(days_since_first_pt_sample, wt = !is.na(days_since_first_pt_sample))


# Group the data by short.ID (patient) and find the date of their first patient sample
first_kind_sample_dates <- cdiff.descriptive3 %>%
  group_by(Short.ID, sample.env) %>%
  summarize(first_kind_sample_date = min(as.Date(Collection.date)))

# Join the original data with the first_pt_sample_dates data to calculate days since the first positive sample
cdiff.descriptive6 <- cdiff.descriptive3 %>%
  left_join(first_kind_sample_dates, by = c("Short.ID","sample.env")) %>%
  mutate(days_since_first_kind_sample = as.numeric(as.Date(Collection.date) - first_kind_sample_date))

# Filter rows for patient samples (sample.env == "patient") and positive samples (C.diff == 1)

# Group data by Short.ID and arrange samples by Collection.date
all_positive_samples <- cdiff.descriptive6 %>%
  filter(C.diff == 1) %>%
  group_by(Short.ID) %>%
  arrange(as.Date(Collection.date) ) %>%
  slice_head(n = 1)

frequency_table <- all_positive_samples %>%
  ungroup() %>% 
  count(days_since_first_kind_sample, wt = !is.na(days_since_first_kind_sample))



colnames(patient_positive_samples)[22] <- "first_positive_date"
colnames(patient_positive_samples)[23] <- "days_since_first_positive"


all_first_positive <- rbind(first_positive_env_samples)


plot.2.4 <- ggplot(all_first_positive, aes(x = days_since_first_positive)) +
  geom_bar(aes(fill = sample.env), position = position_dodge2(width = 0.5, preserve = "single")) + 
  xlab("Days") +
  ylab("Count") +
  theme_cowplot(18) +
  scale_fill_manual(values =  c("#45b7c2", "#024b7a"),
                    name = "Sample Location") +
  theme(legend.position = "top") +
  background_grid()  +
  # xlim(-1,31)+ 
  theme(legend.position = "none",
        aspect.ratio = 0.6)

plot.2.4
