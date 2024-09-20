source("code/GlobalFunctions.R")
source("code/LoadData.R")


summarize_hospital_data <- function(data, hospital_id) {
  # Filter the data for the specified hospital
  hospital_data <- data %>% 
    filter(Hospital == hospital_id)
  
  # Create a list of unique PATIDs with C.diff == 1
  pat_w_cdiff <- hospital_data %>% 
    filter(C.diff == 1) %>% 
    distinct(PATID) %>%
    pull(PATID)
  
  # Convert PATID to numeric (removing "PAT-" prefix and "A" or other suffixes)
  hospital_data <- hospital_data %>%
    mutate(PATID_num = as.numeric(gsub("PAT-|[A-Za-z]", "", PATID)))
  
  # Create subsets based on PATID and whether they are in pat_w_cdiff
  pat_sampling_dat_cdiff <- hospital_data %>% filter(PATID_num < 200 & PATID %in% pat_w_cdiff)
  room_sampling_dat_cdiff <- hospital_data %>% filter(PATID_num >= 200 & PATID_num <= 900 & PATID %in% pat_w_cdiff)
  empty_sampling_dat_cdiff <- hospital_data %>% filter(PATID_num > 900 & PATID %in% pat_w_cdiff)
  
  pat_sampling_dat_nocdiff <- hospital_data %>% filter(PATID_num < 200 & !(PATID %in% pat_w_cdiff))
  room_sampling_dat_nocdiff <- hospital_data %>% filter(PATID_num >= 200 & PATID_num <= 900 & !(PATID %in% pat_w_cdiff))
  empty_sampling_dat_nocdiff <- hospital_data %>% filter(PATID_num > 900 & !(PATID %in% pat_w_cdiff))
  
  # Function to summarize data
  summarize_data <- function(data) {
    data %>%
      summarise(
        unique_patients = n_distinct(PATID),
        patient_samples = sum(sample.env == "Patient"),
        room_env_samples = sum(sample.env == "Room Env."),
        hcp_hands_samples = sum(sample.env == "HCW Hand"),
        total_samples = n()
      )
  }
  
  # Summarize each subset
  pat_sampling_summary_cdiff <- summarize_data(pat_sampling_dat_cdiff) %>% mutate(PATID_range = "<200", Cdiff_status = "Yes")
  room_sampling_summary_cdiff <- summarize_data(room_sampling_dat_cdiff) %>% mutate(PATID_range = "200-900", Cdiff_status = "Yes")
  empty_sampling_summary_cdiff <- summarize_data(empty_sampling_dat_cdiff) %>% mutate(PATID_range = ">900", Cdiff_status = "Yes")
  
  pat_sampling_summary_nocdiff <- summarize_data(pat_sampling_dat_nocdiff) %>% mutate(PATID_range = "<200", Cdiff_status = "No")
  room_sampling_summary_nocdiff <- summarize_data(room_sampling_dat_nocdiff) %>% mutate(PATID_range = "200-900", Cdiff_status = "No")
  empty_sampling_summary_nocdiff <- summarize_data(empty_sampling_dat_nocdiff) %>% mutate(PATID_range = ">900", Cdiff_status = "No")
  
  # Combine the summaries into a single data frame
  summary_combined <- bind_rows(
    pat_sampling_summary_cdiff,
    room_sampling_summary_cdiff,
    empty_sampling_summary_cdiff,
    pat_sampling_summary_nocdiff,
    room_sampling_summary_nocdiff,
    empty_sampling_summary_nocdiff
  )
  
  # Return the combined summary
  return(summary_combined)
}

# Example usage for Hospital A and B
hospital_A_summary <- summarize_hospital_data(gm.metadata, "A")
hospital_B_summary <- summarize_hospital_data(gm.metadata, "B")

# Print the summaries
print(hospital_A_summary)
print(hospital_B_summary)



summarize_hospital_shared_data <- function(data, hospital_id) {
  # Filter the data for the specified hospital
  hospital_data <- data %>% 
    filter(Hospital == hospital_id, 
           sample.env == "Shared Env.")
  
  sumary <- hospital_data %>%
    summarise(
      total_samples = n()
    )
  return(sumary)
}
# Example usage for Hospital A and B
hospital_A_shared <- summarize_hospital_shared_data(gm.metadata, "A")
hospital_B_shared <- summarize_hospital_shared_data(gm.metadata, "B")

