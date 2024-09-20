source("code/GlobalFunctions.R")
source("code/LoadData.R")


# Subset the data
room_data <- gm.metadata[gm.metadata$sample.env == "Room Env.", ]
patient_data <- gm.metadata[gm.metadata$sample.env == "Patient", ]
hcp_data <- gm.metadata[gm.metadata$sample.env == "HCW Hand", ]

# Calculate the number of occurrences where Cdiff.Toxin is 1
room_cdiff_1 <- sum(room_data$C.diff == 1, na.rm = TRUE)
patient_cdiff_1 <- sum(patient_data$C.diff == 1, na.rm = TRUE)
hcp_cdiff_1 <- sum(hcp_data$C.diff == 1, na.rm = TRUE)

# Calculate the odds ratio
odds_ratio <- room_cdiff_1 / patient_cdiff_1
odds_ratio2 <- patient_cdiff_1 / room_cdiff_1


# Subset the data to include only relevant columns
subset_data <- gm.metadata %>%
  filter(as.numeric(gsub("[^0-9]", "", PATID)) <= 200)%>%
  select(C.diff, sample.env) %>%
  filter(!is.na(C.diff)) # Remove rows with NA values in Cdiff.Toxin

# Convert sample.env to factor
subset_data$sample.env <- as.factor(subset_data$sample.env)

# Change the reference level of sample.env to "HCW hand"
subset_data$sample.env <- relevel(subset_data$sample.env, ref = "HCW Hand")

# Fit logistic regression model
logit_model <- glm(C.diff ~ sample.env, family = "binomial", data = subset_data)

# Print summary of the model
summary(logit_model)
# Get coefficient estimates and standard errors from the summary of the model
coefficients <- coef(summary(logit_model))
std_errors <- coefficients[, "Std. Error"]

# Exponentiate coefficients
exp_coef <- exp(coefficients[, "Estimate"])

# Calculate z-value for desired level of confidence (e.g., 95% confidence interval)
z_value <- qnorm(0.975)  # For a two-tailed test, dividing alpha by 2

# Calculate lower and upper bounds of the confidence interval for each exponentiated coefficient
lower_bound <- exp_coef - z_value * std_errors
upper_bound <- exp_coef + z_value * std_errors

# Combine lower and upper bounds into a matrix
CI <- cbind(lower_bound, upper_bound)

# Print confidence intervals
print(CI)

subset_data <- gm.metadata %>%
  filter(as.numeric(gsub("[^0-9]", "", PATID)) <= 200) %>%
  filter(Short.ID != "B003") %>%
  group_by(PATID) %>%
  mutate(
    positive_C.diff = ifelse(any(C.diff == 1 & sample.env == "Patient"), "pos.pat", "neg.pat"),
    positive_room = ifelse(any(C.diff == 1 & sample.env == "Room Env."), "pos.room", "neg.room"),
    positive_HCP_hand = ifelse(any(C.diff == 1 & sample.env == "HCW Hand"), "pos.hcw", "neg.hcw"),
    positive_combined = paste0(positive_C.diff, ".", positive_room, ".", positive_HCP_hand)
  ) %>%
  ungroup()


unique_df <- subset_data %>%
  distinct(PATID, .keep_all = TRUE)

odds.table <- table(unique_df$positive_C.diff, unique_df$positive_HCP_hand)

odds.table.df <- data.frame(neg.HCP = odds.table[,1], pos.HCP = odds.table[,2],
                                       row.names = c("neg.patient","pos.patient"))

#mosaicplot(odds.table)
fisher.test(odds.table)

odds.table <- table(unique_df$positive_C.diff, unique_df$positive_HCP_hand)

odds.table.df <- data.frame(neg.room = odds.table[,1], pos.room = odds.table[,2],
                            row.names = c("neg.patient","pos.patient"))

#mosaicplot(odds.table)
fisher.test(odds.table)

