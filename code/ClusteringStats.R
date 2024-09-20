source("code/ClusteringData.R")
## Probability of being in a cluster number of samples ---------------------

# Calculate the duration of time on the ward (in days) for each Short.ID
summary.cdiff<- cdiff.descriptive3 %>%
  group_by(Short.ID) %>%
  summarize(
    Time_on_Ward = as.numeric(max(Disch) - min(Admit)) + 1,
    Num_Positive_Samples = sum(C.diff == 1),
    Num_Total_Samples = n()
  )


all.cluster.shortIDs


summary.cdiff2 <- summary.cdiff %>%
  mutate(Clustering = ifelse(Short.ID %in% all.cluster.shortIDs, 1, 0))


mod1 = glm(Clustering~Num_Positive_Samples,data=summary.cdiff2,family="binomial")

mod2=glm(Clustering~Num_Positive_Samples + Time_on_Ward,data=summary.cdiff2,family="binomial")

mod3=glm(Clustering~Num_Positive_Samples*Time_on_Ward,data=summary.cdiff2,family="binomial")

value_of_interest <- 2
coefficients <- coef(mod2)
# Calculate the odds
intercept <- coefficients[1]  # Intercept coefficient
coefficient_num_positive_samples <- coefficients[2]  # Coefficient for Num_Positive_Samples
odds <- exp(intercept + coefficient_num_positive_samples * value_of_interest)

# Print the odds
print(odds)

## Logistic regression

## Cluster timing -------

timing.clusters <- data.frame(
  Cluster = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5, 6, 7),
  Pat1 = c("A051", "A056", "A273", "B211", "A239", "A012", "A250", "A046", "A004", "A058", "A279", "B009", "B012"),
  Pat2 = c("A046", "A046", "A046", "A034", "A034", "A034", "A034", "A013", "A013", "A007", "A010", "B008", "B018"),
  Pat1Loc = c("P", "H", "H", "H", "E", "E", "E", "H", "H", "E", "H", "P", "E"),
  time.diff =rep(NA, 13)
)

for(i in 1:nrow(timing.clusters)) {
  node.dat <- dat %>% 
    filter(cluster == timing.clusters$Cluster[i]) %>% 
    select(from, short.ID) %>% 
    distinct(from, short.ID) %>% 
    mutate(sample.loc = from) 
  
  short.name.1 <- unique(node.dat$from)
  
  ## Make the dataset we need
  cdiff.descriptive.1 <- cdiff.descriptive %>% 
    # filter(str_detect(Collection.Sample.ID, "VA")) %>% 
    filter(Collection.Sample.ID %in% short.name.1, negate = TRUE)
  
  cdiff.descriptive.1 <- cdiff.descriptive.1 %>% 
    filter(Short.ID %in% c(timing.clusters$Pat1[i],timing.clusters$Pat2[i]))
  
  Study.Day.1 <- cdiff.descriptive.1 %>% 
    filter(Short.ID == timing.clusters$Pat1[i]) %>% 
    select(Study.Day) %>% 
    slice_head(n = 1)
  
  cdiff.descriptive.2 <- cdiff.descriptive.1 %>% 
    filter(Short.ID %in% c(timing.clusters$Pat2[i]))
  
  closest.day <- cdiff.descriptive.2$Study.Day[which.min(abs(cdiff.descriptive.2$Study.Day - as.numeric(Study.Day.1)))]
  
  timing.clusters$time.diff[i] <- closest.day - as.numeric(Study.Day.1)
  
}

timing.clusters %>%
  group_by(Pat1Loc) %>% 
  summarise(mean = mean(abs(time.diff)), 
            median = median(abs(time.diff))
  ) 

