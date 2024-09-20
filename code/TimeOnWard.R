source("code/GlobalFunctions.R")
source("code/LoadData.R")

## Stats on time on the ward  ---------------------------------------------
time.on.ward.calcs <- cdiff.descriptive %>% 
  distinct(Short.ID, Admit, Disch) %>% 
  mutate(time.on.ward = Disch -Admit) %>%
  na.omit() %>%
  summarise(med.time = median(time.on.ward), mean.time = mean(time.on.ward), 
            st.dev.time = sd(time.on.ward))
sample.mean <- as.numeric(time.on.ward.calcs$mean.time)
sample.sd <-time.on.ward.calcs$st.dev.time


## Merge with color shape 
merged.cdiff.descriptive = merge(cdiff.descriptive, 
                                 global.color.shape, 
                                 by.x = c("Short.ID"),
                                 by.y = c("Room.Loc"),
                                 all = TRUE)

## list of unique patients 
pat.id<-unique(dist.lower.meta.tall.all.dup$Short.ID1)
sorted.id <- c(sort(unique(cdiff.metadata$Short.ID)), "Reference")


