source("code/GlobalFunctions.R")
source("code/LoadData.R")

## Stats on toxin carriage  ------------------------------------------------

toxin.by.cdiff <- cdiff.descriptive2 %>% 
  group_by(C.diff) %>% 
  count(ToxinA)

tox.dat <- cdiff.descriptive2 %>% 
  filter(C.diff == 1) %>% 
  group_by(Short.ID, sample.env) %>% 
  #distinct(Short.ID, Cdiff.Toxin) %>% 
  count(ToxinA)

all.pats.samples.count <- cdiff.descriptive2 %>% 
  filter(sample.env == "Patient") %>%
  group_by(Short.ID) %>% 
  count(Short.ID)
summary(all.pats.samples.count)



toxin.people <- tox.dat %>% filter(sample.env == "Patient") 
tox.neg.people <-tox.dat %>% filter(sample.env == "Patient" & is.na(ToxinA)) 
tox.pos.people <- tox.dat %>% filter(sample.env == "Patient" & ToxinA == "Positive") 

toxin.both <- tox.pos.people$Short.ID[which(tox.pos.people$Short.ID %in% tox.neg.people$Short.ID)]


tox.pos.people <- toxin.people %>% filter(Short.ID %!in% c(toxin.both, unique(tox.neg.people$Short.ID)))
tox.neg.people <- toxin.people %>% filter(Short.ID %!in% c(toxin.both, unique(tox.pos.people$Short.ID)))
toxin.both.people <- toxin.people %>% filter(Short.ID %in% toxin.both)


toxin.not.people <- tox.dat %>% filter(sample.env != "Patient") %>% filter(Short.ID != "ASN") 
length(unique(toxin.not.people$Short.ID))


toxin.pos.not.people <- toxin.not.people %>% filter(ToxinA == "Positive") 
toxin.neg.not.people <- toxin.not.people %>% filter(is.na(ToxinA))

toxin.both.not <- toxin.neg.not.people$Short.ID[which(toxin.neg.not.people$Short.ID %in% toxin.pos.not.people$Short.ID)]

tox.pos.not.people <- toxin.not.people %>% filter(Short.ID %!in% c(toxin.both, unique(toxin.neg.not.people$Short.ID)))
tox.neg.not.people <- toxin.not.people %>% filter(Short.ID %!in% c(toxin.both, unique(toxin.pos.not.people$Short.ID)))
toxin.both.not.people <- toxin.not.people %>% filter(Short.ID %in% toxin.both.not)

#### ~~~~ Toxin status ------

num.samps.pos <- cdiff.descriptive %>% 
  filter(Short.ID %in% tox.pos.people$Short.ID, 
         sample.loc == "PAT") %>%
  group_by(Short.ID) %>%
  summarise(num_samples = n()) %>% 
  summarise(med.samps = median(num_samples), mean.samps = mean(num_samples), 
            st.dev.samps = sd(num_samples), 
            iqr.samps = IQR(num_samples))

num.samps.neg <- cdiff.descriptive %>% 
  filter(Short.ID %in% tox.neg.people$Short.ID, 
         sample.loc == "PAT") %>%
  group_by(Short.ID) %>%
  summarise(num_samples = n()) %>% 
  summarise(med.samps = median(num_samples), mean.samps = mean(num_samples), 
            st.dev.samps = sd(num_samples), 
            iqr.samps = IQR(num_samples))

num.samps.both <- cdiff.descriptive %>% 
  filter(Short.ID %in% toxin.both.people$Short.ID, 
         sample.loc == "PAT") %>%
  group_by(Short.ID) %>%
  summarise(num_samples = n())%>% 
  summarise(med.samps = median(num_samples), mean.samps = mean(num_samples), 
            st.dev.samps = sd(num_samples), 
            iqr.samps = IQR(num_samples))

num.samps.env.pos <- cdiff.descriptive %>% 
  filter(Short.ID %in% tox.pos.people$Short.ID, 
         sample.env == "Room Env.") %>%
  group_by(Short.ID) %>%
  summarise(num_samples = n()) %>% 
  summarise(med.samps = median(num_samples), mean.samps = mean(num_samples), 
            st.dev.samps = sd(num_samples), 
            iqr.samps = IQR(num_samples))

num.samps.env.neg <- cdiff.descriptive %>% 
  filter(Short.ID %in% tox.neg.people$Short.ID, 
         sample.env == "Room Env.") %>%
  group_by(Short.ID) %>%
  summarise(num_samples = n()) %>% 
  summarise(med.samps = median(num_samples), mean.samps = mean(num_samples), 
            st.dev.samps = sd(num_samples), 
            iqr.samps = IQR(num_samples))

num.samps.env.both <- cdiff.descriptive %>% 
  filter(Short.ID %in% toxin.both.people$Short.ID, 
         sample.env == "Room Env.") %>%
  group_by(Short.ID) %>%
  summarise(num_samples = n())%>% 
  summarise(med.samps = median(num_samples), mean.samps = mean(num_samples), 
            st.dev.samps = sd(num_samples), 
            iqr.samps = IQR(num_samples))
