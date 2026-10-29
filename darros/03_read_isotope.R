library(tidyverse)

iso<-read.csv('data/darros-SIA/darros_isotope_results.csv') %>% 
  select(-X, -X.1, -X.2) %>% 
  janitor::clean_names() %>% 
  filter(sample_id != 'S1303TC rr') %>% # replicate?
  mutate(sample_id = ifelse(sample_id == 'S040302ZPHG', 'S0403ZPHG', sample_id),
         sample_id = ifelse(sample_id == 'S0401DETCDC', 'S0401DETDC', sample_id)) # fix label transcribe error (CDC)

meta<-read.csv('data/darros-SIA/darros_isotope_samples.csv') %>% 
  select(unique_sample_id:sampling_date) %>% 
  mutate(sample_id = str_replace_all(tube_id, '-', ''))

iso <- iso %>% left_join(meta)

write.csv(iso, file = 'data/darros-SIA/basal_SIA_2025.csv', row.names=FALSE)
