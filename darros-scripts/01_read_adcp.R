library(tidyverse)
library(janitor)

adcp<-rbind(
  read.csv('data/darros-ADCP/Eco358_202511180834.csv', skip=28) %>% mutate(site = 'Manta'),
  read.csv('data/darros-ADCP/Eco773_202511180838.csv', skip=28) %>% mutate(site = 'Channel'))

heads<-colnames(adcp %>% clean_names())
units<-adcp[1,] %>% str_replace_all('/', '_')

# drop units
adcp<-adcp %>% filter(!Depth == 'm')

# add colnames with units
colnames(adcp)<-paste(heads, units, sep='_')
colnames(adcp)[11]<-'site'

# add date and time
adcp<-adcp %>% mutate(time = strptime(.[[1]], format="%d/%m/%Y %H:%M:%S"),
                       date = dmy(str_split_fixed(.[[1]], ' ', 2)[,1]),
                       time24 = hms(str_split_fixed(.[[1]], ' ', 2)[,2])) %>% 
                select(-1) %>% 
                mutate(across(c(upper_speed_m_s:tilt_deg), as.numeric)) %>% 
              # drop unit on land
                filter(depth_m > 16)


## estimate degree cooling



