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

write.csv(iso, file = 'data/darros-SIA/darros_basal_SIA_2025.csv', row.names=FALSE)


ggplot(iso, aes(tax, x13c, fill=decalcified)) + geom_boxplot()
ggplot(iso, aes(tax, x15n, fill=decalcified)) + geom_boxplot()

dc_diff<-iso %>% 
  group_by(tax, unique_sample_id, decalcified) %>% 
  summarise(x13c = mean(x13c)) %>% ungroup() %>% 
  group_by(unique_sample_id) %>% mutate(n = length(x13c)) %>% 
  filter(n == 2) %>% 
  pivot_wider(id_cols=-c(n), names_from = decalcified, values_from = x13c) %>% 
  mutate(diff = yes - no, prop = diff / no * 100) %>% 
  filter(!is.na(diff))

dc_diff %>% group_by(tax) %>% summarise(mean(prop))

ggplot(dc_diff, aes(unique_sample_id, no)) + geom_point() +
  geom_point(aes(unique_sample_id, yes), col='red') +
  facet_wrap(~tax, scales='free_x') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90)) +
  labs(x = '', y='13C')

ggplot(iso, aes(site, x13c)) + geom_point() + facet_wrap(~tax, scales='free_x')
ggplot(iso, aes(site, x15n)) + geom_point() + facet_wrap(~tax, scales='free_x')

iso_mean<-iso %>% group_by(tax) %>% 
  summarise(x13c = mean(x13c), x15n = mean(x15n))

pdf(file = 'darros/basal_iso_niche.pdf', height=5, width=8)
ggplot(iso %>% filter(x15n > 0), aes(x13c, x15n, fill=tax, col=tax, group=tax)) + 
  geom_point(alpha=0.5, size=1.2) + 
    stat_ellipse(geom='polygon', type='norm', level=0.4, alpha=0.5) +
  ggrepel::geom_text_repel(data = iso_mean, aes(label = tax, col=tax),
                           show.legend = FALSE, size = 3) +
  theme_bw() + guides(fill='none', colour='none')
dev.off()
