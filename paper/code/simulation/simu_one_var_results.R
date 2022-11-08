library(tidyverse)
theme_set(theme_bw())
df <- list()
for (i in 1:1000) {
  load(paste0('paper/results/simulation_one_var/simu_', i, '.RData'))
  df[[i]] <- simu_results
}
sim_results <- dplyr::bind_rows(df)
by_iter_summary <- sim_results %>%
  group_by(iter, model_unique, n_obs) %>%
  summarise(correct = mean(correct),
            ess = mean(ess), 
            theta = mean(theta), 
            ARI = mean(ARI)) %>%
  mutate(model_unique = as.factor(model_unique)) %>%
  mutate(`Number of observations per location` = as.factor(n_obs))
type_df <- data.frame('model_unique' = factor(1:5),
                      'model_label' = factor(c('IS','Gibbs',
                                               'P[I]', 'Independence', 
                                               'Independence~No~MRF'),
                                             levels = c('IS','Gibbs',
                                                        'P[I]', 'Independence', 
                                                        'Independence~No~MRF')))

type_df <- data.frame('model_unique' = factor(1:5),
                      'model_label' = factor(c('B) IS', 'A) Liang et al.',
                                               'P[I]', 'Independent~MRF', 
                                               'C) Independent'),
                                             levels = c('A) Liang et al.', 'B) IS',
                                                        'P[I]', 'Independent~MRF', 
                                                        'C) Independent')))


ggplot(mapping = aes(x = iter, y = correct, group = as.factor(model_label), 
                     color = model_label, linetype = model_label,
                     shape = model_label)) +
  geom_point(data = by_iter_summary %>% left_join(type_df) %>%filter(iter %%5 == 0,
                                                                     model_unique %in% c(1,2,5))) +
  geom_line(data = by_iter_summary %>% left_join(type_df) %>% filter(model_unique %in% c(1,2,5))) + 
  facet_wrap(~`Number of observations per location`, labeller = label_both, scales = 'free_y') + 
  labs(x = 'EM iteration', y = 'Percent of locations\nassigned in true cluster',
       color = 'Estimation Strategy',
       linetype = 'Estimation Strategy',  shape = 'Estimation Strategy') + 
  # scale_color_discrete(labels = scales::label_parse())+
  # scale_linetype_discrete(labels = scales::label_parse())+
  # scale_shape_discrete(labels = scales::label_parse())+
  theme(legend.text.align = 0, legend.position = 'bottom',text = element_text(size = 15))
ggsave(filename = 'paper/images/simu_one_var_percent_correct.png', height = 4, width =  10)

# type_df2 <- data.frame('model_unique' = factor(1:5),
#                       'model_label' = factor(c('Importance~Sampling','Gibbs',
#                                                'P[I]', 'Independence', 
#                                                'Independence~No~MRF'),
#                                              levels = c('Importance~Sampling','Gibbs',
#                                                         'P[I]', 'Independence', 
#                                                         'Independence~No~MRF')))
ggplot(mapping = aes(x = iter, y = correct, group = as.factor(model_label), 
                     color = model_label, linetype = model_label,
                     shape = model_label)) +
  geom_point(data = by_iter_summary %>% left_join(type_df) %>%
               filter(iter %%5 == 0, model_unique %in% c(1,2,4)),
             size = 2) +
  geom_line(data = by_iter_summary %>% left_join(type_df) %>% 
              filter(model_unique %in% c(1,2,4)), size = 1) + 
  facet_wrap(~`Number of observations per location`, labeller = label_both, scales = 'free_y') + 
  labs(x = 'EM iteration', y = 'Percent of locations\nassigned in true cluster',
       color = 'Estimation Strategy',
       linetype = 'Estimation Strategy',  shape = 'Estimation Strategy') + 
  # scale_color_discrete(labels = scales::label_parse())+
  # scale_linetype_discrete(labels = scales::label_parse())+
  # scale_shape_discrete(labels = scales::label_parse())+
  theme(legend.text.align = 0, legend.position = 'bottom',text = element_text(size = 18)) + 
  geom_vline(xintercept = 5, linetype = 'dotted')
ggsave(filename = 'paper/images/simu_one_var_percent_correct2.png', height = 4, width =  10)

ggplot(mapping = aes(x = iter, y = ess, group = as.factor(`Number of observations per location`), 
                     color = as.factor(`Number of observations per location`), 
                     linetype = as.factor(`Number of observations per location`),
                     shape = as.factor(`Number of observations per location`))) +
  geom_point(data = type_df %>% left_join(by_iter_summary) %>%#filter(iter %%5 == 0) %>%
               filter(model_unique == 1)%>%
               mutate(ess = ifelse(iter <= 5, NA, ess))) +
  geom_line(data = type_df %>% left_join(by_iter_summary)  %>%
              filter(model_unique == 1) %>%
              mutate(ess = ifelse(iter <= 5, NA, ess))) + 
  labs(x = 'EM iteration', y = 'Average Effective Sample\nSize of 30 MC iterations',
       color = 'Observations\nper location',
       linetype = 'Observations\nper location',  shape = 'Observations\nper location') + 
  scale_color_discrete(labels = scales::label_parse())+
  scale_linetype_discrete(labels = scales::label_parse())+
  scale_shape_discrete(labels = scales::label_parse())+
  scale_y_continuous(limits = c(0, 20)) + 
  theme(legend.text.align = 0, legend.position = 'bottom',text = element_text(size = 15))
ggsave(filename = 'paper/images/simu_one_var_ess.png', height = 4, width =  5)


ggplot(mapping = aes(x = iter, y = ARI, group = as.factor(model_label), 
                     color = model_label, linetype = model_label,
                     shape = model_label)) +
  geom_point(data = by_iter_summary %>% left_join(type_df) %>%filter(iter %%5 == 0,
                                                                     model_unique %in% c(1,2,5))) +
  geom_line(data = by_iter_summary %>% left_join(type_df) %>% filter(model_unique %in% c(1,2,5))) + 
  facet_wrap(~`Number of observations per location`, labeller = label_both, scales = 'free_y') + 
  labs(x = 'EM iteration', y = 'Adjusted Rand Index',
       color = 'Estimation Strategy',
       linetype = 'Estimation Strategy',  shape = 'Estimation Strategy') + 
  # scale_color_discrete(labels = scales::label_parse())+
  # scale_linetype_discrete(labels = scales::label_parse())+
  # scale_shape_discrete(labels = scales::label_parse())+
  theme(legend.text.align = 0, legend.position = 'bottom',text = element_text(size = 15))
ggsave(filename = 'paper/images/simu_one_var_ARI.png', height = 4, width =  10)

ggplot(mapping = aes(x = iter, y = theta, group = as.factor(model_label), 
                     color = model_label, linetype = model_label,
                     shape = model_label)) +
  geom_point(data = by_iter_summary %>% left_join(type_df) %>%filter(iter %%5 == 0)) +
  geom_line(data = by_iter_summary %>% left_join(type_df)) + 
  facet_wrap(~`Number of observations per location`, labeller = label_both, scales = 'free_y') + 
  labs(x = 'EM iteration', y = 'Percent of locations\nassigned in true cluster',
       color = 'Estimation Strategy',
       linetype = 'Estimation Strategy',  shape = 'Estimation Strategy') + 
  # scale_color_discrete(labels = scales::label_parse())+
  # scale_linetype_discrete(labels = scales::label_parse())+
  # scale_shape_discrete(labels = scales::label_parse())+
  theme(legend.text.align = 0, legend.position = 'bottom',text = element_text(size = 15))
