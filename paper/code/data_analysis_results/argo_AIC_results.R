library(tidyverse)
theme_set(theme_bw())
G <- 5
n_pcs <- c(5, 8:13)
test_pcs <- 5:12
data_df <- list()
for (q in n_pcs) {
  for (q2 in test_pcs) {
    load(paste0('paper/results/oxy_none_NA_', G, '_', q, '_12_15_15/AIC_', q2, '.RData'))
    data_df <- c(data_df, list(data.frame('orig_pc' = q, 'test_pc' = q2, 'AIC' = model$AIC, 'BIC' = model$BIC)))
  }
}

data_df <- dplyr::bind_rows(data_df)
ggplot(data = data_df, aes(x = test_pc, y = AIC, color = factor(orig_pc))) + 
  geom_point() + 
  geom_line() + 
  labs(x = 'Number of response PCs', 
       color = 'Number of\npredictor\nPCs', title = paste0('AIC results with ', G, ' Clusters'))
ggsave(paste0('paper/results/oxy_none_NA_', G, '_12_12_15_15/AIC_plot.png'), 
       width = 7, height = 4)

ggplot(data = data_df, aes(x = test_pc, y = AIC, color = factor(orig_pc))) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~factor(orig_pc), scales = 'free_y', ncol = 2) +
  labs(x = 'Number of response PCs', 
       color = 'Number of\npredictor\nPCs', title = paste0('AIC results with ', G, ' Clusters'))
ggsave(paste0('paper/results/oxy_none_NA_', G, '_12_12_15_15/AIC_scale_plot.png'), 
       width = 5, height = 6)
