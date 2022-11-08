library(fda)
library(orthogonalsplinebasis)
library(Matrix)
library(tidyverse)
grid <- expand.grid('P'= seq(5, 200, by = 15), 'n' = seq(5, 400, by = 25))

n_entries <- rep(0, nrow(grid))
n_bspline <- rep(0, nrow(grid))
n_obspline <- rep(0, nrow(grid))
for (i in 1:nrow(grid)) {
  P <- grid[i,1]
  n  <- grid[i,2]
  
  knots <- seq(0, 1, length.out= P-2)
  
#  p <- runif(n)
  p <- seq(0, 1, length.out = n)
  p <- sort(p)
  
  B <- as(eval.basis(create.bspline.basis(knots), evalarg = p), 'sparseMatrix')
  B_tilde <- as(orthogonalsplinebasis::evaluate(OrthogonalSplineBasis(expand.knots(knots)), x =  p),
                'sparseMatrix')
  n_entries[i] <- prod(dim(B))
  n_bspline[i] <- length(B@i)
  n_obspline[i] <- length(B_tilde@i)
  if (i%% 100 == 0){
    print(i)
  }
}

df <- data.frame(grid, n_entries, n_bspline, n_obspline)
df_long <- tidyr::pivot_longer(df, cols = c(n_bspline, n_obspline), names_to = 'type', values_to = 'entries')

type_df <- data.frame(type = c('n_bspline', 'n_obspline'),
                      type_label = c('B-splines', 'Orthogonal B-splines'))

ggplot(df_long, aes(x = P, y = n, fill = entries))+
  geom_raster() + 
  facet_wrap(~type)
ggplot(df_long, aes(x = P, y = entries, group = n, color = n))+
  facet_wrap(~type) + 
  geom_line() 

ggplot(df_long %>%left_join(type_df), aes(x = n, y = entries, group = P, color = P))+
  facet_wrap(~type_label) + 
  geom_line() +
  theme_bw() +
  labs(x = 'Number of Measurements', y = 'Number of Non-Zero Entries', 
       color = 'Number of\nBasis Functions')+
  theme(legend.position = 'bottom')+
  scale_color_viridis_c(direction = -1)
ggsave(filename = 'paper/images/bspline_comparison.png', height = 3, width = 8)

ggplot(df_long %>%left_join(type_df), aes(x = n, y = entries/n_entries * 100, group = P, color = P))+
  facet_wrap(~type_label) + 
  geom_line() +
  theme_bw() +
  labs(x = 'Number of Measurements', y = 'Number of Non-Zero Entries', 
       color = 'Number of\nBasis Functions')+
  theme(legend.position = 'bottom')+
  scale_color_viridis_c(direction = -1)




