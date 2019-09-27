Nsim <- 1000
theta25p_impact <- c()
theta75p_impact <- c()
for (i in 1:Nsim) {
  
  n_env <- new.env()
  sys.source('brd_seir_sim.R', n_env)
  theta25p_impact[i] <- n_env$avg_cost_25perc
  theta75p_impact[i] <- n_env$avg_cost_75perc
  cat('Impact simulation', i, 'of', Nsim, 'complete\n')
  
}

brd_impact_dt <- data.table(c25 = theta25p_impact,
                            c75 = theta75p_impact,
                            impact = theta25p_impact - theta75p_impact)
fwrite(brd_impact_dt, 'brd_impact.csv')

#### ####
remove(list = objects())
library(data.table)
library(ggplot2)
library(latex2exp)

impact_dt <- fread('brd_impact.csv')

impact_rvs <- rnorm(1000,
                     mean(impact_dt$c25) - mean(impact_dt$c75),
                     sqrt((var(impact_dt$c25) + var(impact_dt$c75))/nrow(impact_dt)))
impact_dist_dt <- data.table(rv = impact_rvs)

ggplot(impact_dist_dt, aes(rv)) +
  geom_density(alpha = 0.4, size = 0.1, color = '#978B88', fill = '#978B88')

impact_conv_vec <- c()
impact_var_vec <- c()
impact_upper <- c()
impact_lower <- c()
crit <- qt(1-.05/2, nrow(impact_dt))
for (i in 1:nrow(impact_dt)) {
  
  impact_conv_vec[i] <- mean(impact_dt$impact[1:i])
  impact_var_vec[i] <- var(impact_dt$impact[1:i])
  impact_upper[i] <- impact_conv_vec[i] + crit*sqrt(impact_var_vec[i])
  impact_lower[i] <- impact_conv_vec[i] - crit*sqrt(impact_var_vec[i])
  
}

conv_dt <- data.table(conv = impact_conv_vec,
                      upper = impact_upper,
                      lower = impact_lower,
                      itr = 1:nrow(impact_dt))

ggplot(conv_dt, aes(itr)) +
  geom_line(aes(y = conv, colour = 'mean'), alpha = 0.75) +
  geom_line(aes(y = upper, colour = 'upper')) +
  geom_line(aes(y = lower, colour = 'lower')) +
  scale_color_manual(name = '',
                     values = c('mean' = 'black', 'upper' = 'gray', 'lower' = 'gray')) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray', alpha = 0.4) +
  ylab(TeX('$\\mu$')) +
  xlab(TeX('Iteration')) +
  ggtitle(TeX('Expected Cost Impact Convergence')) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none')
  
                     
conv_dt$upper[nrow(conv_dt)]
conv_dt$lower[nrow(conv_dt)]



