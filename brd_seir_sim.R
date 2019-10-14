# SEIR BRD Simulation using Genetic Marker Selection Bayesian Methods

# susceptible, exposed, infectious, recovered
# N total population = sum(s,e,i,r)

# parameters informed by Snowder, Van Vleck, Cundiff et al. (2006)

# beta: infectious rate, transmitting disease between S and E, .17 feedlot incidence rate

# sigma: incubation rate, rate of latent subjects becoming infectious 
#        average duration of incubation is 1/sigma
#        characterized as 1-theta

# theta: associated brd resistance rate from genetic marker selection
#        random variable with conditional distribution on m, marker selection
#        prior on theta assumed beta
#        likelihood on m assumed binomial for ms successes in cow-calf population
#        likelihood parameter theta implies probability of cow-calf operation adopting genetic marker selection

# gamma: recovery rate 1/D, determined by avg duration D of infection

# delta: death rate .039

remove(list = objects())
library(deSolve)
library(data.table)
library(ggplot2)
library(latex2exp)

##### global fixed params #####

# vary population and pen size to evaluate outcomes under different operation sizes
Npop <- 150
pen_size <- 1
#Npop <- 7963 # population, 2018 avg feedlot head capacity
#pen_size <- 53 # 2018 avg feedlot pens
head_pen <- Npop / pen_size
delta_val <- mean(c(.039, .01))
gamma_val <- 1/61.81
beta_val <- 0.17*0.25*head_pen
fixed_params <- c(beta = beta_val, delta = delta_val, gamma = gamma_val)
t <- seq(0, 182, 1) # representative feedlot days on feed

# initialize state variables and time bounds

e_init <- 53 # assume one calf exposed to BRD per pen
i_init <- 0.128*Npop # 0.128 BRD incidence rate pre-feedlot introduction
r_init <- 0
s_init <- Npop - sum(e_init, i_init, r_init)
init_vals <- c(S = s_init, E = e_init, I = i_init, R = r_init) / Npop

##### Model specification #####

m_draw <- function() {
  
  m_vec <- c()
  for (i in 1:pen_size) {
    
    m_vec[i] <- runif(1, 0, head_pen)
    
  }
  
  M <- sum(m_vec)
  return(list(M = sum(m_vec), draws = m_vec))
  
}

theta_draw <- function(s) {

  # posterior hyperparameters from pseudo obs, set up and draw
  # psuedo (a) avg number of successes, implied genetic marker selection calves
  # psuedo (b) number of failures, implied non-genetic marker selection calves
  
  m_rv <- m_draw()
  a_param <- m_rv[['M']] + mean(m_rv[['draws']])
  b_param <- sum(Npop - m_rv[['draws']]) + Npop - mean(m_rv[['draws']])
  theta_vec <- rbeta(s, a_param, b_param)
  
  return(list(draws = theta_vec, a = a_param, b = b_param, m = m_rv[['draws']]))
  
}

brd_seir <- function(time, state, parameters) {
  
  S <- state[['S']]
  E <- state[['E']]
  I <- state[['I']]
  R <- state[['R']]
  
  with(as.list(c(state, parameters)), {
    
    # specify differential equations
    dsdt <- -beta*S*I
    dedt <- beta*S*I - (1-theta)*E - theta*E
    didt <- (1-theta)*E - gamma*I - delta*I
    drdt <- gamma*I + theta*E
    
    list(c(dsdt, dedt, didt, drdt), theta = theta)
    
  })
  
}

##### simulation #####

# brd sim
Nsim <- 500
theta_vec <- c()
# m_prop_vec <- c()
brd_sim <- list()
for (i in 1:Nsim) {
  
  # post_draw <- theta_draw(1)
  # m_prop_vec[i] <- post_draw[['m']] / Npop
  # theta <- post_draw[['draws']]
  theta <- theta_draw(1)[['draws']]
  theta_vec[i] <- theta
  sim <- data.table(ode(init_vals, t, brd_seir, fixed_params))
  brd_sim <- append(brd_sim, list(sim))
  
  if (i == Nsim) {
    
    print('SEIR simulation complete')
    
  }
  
}

# posterior evaluation with density plot

# theta_dens <- data.table(theta = theta_vec)
# 
# ggplot(theta_dens, aes(theta)) +
#   geom_density(alpha = 0.3, size = 0.2, color = 'coral', fill = 'coral')

#### feedlost infectious cost impact ####

# seperate 25th and 75th theta percentiles for simulation evaluation

theta_percs <- function(ineq, q) {
  
  if (ineq == 'leq') {
    
    theta_perc_idx <- which(theta_vec <= quantile(theta_vec, c(q)))
    brd_sim_perc <- list()
    
    for (i in theta_perc_idx) {
      
      brd_sim_perc <- append(brd_sim_perc, list(brd_sim[[i]]))
      
    }
    
  }
  
  if(ineq == 'geq') {
    
    theta_perc_idx <- which(theta_vec >= quantile(theta_vec, c(q)))
    brd_sim_perc <- list()
    
    for (i in theta_perc_idx) {
      
      brd_sim_perc <- append(brd_sim_perc, list(brd_sim[[i]]))
      
    }
    
  }
  
  return(brd_sim_perc)
  
}

brd_sim25perc <- theta_percs(ineq = 'leq', q = .25)
brd_sim75perc <- theta_percs(ineq = 'geq', q = .75)

# compare feedlot cost impacts based on brd resistance rate percentiles
# brd costs occur during period spanning infectious population
# estimated $13.9 loss per animal associated with lower gains and treatment costs
# 31 profit per head based on previous monte carlo study
cpd <- 13.9 / max(t)
plph <- 31.65

cost_vec_25perc <- c()
cost_vec_75perc <- c()
for (i in 1:length(brd_sim25perc)) {
  
  cost_vec_25perc[i] <- (sum(cpd*brd_sim25perc[[i]]$I*Npop) +
                           plph*(1-max(brd_sim25perc[[i]]$R))*Npop)
  
  cost_vec_75perc[i] <- (sum(cpd*brd_sim75perc[[i]]$I*Npop) +
                           plph*(1-max(brd_sim75perc[[i]]$R))*Npop)
  
}

avg_cost_25perc <- mean(cost_vec_25perc)
avg_cost_75perc <- mean(cost_vec_75perc)
# impact <- avg_cost_25perc - avg_cost_75perc
# cat(avg_cost_25perc, '-', avg_cost_75perc, '=', impact)
 
#### epi transition plot ####

# theta 25th percentile plot

plot_fun <- function(ineq, sim, y) {
  
  t <- seq(0, 182, 1)
  
  if (ineq == 'leq') {

    sub <- TeX('$\\bar{\\theta}_{25^{th} percentile}$ Associated BRD Resistance Rate')
    
  }
  
  if (ineq == 'geq') {

    sub <- TeX('$\\bar{\\theta}_{\\geq 75^{th} percentile}$ Associated BRD Resistance Rate')
    
  }
  
  plot <- ggplot(sim, aes(t)) +
    geom_line(aes(y = y$S, colour = 'Susceptible')) +
    geom_line(aes(y = y$E, colour = 'Exposed')) +
    geom_line(aes(y = y$I, colour = 'Infectious')) +
    geom_line(aes(y = y$R, colour = 'Recovered')) +
    ylim(-.05, 1) +
    labs(title = 'BRD Transition plot',
         subtitle = sub,
         x = TeX('$t$'),
         y = TeX('Population Proportion')) +
    scale_color_manual(name = '',
                       values = c('Susceptible' = 'blue', 'Exposed' = 'orange',
                                  'Infectious' = 'red', 'Recovered' = 'green'),
                       limits = c('Susceptible', 'Exposed', 'Infectious', 'Recovered')) +
    theme(plot.title = element_text(family = 'serif', hjust = 0.5),
          plot.subtitle = element_text(family = 'serif', hjust = 0.5),
          axis.title = element_text(family = 'serif'),
          legend.text = element_text(family = 'serif'))
  
  return(plot)
  
}

theta <- mean(theta_vec[theta_vec <= quantile(theta_vec, c(.25))])
sim_25th <- data.table(ode(init_vals, t, brd_seir, fixed_params))
y_25th <- sim_25th[, .(S, E, I, R)]
plot_25th <- plot_fun('leq', sim_25th, y_25th)

theta <- mean(theta_vec[theta_vec >= quantile(theta_vec, c(.75))])
sim_75th <- data.table(ode(init_vals, t, brd_seir, fixed_params))
y_75th <- sim_75th[, .(S, E, I, R)]
plot_75th <- plot_fun('geq', sim_75th, y_75th)
