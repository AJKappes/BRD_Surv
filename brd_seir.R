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

##### global fixed params #####

beta_val <- 0.17
delta_val <- mean(c(.039, .01))
gamma_val <- 1/61.81
fixed_params <- c(beta = beta_val, delta = delta_val, gamma = gamma_val)

# initialize state variables and time bounds
# vary population and pen size to evaluate outcomes under different operation sizes
Npop <- 7963 # population, 2018 avg feedlot head capacity
pen_size <- 53 # 2018 avg feedlot pens
head_pen <- Npop / pen_size
t <- seq(0, 182, 1) # representative feedlot days on feed

e_init <- 53 # assume one calf exposed to BRD per pen
i_init <- 0.128*Npop # 0.128 BRD incidence rate pre-feedlot introduction
r_init <- 0
s_init <- Npop - sum(e_init, i_init, r_init)
init_vals <- c(S = s_init, E = e_init, I = i_init, R = r_init) / Npop

##### Model specification #####

# initialize proportion of genetic marker selected in each pen
m_vec <- c()
for (i in 1:pen_size) {
  
  u <- runif(1)
  m_vec[i] <- u*head_pen

}

# equate summation of independent genetic marker selection population for posterior
M <- sum(m_vec)

theta_draw <- function(s) {

  # posterior hyperparameters from pseudo obs, set up and draw
  # psuedo (a) avg number of successes, implied genetic marker selection calves
  # psuedo (b) number of failures, implied non-genetic marker selection calves
  a_param <- M + mean(m_vec)
  b_param <- sum(Npop - m_vec) + Npop - mean(m_vec)
  theta_vec <- rbeta(s, a_param, b_param)
  
  return(list(draw = theta_vec, a = a_param, b = b_param))
  
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

# posterior evaluation
Nden <- 1000
theta_dens <- data.table(theta = theta_draw(Nden)[['draw']])
ggplot(theta_dens, aes(theta)) +
  geom_density(alpha = 0.3, size = 0.2, color = 'coral', fill = 'coral')

# brd sim
Nsim <- 2
theta_vec <- c()
brd_sim <- list()
for (i in 1:Nsim) {
  
  theta <- theta_draw(1)[['draw']]
  theta_vec[i] <- theta
  sim <- data.table(ode(init_vals, seq(0, 10, 1), brd_seir, fixed_params))
  brd_sim <- append(brd_sim, list(sim))
  
  if (i == Nsim) {
    
    print('Simulation complete')
    
  }
  
}


