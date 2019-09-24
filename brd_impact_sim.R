total_sim <- 3
impact_vec <- c()

for (i in 1:total_sim) {
  
  n_env <- new.env()
  sys.source('brd_seir_sim.R', n_env)
  impact_vec[i] <- n_env$brd_impact
  cat('Impact Simulation', i, 'of', total_sim, 'complete\n')
  
}




