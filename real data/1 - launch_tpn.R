# Launch script for the real-data application with the Toroidal Projected
# Normal model. The script runs multiple seeds and both covariance structures.

# Prefix used in the output file names.
name_sim_general <- "REAL"

# MCMC settings shared by all real-data runs.
m_mcmc <- 10
iter_mcmc <- 50000
thin_mcmc <- 30
burnin_mcmc <- 20000

# Adaptive MCMC settings.
batch_mcmc <- 20
a_mcmc <- 10000 / 2
b_mcmc <- 12000 / 2
alpha_target <- 0.4
par_adapt_mcmc <- 5
nu_app <- 2


# The args vector selects:
#   1. seed
#   2. initialization setting
#   3. subset/data setting
#   4. ESS switch
#   5. ESS type
#   6. number of covariance proposal checks
#   7. independence setting
#   8. MCMC multiplier
#   9. mixture switch
#  10. maximum number of mixture components
for(iseed in 1:5)
{
  for(iind in 1:2)
  {
    # Reset the output-name prefix before each run; the sourced script appends
    # the model and scenario information.
    name_sim <- name_sim_general
    args <- c(iseed, 1,1,1,1,40,iind,1,2,1)
    source("real data/tpn.R", , echo = TRUE)
  }
}
