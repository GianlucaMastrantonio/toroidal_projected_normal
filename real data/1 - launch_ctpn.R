name_sim_general <- "REAL"
m_mcmc <- 10
iter_mcmc <- 50000
thin_mcmc <- 30
burnin_mcmc <- 20000
batch_mcmc <- 20
a_mcmc <- 10000 / 2
b_mcmc <- 12000 / 2
alpha_target <- 0.4
par_adapt_mcmc <- 5
nu_app <- 2



for(iseed in 1:5)
{
  for(iind in 1:2)
  {
    name_sim <- name_sim_general
    args <- c(iseed, 1,1,1,1,40,iind,1,2,1)
    source("real data/ctpn.R", , echo = TRUE)
  }
}
