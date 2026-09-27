name_sim_general <- "Simulation"
m_mcmc <- 10
iter_mcmc <- 50000
thin_mcmc <- 30
burnin_mcmc <- 20000
batch_mcmc <- 20
a_mcmc <- 10000 / 2
b_mcmc <- 12000 / 2
alpha_target <- 0.4
par_adapt_mcmc <- 5000
nu_app <- 2



for (id_for in 1:4)
{
  for (ik_for in 4:4)
  {
    for (isigma_for in 1:2)
    {
      for (ichain_for in 1:5)
      {
        for (isigma2_for in 1:1)
        {
          for (idata_for in 1:25)
          {
            for (iin_for in 1:2)
            {
              for (iinit_for in 2:2)
              {
                for (iess_for in 1:1)
                {
                  for (typeess_for in 1:1)
                  {
                    for (ntry_for in 40:40)
                    {
                      name_sim <- name_sim_general
                      args <- c(id_for, ik_for, isigma_for, ichain_for, isigma2_for, idata_for, iin_for, iinit_for, iess_for, typeess_for, ntry_for)
                      source("simulations/ctpn.R", echo = TRUE)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
