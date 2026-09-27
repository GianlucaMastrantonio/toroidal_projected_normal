# toroidal_projected_normal

This repository contains the code and data used for the analyses in the paper *An interpretable family of projected normal distributions and a related copula model for Bayesian analysis of hypertoroidal data*.

The repository has two main purposes:

1. It contains the local R package `toroidalPNcopula`, which provides the main MCMC functions.
2. It contains the scripts needed to reproduce the simulation study and real-data analysis from the paper.

## Install the Package

Before using the functions or running the reproduction scripts, install the local package from the root of this repository:

```r
install.packages("toroidalPNcopula", repos = NULL, type = "source")
```

Then load it with:

```r
library(toroidalPNcopula)
```

During package development, the package can also be loaded without installing it:

```r
devtools::load_all("toroidalPNcopula")
```

## Using the Package

The package exposes two main MCMC functions:

- `mcmc_tpn()` fits the Toroidal Projected Normal model.
- `mcmc_cwc()` fits the related wrapped Cauchy copula model.

Both functions expect the data as a matrix `theta`, where rows are observations, columns are angular variables, and all angles are in radians.

### Required Fixed Options

Two arguments are kept in the function interface for development reasons. For the analyses in this repository, and for ordinary use of the current implementation, they must always be set as:

```r
do_only_ESS = TRUE
type_ess = 1
```

The launch scripts already set these values.

### Main Argument Groups

The MCMC length is controlled by:

- `burnin`: number of burn-in iterations.
- `thin`: thinning interval.
- `iterations`: total number of MCMC iterations.

The prior arguments common to both models are:

- `prior_mu_mean`: prior mean for the angular location parameter `mu`.
- `prior_mu_var`: prior variance for `mu`.
- `prior_sigma_nu`: degrees of freedom of the inverse-Wishart prior on the covariance matrix.
- `prior_sigma_psi`: scale matrix of the inverse-Wishart prior.

For `mcmc_tpn()`, the model-specific prior arguments are:

- `prior_kappa_mean`: prior mean for `kappa`.
- `prior_kappa_var`: prior variance for `kappa`.

For `mcmc_cwc()`, the model-specific prior arguments are:

- `prior_rho_a`: first beta-prior parameter for `rho`.
- `prior_rho_b`: second beta-prior parameter for `rho`.

The initial values are:

- `mu_init`: initial value for `mu`.
- `sigma_init`: initial covariance matrix.
- `r_init`: initial latent radial variables.
- `kappa_init`: initial value for `kappa`, used only by `mcmc_tpn()`.
- `rho_init`: initial value for `rho`, used only by `mcmc_cwc()`.

The adaptive sampler arguments are:

- `adapt_batch`: number of iterations between adaptation updates.
- `adapt_a`, `adapt_b`: adaptation-rate parameters.
- `adapt_alpha_target`: target acceptance probability.
- `sd_mu_scal`: initial proposal scale for `mu`.
- `sd_rho_scal`: initial proposal scale for `rho`, used only by `mcmc_cwc()`.
- `par_sigma_adapt`: adaptation parameter for covariance updates.
- `n_test_sigma`: number of covariance proposals checked in each covariance update.

Missing values can be passed through `na_index`, a list of missing-value indices by dimension. If there are no missing values, use:

```r
na_index = list(NA)
```

### Function Call Templates

For the Toroidal Projected Normal model:

```r
fit_tpn <- mcmc_tpn(
  theta = theta,
  burnin = burnin,
  thin = thin,
  iterations = iterations,
  prior_mu_mean = prior_mu_mean,
  prior_mu_var = prior_mu_var,
  prior_kappa_mean = prior_kappa_mean,
  prior_kappa_var = prior_kappa_var,
  prior_sigma_nu = prior_sigma_nu,
  prior_sigma_psi = prior_sigma_psi,
  mu_init = mu_init,
  kappa_init = kappa_init,
  sigma_init = sigma_init,
  r_init = r_init,
  adapt_batch = adapt_batch,
  adapt_a = adapt_a,
  adapt_b = adapt_b,
  adapt_alpha_target = adapt_alpha_target,
  sd_mu_scal = sd_mu_scal,
  par_sigma_adapt = par_sigma_adapt,
  na_index = list(NA),
  do_only_ESS = TRUE,
  n_test_sigma = n_test_sigma,
  type_ess = 1
)
```

For the wrapped Cauchy copula model:

```r
fit_cwc <- mcmc_cwc(
  theta = theta,
  burnin = burnin,
  thin = thin,
  iterations = iterations,
  prior_mu_mean = prior_mu_mean,
  prior_mu_var = prior_mu_var,
  prior_rho_a = prior_rho_a,
  prior_rho_b = prior_rho_b,
  prior_sigma_nu = prior_sigma_nu,
  prior_sigma_psi = prior_sigma_psi,
  mu_init = mu_init,
  rho_init = rho_init,
  sigma_init = sigma_init,
  r_init = r_init,
  adapt_batch = adapt_batch,
  adapt_a = adapt_a,
  adapt_b = adapt_b,
  adapt_alpha_target = adapt_alpha_target,
  sd_mu_scal = sd_mu_scal,
  sd_rho_scal = sd_rho_scal,
  par_sigma_adapt = par_sigma_adapt,
  na_index = list(NA),
  do_only_ESS = TRUE,
  n_test_sigma = n_test_sigma,
  type_ess = 1
)
```

The returned object contains posterior samples of the model parameters, latent variables, missing-value imputations when requested, and sampler diagnostics.

## Reproducing the Paper Results

The scripts in this repository reproduce the simulation study and the real-data application from the paper. These scripts call the package functions internally, set the required MCMC options, and save outputs in the corresponding output folders.

Run all commands from the root of the repository.

### Simulation Study

Posterior samples for all simulated datasets and chains can be obtained with the launch files in `simulations/`.

For the Toroidal Projected Normal model:

```r
source("simulations/1 - launch_tpn.R")
```

For the copula model:

```r
source("simulations/1 - launch_ctpn.R")
```

Posterior samples are saved in:

```text
simulations/output/
```

After the posterior samples have been generated, compute the model diagnostics with:

```r
source("simulations/2 - post_analisi.R")
```

Diagnostic outputs are saved in:

```text
simulations/output diagnostic/
```

### Real-Data Application

Posterior samples for the real-data application can be obtained with the launch files in `real data/`.

For the Toroidal Projected Normal model:

```r
source("real data/1 - launch_tpn.R")
```

For the copula model:

```r
source("real data/1 - launch_ctpn.R")
```

Posterior samples are saved in:

```text
real data/output/
```

After the posterior samples have been generated, compute the model diagnostics with:

```r
source("real data/2 - post_analisi.R")
```

Diagnostic outputs are saved in:

```text
real data/output diagnostic/
```

### Figures and Tables

After running the simulation and real-data scripts, figures and tables can be generated from the scripts in `plots and tables/`.

To generate correlation plots:

```r
source("plots and tables/Correlations.R")
```

To generate plots for the real-data application:

```r
source("plots and tables/PlotReal.R")
```

## Notes

For questions or problems, contact gianluca.mastrantonio@polito.it.

## License

See [LICENSE](LICENSE) for details.
