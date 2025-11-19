# IRT-M Continuous Sampler Simulations

This directory contains simulation scripts used to test and validate the continuous extension of the IRT-M sampler currently under development.

The simulations are intended to stress-test behavior under different latent-trait and item-parameter distributions before finalizing code for integration into the package.

## Contents

### Simulation Scripts
- **sim1-normal.R**  
  Latent traits drawn from a standard normal distribution; baseline behavior check for the continuous sampler.

- **sim2-skew.R**  
  Tests performance under skewed latent distributions (e.g., log-normal, shifted gamma).

- **sim3-heavytail.R**  
  Heavy-tailed latent distributions (e.g., t-distribution) to evaluate robustness under extreme values.

- **sim4-bounded.R**  
  Latent traits constrained to finite intervals; checks sampler stability under bounded support.

### Exploratory Notes
- **exploratory-notes.Rmd**  
  Running commentary, diagnostics, plots, and early analysis from simulation runs.

## How to Run These Simulations

Until the updated sampler is fully merged and installed as a package, load the local development version directly with:

```r
devtools::load_all("path/to/IRTM-package")
