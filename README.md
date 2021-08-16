## Sapwood rating curves

This software package fits a model to sapwood data from scots pine (Pinus sylvestris L.) using heartwood data and tree ring width using a nonlinear model described in Edvardsson et al. (2021). Three models are implemented:

`sapwood_fit_plw()` - Parabolic-linear model with constant variance and lognormal residuals that uses tree ring width, described as Model 1 in the article

`sapwood_fit_pl()` - Parabolic-linear model with constant variance and lognormal residuals that does not use tree ring width, described as Model 2 in the article

`sapwood_fit_l()` - Linear model with constant variance and lognormal residuals, described as Model 3 in the article

### Installation
```
# Install development version from GitHub
devtools::install_github("karirogg/sapwood_scots_pine")
```

### Getting started

### References
Edvardsson et. al.
