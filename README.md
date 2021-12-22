## SapwoodModel

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
Get to know the package by trying to fit your first model to the data. Fetch some data by calling
```
data(smaland)
```
and run the model with:
```
model <- sapwood_fit_pl(S~H, smaland)
```
Note that `sapwood_fit_pl` can be replaced with `sapwood_fit_l` or `sapwood_fit_plw` (see `?sapwood_fit` for usage documentation).

From here, you can see some information about the fit by running
```
summary(model)
plot(model)
predict(model, 220:300)
model$parameter_CI
```

and more (see `?plot.sapwood_fit`, `?predict.sapwood_fit` etc.)


### References
Edvardsson et. al.
