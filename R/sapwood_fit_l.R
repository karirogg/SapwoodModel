source("sapwood_utils.R")

#' @export
sapwood_fit_l <- function(formula,
                          dat,
                          alpha=0.05,
                          mu_theta_1 = log(0.1),
                          sd_theta_1 = (log(0.1)+4*log(10))/2,
                          mu_theta_2 = log(0.1),
                          sd_theta_2 = (log(0.1)+4*log(10))/2,
                          H_0 = 100) {
    stopifnot('formula' %in% class(formula))
    stopifnot('data.frame' %in% class(dat) | 'tibble' %in% class(dat))
    if(length(all.vars(formula)) != 2)
        stop(paste0("Formula ", formula, " not valid. Please specify a formula on the form S~H"))
    return(sapwood_fit_raw(formula=formula,
                           dat = dat,
                           fit_function = log_linear_estimate,
                           parameters = c('theta_0', 'theta_1'),
                           alpha = alpha,
                           transformation = log,
                           inverse_transformation = exp,
                           sigma_function = sigma_i_function,
                           type="linear",
                           H_0 = H_0))
}
