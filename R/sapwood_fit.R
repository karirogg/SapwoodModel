##' @name sapwood_fit
##' @rdname sapwood_fit
##'
##' @title A model for sapwood rings in scots pine
##'
##' @description The function takes in heartwood and sapwood (and possibly tree ring width) data from a dataset and returns a fit and some information about the fit, such as prediction intervals, figures, confidence intervals for parameters etc. See Edvardsson et al.
##'
#' @param formula Formula for the fit. If using \code{sapwood_fit_pl} or \code{sapwood_fit_l},
#' the formula should be on the form \code{S~H}, where S is the name of the column containing the
#' number of sapwood rings and H is the name of the column containing the number of heartwood rings.
#' If using \code{sapwood_fit_plw}, the formula should be on the form \code{S~H+W}, where S and H are as before,
#' and W is the name of the column including data on mean tree ring width.
#' @param dat Dataset to be fitted to. Column names should match formula arguments
#' @param alpha Confidence of the fit, for prediction, confidence for median and parameter confidence intervals. defaults to 0.05 (which corresponds to 95\% confidence)
#' @param mu_theta_1 Regularization mu for theta_1 (see Edvardsson et al. 2021)
#' @param sd_theta_1 Regularization sigma for theta_1 (see Edvardsson et al. 2021)
#' @param mu_theta_2 Regularization mu for theta_2 (see Edvardsson et al. 2021)
#' @param sd_theta_2 Regularization sigma for theta_2 (see Edvardsson et al. 2021)
#' @param H_0 Cutoff point for parabolic-linear model. Only relevant for \code{sapwood_fit_pl} and \code{sapwood_fit_plw}.
#' @return The functions return objects of class "sapwood_fit". An object of class "plm" is a list containing the following components:\cr
#' \item{\code{parameter_CI}}{Confidence intervals for the parameters of the model.}
#' \item{\code{predictions}}{Predictions for the model for H between 0 and 500, if \code{sapwood_fit_l} or \code{sapwood_fit_pl} are used, prediction and confidence intervals for the fit is included. If \code{sapwood_fit_plw} is used, only median of the prediction is returned.}
#' \item{\code{residuals}}{Standardized residuals of the model. See Edvardsson et al. 2021.}
#' \item{\code{AIC}}{AIC of the model}
#' \item{\code{formula}}{The formula for the model, as specified in the input.}
#' \item{\code{type}}{Type of the model, "parabolic_linear_W", "parabolic_linear" or "linear" (Models 1,2 and 3, respectively).}
#' \item{\code{alpha}}{Confidence of the model}
#'
#' @examples
#' data(dat_TA)
#' fit <- sapwood_fit_pl(S~H, dat_TA)
#' plot(fit)
#' plot(fit, type="residual")
#' plot(fit, type="qq")
#' residuals(fit)
#' AIC(fit)
#'
NULL

#' @rdname sapwood_fit
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

#' @rdname sapwood_fit
#' @export
sapwood_fit_pl <- function(formula,
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
        stop(paste0("Formula ", deparse(formula), " not valid. Please specify a formula on the form S~H"))
    return(sapwood_fit_raw(formula=formula,
                           dat = dat,
                           fit_function = log_parabolic_estimate,
                           parameters = c('theta_0', 'theta_1', 'theta_2'),
                           alpha = alpha,
                           transformation = log,
                           inverse_transformation = exp,
                           sigma_function = sigma_i_function,
                           regularization_term = parabolic_regularization_term,
                           theta_1_index = 2,
                           mu_theta_1 = mu_theta_1,
                           sd_theta_1 = sd_theta_1,
                           theta_2_index = 3,
                           mu_theta_2 = mu_theta_2,
                           sd_theta_2 = sd_theta_2,
                           type="parabolic_linear",
                           H_0 = H_0))
}

#' @rdname sapwood_fit
#' @export
sapwood_fit_plw <- function(formula,
                            dat,
                            alpha=0.05,
                            mu_theta_1 = log(0.1),
                            sd_theta_1 = (log(0.1)+4*log(10))/2,
                            mu_theta_2 = log(0.1),
                            sd_theta_2 = (log(0.1)+4*log(10))/2,
                            H_0 = 100) {
    stopifnot('formula' %in% class(formula))
    stopifnot('data.frame' %in% class(dat) | 'tibble' %in% class(dat))
    if(length(all.vars(formula)) != 3)
        stop(paste0("Formula ", formula, " not valid. Please specify a formula on the form S~H+W"))
    return(sapwood_fit_raw(formula=formula,
                           dat = dat,
                           fit_function = log_parabolic_estimate_W,
                           parameters = c('theta_0', 'theta_1', 'theta_2', 'beta_3'),
                           alpha = alpha,
                           transformation = log,
                           inverse_transformation = exp,
                           sigma_function = sigma_i_function,
                           regularization_term = parabolic_regularization_term,
                           theta_1_index = 2,
                           mu_theta_1 = mu_theta_1,
                           sd_theta_1 = sd_theta_1,
                           theta_2_index = 3,
                           mu_theta_2 = mu_theta_2,
                           sd_theta_2 = sd_theta_2,
                           type="parabolic_linear_W",
                           H_0 = H_0))
}
