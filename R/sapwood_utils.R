library(dplyr)
library(ggplot2)
library(here)
library(Hmisc)
library(stringr)
library(tidyr)

# Fit function for main model
log_parabolic_estimate <- function(theta, H, H_0, Z) {
    theta_0 <- theta[1][[1]]
    theta_1 <- theta[2][[1]]
    theta_2 <- theta[3][[1]]

    log(exp(theta_0) + exp(theta_2)*H+exp(theta_1)*((H <= H_0)*(H-H^2/(2*H_0))+(H>H_0)*1/2*H_0))
}

# Fit function for linear model with lognormal residuals
log_linear_estimate <- function(theta, H, H_0, Z) {
    theta_0 <- theta[1][[1]]
    theta_1 <- theta[2][[1]]

    log(exp(theta_0) + exp(theta_1)*H)
}

# Fit function for main model with the addition of the Z parameter
log_parabolic_estimate_W <- function(theta, H, H_0, Z) {
    theta_0 <- theta[1][[1]]
    theta_1 <- theta[2][[1]]
    theta_2 <- theta[3][[1]]
    beta_3 <- theta[4][[1]]

    log(exp(theta_0) + exp(theta_2)*H+exp(theta_1)*((H <= H_0)*(H-H^2/(2*H_0))+(H>H_0)*1/2*H_0) + beta_3*Z)
}

# Make sure data has constant variance on real scale, see article
sigma_i_function <- function(tau, mu) {
    sqrt(log(1/2+1/2*sqrt(1+4*exp(tau-2*mu))))
}

# Regularization of parameters to increase quality of confidence intervals
parabolic_regularization_term <- function(parameters, theta_1_index, mu_theta_1, sd_theta_1, theta_2_index, mu_theta_2, sd_theta_2) {
    -(parameters[theta_1_index][[1]]-mu_theta_1)^2/(2*sd_theta_1)-
        (parameters[theta_2_index][[1]]-mu_theta_2)^2/(2*sd_theta_2)
}

# Log likelihood for lognormal distribution
calculate_log_likelihood <- function(parameters, x, y, Z, fit_function, H_0, sigma_function) {
    tau <- parameters[length(parameters)][[1]]

    mu <- fit_function(theta=parameters, H=x, H_0=H_0, Z=Z)
    sigma <- sigma_function(tau, mu)

    # adding y in the end for likelihood function of the lognormal distribution
    -sum((y-mu)^2/(2*sigma^2) + log(sqrt(2*pi)*sigma) + y)
}

calculate_negative_log_likelihood <- function(parameters, x, y, Z, fit_function, H_0, sigma_function, regularization_term=function(x){0}, ...) {
    -(calculate_log_likelihood(parameters=unlist(parameters),
                               x=x,
                               y=y,
                               Z=Z,
                               fit_function=fit_function,
                               H_0 = H_0,
                               sigma_function = sigma_function) +
          regularization_term(parameters,...))
}

# Not used
log_likelihood_vectorized <- function(theta_matrix, ...) {
    values = c()
    for(i in 1:nrow(theta_matrix)) {
        values[i] = calculate_log_likelihood(unlist(theta_matrix[i,]), ...)
    }
    values
}


# Plot QQ plot with confidence intervals
#' @importFrom ggplot2 ggplot aes geom_point geom_line xlab ylab coord_cartesian
#' @importFrom stats qbeta qnorm
normalqqplot <- function(x){
    alpha <- 0.05
    n <- length(x)
    pplot <- (1:n)/(n + 1)
    plow <- qbeta(alpha/2,(1:n),n-(1:n)+1)
    pupp <- qbeta(1-alpha/2,(1:n),n-(1:n)+1)
    qqnorm <- qnorm(pplot)
    qx_low <- qnorm(plow)
    qx_upp <- qnorm(pupp)
    e_sort <- sort(x)
    index_e = order(x)
    p <- ggplot(data=NULL,aes(x=qqnorm, y=e_sort)) +
        geom_point() +
        geom_line(y = qx_low, lty=2) +
        geom_line(y = qx_upp, lty=2) +
        geom_line(y = qqnorm) +
        xlab("Standard normal quantiles") +
        ylab("Standardized residuals") +
        coord_cartesian(xlim=c(-2.3,2.3), ylim=c(-2.3,2.3))
    p
}

# Fit data with model, return plots, predictions, parameter CI and model performance
#' @importFrom stats as.formula quantile rnorm lm nlminb
#' @importFrom dplyr %>% rename mutate filter tibble bind_rows as_tibble group_by summarise n inner_join arrange select left_join if_else
#' @importFrom utils head
#' @importFrom stringr str_detect
sapwood_fit_raw <- function(formula = as.formula("S~H"),
                            type = "parabolic_linear",
                            dat,
                            fit_function,
                            parameters,
                            alpha = 0.05,
                            regularization_term = function(parameters) {0},
                            transformation = function(x){x},
                            inverse_transformation = function(x){x},
                            H_0,
                            sigma_function = function(tau,mu){exp(tau)}, ...) {
    formula_args <- all.vars(formula)


    dat <- dat[,formula_args]
    # Transform data accordingly and rename columns to needs
    dat <- dat %>%
        rename(H = !!formula_args[2][[1]],
               S = !!formula_args[1][[1]]) %>%
        mutate(transformed_S = transformation(S)) %>%
        filter(H != 0)

    if(length(formula_args) == 3) {
        dat <- dat %>% rename(W = !!formula_args[3][[1]])
    }

    start_list <- list()

    if("W" %in% colnames(dat)) {
        # Add transformed mean TRW (if included)
        Z_model <- lm(W^(-1/2) ~ H, data=dat)
        dat <- dat %>% filter(W != 0) %>% mutate(Z = residuals(Z_model))
    }

    # Initializing parameters for optimization
    for(parameter in parameters) start_list[parameter] = 0

    # exp(tau) is the variance of the data on real scale (constant)
    start_list['tau'] = 0

    Z_input <- rep(0,nrow(dat))
    if("Z" %in% colnames(dat)) Z_input <- dat$Z

    # Fitting maximum likelihood estimate for log likelihood function
    # with regularization term
    fit <- nlminb(start = unlist(start_list),
                  objective = calculate_negative_log_likelihood,
                  x = dat$H,
                  y = transformation(dat$S),
                  Z = Z_input,
                  fit_function = fit_function,
                  sigma_function = sigma_function,
                  regularization_term = regularization_term,
                  H_0 = H_0,
                  ...)

    kmax <- 500

    tau = fit$par['tau'][[1]]

    best_fit_function <- function(H,Z) {
        fit_function(unlist(fit$par), H, H_0, Z)
    }

    # Bootstrap for confidence interval
    # Resampling from standardized residuals, prediction and confidence intervals calculated pointwise.
    residuals <- (transformation(dat$S)-best_fit_function(dat$H,Z_input))/sigma_function(tau, best_fit_function(dat$H,Z_input))

    n_samples <- 1000
    sample_size <- nrow(dat)
    n_prediction_samples <- 100

    # Contains the bootstrap estimates in every point to create confidence interval for the fit
    bootstrap_medians <- tibble()
    # Contains the bootstrap estimates for the parameters in every resample to create CI for parameters
    bootstrap_parameters <- tibble()

    if(type != "parabolic_linear_W") {
        # Contains the bootstrap estimates in every point to create prediction interval for the fit
        bootstrap_predictions <- rep(0,(kmax+1)*n_samples*n_prediction_samples)
    }

    mle = fit$par
    if(length(parameters) > 2) mle["beta_1+beta_2"] = exp(fit$par["theta_1"][[1]]) + exp(fit$par["theta_2"][[1]])
    mle["sigma"] = exp(fit$par["tau"][[1]]/2)

    bootstrap_prediction_H <- rep(0:kmax, each=n_prediction_samples)

    for(i in 1:n_samples) {
        if(i %% floor(n_samples/100) == 0) {
            cat("|")
        }

        # Sampling
        sample_indices <- sample(x=1:sample_size, size=sample_size, replace=T)

        # Standardized residuals calculated to residuals
        residuals_bootstrap <- residuals[sample_indices]*sigma_function(tau, best_fit_function(dat$H, Z_input))

        # Fetch data to perform fit on by adding the residuals to mu_i
        S_bootstrap <- best_fit_function(dat$H, Z_input)+residuals_bootstrap

        start_list = fit$par

        # Fitting maximum likelihood estimate for log likelihood function
        # with regularization term
        bootstrap_fit <- nlminb(start = unlist(start_list),
                                objective = calculate_negative_log_likelihood,
                                x = dat$H,
                                y = S_bootstrap,
                                Z = Z_input,
                                fit_function = fit_function,
                                sigma_function = sigma_function,
                                regularization_term = regularization_term,
                                H_0 = H_0,
                                ...)

        bootstrap_medians <-  bootstrap_medians %>%
            bind_rows(tibble(H=0:kmax,
                             value=fit_function(head(unlist(bootstrap_fit$par),-1),
                                                0:kmax,
                                                H_0,
                                                0)))

        parameter_list <- bootstrap_fit$par
        if(length(parameters) > 2)
            parameter_list["beta_1+beta_2"] = exp(bootstrap_fit$par["theta_1"][[1]]) +
                exp(bootstrap_fit$par["theta_2"][[1]])
        parameter_list["sigma"] = exp(bootstrap_fit$par["tau"][[1]]/2)
        parameter_tibble <- as_tibble(parameter_list) %>% mutate(term=names(parameter_list)) %>%
            mutate(replication = i)

        bootstrap_parameters <- bootstrap_parameters %>% bind_rows(parameter_tibble)
    }
    cat("\n")

    prediction_intervals <- tibble()

    if(type != "parabolic_linear_W") {
        confidence_intervals_mu <- bootstrap_medians %>%
            group_by(H) %>%
            summarise(conf.lower = round(inverse_transformation(quantile(value, alpha/2))),
                      conf.upper = round(inverse_transformation(quantile(value, 1-alpha/2))))
        for(i in 0:kmax) {
            H_medians <- bootstrap_medians %>% filter(H == i)
            bootstrap_predictions <- tibble(H=rep(H_medians$H, each=100), value=rep(H_medians$value, each=100)) %>%
                mutate(log_prediction = rnorm(n(), value, sigma_function(tau, best_fit_function(H))),
                       prediction = round(exp(log_prediction))) %>%
                group_by(H) %>%
                summarise(pred.lower = quantile(prediction, alpha/2),
                          pred.upper = quantile(prediction, 1-alpha/2))
            prediction_intervals <- prediction_intervals %>% bind_rows(bootstrap_predictions)
        }

        prediction_intervals <- prediction_intervals %>%
            mutate(median = inverse_transformation(best_fit_function(H, 0))) %>%
            inner_join(confidence_intervals_mu, by="H") %>%
            arrange(H) %>%
            select(H, median, pred.lower, pred.upper, conf.lower, conf.upper)
    } else {
        prediction_intervals <- dat %>% mutate(median = best_fit_function(H, Z_input)) %>%
                                        select(H,Z, median) %>%
                                        arrange(H,Z)
    }

    parameter_medians <- as_tibble(mle) %>% mutate(term = names(mle)) %>% rename(median=value)
    parameter_CI <- bootstrap_parameters %>%
        group_by(term) %>%
        summarise(conf.lower = quantile(value, alpha/2),
                  conf.upper = quantile(value, 1-alpha/2)) %>%
        ungroup %>%
        left_join(parameter_medians, by="term")

    # AIC calculations (length(parameters)+1 used because tau is not included in the parameters vector)
    model_AIC <-  2*(length(parameters)+1) -
        2*(calculate_log_likelihood(parameters = unlist(fit$par),
                                    x = dat$H,
                                    y = dat$transformed_S,
                                    Z = Z_input,
                                    fit_function = fit_function,
                                    H_0 = H_0,
                                    sigma_function = sigma_function))

    # Combining parameter CI, including median
    # exp(theta_0), exp(theta_1) and exp(theta_2) are easier to interpret so CI are based on them
    # as well as sqrt(exp(tau)) (standard deviation of the fit)
    transformed_parameter_CI <- parameter_CI %>%
        filter(str_detect(term, c("theta"))) %>%
        mutate(median = exp(median),
               conf.lower = exp(conf.lower),
               conf.upper = exp(conf.upper),
               term = paste0("b", substr(term, 3, nchar(term))))

    transformed_parameter_CI <- bind_rows(transformed_parameter_CI,
                                          parameter_CI) %>%
        arrange(parameter) %>%
        select(term, conf.lower, median, conf.upper) %>%
        filter(term %in% c("beta_0", "beta_1", "beta_2", "beta_3", "beta_1+beta_2", "sigma"))

    out <- list()

    attr(out, "class") <- "sapwood_fit"

    out$dat <- dat
    out$parameter_CI <- transformed_parameter_CI
    out$predictions <- prediction_intervals
    out$residuals <- residuals
    out$AIC <- model_AIC
    out$formula <- formula
    out$type <- type
    out$alpha <- alpha
    out$data <- dat
    out$sigma_function <- sigma_function
    out$best_fit_function <- best_fit_function
    out$fit_function <- fit_function
    out$bootstrap_parameters <- bootstrap_parameters
    out$mle <- mle
    out$H_0 <- H_0
    out$n_samples <- n_samples
    out$n_prediction_samples <- n_prediction_samples
    if(type == "parabolic_linear_W") out$Z_model = Z_model

    out
}

#' @export
AIC.sapwood_fit <- function(object,...) {
    object$AIC
}

#' Prediction for sapwood rings
#'
#' Obtain predictions for a model of type 'sapwood_fit',
#' including prediction and confidence intervals.
#' Optionally, one can predict how many sapwood rings were in a sample
#' that contains some remaining rings but not all.
#'
#' @param object an object of class "sapwood_fit", a result from a call to \code{sapwood_fit_l}, \code{sapwood_fit_pl} or \code{sapwood_fit_plw}.
#' @param newdata an optional data frame/tibble/vector in which to look for variables with which to predict. If omitted, the fitted values are used. If a column \code{remaining} is in the \code{newdata} tibble/data frame, it will be treated as remaining sapwood rings and the prediction interval will be based on that
#' @param confidence Confidence level in prediction (1-alpha). Defaults to 0.95.
#' @export
#' @importFrom dplyr mutate filter %>% select bind_rows left_join
#' @importFrom stats median
predict.sapwood_fit <- function(object, newdata=NULL, confidence=0.95) {
    if(is.null(newdata) & confidence == 1-object$alpha) {
        return(mutate(object$predictions,median = round(median)))
    }

    if(!("data.frame" %in% class(newdata))) return(filter(object$predictions, object$predictions$H %in% newdata))
    if(!("remaining" %in% colnames(newdata))) {
        newdata <- newdata %>% mutate(remaining = 0)
    }

    predictions_no_remaining <- newdata %>% filter(remaining == 0 | is.na(remaining)) %>%
        select(-remaining) %>%
        predict_util(object, ., confidence) %>%
        mutate(remaining = 0)

    predictions_remaining <- newdata %>% filter(remaining != 0 & !is.na(remaining)) %>%
        predict_util(object, ., confidence)

    predictions <- bind_rows(predictions_no_remaining,
                             predictions_remaining)

    newdata %>% left_join(predictions)
}

#' @importFrom dplyr %>% filter tibble slice n mutate group_by summarise ungroup select if_else rename
#' @importFrom tidyr pivot_wider
#' @importFrom stats median rnorm quantile
predict_util <- function(object, newdata=NULL, confidence=0.95) {
    if((object$type == "parabolic_linear" | object$type == "linear") & !("remaining" %in% colnames(newdata)) & confidence == 1-object$alpha) {
        return(filter(object$predictions, object$predictions$H %in% newdata$H))
    }
    else {
        alpha <- 1-confidence

        n_samples <- 1000
        sample_size <- nrow(object$dat)
        n_prediction_samples <- 100

        bootstrap_Z <- 0
        bootstrap_W <- 0
        if(object$type == "parabolic_linear_W") {
            len <- nrow(object$bootstrap_parameters)/length(unique(object$bootstrap_parameters$term))
            bootstrap_Z <- rep(newdata$W-predict(object$Z_model, newdata=tibble(H=newdata$H)), len)
            bootstrap_W <- rep(newdata$W, len)
        }

        predictions <- object$bootstrap_parameters %>% filter(term != "beta_1+beta_2") %>%
            pivot_wider(names_from=term, values_from=value) %>%
            slice(rep(1:n(),each=nrow(newdata))) %>%
            mutate(H = rep_len(newdata$H, n()),
                   W = bootstrap_W,
                   Z = bootstrap_Z,
                   remaining = 0,
                   median = 0)

        if(!("beta_3" %in% colnames(predictions))) predictions <- predictions %>% mutate(beta_3 = 0)
        if(!("theta_2" %in% colnames(predictions))) predictions <- predictions %>% mutate(theta_2 = 0)
        if("remaining" %in% colnames(newdata)) predictions <- predictions %>% mutate(remaining = rep_len(newdata$remaining, n()))

        for(i in 1:nrow(predictions)) {
            predictions$median[i] = object$fit_function(c(predictions$theta_0[i],
                                                          predictions$theta_1[i],
                                                          predictions$theta_2[i],
                                                          predictions$beta_3[i]),
                                                        H = predictions$H[i],
                                                        Z = predictions$Z[i],
                                                        H_0 = object$H_0)
        }

        predictions <- predictions %>%
            slice(rep(1:n(), times=n_prediction_samples)) %>%
            mutate(sigma_i = object$sigma_function(object$mle['tau'][[1]], median),
                   log_prediction = rnorm(n(), median, sigma_i),
                   prediction = round(exp(log_prediction))) %>%
            group_by(H,Z,W, remaining) %>%
            filter(prediction >= remaining) %>%
            summarise(pred.lower = quantile(prediction, alpha/2),
                      pred.upper = quantile(prediction, 1-alpha/2),
                      conf.lower = round(exp(quantile(median,alpha/2))),
                      conf.upper = round(exp(quantile(median,1-alpha/2))),
                      median = round(quantile(prediction,0.5))) %>%
            ungroup()

        if(!("remaining" %in% colnames(newdata))) {
            predictions <- predictions %>% mutate(median = round(exp(object$best_fit_function(H,Z))))
        }

        predictions <- predictions %>% select(H,W,remaining, median,pred.lower,pred.upper,conf.lower,conf.upper)

        if(object$type != "parabolic_linear_W") predictions <- predictions %>% select(c(-W))
        if(!("remaining" %in% colnames(newdata))) predictions <- predictions %>% select(c(-remaining))
        else predictions <- predictions %>% select(c(-conf.lower,-conf.upper))

        return(predictions)
    }
}

#' @export
residuals.sapwood_fit <- function(object,...) {
    object$residuals
}

#' @export
print.sapwood_fit <- function(x,...) {
    cat(paste0(" ",class(x), " - Call:\n\n"),
        if_else(x$type == "parabolic_linear", "Linear-parabolic model (no mean TRW), Model 2\n\n",
                if_else(x$type == "linear", "Linear model (no mean TRW), Model 3\n\n",
                        "Linear-parabolic model, using mean TRW, Model 1\n\n")),
        paste(deparse(x$formula), collapse = "\n"))
}

#' Summary method for sapwood fits
#'
#' \code{summary} method for class "sapwood_fit"
#'
#' @param object an object of class "sapwood_fit", a result from a call to \code{sapwood_fit_l}, \code{sapwood_fit_pl} or \code{sapwood_fit_plw}.
#' @export
summary.sapwood_fit <- function(object) {
    cat(paste0(" ",class(object), " - Call:\n\n"),
        if_else(object$type == "parabolic_linear", "Linear-parabolic model (no mean TRW), Model 2\n\n",
                if_else(object$type == "linear", "Linear model (no mean TRW), Model 3\n\n",
                        "Linear-parabolic model, using mean TRW, Model 1\n\n")),
        paste(deparse(object$formula), collapse = "\n"),
        "\n\n",
        "AIC:",
        AIC.sapwood_fit(object),
        "\n\n",
        "Coefficients:\n")
    print(object$parameter_CI)
}

#' Plot method for sapwood models
#'
#' Visualize sapwood model.
#'
#' @param object an object of class "sapwood_fit", a result from a call to \code{sapwood_fit_l}, \code{sapwood_fit_pl} or \code{sapwood_fit_plw}.
#' @param type Type of plot. Possible types are:
#' \itemize{
#' \item{"fit"}{Plots sapwood rings versus heartwood rings and the best fit, including prediction and/or confidence intervals for the median (if wanted). Not available for an object of type "parabolic_linear_W' (result from \code{sapwood_fit_plw})}
#' \item{"residual"}{Plots standardized residuals on log scale of the model as well as dotted lines for z_0.025 and z_0.975}
#' \item{"qq"}{Plots a QQ plot of the standardized residuals}
#' }
#' @param xlim,ylim Limits of the axes of the plot. Given as a vector similar to c(0,100). \code{xlim} defaults to the range from 0 to the maximum value of heartwood in the data (with some padding). \code{ylim} defaults to 0 to the maximum value of predicted data with some padding.
#' @param prediction If TRUE, prediction bands are plotted on the plot. Only used for type = "fit".
#' @param confidence If TRUE, confidence bands to the median are plotted to the plot. Only used for type = "fit".
#'
#' @examples
#' data(smaland)
#' fit <- sapwood_fit_pl(S~H, smaland)
#' plot(fit, xlim=c(0,220), ylim=c(0,150))
#' plot(fit)
#' plot(fit, type="residual")
#' plot(fit, type="qq")
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_line xlab ylab coord_cartesian theme_set theme scale_x_continuous scale_y_continuous theme_classic geom_hline
#' @importFrom dplyr %>% tibble
#' @importFrom stats qnorm
plot.sapwood_fit <- function(object, type='fit', xlim=NULL, ylim=NULL, prediction=T, confidence=F) {
    if(object$type == "parabolic_linear_W")
        stop("Not possible to plot model, use type='parabolic_linear' or type='linear' instead.")
    if(is.null(xlim)) xlim = c(0,max(object$data$H)+20)
    if(is.null(ylim)) ylim = c(0,max(object$predictions$pred.upper[xlim])+20)
    theme_set(theme_classic(base_size = 12) + theme(legend.position = "none"))
    if(type == "fit") {
        p <- object$data %>% ggplot(aes(x=H, y=S)) +
            geom_point() +
            geom_line(data=object$predictions, aes(x=H, y=median), color = "black") +
            coord_cartesian(xlim=xlim, ylim = ylim) +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0)) +
            xlab("Number of heartwood rings") +
            ylab("Number of sapwood rings")
        if(prediction)
            p <- p + geom_line(data=object$predictions, aes(x=H, y=pred.upper), lty=2) +
                geom_line(data=object$predictions, aes(x=H, y=pred.lower), lty=2)
        if(confidence)
            p <- p + geom_line(data=object$predictions, aes(x=H, y=conf.upper), lty=2) +
                geom_line(data=object$predictions, aes(x=H, y=conf.lower), lty=2)

        return(p)
    } else if(type == "residual") {
        p <- tibble(H = object$data$H, residuals = object$residuals) %>%
            ggplot(aes(x = H, y = residuals)) +
            xlab("Number of heartwood rings") +
            ylab("Standardized residuals") +
            coord_cartesian(xlim = xlim, ylim=c(-3.2,3.2)) +
            geom_point() +
            geom_hline(yintercept = 0) +
            geom_hline(yintercept = qnorm(object$alpha/2), lty=2) +
            geom_hline(yintercept = qnorm(1-object$alpha/2), lty=2)
        return(p)
    } else if(type == "qq") {
        return(normalqqplot(object$residuals))
    }

    stop(paste0("Type '", type, "' not available"))
}
