#' Validate parameter lengths for decay model functions
#'
#' @param ... Named arguments to validate
#' @param required_lengths List specifying which parameters should have the same length
#' @param function_name Name of the function for error messages (optional)
#' @return NULL if validation passes, throws error if validation fails
#' @keywords internal
#' @noRd
validate_parameter_lengths <- function(
  ...,
  required_lengths = NULL,
  function_name = NULL
) {
  args <- list(...)

  # If required_lengths is not specified, assume all parameters should have the same length
  if (is.null(required_lengths)) {
    param_lengths <- sapply(args, length)

    if (length(unique(param_lengths)) != 1) {
      param_names <- names(args)
      length_info <- paste(
        paste0(param_names, "=", param_lengths),
        collapse = ", "
      )

      error_msg <- paste0(
        "All parameters must have the same length. ",
        "Lengths: ",
        length_info
      )

      if (!is.null(function_name)) {
        error_msg <- paste0("In ", function_name, ": ", error_msg)
      }

      stop(error_msg)
    }
  } else {
    # Validate specific groups of parameters
    for (group_name in names(required_lengths)) {
      group_params <- required_lengths[[group_name]]

      # Check if all required parameters are present
      missing_params <- setdiff(group_params, names(args))
      if (length(missing_params) > 0) {
        stop(
          "Missing required parameters: ",
          paste(missing_params, collapse = ", ")
        )
      }

      # Get lengths for this group
      group_lengths <- sapply(args[group_params], length)

      if (length(unique(group_lengths)) != 1) {
        length_info <- paste(
          paste0(group_params, "=", group_lengths),
          collapse = ", "
        )

        error_msg <- paste0(
          group_name,
          " must have the same length. ",
          "Lengths: ",
          length_info
        )

        if (!is.null(function_name)) {
          error_msg <- paste0("In ", function_name, ": ", error_msg)
        }

        stop(error_msg)
      }
    }
  }

  return(invisible(NULL))
}

#' DecayModel class
#'
#' @description A class to represent a decay model
#' @importFrom R6 R6Class
#' @export
DecayModel <- R6::R6Class(
  "DecayModel",
  public = list(
    #' @field name Model name
    name = NULL,

    #' @field decay_function The decay function
    decay_function = NULL,

    #' @field param_config Parameter configuration
    param_config = NULL,

    #' @field optimization_results Results from the last optimization
    optimization_results = NULL,

    #' @field nlrob_samples List of nlrob samples
    nlrob_samples = NULL,

    #' @field ve_summary Summary of VE samples
    ve_summary = NULL,

    #' @description Initialize a specific decay model
    #' @param name Model name
    #' @param decay_function Function that implements the decay model
    #' @param param_config Parameter configuration list
    initialize = function(name, decay_function, param_config) {
      self$name <- name
      self$decay_function <- decay_function
      self$param_config <- param_config
      self$optimization_results <- NULL
    },

    #' @description Evaluate the model with given parameters
    #' @param time_steps Time steps to evaluate
    #' @param params Named vector of parameters
    #' @return Vector of decay values
    evaluate = function(time_steps, params) {
      param_list <- as.list(params)
      do.call(
        self$decay_function,
        c(list(time_steps = time_steps), param_list)
      )
    },

    #' @description Create formula for nls fitting
    #' @return Formula suitable for nls
    create_nls_formula = function() {
      decay_func <- self$decay_function

      # Create the formula based on the model name and parameters
      if (self$name == "exponential") {
        formula <- estimates ~ decay_func(w_start, VE0, decay_rate)
      } else if (self$name == "logistic") {
        formula <- estimates ~ decay_func(w_start, VE0, decay_rate, constant)
      } else {
        # For custom models, we need to construct the formula dynamically
        param_names <- self$param_config$names
        args <- paste(param_names, collapse = ", ")
        formula_str <- paste0("estimates ~ decay_func(w_start, ", args, ")")
        formula <- as.formula(formula_str)
      }

      return(formula)
    },

    #' @description Fit the model using 3-step strategy (nlxb -> nls -> bootstrap)
    #' @param estimate Log estimates (in the HR scale)
    #' @param ci_low Lower confidence interval bounds (in the HR scale)
    #' @param ci_up Upper confidence interval bounds (in the HR scale)
    #' @param periods_starts Vector of periods starts
    #' @param periods_ends Vector of periods ends
    #' @param seed Seed for reproducibility (default is 123)
    #' @param maxiteritations Maximum number of iterations for nls (default is 500)
    #' @param type_of_fit Type of fit (default is "midpoint")
    #' @importFrom stats nls AIC BIC coef
    #' @importFrom nlsr nlxb
    #' @importFrom robustbase nlrob
    #' @return List with fitting results
    fit_to_determine_bic = function(
      estimate,
      ci_low,
      ci_up,
      periods_starts,
      periods_ends,
      seed = 123,
      maxiteritations = 50,
      type_of_fit = "midpoint"
    ) {
      match.arg(type_of_fit, c("midpoint", "average"))
      # Validate parameter lengths
      validate_parameter_lengths(
        estimate = estimate,
        ci_low = ci_low,
        ci_up = ci_up,
        periods_starts = periods_starts,
        periods_ends = periods_ends,
        function_name = "fit_to_determine_bic"
      )

      # Get period-averaged decay functions
      decay_functions_period <- decay_functions_period_avg()

      # Get the appropriate decay function
      if (self$name == "exponential") {
        decay_func_period <- decay_functions_period[[paste0(
          "exponential_",
          type_of_fit
        )]]
      } else if (self$name == "logistic") {
        decay_func_period <- decay_functions_period[[paste0(
          "logistic_",
          type_of_fit
        )]]
      } else {
        stop("Unknown model name: ", self$name)
      }

      # Prepare data
      log_estimate <- log(estimate)
      se <- get_se_from_ci(log(ci_low), log(ci_up))
      timing <- build_time_objects(periods_starts, periods_ends)$timing
      n_weeks <- build_time_objects(periods_starts, periods_ends)$n_weeks

      # Create formula for fitting - we need to use a different approach for nls
      # We'll create the formula as a string and convert it
      if (self$name == "exponential") {
        formula_str <- "log_estimate ~ decay_func_period(w_start, nweeks, VE0, decay_rate)"
      } else if (self$name == "logistic") {
        formula_str <- "log_estimate ~ decay_func_period(w_start, nweeks, VE0, decay_rate, constant)"
      }

      # Create the formula with the decay function in the environment
      formula <- as.formula(formula_str)

      decay_func_period <<- decay_func_period

      nlrob_fit <- nlrob(
        as.formula(formula_str),
        data = data.frame(
          log_estimate = log_estimate,
          w_start = timing$diff_start,
          nweeks = n_weeks
        ),
        start = self$param_config$start,
        lower = self$param_config$lower,
        upper = self$param_config$upper,
        control = list(maxiter = maxiteritations, warnOnly = TRUE),
        # weights = 1 / se^2,
        maxit = maxiteritations,
        method = "M",
        algorithm = "port"
      )

      has_converged <- nlrob_fit$status == "converged"

      # Create results data frame (similar to optim results for compatibility)
      result_df <- as.data.frame(t(coef(nlrob_fit)))
      names(result_df) <- self$param_config$names
      result_df$aic <- AIC(nlrob_fit)
      result_df$bic <- BIC(nlrob_fit)
      result_df$model <- self$name

      # Store the results
      self$optimization_results <- list(
        params = result_df,
        nls_fit = nlrob_fit,
        coefficients = coef(nlrob_fit),
        AIC = AIC(nlrob_fit),
        BIC = BIC(nlrob_fit),
        summary = summary(nlrob_fit),
        decay_func = self$decay_function
      )

      return(self$optimization_results)
    },

    #' @description Generate VE based on boot fit
    #' @param estimate Log estimates (in the HR scale)
    #' @param ci_low Lower confidence interval bounds (in the HR scale)
    #' @param ci_up Upper confidence interval bounds (in the HR scale)
    #' @param periods_starts Vector of periods starts
    #' @param periods_ends Vector of periods ends
    #' @param n_weekly_ve_to_generate Number of weekly VE to generate (default is 1000)
    #' @param type_of_fit Type of fit (default is "midpoint")
    #' @param seed Seed for reproducibility (default is 123)
    #' @return Matrix of VE samples
    #' @importFrom MASS mvrnorm
    generate_VE_based_on_resampling = function(
      estimate,
      ci_low,
      ci_up,
      periods_starts,
      periods_ends,
      n_weekly_ve_to_generate = 1000,
      type_of_fit = "midpoint",
      seed = 123
    ) {
      set.seed(seed)
      match.arg(type_of_fit, c("midpoint", "average"))
      # Validate parameter lengths
      validate_parameter_lengths(
        periods_starts = periods_starts,
        periods_ends = periods_ends,
        function_name = "generate_VE_based_on_resampling"
      )

      if (is.null(self$optimization_results)) {
        stop(
          "No optimization results available. Run fit_to_determine_bic() first."
        )
      }

      generate_8week_sample <- function(log_estimate, sd) {
        rnorm(n = length(log_estimate), mean = log_estimate, sd = sd)
      }

      log_estimate <- log(estimate)
      se <- get_se_from_ci(log(ci_low), log(ci_up))
      timing <- build_time_objects(periods_starts, periods_ends)$timing
      n_weeks <- build_time_objects(periods_starts, periods_ends)$n_weeks

      eight_weeks_samples <- lapply(1:n_weekly_ve_to_generate, function(i) {
        generate_8week_sample(log_estimate, se)
      })

      fit_1_sample_nlrob <- function(
        log_estimate,
        se,
        timing,
        n_weeks,
        model,
        decay_func_period,
        correct_se_8weeks = FALSE
      ) {
        if (model$name == "exponential" & type_of_fit == "midpoint") {
          formula_str <- "log_estimate ~ decay_functions_period_avg()$exponential_midpoint(w_start, nweeks, VE0, decay_rate)"
        } else if (model$name == "logistic" & type_of_fit == "midpoint") {
          formula_str <- "log_estimate ~ decay_functions_period_avg()$logistic_midpoint(w_start, nweeks, VE0, decay_rate, constant)"
        } else if (model$name == "exponential" & type_of_fit == "average") {
          formula_str <- "log_estimate ~ decay_functions_period_avg()$exponential_average(w_start, nweeks, VE0, decay_rate)"
        } else if (model$name == "logistic" & type_of_fit == "average") {
          formula_str <- "log_estimate ~ decay_functions_period_avg()$logistic_average(w_start, nweeks, VE0, decay_rate, constant)"
        }

        nlrob_fit <- suppressWarnings(nlrob(
          as.formula(formula_str),
          data = data.frame(
            log_estimate = log_estimate,
            w_start = timing$diff_start,
            nweeks = n_weeks
          ),
          start = model$param_config$start,
          lower = model$param_config$lower,
          upper = model$param_config$upper,
          control = list(maxiter = 50, warnOnly = TRUE),
          # weights = 1 / se^2,
          maxit = 50,
          method = "M",
          algorithm = "port"
        ))

        weekly_VE <- do.call(
          model$decay_function,
          c(list(time_steps = 1:n_weeks), as.list(coef(nlrob_fit)))
        )

        sigma <- summary(nlrob_fit)$Scale
        if (correct_se_8weeks) {
          sigma <- sigma * sqrt(8)
        }

        weekly_VE_as_logHR <- log(1 - weekly_VE)
        weekly_VE_as_logHR_with_noise <- rnorm(
          length(weekly_VE_as_logHR),
          mean = weekly_VE_as_logHR,
          sd = sigma
        )
        weekly_VE_with_noise <- 1 - exp(weekly_VE_as_logHR_with_noise)

        weekly_VE_with_noise
      }

      fit_1_sample_nlrob_robust <- purrr::possibly(fit_1_sample_nlrob)

      nlrob_samples <- purrr::map(1:n_weekly_ve_to_generate, function(i) {
        log_estimate <- eight_weeks_samples[[i]]
        set.seed(seed + i)
        fit_1_sample_nlrob_robust(
          log_estimate = log_estimate,
          se = se,
          timing = timing,
          n_weeks = n_weeks,
          model = self,
          decay_func_period = decay_func_period,
          correct_se_8weeks = FALSE
        )
      })

      self$nlrob_samples <- nlrob_samples

      return(self$nlrob_samples)
    },

    #' @description Generate VE summary statistics from bootstrap samples
    #' @param periods_starts Vector of periods starts
    #' @param periods_ends Vector of periods ends
    #' @param probs Vector of quantiles to compute (default is c(0.025, 0.5, 0.975))
    #' @return Data frame with summary statistics
    #' @importFrom dplyr group_by summarise
    generate_VE_summary = function(
      periods_starts,
      periods_ends,
      probs = c(0.025, 0.5, 0.975)
    ) {
      # Validate parameter lengths
      validate_parameter_lengths(
        periods_starts = periods_starts,
        periods_ends = periods_ends,
        function_name = "generate_VE_summary"
      )
      if (is.null(self$nlrob_samples)) {
        stop(
          "No nlrob samples available. Run generate_VE_based_on_resampling() first."
        )
      }
      # Generate VE samples
      ve_samples <- self$nlrob_samples

      n_weeks <- build_time_objects(
        periods_starts,
        periods_ends
      )$n_weeks

      ve_samples <- lapply(seq_len(length(ve_samples)), function(x) {
        data.frame(
          period_start = min(periods_starts) + 0:(n_weeks - 1) * 7,
          VE = ve_samples[[x]],
          sample = x
        )
      })
      ve_samples <- do.call(rbind, ve_samples)

      # Calculate summary statistics
      ve_summary <- ve_samples %>%
        dplyr::group_by(period_start) %>%
        dplyr::summarise(
          VE_median = median(VE),
          VE_mean = mean(VE),
          VE_sd = sd(VE),
          VE_quantile_025 = quantile(VE, 0.025),
          VE_quantile_975 = quantile(VE, 0.975),
          .groups = "drop"
        )

      self$ve_summary <- ve_summary

      return(self$ve_summary)
    },

    #' @description Plot VE with uncertainty bands from bootstrap
    #' @param periods_starts Vector of periods starts
    #' @param periods_ends Vector of periods ends
    #' @param title Plot title (default is NULL)
    #' @param color Color for the line and ribbon (default is "blue")
    #' @param alpha Transparency for the ribbon (default is 0.3)
    #' @param show_points Whether to show points on the line (default is TRUE)
    #' @param probs Vector of quantiles for confidence intervals (default is c(0.025, 0.975))
    #' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point theme_minimal labs
    #' @return ggplot object
    plot_VE_with_uncertainty = function(
      periods_starts,
      periods_ends,
      title = NULL,
      color = "blue",
      alpha = 0.3,
      show_points = TRUE,
      probs = c(0.025, 0.975)
    ) {
      # Validate parameter lengths
      validate_parameter_lengths(
        periods_starts = periods_starts,
        periods_ends = periods_ends,
        function_name = "plot_VE_with_uncertainty"
      )

      if (is.null(self$optimization_results)) {
        stop("No optimization results available. Run fit_3step() first.")
      }

      # Generate summary statistics
      ve_summary <- self$generate_VE_summary(
        periods_starts,
        periods_ends,
        probs
      )

      # Create plot
      p <- ggplot(ve_summary, aes(x = period_start, y = VE_median)) +
        geom_ribbon(
          aes(ymin = VE_quantile_025, ymax = VE_quantile_975),
          alpha = alpha,
          fill = color
        ) +
        geom_line(color = color, linewidth = 1) +
        theme_minimal() +
        labs(
          x = "Period Start",
          y = "Vaccine Effectiveness",
          title = ifelse(
            is.null(title),
            paste(self$name, "model"),
            title
          ),
          subtitle = paste(
            "Ribbon shows",
            paste0(round((probs[length(probs)] - probs[1]) * 100), "%"),
            "prediction intervals"
          )
        )

      if (show_points) {
        p <- p + geom_point(color = color, size = 2)
      }

      return(p)
    },

    #' @description Plot VE with uncertainty bands from bootstrap
    #' @param estimate Log estimates (in the HR scale)
    #' @param ci_low Lower confidence interval bounds (in the HR scale)
    #' @param ci_up Upper confidence interval bounds (in the HR scale)
    #' @param periods_starts Vector of periods starts
    #' @param periods_ends Vector of periods ends
    #' @param title Plot title (default is NULL)
    #' @param color Color for the line and ribbon (default is "blue")
    #' @param alpha Transparency for the ribbon (default is 0.3)
    #' @param show_points Whether to show points on the line (default is TRUE)
    #' @param probs Vector of quantiles for confidence intervals (default is c(0.025, 0.975))
    #' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point theme_minimal labs
    #' @return ggplot object
    plot_VE_with_uncertainty_plus_raw_data = function(
      estimate,
      ci_low,
      ci_up,
      periods_starts,
      periods_ends,
      title = NULL,
      color = "blue",
      alpha = 0.3,
      show_points = TRUE,
      probs = c(0.025, 0.975)
    ) {
      # Validate parameter lengths
      validate_parameter_lengths(
        periods_starts = periods_starts,
        periods_ends = periods_ends,
        function_name = "plot_VE_with_uncertainty"
      )

      if (is.null(self$optimization_results)) {
        stop("No optimization results available. Run fit_3step() first.")
      }

      # Generate summary statistics
      ve_summary <- self$generate_VE_summary(
        periods_starts,
        periods_ends,
        probs
      )

      raw_ve <- data.frame(
        VE = 1 - estimate,
        VE_low = 1 - ci_up,
        VE_high = 1 - ci_low,
        period_start = periods_starts + 28
      )

      # Create plot
      p <- ggplot() +
        geom_ribbon(
          data = ve_summary,
          aes(
            x = period_start,
            y = VE_median,
            ymin = VE_quantile_025,
            ymax = VE_quantile_975
          ),
          alpha = alpha,
          fill = color
        ) +
        geom_line(
          data = ve_summary,
          aes(x = period_start, y = VE_median),
          color = color,
          linewidth = 1
        ) +
        geom_point(data = raw_ve, aes(x = period_start, y = VE)) +
        geom_errorbar(
          data = raw_ve,
          aes(x = period_start, ymin = VE_low, ymax = VE_high),
          width = 0.1
        ) +
        theme_minimal() +
        labs(
          x = "Period Start",
          y = "Vaccine Effectiveness",
          title = ifelse(
            is.null(title),
            paste(self$name, "model"),
            title
          ),
          subtitle = paste(
            "Ribbon shows",
            paste0(round((probs[length(probs)] - probs[1]) * 100), "%"),
            "prediction intervals\n",
            "Points and error bars show raw data"
          )
        )

      if (show_points) {
        p <- p + geom_point(color = color, size = 2)
      }

      return(p)
    },

    #' @description Plot individual bootstrap samples
    #' @param periods_starts Vector of periods starts
    #' @param periods_ends Vector of periods ends
    #' @param n_samples Number of bootstrap samples to show (default is 5)
    #' @param title Plot title (default is NULL)
    #' @param alpha Transparency for the lines (default is 0.7)
    #' @importFrom ggplot2 ggplot aes geom_line theme_minimal labs scale_color_brewer
    #' @importFrom dplyr filter
    #' @return ggplot object
    plot_bootstrap_samples = function(
      periods_starts,
      periods_ends,
      n_samples = 5,
      title = NULL,
      alpha = 0.7
    ) {
      # Validate parameter lengths
      validate_parameter_lengths(
        periods_starts = periods_starts,
        periods_ends = periods_ends,
        function_name = "plot_bootstrap_samples"
      )

      if (is.null(self$optimization_results)) {
        stop("No optimization results available. Run fit_3step() first.")
      }

      # Generate VE samples
      ve_samples <- self$generate_VE_based_on_boot_fit(
        periods_starts,
        periods_ends
      )

      # Filter to show only n_samples
      ve_samples_subset <- ve_samples %>%
        filter(boot_sample %in% 1:n_samples)

      best_params <- self$optimization_results$boot_fit$coefboot[4, ]
      pred <- decay_functions()$logistic(
        1:56,
        best_params[1],
        best_params[2],
        best_params[3]
      )
      ve <- 1 - exp(pred)
      plot(ve)

      # Create plot
      p <- ve_samples_subset %>%
        ggplot(aes(x = period_start, y = VE, color = factor(boot_sample))) +
        geom_line(alpha = alpha) +
        theme_minimal() +
        labs(
          x = "Period Start",
          y = "Vaccine Effectiveness",
          title = ifelse(
            is.null(title),
            paste("Individual Bootstrap Samples -", self$name, "Model"),
            title
          ),
          color = "Bootstrap Sample"
        ) +
        scale_color_brewer(palette = "Set1")

      return(p)
    },

    #' @description Get optimal parameters from last optimization
    #' @return Named vector of optimal parameters
    get_optimal_params = function() {
      if (is.null(self$optimization_results)) {
        stop("No optimization results available. Run fit_3step() first.")
      }

      # Extract parameters (excluding aic, bic and model columns)
      params <- self$optimization_results$params[,
        !names(self$optimization_results$params) %in% c("aic", "bic", "model")
      ]
      return(params)
    },

    #' @description Get AIC value from last optimization
    #' @return AIC value
    get_aic = function() {
      if (is.null(self$optimization_results)) {
        stop("No optimization results available. Run fit_3step() first.")
      }
      return(self$optimization_results$params$aic)
    },

    #' @description Get BIC value from last optimization
    #' @return BIC value
    get_bic = function() {
      if (is.null(self$optimization_results)) {
        stop("No optimization results available. Run fit_3step() first.")
      }
      return(self$optimization_results$params$bic)
    },

    #' @description Predict using optimal parameters
    #' @param estimate estimates (in the HR scale)
    #' @param periods_starts Vector of periods starts
    #' @param periods_ends Vector of periods ends
    #' @param aggregate_by_8_weeks Whether to aggregate the predictions by 8 weeks (default is FALSE)
    #' @param criterion Criterion to use for prediction (default is "aic")
    #' @return Vector of predictions
    predict_optimal = function(
      estimate,
      periods_starts,
      periods_ends,
      aggregate_by_8_weeks = FALSE,
      criterion = "aic"
    ) {
      # Validate parameter lengths
      validate_parameter_lengths(
        estimate = estimate,
        periods_starts = periods_starts,
        periods_ends = periods_ends,
        function_name = "predict_optimal"
      )

      optimal_params <- self$get_optimal_params()
      return(self$predict(
        optimal_params,
        estimate,
        periods_starts,
        periods_ends,
        aggregate_by_8_weeks,
        criterion
      ))
    },

    #' @description Predict decay model for given parameters
    #' @param params Named vector of parameters
    #' @param estimate estimates (in the HR scale)
    #' @param periods_starts Vector of periods starts
    #' @param periods_ends Vector of periods ends
    #' @param aggregate_by_8_weeks Whether to aggregate the predictions by 8 weeks (default is FALSE)
    #' @param criterion Criterion to use for prediction (default is "aic")
    #' @return Vector of predictions
    predict = function(
      params,
      estimate,
      periods_starts,
      periods_ends,
      aggregate_by_8_weeks = FALSE,
      criterion = "aic"
    ) {
      # Validate parameter lengths
      validate_parameter_lengths(
        estimate = estimate,
        periods_starts = periods_starts,
        periods_ends = periods_ends,
        function_name = "predict"
      )

      timing <- build_time_objects(periods_starts, periods_ends)$timing
      n_weeks <- build_time_objects(periods_starts, periods_ends)$n_weeks

      time_steps <- 1:n_weeks
      param_list <- as.list(params)

      log_estimate <- log(estimate)

      ve_t <- do.call(
        self$decay_function,
        c(list(time_steps = time_steps), param_list)
      )
      log_estimate_fitted <- log(1 - ve_t)

      if (!aggregate_by_8_weeks) {
        return(log_estimate_fitted)
      }

      predictions <- sapply(1:length(log_estimate), function(x) {
        start_week <- timing$diff_start[x] + 1
        end_week <- min(start_week + 7, n_weeks)
        mean(log_estimate_fitted[start_week:end_week])
      })

      return(predictions)
    },

    #' @description Plot using optimal parameters
    #' @param estimate Vector of estimates (in the HR scale)
    #' @param ci_low Lower confidence interval (in the HR scale)
    #' @param ci_up Upper confidence interval (in the HR scale)
    #' @param periods_starts Vector of periods starts
    #' @param periods_ends Vector of periods ends
    #' @param title Plot title
    #' @param color Color of the line
    #' @param age_group Age group
    #' @param site Site
    #' @param season Season
    #' @param outcome Outcome
    #' @importFrom ggplot2 ggplot aes geom_errorbar geom_line geom_point theme_minimal labs scale_color_brewer facet_wrap
    #' @return ggplot object
    plot_optimal = function(
      estimate,
      ci_low,
      ci_up,
      periods_starts,
      periods_ends,
      title = NULL,
      color = "red",
      age_group,
      site,
      season,
      outcome
    ) {
      # Validate parameter lengths
      validate_parameter_lengths(
        estimate = estimate,
        ci_low = ci_low,
        ci_up = ci_up,
        periods_starts = periods_starts,
        periods_ends = periods_ends,
        function_name = "plot_optimal"
      )

      optimal_params <- self$get_optimal_params()
      return(self$plot(
        optimal_params,
        estimate,
        ci_low,
        ci_up,
        periods_starts,
        periods_ends,
        title,
        color,
        age_group,
        site,
        season,
        outcome
      ))
    },

    #' @description Plot decay model
    #' @param params Named vector of parameters
    #' @param estimate Vector of estimates (in the HR scale)
    #' @param ci_low Lower confidence interval (in the HR scale)
    #' @param ci_up Upper confidence interval (in the HR scale)
    #' @param periods_starts Vector of periods starts
    #' @param periods_ends Vector of periods ends
    #' @param title Plot title
    #' @param color Color of the line (default is "red")
    #' @param age_group Age group
    #' @param site Site
    #' @param season Season
    #' @param outcome Outcome
    #' @param criterion Criterion to use for prediction (default is "aic")
    #' @importFrom ggplot2 ggplot aes geom_errorbar geom_line geom_point theme_minimal labs scale_color_brewer facet_wrap
    #' @return ggplot object
    plot = function(
      params,
      estimate,
      ci_low,
      ci_up,
      periods_starts,
      periods_ends,
      title = NULL,
      color = "red",
      age_group,
      site,
      season,
      outcome
    ) {
      # Validate parameter lengths
      validate_parameter_lengths(
        estimate = estimate,
        ci_low = ci_low,
        ci_up = ci_up,
        periods_starts = periods_starts,
        periods_ends = periods_ends,
        function_name = "plot"
      )

      timing <- build_time_objects(periods_starts, periods_ends)$timing
      n_weeks <- build_time_objects(periods_starts, periods_ends)$n_weeks

      # Generate fitted curve
      time_steps <- 1:n_weeks
      param_list <- as.list(params)

      ve_fit <- do.call(
        self$decay_function,
        c(list(time_steps = time_steps), param_list)
      )

      fit_data <- data.frame(
        x = as.Date(timing$date_min)[1] + time_steps * 7,
        y = ve_fit
      )

      # Generate predictions
      predictions <- self$predict(
        params,
        estimate,
        periods_starts,
        periods_ends,
        aggregate_by_8_weeks = TRUE,
        criterion = criterion
      )
      predictions_ve <- 1 - exp(predictions)

      pred_data <- data.frame(
        x = as.Date(timing$date_min) + (8 * 7 / 2),
        y = predictions_ve
      )

      data_ve <- data.frame(
        date_min = as.Date(timing$date_min) + (8 * 7 / 2),
        ve = 1 - estimate,
        ve_low = 1 - ci_up,
        ve_high = 1 - ci_low
      )

      if (is.null(title)) {
        title <- paste(self$name, "Decay Model")
      }

      ggplot(data_ve, aes(x = date_min, y = ve)) +
        geom_errorbar(
          aes(ymin = ve_low, ymax = ve_high),
          alpha = 0.75,
          linewidth = 1.2
        ) +
        geom_point(
          size = 1.5
        ) +
        geom_line(
          data = fit_data,
          aes(x = x, y = y),
          color = color,
          linewidth = 2
        ) +
        geom_point(
          data = pred_data,
          aes(x = x, y = y),
          color = color,
          size = 2
        ) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
        theme_minimal() +
        labs(
          x = "Date",
          y = "Vaccine Effectiveness",
          title = title,
          subtitle = paste(
            "Age group:",
            age_group,
            "Site:",
            site,
            "Season:",
            season,
            "Outcome:",
            outcome
          )
        )
    },

    #' @description Print model information
    print = function() {
      cat("DecayModel:", self$name, "\n")
      cat(
        "Parameters:",
        paste(self$param_config$names, collapse = ", "),
        "\n"
      )
      if (!is.null(self$optimization_results)) {
        cat("AIC:", self$get_aic(), "\n")
        cat("BIC:", self$get_bic(), "\n")
        cat(
          "Optimal parameters:",
          paste(
            names(self$get_optimal_params()),
            "=",
            round(self$get_optimal_params(), 4),
            collapse = ", "
          ),
          "\n"
        )
      } else {
        cat("No optimization results available\n")
      }
    },

    #' @description Build final dataset with weekly and 8-weeks predictions
    #' @param estimate Vector of estimates (in the HR scale)
    #' @param ci_low Vector of lower confidence interval bounds (in the HR scale)
    #' @param ci_up Vector of upper confidence interval bounds (in the HR scale)
    #' @param periods_starts Vector of periods starts
    #' @param periods_ends Vector of periods ends
    #' @param include_weekly Whether to include weekly predictions (default is TRUE)
    #' @param include_8_weeks Whether to include 8-weeks predictions (default is TRUE)
    #' @param criterion Criterion to use for prediction (default is "aic")
    #' @return List containing weekly_data and data_8_weeks data frames
    build_prediction_dataset = function(
      estimate,
      ci_low,
      ci_up,
      periods_starts,
      periods_ends,
      include_weekly = TRUE,
      include_8_weeks = TRUE,
      criterion = "aic"
    ) {
      # Validate parameter lengths
      validate_parameter_lengths(
        estimate = estimate,
        ci_low = ci_low,
        ci_up = ci_up,
        periods_starts = periods_starts,
        periods_ends = periods_ends,
        function_name = "build_prediction_dataset"
      )

      if (is.null(self$optimization_results)) {
        stop("No optimization results available. Run fit_3step() first.")
      }

      result <- list()

      # Get timing information
      time_parameters <- build_time_objects(periods_starts, periods_ends)
      n_weeks <- time_parameters$n_weeks

      if (include_weekly) {
        # Get weekly predictions
        weekly_predictions <- self$predict_optimal(
          estimate = estimate,
          periods_starts = periods_starts,
          periods_ends = periods_ends,
          aggregate_by_8_weeks = FALSE,
          criterion = criterion
        )

        # Convert to VE scale
        weekly_predictions_ve <- 1 - exp(weekly_predictions)

        # Create weekly data frame
        weekly_data <- data.frame(
          week = periods_starts[1] + (0:(n_weeks - 1) * 7),
          predicted_ve = weekly_predictions_ve,
          model_used = self$name
        )

        result$weekly_data <- weekly_data
      }

      if (include_8_weeks) {
        # Get 8-weeks predictions
        predictions_8_weeks <- self$predict_optimal(
          estimate = estimate,
          periods_starts = periods_starts,
          periods_ends = periods_ends,
          aggregate_by_8_weeks = TRUE,
          criterion = criterion
        )

        # Convert to VE scale
        predictions_ve_8_weeks <- 1 - exp(predictions_8_weeks)

        # Create 8-weeks data frame
        data_8_weeks <- data.frame(
          period_start = periods_starts,
          period_end = periods_ends,
          monthly_ve = 1 - estimate,
          monthly_ve_low = 1 - ci_up,
          monthly_ve_high = 1 - ci_low,
          predicted_ve = predictions_ve_8_weeks,
          model_used = self$name
        )

        result$data_8_weeks <- data_8_weeks
      }

      return(result)
    }
  )
)

# Create a custom model
# custom_model <- DecayModel$new(
#   name = "custom_decay",
#   decay_function = function(time_steps, VE0, rate, ...) {
#     VE0 * (1 - rate * time_steps)
#   },
#   param_config = list(
#     start = c(VE0 = 0.8, rate = 0.01),
#     lower = c(0.01, 0.001),
#     upper = c(0.99, 0.1),
#     names = c("VE0", "rate")
#   )
# )

#' Compare multiple decay models
#'
#' @param model_names Vector of model names to use (default is c("exponential", "logistic"))
#' @param estimate Vector of estimates (in the HR scale)
#' @param ci_low Vector of lower confidence interval bounds (in the HR scale)
#' @param ci_up Vector of upper confidence interval bounds (in the HR scale)
#' @param periods_starts Vector of periods starts
#' @param periods_ends Vector of periods ends
#' @param decay_funs List of decay functions to use (default is list(exponential = exponential_decay(), logistic = logistic_decay(), logistic_simple = logistic_simple_decay()))
#' @param param_configs List of parameter configurations to use (default is list(exponential = list(start = c(VE0 = 0.7, decay_rate = 0.1), lower = c(0.01, 0.001), upper = c(0.99, 1.0), names = c("VE0", "decay_rate")), logistic = list(start = c(VE0 = 0.7, decay_rate = 0.1, constant = 0.5), lower = c(0.01, 0.001, 0.01), upper = c(0.99, 1.0, 2.0), names = c("VE0", "decay_rate", "constant")), logistic_simple = list(start = c(VE0 = 0.7, decay_rate = 0.1), lower = c(0.01, 0.001), upper = c(0.99, 1.0), names = c("VE0", "decay_rate"))))
#' @param nboot Number of bootstrap samples (default is 100)
#' @param seed Seed to use for reproducibility (default is 123)
#' @importFrom purrr map list_rbind
#' @return List containing models, results, comparison, and summary
#' @export
#' @rdname compare_decay_models
compare_decay_models <- function(
  model_names = c("exponential", "logistic", "logistic_simple"),
  estimate,
  ci_low,
  ci_up,
  periods_starts,
  periods_ends,
  decay_funs = list(
    exponential = exponential_decay(),
    logistic = logistic_decay(),
    logistic_simple = logistic_simple_decay()
  ),
  param_configs = list(
    exponential = list(
      start = c(VE0 = 0.7, decay_rate = 0.1),
      lower = c(0.01, 0.001),
      upper = c(0.99, 1.0),
      names = c("VE0", "decay_rate")
    ),
    logistic = list(
      start = c(VE0 = 0.7, decay_rate = 0.1, constant = 0.5),
      lower = c(0.01, 0.001, 0.01),
      upper = c(0.99, 1.0, 2.0),
      names = c("VE0", "decay_rate", "constant")
    ),
    logistic_simple = list(
      start = c(VE0 = 0.7, decay_rate = 0.1),
      lower = c(0.01, 0.001),
      upper = c(0.99, 1.0),
      names = c("VE0", "decay_rate")
    )
  ),
  seed = 123,
  type_of_fit = "midpoint"
) {
  # Validate parameter lengths
  validate_parameter_lengths(
    estimate = estimate,
    ci_low = ci_low,
    ci_up = ci_up,
    periods_starts = periods_starts,
    periods_ends = periods_ends,
    function_name = "compare_decay_models"
  )

  # Create model instances
  models <- list()
  results <- list()

  for (model_name in model_names) {
    if (!model_name %in% names(decay_funs)) {
      warning(
        "Model '",
        model_name,
        "' not found in decay_functions. Skipping."
      )
      next
    }

    if (!model_name %in% names(param_configs)) {
      warning(
        "Parameter configuration for '",
        model_name,
        "' not found. Skipping."
      )
      next
    }

    # Create model instance
    models[[model_name]] <- DecayModel$new(
      name = model_name,
      decay_function = decay_funs[[model_name]],
      param_config = param_configs[[model_name]]
    )

    # Optimize the model
    results[[model_name]] <- models[[model_name]]$fit_to_determine_bic(
      estimate = estimate,
      ci_low = ci_low,
      ci_up = ci_up,
      periods_starts = periods_starts,
      periods_ends = periods_ends,
      seed = seed,
      type_of_fit = type_of_fit
    )
  }

  # Combine results for comparison
  all_results <- list_rbind(map(results, "params"))

  return(list(
    models = models,
    results = results,
    comparison = all_results,
    summary = all_results[, c("model", "aic", "bic")]
  ))
}

#' Plot model comparison
#'
#' @param comparison_result List containing models, results, comparison, and summary
#' @param estimate Vector of estimates (in the HR scale)
#' @param ci_low Vector of lower confidence interval bounds (in the HR scale)
#' @param ci_up Vector of upper confidence interval bounds (in the HR scale)
#' @param periods_starts Vector of periods starts
#' @param periods_ends Vector of periods ends
#' @param colors Vector of colors to use for the models (default is c("red", "blue", "green", "purple", "orange"))
#' @param criterion Criterion to use for model selection (default is "bic")
#' @param age_group Age group
#' @param site Site
#' @param season Season
#' @param outcome Outcome
#' @importFrom ggplot2 ggplot aes geom_errorbar geom_line geom_point theme_minimal labs scale_color_brewer facet_wrap
#' @return ggplot object
#' @export
#' @rdname compare_decay_models
plot_model_comparison <- function(
  comparison_result,
  estimate,
  ci_low,
  ci_up,
  periods_starts,
  periods_ends,
  colors = c("red", "blue", "green", "purple", "orange"),
  criterion = "bic",
  age_group,
  site,
  season,
  outcome
) {
  # Validate parameter lengths
  validate_parameter_lengths(
    estimate = estimate,
    ci_low = ci_low,
    ci_up = ci_up,
    periods_starts = periods_starts,
    periods_ends = periods_ends,
    function_name = "plot_model_comparison"
  )

  # Get timing information
  timing <- build_time_objects(periods_starts, periods_ends)$timing
  n_weeks <- build_time_objects(periods_starts, periods_ends)$n_weeks

  # Prepare data for plotting
  comparison_data <- data.frame()
  prediction_data <- data.frame()

  for (i in seq_along(comparison_result$models)) {
    model_name <- names(comparison_result$models)[i]
    model <- comparison_result$models[[model_name]]
    color <- colors[((i - 1) %% length(colors)) + 1]

    # Get optimal parameters
    optimal_params <- model$get_optimal_params()

    # Generate fitted curve
    time_steps <- 1:n_weeks
    param_list <- as.list(optimal_params)
    ve_fit <- do.call(
      model$decay_function,
      c(list(time_steps = time_steps), param_list)
    )

    model_data <- data.frame(
      x = as.Date(timing$date_min)[1] + time_steps * 7,
      y = ve_fit,
      model = model_name,
      type = "fitted_curve",
      color = color
    )

    # Generate predictions for observation periods
    predictions <- model$predict_optimal(
      estimate,
      periods_starts,
      periods_ends,
      aggregate_by_8_weeks = TRUE,
      criterion = criterion
    )
    predictions_ve <- 1 - exp(predictions)

    pred_data <- data.frame(
      x = as.Date(timing$date_min) + (8 * 7 / 2),
      y = predictions_ve,
      model = model_name,
      type = "predictions",
      color = color
    )

    comparison_data <- rbind(comparison_data, model_data, pred_data)
  }

  # Create data frame for observed values
  data_ve <- data.frame(
    date_min = as.Date(timing$date_min) + (8 * 7 / 2),
    ve = 1 - estimate,
    ve_low = 1 - ci_up,
    ve_high = 1 - ci_low
  )

  # Create the plot
  p <- ggplot(data_ve, aes(x = date_min, y = ve)) +
    geom_errorbar(
      aes(ymin = ve_low, ymax = ve_high),
      alpha = 0.5
    ) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(
      x = "Date",
      y = "Vaccine Effectiveness",
      title = paste(
        "Comparison of Decay Models for",
        age_group,
        site,
        season,
        outcome
      ),
      subtitle = "Lines show fitted curves, colored points show weekly model predictions, colored squares show 8-weeks model predictions.\nVaccine effectiveness from Cox models are shown in black.\n8-weeks values are shown at the middle of the 8-weeks period.",
      caption = paste(
        "Best fit is provided by the",
        comparison_result$summary$model[which.min(
          comparison_result$summary[[criterion]]
        )],
        "model (",
        criterion,
        " =",
        round(min(comparison_result$summary[[criterion]]), 4),
        "vs",
        round(max(comparison_result$summary[[criterion]]), 4),
        "for worst model)."
      )
    )

  # Add fitted curves
  for (i in seq_along(comparison_result$models)) {
    model_name <- names(comparison_result$models)[i]
    color <- colors[((i - 1) %% length(colors)) + 1]

    curve_data <- comparison_data[
      comparison_data$model == model_name &
        comparison_data$type == "fitted_curve",
    ]
    pred_data <- comparison_data[
      comparison_data$model == model_name &
        comparison_data$type == "predictions",
    ]

    p <- p +
      geom_line(
        data = curve_data,
        aes(x = x, y = y),
        color = color,
        linewidth = 1,
        alpha = 0.5
      ) +
      geom_point(
        data = curve_data,
        aes(x = x, y = y),
        color = color,
        size = 1,
        alpha = 0.8
      ) +
      geom_point(
        data = pred_data,
        aes(x = x, y = y),
        color = color,
        size = 4,
        alpha = 0.8,
        shape = 15
      )
  }
  p <- p +
    scale_color_brewer(type = "qual", palette = "Set1") +
    facet_wrap(~model)

  return(p)
}

#' Get the best model based on lowest AIC
#'
#' @param comparison_result List containing models, results, comparison, and summary
#' @param criterion Criterion to use for best model selection (default is "aic")
#' @return List containing best model name, instance, AIC, and all results
#' @export
#' @rdname compare_decay_models
get_best_model <- function(comparison_result, criterion = "aic") {
  if (!is.list(comparison_result) || !"summary" %in% names(comparison_result)) {
    stop("comparison_result must be the output of compare_decay_models()")
  }

  if (!criterion %in% c("aic", "bic")) {
    stop("criterion must be either 'aic' or 'bic'")
  }

  # Find the model with the lowest AIC
  best_model_name <- comparison_result$summary$model[which.min(
    comparison_result$summary[[criterion]]
  )]
  best_criterion <- min(comparison_result$summary[[criterion]])

  # Get the best model instance
  best_model <- comparison_result$models[[best_model_name]]

  return(list(
    model_name = best_model_name,
    model = best_model,
    criterion = best_criterion,
    all_results = comparison_result$summary
  ))
}

#' Get the best model name only
#'
#' @param comparison_result List containing models, results, comparison, and summary
#' @param criterion Criterion to use for best model selection (default is "aic")
#' @return Best model name
#' @export
#' @rdname compare_decay_models
get_best_model_name <- function(comparison_result, criterion = "aic") {
  if (!is.list(comparison_result) || !"summary" %in% names(comparison_result)) {
    stop("comparison_result must be the output of compare_decay_models()")
  }

  return(comparison_result$summary$model[which.min(
    comparison_result$summary[[criterion]]
  )])
}

#' Get the best model instance only
#'
#' @param comparison_result List containing models, results, comparison, and summary
#' @param criterion Criterion to use for best model selection (default is "aic")
#' @return Best model instance
#' @export
#' @rdname compare_decay_models
get_best_model_instance <- function(comparison_result, criterion = "aic") {
  if (!is.list(comparison_result) || !"summary" %in% names(comparison_result)) {
    stop("comparison_result must be the output of compare_decay_models()")
  }

  best_model_name <- comparison_result$summary$model[which.min(
    comparison_result$summary[[criterion]]
  )]
  return(comparison_result$models[[best_model_name]])
}
