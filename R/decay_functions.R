#' List of decay functions
#'
#' @return List of decay functions
#' @export
#' @rdname decay_functions
decay_functions <- function() {
  list(
    exponential = function(time_steps, VE0, decay_rate, ...) {
      VE0 * exp(-decay_rate * time_steps)
    },
    logistic = function(time_steps, VE0, decay_rate, constant, ...) {
      VE0 / (1 + constant * exp(decay_rate * time_steps))
    },
    logistic_simple = function(time_steps, VE0, decay_rate, ...) {
      VE0 / (1 + exp(decay_rate * time_steps))
    }
  )
}

#' Logistic decay function
#'
#' @return Logistic decay function
#' @export
#' @rdname decay_functions
logistic_decay <- function() {
  decay_functions()$logistic
}

#' Simplified logistic decay function (without constant parameter)
#'
#' @return Simplified logistic decay function
#' @export
#' @rdname decay_functions
logistic_simple_decay <- function() {
  decay_functions()$logistic_simple
}

#' Exponential decay function
#'
#' @return Exponential decay function
#' @export
#' @rdname decay_functions
exponential_decay <- function() {
  decay_functions()$exponential
}

#' List of period-averaged decay functions (log(1-VE) scale)
#'
#' These functions return the mean log(1-VE) over 8-week periods, as used in the toy code.
#' @return Vector of averaged log(1-VE) values for each period
#' @export
#' @rdname decay_functions_period_avg
decay_functions_period_avg <- function() {
  # Get the base decay functions
  base_decay_funs <- decay_functions()

  list(
    exponential_average = function(w_start, nweeks, VE0, decay_rate, ...) {
      time <- 1:max(nweeks)
      # Use the base exponential decay function
      ve_t <- base_decay_funs$exponential(time, VE0, decay_rate)
      weekly_est <- log(1 - ve_t)
      sapply(seq_along(w_start), function(x) {
        mean(weekly_est[(w_start[x] + 1):(w_start[x] + 8)])
      })
    },
    exponential_midpoint = function(w_start, nweeks, VE0, decay_rate, ...) {
      sapply(seq_along(w_start), function(x) {
        log(
          1 -
            base_decay_funs$exponential(
              w_start[x] + 4.5,
              VE0,
              decay_rate
            )
        )
      })
    },
    logistic_average = function(
      w_start,
      nweeks,
      VE0,
      decay_rate,
      constant,
      ...
    ) {
      time <- 1:max(nweeks)
      # Use the base logistic decay function
      ve_t <- base_decay_funs$logistic(time, VE0, decay_rate, constant)
      weekly_est <- log(1 - ve_t)
      sapply(seq_along(w_start), function(x) {
        mean(weekly_est[(w_start[x] + 1):(w_start[x] + 8)])
      })
    },

    logistic_midpoint = function(
      w_start,
      nweeks,
      VE0,
      decay_rate,
      constant,
      ...
    ) {
      sapply(seq_along(w_start), function(x) {
        log(
          1 -
            base_decay_funs$logistic(
              w_start[x] + 4.5,
              VE0,
              decay_rate,
              constant
            )
        )
      })
    }
  )
}

#' Logistic period-averaged decay function (log(1-VE) scale)
#' @export
#' @rdname decay_functions_period_avg
logistic_decay_period_avg <- function() {
  decay_functions_period_avg()$logistic
}

#' Simplified logistic period-averaged decay function (log(1-VE) scale)
#' @export
#' @rdname decay_functions_period_avg
logistic_simple_decay_period_avg <- function() {
  decay_functions_period_avg()$logistic_simple
}

#' Exponential period-averaged decay function (log(1-VE) scale)
#' @export
#' @rdname decay_functions_period_avg
exponential_decay_period_avg <- function() {
  decay_functions_period_avg()$exponential
}
