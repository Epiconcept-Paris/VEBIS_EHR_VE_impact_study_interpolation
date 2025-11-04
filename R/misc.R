#' Get standard error from confidence interval
#'
#' @param ci_low Lower confidence interval
#' @param ci_high Upper confidence interval
#' @return Standard error
#' @export
get_se_from_ci <- function(ci_low, ci_high) {
  (ci_high - ci_low) / 3.92
}

#' Build time objects for interpolation
#'
#' @param periods_starts Vector of periods starts
#' @param periods_ends Vector of periods ends
#' @return A list with the timing and the number of weeks
#' @export
build_time_objects <- function(periods_starts, periods_ends) {
  timing <- data.frame(
    date_min = periods_starts,
    date_max = periods_ends,
    diff_start = c(
      0,
      floor(as.numeric(difftime(
        periods_starts[-1],
        periods_starts[1],
        units = "weeks"
      )))
    )
  )

  n_weeks <- round(as.numeric(difftime(
    max(timing$date_max),
    min(timing$date_min),
    units = "weeks"
  )))

  return(list(timing = timing, n_weeks = n_weeks))
}
