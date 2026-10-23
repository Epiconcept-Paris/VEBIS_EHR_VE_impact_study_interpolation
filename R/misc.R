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

#' Add VE line and ribbons to a plot
#'
#' With `show_clean = TRUE`, the line is the noise-free interpolated VE, the
#' inner ribbon is its bootstrap interval and the outer ribbon is the interval
#' of the samples with measurement error added. With `show_clean = FALSE`, the
#' line is the median of the noisy samples and a single ribbon is drawn.
#'
#' @param p ggplot object to add the layers to
#' @param ve_summary Output of `generate_VE_summary()`
#' @param color Color for the line and ribbons
#' @param alpha Transparency for the ribbons
#' @param show_clean Whether to separate the noise-free VE from the
#'   measurement error envelope
#' @param show_points Whether to show points on the line
#' @return ggplot object
#' @importFrom ggplot2 aes geom_ribbon geom_line geom_point
#' @noRd
add_ve_layers <- function(p, ve_summary, color, alpha, show_clean, show_points) {
  if (show_clean && !"VE_clean_median" %in% names(ve_summary)) {
    stop("No noise-free VE available. Re-run generate_VE_based_on_resampling().")
  }
  ve_summary$VE_line <- if (show_clean) {
    ve_summary$VE_clean_median
  } else {
    ve_summary$VE_median
  }

  p <- p +
    geom_ribbon(
      data = ve_summary,
      aes(
        x = period_start,
        ymin = VE_quantile_025,
        ymax = VE_quantile_975
      ),
      alpha = alpha,
      fill = color
    )

  if (show_clean) {
    p <- p +
      geom_ribbon(
        data = ve_summary,
        aes(
          x = period_start,
          ymin = VE_clean_quantile_025,
          ymax = VE_clean_quantile_975
        ),
        alpha = alpha,
        fill = color
      )
  }

  p <- p +
    geom_line(
      data = ve_summary,
      aes(x = period_start, y = VE_line),
      color = color,
      linewidth = 1
    )

  if (show_points) {
    p <- p +
      geom_point(
        data = ve_summary,
        aes(x = period_start, y = VE_line),
        color = color,
        size = 2
      )
  }

  p
}

#' Subtitle for VE plots with uncertainty
#'
#' @param probs Vector of quantiles
#' @param show_clean Whether the noise-free VE is drawn separately
#' @return Character string
#' @noRd
ve_plot_subtitle <- function(probs, show_clean) {
  level <- paste0(round((probs[length(probs)] - probs[1]) * 100), "%")
  if (show_clean) {
    paste0(
      "Line: interpolated VE. Dark ribbon: ", level, " interval of the fit.\n",
      "Light ribbon: ", level, " prediction interval incl. measurement error"
    )
  } else {
    paste("Ribbon shows", level, "prediction intervals")
  }
}
