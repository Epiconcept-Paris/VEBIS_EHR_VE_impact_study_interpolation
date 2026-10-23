fit_exponential_model <- function() {
  data("sample_8_weeks_data", package = "VEinterpolation", envir = environment())
  d <- sample_8_weeks_data

  model <- DecayModel$new(
    name = "exponential",
    decay_function = exponential_decay(),
    param_config = list(
      start = c(VE0 = 0.7, decay_rate = 0.1),
      lower = c(0.01, 0.001),
      upper = c(0.99, 1.0),
      names = c("VE0", "decay_rate")
    )
  )
  model$fit_to_determine_bic(
    estimate = d$estimate,
    ci_low = d$CIlow,
    ci_up = d$CIhigh,
    periods_starts = d$date_min,
    periods_ends = d$date_max,
    seed = 666
  )
  model$generate_VE_based_on_resampling(
    estimate = d$estimate,
    ci_low = d$CIlow,
    ci_up = d$CIhigh,
    periods_starts = d$date_min,
    periods_ends = d$date_max,
    n_weekly_ve_to_generate = 50
  )
  list(model = model, data = d)
}

test_that("noise-free VE curves are never negative", {
  fit <- fit_exponential_model()
  clean <- purrr::compact(fit$model$nlrob_samples_clean)

  expect_gt(length(clean), 0)
  for (curve in clean) {
    expect_true(all(curve >= 0))
    expect_true(all(curve <= 1))
  }
})

test_that("noise-free VE curves decrease over time for the exponential model", {
  fit <- fit_exponential_model()
  clean <- purrr::compact(fit$model$nlrob_samples_clean)

  for (curve in clean) {
    expect_true(all(diff(curve) <= 0))
  }
})

test_that("clean and noisy samples are aligned and the noise is what differs", {
  fit <- fit_exponential_model()

  expect_length(fit$model$nlrob_samples_clean, length(fit$model$nlrob_samples))
  expect_equal(
    lengths(fit$model$nlrob_samples_clean),
    lengths(fit$model$nlrob_samples)
  )

  ok <- !vapply(fit$model$nlrob_samples, is.null, logical(1))
  clean <- fit$model$nlrob_samples_clean[ok][[1]]
  noisy <- fit$model$nlrob_samples[ok][[1]]
  expect_false(isTRUE(all.equal(clean, noisy)))
})

test_that("noisy VE is bounded above by 1 (negative values only come from the noise term)", {
  fit <- fit_exponential_model()
  clean <- purrr::compact(fit$model$nlrob_samples_clean)
  noisy <- purrr::compact(fit$model$nlrob_samples)

  # VE = 1 - exp(log(1 - VE_clean) + eps): VE < 0 requires eps > -log(1 - VE_clean)
  expect_true(all(unlist(noisy) <= 1))
})

test_that("generate_VE_summary exposes clean columns inside the noisy interval on average", {
  fit <- fit_exponential_model()
  d <- fit$data
  summary <- fit$model$generate_VE_summary(d$date_min, d$date_max)

  expect_true(all(
    c("VE_clean_median", "VE_clean_quantile_025", "VE_clean_quantile_975") %in%
      names(summary)
  ))
  expect_true(all(summary$VE_clean_quantile_025 >= 0))
  expect_true(all(summary$VE_clean_quantile_025 <= summary$VE_clean_median))
  expect_true(all(summary$VE_clean_median <= summary$VE_clean_quantile_975))
  # Measurement error widens the interval
  expect_gt(
    mean(summary$VE_quantile_975 - summary$VE_quantile_025),
    mean(summary$VE_clean_quantile_975 - summary$VE_clean_quantile_025)
  )
})

test_that("plots work with show_clean TRUE and FALSE", {
  fit <- fit_exponential_model()
  d <- fit$data

  p_clean <- fit$model$plot_VE_with_uncertainty(d$date_min, d$date_max)
  p_old <- fit$model$plot_VE_with_uncertainty(
    d$date_min,
    d$date_max,
    show_clean = FALSE
  )
  p_raw <- fit$model$plot_VE_with_uncertainty_plus_raw_data(
    d$estimate, d$CIlow, d$CIhigh, d$date_min, d$date_max
  )

  expect_s3_class(p_clean, "ggplot")
  expect_s3_class(p_old, "ggplot")
  expect_s3_class(p_raw, "ggplot")
  # line + 2 ribbons + points vs line + 1 ribbon + points
  expect_length(p_clean$layers, 4)
  expect_length(p_old$layers, 3)
  expect_no_error(ggplot2::ggplot_build(p_clean))
})

test_that("show_clean = TRUE errors if no clean samples were generated", {
  fit <- fit_exponential_model()
  d <- fit$data
  summary <- fit$model$generate_VE_summary(d$date_min, d$date_max)
  summary$VE_clean_median <- NULL

  expect_error(
    add_ve_layers(
      ggplot2::ggplot(), summary, "blue", 0.3,
      show_clean = TRUE, show_points = FALSE
    ),
    "noise-free"
  )
})
