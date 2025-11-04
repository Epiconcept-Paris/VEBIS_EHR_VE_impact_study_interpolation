# Unit tests for DecayModel R6 class

# Load test data
data("sample_8_weeks_data")

# Test data setup
test_data <- sample_8_weeks_data
test_estimate <- test_data$estimate
test_ci_low <- test_data$CIlow
test_ci_up <- test_data$CIhigh
test_periods_starts <- test_data$date_min
test_periods_ends <- test_data$date_max

# Helper function to create a test model
create_test_model <- function(model_name = "exponential") {
  decay_funs <- list(
    exponential = exponential_decay(),
    logistic = logistic_decay(),
    logistic_simple = logistic_simple_decay()
  )

  param_configs <- list(
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
  )

  DecayModel$new(
    name = model_name,
    decay_function = decay_funs[[model_name]],
    param_config = param_configs[[model_name]]
  )
}

# Test initialization
test_that("DecayModel initializes correctly", {
  model <- create_test_model("exponential")

  expect_s3_class(model, "DecayModel")
  expect_equal(model$name, "exponential")
  expect_type(model$decay_function, "closure")
  expect_type(model$param_config, "list")
  expect_null(model$optimization_results)

  # Check parameter configuration structure
  expect_named(model$param_config, c("start", "lower", "upper", "names"))
  expect_equal(length(model$param_config$start), 2)
  expect_equal(length(model$param_config$lower), 2)
  expect_equal(length(model$param_config$upper), 2)
  expect_equal(length(model$param_config$names), 2)
})

test_that("DecayModel initializes with logistic model", {
  model <- create_test_model("logistic")

  expect_equal(model$name, "logistic")
  expect_equal(length(model$param_config$start), 3)
  expect_equal(length(model$param_config$lower), 3)
  expect_equal(length(model$param_config$upper), 3)
  expect_equal(length(model$param_config$names), 3)
})

test_that("DecayModel initializes with logistic_simple model", {
  model <- create_test_model("logistic_simple")

  expect_equal(model$name, "logistic_simple")
  expect_equal(length(model$param_config$start), 2)
  expect_equal(length(model$param_config$lower), 2)
  expect_equal(length(model$param_config$upper), 2)
  expect_equal(length(model$param_config$names), 2)
})

# Test evaluate method
test_that("evaluate method works correctly", {
  model <- create_test_model("exponential")

  time_steps <- 1:10
  params <- c(VE0 = 0.8, decay_rate = 0.1)

  result <- model$evaluate(time_steps, params)

  expect_equal(length(result), 10)
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
  expect_true(result[1] > result[10]) # Decay should decrease over time
})

test_that("evaluate method works with logistic model", {
  model <- create_test_model("logistic")

  time_steps <- 1:10
  params <- c(VE0 = 0.8, decay_rate = 0.1, constant = 0.5)

  result <- model$evaluate(time_steps, params)

  expect_equal(length(result), 10)
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})

test_that("evaluate method works with logistic_simple model", {
  model <- create_test_model("logistic_simple")

  time_steps <- 1:10
  params <- c(VE0 = 0.8, decay_rate = 0.1)

  result <- model$evaluate(time_steps, params)

  expect_equal(length(result), 10)
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})

# Test create_nls_formula method
test_that("create_nls_formula works for exponential model", {
  model <- create_test_model("exponential")

  formula <- model$create_nls_formula()

  expect_s3_class(formula, "formula")
  expect_equal(as.character(formula)[2], "estimates")
  expect_true(grepl("decay_func", as.character(formula)[3]))
})

test_that("create_nls_formula works for logistic model", {
  model <- create_test_model("logistic")

  formula <- model$create_nls_formula()

  expect_s3_class(formula, "formula")
  expect_equal(as.character(formula)[2], "estimates")
  expect_true(grepl("decay_func", as.character(formula)[3]))
})

# Test fit_3step method
test_that("fit_3step works with exponential model", {
  model <- create_test_model("exponential")

  # Use minimal bootstrap for testing
  result <- model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10, # Small number for testing
    seed = 123
  )

  expect_type(result, "list")
  expect_named(
    result,
    c(
      "params",
      "nlxb_fit",
      "nls_fit",
      "boot_fit",
      "covariance_matrix",
      "coefficients",
      "AIC",
      "BIC",
      "summary",
      "decay_func"
    )
  )
  expect_s3_class(result$params, "data.frame")
  expect_true(nrow(result$params) == 1)
  expect_true("aic" %in% names(result$params))
  expect_true("bic" %in% names(result$params))
  expect_true("model" %in% names(result$params))
})

test_that("fit_3step works with logistic model", {
  model <- create_test_model("logistic")

  result <- model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  expect_type(result, "list")
  expect_s3_class(result$params, "data.frame")
  expect_true(nrow(result$params) == 1)
})

test_that("fit_3step works with logistic_simple model", {
  model <- create_test_model("logistic_simple")

  result <- model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  expect_type(result, "list")
  expect_s3_class(result$params, "data.frame")
  expect_true(nrow(result$params) == 1)
})

test_that("fit_3step fails with unknown model name", {
  # Create a model with unknown name
  model <- DecayModel$new(
    name = "unknown_model",
    decay_function = exponential_decay(),
    param_config = list(
      start = c(VE0 = 0.7, decay_rate = 0.1),
      lower = c(0.01, 0.001),
      upper = c(0.99, 1.0),
      names = c("VE0", "decay_rate")
    )
  )

  expect_error(
    model$fit_3step(
      estimate = test_estimate,
      ci_low = test_ci_low,
      ci_up = test_ci_up,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends,
      nboot = 10,
      seed = 123
    ),
    "Unknown model name: unknown_model"
  )
})

# Test get_optimal_params method
test_that("get_optimal_params works after fitting", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  params <- model$get_optimal_params()

  expect_type(params, "list")
  expect_named(params, c("VE0", "decay_rate"))
  expect_true(all(params >= 0))
})

test_that("get_optimal_params fails without fitting", {
  model <- create_test_model("exponential")

  expect_error(
    model$get_optimal_params(),
    "No optimization results available. Run fit_3step\\(\\) first."
  )
})

# Test get_aic and get_bic methods
test_that("get_aic and get_bic work after fitting", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  aic <- model$get_aic()
  bic <- model$get_bic()

  expect_type(aic, "double")
  expect_type(bic, "double")
})

test_that("get_aic and get_bic fail without fitting", {
  model <- create_test_model("exponential")

  expect_error(
    model$get_aic(),
    "No optimization results available. Run fit_3step\\(\\) first."
  )

  expect_error(
    model$get_bic(),
    "No optimization results available. Run fit_3step\\(\\) first."
  )
})

# Test predict method
test_that("predict method works", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  params <- model$get_optimal_params()

  # Test weekly predictions
  weekly_pred <- model$predict(
    params = params,
    estimate = test_estimate,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    aggregate_by_8_weeks = FALSE
  )

  expect_equal(
    length(weekly_pred),
    build_time_objects(test_periods_starts, test_periods_ends)$n_weeks
  )

  # Test 8-weeks predictions
  pred_8_weeks <- model$predict(
    params = params,
    estimate = test_estimate,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    aggregate_by_8_weeks = TRUE
  )

  expect_equal(length(pred_8_weeks), length(test_estimate))
})

# Test predict_optimal method
test_that("predict_optimal method works", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  # Test weekly predictions
  weekly_pred <- model$predict_optimal(
    estimate = test_estimate,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    aggregate_by_8_weeks = FALSE
  )

  expect_equal(
    length(weekly_pred),
    build_time_objects(test_periods_starts, test_periods_ends)$n_weeks
  )

  # Test 8-weeks predictions
  pred_8_weeks <- model$predict_optimal(
    estimate = test_estimate,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    aggregate_by_8_weeks = TRUE
  )

  expect_equal(length(pred_8_weeks), length(test_estimate))
})

test_that("predict_optimal fails without fitting", {
  model <- create_test_model("exponential")

  expect_error(
    model$predict_optimal(
      estimate = test_estimate,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends
    ),
    "No optimization results available. Run fit_3step\\(\\) first."
  )
})

# Test generate_VE_based_on_boot_fit method
test_that("generate_VE_based_on_boot_fit works", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  ve_samples <- model$generate_VE_based_on_boot_fit(
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends
  )

  expect_s3_class(ve_samples, "data.frame")
  expect_named(ve_samples, c("period_start", "period_end", "VE", "ve_sample"))
  expect_true(nrow(ve_samples) > 0)
  expect_true(all(ve_samples$VE >= 0))
  expect_true(all(ve_samples$VE <= 1))
})

test_that("generate_VE_based_on_boot_fit fails without fitting", {
  model <- create_test_model("exponential")

  expect_error(
    model$generate_VE_based_on_boot_fit(
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends
    ),
    "No optimization results available. Run fit_3step\\(\\) first."
  )
})

# Test generate_VE_summary method
test_that("generate_VE_summary works", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  ve_summary <- model$generate_VE_summary(
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends
  )

  expect_s3_class(ve_summary, "data.frame")
  expect_named(
    ve_summary,
    c("period_start", "VE_median", "VE_mean", "VE_sd", "ci_low", "ci_up")
  )
  expect_true(nrow(ve_summary) > 0)
  expect_true(all(ve_summary$VE_median >= 0))
  expect_true(all(ve_summary$VE_median <= 1))
})

test_that("generate_VE_summary works with custom quantiles", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  ve_summary <- model$generate_VE_summary(
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    probs = c(0.1, 0.5, 0.9)
  )

  expect_s3_class(ve_summary, "data.frame")
  expect_named(
    ve_summary,
    c("period_start", "VE_median", "VE_mean", "VE_sd", "ci_low", "ci_up")
  )
})

test_that("generate_VE_summary fails without fitting", {
  model <- create_test_model("exponential")

  expect_error(
    model$generate_VE_summary(
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends
    ),
    "No optimization results available. Run fit_3step\\(\\) first."
  )
})

# Test plot methods
test_that("plot_VE_with_uncertainty works", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  p <- model$plot_VE_with_uncertainty(
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_VE_with_uncertainty works with custom parameters", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  p <- model$plot_VE_with_uncertainty(
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    title = "Custom Title",
    color = "red",
    alpha = 0.5,
    show_points = FALSE,
    probs = c(0.1, 0.9)
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_VE_with_uncertainty fails without fitting", {
  model <- create_test_model("exponential")

  expect_error(
    model$plot_VE_with_uncertainty(
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends
    ),
    "No optimization results available. Run fit_3step\\(\\) first."
  )
})


# Test plot and plot_optimal methods
test_that("plot method works", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  params <- model$get_optimal_params()

  p <- model$plot(
    params = params,
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    title = "Test Plot",
    color = "blue",
    age_group = "65-79",
    site = "Denmark",
    season = "2023-2024",
    outcome = "death"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_optimal method works", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  p <- model$plot_optimal(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    title = "Test Plot",
    color = "green",
    age_group = "65-79",
    site = "Denmark",
    season = "2023-2024",
    outcome = "death"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_optimal fails without fitting", {
  model <- create_test_model("exponential")

  expect_error(
    model$plot_optimal(
      estimate = test_estimate,
      ci_low = test_ci_low,
      ci_up = test_ci_up,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends
    ),
    "No optimization results available. Run fit_3step\\(\\) first."
  )
})

# Test build_prediction_dataset method
test_that("build_prediction_dataset works with both weekly and 8-weeks", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  result <- model$build_prediction_dataset(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    include_weekly = TRUE,
    include_8_weeks = TRUE
  )

  expect_type(result, "list")
  expect_named(result, c("weekly_data", "data_8_weeks"))
  expect_s3_class(result$weekly_data, "data.frame")
  expect_s3_class(result$data_8_weeks, "data.frame")
  expect_true(nrow(result$weekly_data) > 0)
  expect_true(nrow(result$data_8_weeks) > 0)
})

test_that("build_prediction_dataset works with weekly only", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  result <- model$build_prediction_dataset(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    include_weekly = TRUE,
    include_8_weeks = FALSE
  )

  expect_type(result, "list")
  expect_named(result, "weekly_data")
  expect_s3_class(result$weekly_data, "data.frame")
  expect_true(nrow(result$weekly_data) > 0)
})

test_that("build_prediction_dataset works with 8-weeks only", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  result <- model$build_prediction_dataset(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    include_weekly = FALSE,
    include_8_weeks = TRUE
  )

  expect_type(result, "list")
  expect_named(result, "data_8_weeks")
  expect_s3_class(result$data_8_weeks, "data.frame")
  expect_true(nrow(result$data_8_weeks) > 0)
})

test_that("build_prediction_dataset fails without fitting", {
  model <- create_test_model("exponential")

  expect_error(
    model$build_prediction_dataset(
      estimate = test_estimate,
      ci_low = test_ci_low,
      ci_up = test_ci_up,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends
    ),
    "No optimization results available. Run fit_3step\\(\\) first."
  )
})

# Test print method
test_that("print method works", {
  model <- create_test_model("exponential")

  # Test print without fitting
  expect_output(print(model), "DecayModel: exponential")
  expect_output(print(model), "Parameters: VE0, decay_rate")
  expect_output(print(model), "No optimization results available")

  # Fit the model and test print with results
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  expect_output(print(model), "DecayModel: exponential")
  expect_output(print(model), "Parameters: VE0, decay_rate")
  expect_output(print(model), "AIC:")
  expect_output(print(model), "BIC:")
  expect_output(print(model), "Optimal parameters:")
})

# Test edge cases and error handling
test_that("fit_3step handles edge cases", {
  model <- create_test_model("exponential")

  # Test with very small bootstrap sample
  result <- model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 1,
    seed = 123
  )

  expect_type(result, "list")
  expect_s3_class(result$params, "data.frame")
})

test_that("fit_3step handles different seeds", {
  model1 <- create_test_model("exponential")
  model2 <- create_test_model("exponential")

  # Fit with different seeds
  result1 <- model1$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  result2 <- model2$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 456
  )

  expect_type(result1, "list")
  expect_type(result2, "list")
  expect_s3_class(result1$params, "data.frame")
  expect_s3_class(result2$params, "data.frame")
})

# Test data validation
test_that("fit_3step validates input data", {
  model <- create_test_model("exponential")

  # Test with mismatched lengths
  expect_error(
    model$fit_3step(
      estimate = test_estimate[1:3],
      ci_low = test_ci_low,
      ci_up = test_ci_up,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends,
      nboot = 10,
      seed = 123
    )
  )
})

# Test custom model creation (commented out in original code)
test_that("custom model can be created", {
  custom_model <- DecayModel$new(
    name = "custom_decay",
    decay_function = function(time_steps, VE0, rate, ...) {
      VE0 * (1 - rate * time_steps)
    },
    param_config = list(
      start = c(VE0 = 0.8, rate = 0.01),
      lower = c(0.01, 0.001),
      upper = c(0.99, 0.1),
      names = c("VE0", "rate")
    )
  )

  expect_s3_class(custom_model, "DecayModel")
  expect_equal(custom_model$name, "custom_decay")
  expect_type(custom_model$decay_function, "closure")

  # Test evaluate method with custom model
  time_steps <- 1:5
  params <- c(VE0 = 0.8, rate = 0.01)
  result <- custom_model$evaluate(time_steps, params)

  expect_equal(length(result), 5)
  expect_true(all(result >= 0))
})

# Test comparison functions
test_that("compare_decay_models works with default models", {
  result <- compare_decay_models(
    model_names = c("exponential", "logistic"),
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  expect_type(result, "list")
  expect_named(result, c("models", "results", "comparison", "summary"))
  expect_type(result$models, "list")
  expect_type(result$results, "list")
  expect_s3_class(result$comparison, "data.frame")
  expect_s3_class(result$summary, "data.frame")
  expect_equal(nrow(result$summary), 2)
  expect_true(all(c("exponential", "logistic") %in% result$summary$model))
})

test_that("compare_decay_models works with all three models", {
  result <- compare_decay_models(
    model_names = c("exponential", "logistic", "logistic_simple"),
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  expect_type(result, "list")
  expect_named(result, c("models", "results", "comparison", "summary"))
  expect_equal(nrow(result$summary), 3)
  expect_true(all(
    c("exponential", "logistic", "logistic_simple") %in% result$summary$model
  ))
})

test_that("compare_decay_models handles missing models gracefully", {
  expect_warning(
    result <- compare_decay_models(
      model_names = c("exponential", "nonexistent_model"),
      estimate = test_estimate,
      ci_low = test_ci_low,
      ci_up = test_ci_up,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends,
      nboot = 10,
      seed = 123
    )
  )

  expect_type(result, "list")
  expect_named(result, c("models", "results", "comparison", "summary"))
  expect_equal(nrow(result$summary), 1)
  expect_equal(result$summary$model, "exponential")
})

test_that("plot_model_comparison works", {
  comparison_result <- compare_decay_models(
    model_names = c("exponential", "logistic"),
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  p <- plot_model_comparison(
    comparison_result = comparison_result,
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    colors = c("red", "blue"),
    criterion = "aic",
    age_group = "65-79",
    site = "Denmark",
    season = "2023-2024",
    outcome = "death"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_model_comparison works with different criterion", {
  comparison_result <- compare_decay_models(
    model_names = c("exponential", "logistic"),
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  p <- plot_model_comparison(
    comparison_result = comparison_result,
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    criterion = "bic",
    age_group = "65-79",
    site = "Denmark",
    season = "2023-2024",
    outcome = "death"
  )

  expect_s3_class(p, "ggplot")
})

test_that("get_best_model works with AIC criterion", {
  comparison_result <- compare_decay_models(
    model_names = c("exponential", "logistic"),
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  best_model <- get_best_model(comparison_result, criterion = "aic")

  expect_type(best_model, "list")
  expect_named(best_model, c("model_name", "model", "criterion", "all_results"))
  expect_s3_class(best_model$model, "DecayModel")
  expect_type(best_model$criterion, "double")
  expect_s3_class(best_model$all_results, "data.frame")
})

test_that("get_best_model works with BIC criterion", {
  comparison_result <- compare_decay_models(
    model_names = c("exponential", "logistic"),
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  best_model <- get_best_model(comparison_result, criterion = "bic")

  expect_type(best_model, "list")
  expect_named(best_model, c("model_name", "model", "criterion", "all_results"))
  expect_s3_class(best_model$model, "DecayModel")
  expect_type(best_model$criterion, "double")
})

test_that("get_best_model fails with invalid criterion", {
  comparison_result <- compare_decay_models(
    model_names = c("exponential", "logistic"),
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  expect_error(
    get_best_model(comparison_result, criterion = "invalid"),
    "criterion must be either 'aic' or 'bic'"
  )
})

test_that("get_best_model fails with invalid input", {
  expect_error(
    get_best_model("not_a_comparison_result"),
    "comparison_result must be the output of compare_decay_models\\(\\)"
  )
})

test_that("get_best_model_name works", {
  comparison_result <- compare_decay_models(
    model_names = c("exponential", "logistic"),
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  best_model_name <- get_best_model_name(comparison_result, criterion = "aic")

  expect_type(best_model_name, "character")
  expect_true(best_model_name %in% c("exponential", "logistic"))
})

test_that("get_best_model_instance works", {
  comparison_result <- compare_decay_models(
    model_names = c("exponential", "logistic"),
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  best_model_instance <- get_best_model_instance(
    comparison_result,
    criterion = "aic"
  )

  expect_s3_class(best_model_instance, "DecayModel")
  expect_true(best_model_instance$name %in% c("exponential", "logistic"))
})

test_that("get_best_model_name and get_best_model_instance fail with invalid input", {
  expect_error(
    get_best_model_name("not_a_comparison_result"),
    "comparison_result must be the output of compare_decay_models\\(\\)"
  )

  expect_error(
    get_best_model_instance("not_a_comparison_result"),
    "comparison_result must be the output of compare_decay_models\\(\\)"
  )
})

# Test parameter length validation
test_that("fit_3step validates parameter lengths correctly", {
  model <- create_test_model("exponential")

  # Test with mismatched lengths
  expect_error(
    model$fit_3step(
      estimate = test_estimate[1:5],
      ci_low = test_ci_low[1:6],
      ci_up = test_ci_up[1:5],
      periods_starts = test_periods_starts[1:5],
      periods_ends = test_periods_ends[1:5],
      nboot = 10,
      seed = 123
    ),
    "All parameters must have the same length"
  )

  # Test with all correct lengths
  expect_no_error(
    model$fit_3step(
      estimate = test_estimate,
      ci_low = test_ci_low,
      ci_up = test_ci_up,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends,
      nboot = 10,
      seed = 123
    )
  )
})

test_that("predict_optimal validates parameter lengths correctly", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  # Test with mismatched lengths
  expect_error(
    model$predict_optimal(
      estimate = test_estimate[1:5],
      periods_starts = test_periods_starts[1:6],
      periods_ends = test_periods_ends[1:5]
    )
  )

  # Test with all correct lengths
  expect_no_error(
    model$predict_optimal(
      estimate = test_estimate,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends
    )
  )
})

test_that("predict validates parameter lengths correctly", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  params <- model$get_optimal_params()

  # Test with mismatched lengths
  expect_error(
    model$predict(
      params = params,
      estimate = test_estimate[1:5],
      periods_starts = test_periods_starts[1:6],
      periods_ends = test_periods_ends[1:5]
    )
  )

  # Test with all correct lengths
  expect_no_error(
    model$predict(
      params = params,
      estimate = test_estimate,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends
    )
  )
})

test_that("plot_optimal validates parameter lengths correctly", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  # Test with mismatched lengths
  expect_error(
    model$plot_optimal(
      estimate = test_estimate[1:5],
      ci_low = test_ci_low[1:6],
      ci_up = test_ci_up[1:5],
      periods_starts = test_periods_starts[1:5],
      periods_ends = test_periods_ends[1:5],
      age_group = "65-79",
      site = "Denmark",
      season = "2023-2024",
      outcome = "death"
    ),
    "All parameters must have the same length"
  )

  # Test with all correct lengths
  expect_no_error(
    model$plot_optimal(
      estimate = test_estimate,
      ci_low = test_ci_low,
      ci_up = test_ci_up,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends,
      age_group = "65-79",
      site = "Denmark",
      season = "2023-2024",
      outcome = "death"
    )
  )
})

test_that("plot validates parameter lengths correctly", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  params <- model$get_optimal_params()

  # Test with mismatched lengths
  expect_error(
    model$plot(
      params = params,
      estimate = test_estimate[1:5],
      ci_low = test_ci_low[1:6],
      ci_up = test_ci_up[1:5],
      periods_starts = test_periods_starts[1:5],
      periods_ends = test_periods_ends[1:5],
      age_group = "65-79",
      site = "Denmark",
      season = "2023-2024",
      outcome = "death"
    ),
    "All parameters must have the same length"
  )

  # Test with all correct lengths
  expect_no_error(
    model$plot(
      params = params,
      estimate = test_estimate,
      ci_low = test_ci_low,
      ci_up = test_ci_up,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends,
      age_group = "65-79",
      site = "Denmark",
      season = "2023-2024",
      outcome = "death"
    )
  )
})

test_that("build_prediction_dataset validates parameter lengths correctly", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  # Test with mismatched lengths
  expect_error(
    model$build_prediction_dataset(
      estimate = test_estimate[1:5],
      ci_low = test_ci_low[1:6],
      ci_up = test_ci_up[1:5],
      periods_starts = test_periods_starts[1:5],
      periods_ends = test_periods_ends[1:5]
    ),
    "All parameters must have the same length"
  )

  # Test with all correct lengths
  expect_no_error(
    model$build_prediction_dataset(
      estimate = test_estimate,
      ci_low = test_ci_low,
      ci_up = test_ci_up,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends
    )
  )
})

test_that("generate_VE_based_on_boot_fit validates parameter lengths correctly", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  # Test with mismatched lengths
  expect_error(
    model$generate_VE_based_on_boot_fit(
      periods_starts = test_periods_starts[1:5],
      periods_ends = test_periods_ends[1:6]
    )
  )

  # Test with all correct lengths
  expect_no_error(
    model$generate_VE_based_on_boot_fit(
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends
    )
  )
})

test_that("generate_VE_summary validates parameter lengths correctly", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  # Test with mismatched lengths
  expect_error(
    model$generate_VE_summary(
      periods_starts = test_periods_starts[1:5],
      periods_ends = test_periods_ends[1:6]
    )
  )

  # Test with all correct lengths
  expect_no_error(
    model$generate_VE_summary(
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends
    )
  )
})

test_that("plot_VE_with_uncertainty validates parameter lengths correctly", {
  model <- create_test_model("exponential")

  # Fit the model first
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  # Test with mismatched lengths
  expect_error(
    model$plot_VE_with_uncertainty(
      periods_starts = test_periods_starts[1:5],
      periods_ends = test_periods_ends[1:6]
    )
  )

  # Test with all correct lengths
  expect_no_error(
    model$plot_VE_with_uncertainty(
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends
    )
  )
})


test_that("compare_decay_models validates parameter lengths correctly", {
  # Test with mismatched lengths
  expect_error(
    compare_decay_models(
      model_names = c("exponential"),
      estimate = test_estimate[1:5],
      ci_low = test_ci_low[1:6],
      ci_up = test_ci_up[1:5],
      periods_starts = test_periods_starts[1:5],
      periods_ends = test_periods_ends[1:5],
      nboot = 10,
      seed = 123
    ),
    "All parameters must have the same length"
  )

  # Test with all correct lengths
  expect_no_error(
    compare_decay_models(
      model_names = c("exponential"),
      estimate = test_estimate,
      ci_low = test_ci_low,
      ci_up = test_ci_up,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends,
      nboot = 10,
      seed = 123
    )
  )
})

test_that("plot_model_comparison validates parameter lengths correctly", {
  # Create comparison result first
  comparison_result <- compare_decay_models(
    model_names = c("exponential"),
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  # Test with mismatched lengths
  expect_error(
    plot_model_comparison(
      comparison_result = comparison_result,
      estimate = test_estimate[1:5],
      ci_low = test_ci_low[1:6],
      ci_up = test_ci_up[1:5],
      periods_starts = test_periods_starts[1:5],
      periods_ends = test_periods_ends[1:5],
      age_group = "65-79",
      site = "Denmark",
      season = "2023-2024",
      outcome = "death"
    ),
    "All parameters must have the same length"
  )

  # Test with all correct lengths
  expect_no_error(
    plot_model_comparison(
      comparison_result = comparison_result,
      estimate = test_estimate,
      ci_low = test_ci_low,
      ci_up = test_ci_up,
      periods_starts = test_periods_starts,
      periods_ends = test_periods_ends,
      age_group = "65-79",
      site = "Denmark",
      season = "2023-2024",
      outcome = "death"
    )
  )
})

test_that("parameter length validation provides detailed error messages", {
  model <- create_test_model("exponential")

  # Test detailed error message for fit_3step
  expect_error(
    model$fit_3step(
      estimate = test_estimate[1:5],
      ci_low = test_ci_low[1:6],
      ci_up = test_ci_up[1:7],
      periods_starts = test_periods_starts[1:8],
      periods_ends = test_periods_ends[1:9],
      nboot = 10,
      seed = 123
    ),
    "Lengths: estimate=5, ci_low=6, ci_up=7, periods_starts=8, periods_ends=9"
  )

  # Test detailed error message for predict_optimal
  model$fit_3step(
    estimate = test_estimate,
    ci_low = test_ci_low,
    ci_up = test_ci_up,
    periods_starts = test_periods_starts,
    periods_ends = test_periods_ends,
    nboot = 10,
    seed = 123
  )

  expect_error(
    model$predict_optimal(
      estimate = test_estimate[1:5],
      periods_starts = test_periods_starts[1:6],
      periods_ends = test_periods_ends[1:7]
    ),
    "Lengths: estimate=5, periods_starts=6, periods_ends=7"
  )
})

# Test the generic validate_parameter_lengths function
test_that("validate_parameter_lengths works correctly", {
  # Test with all parameters having the same length
  expect_no_error(
    validate_parameter_lengths(
      estimate = c(1, 2, 3),
      ci_low = c(0.5, 1.5, 2.5),
      ci_up = c(1.5, 2.5, 3.5),
      periods_starts = c(1, 2, 3),
      periods_ends = c(2, 3, 4)
    )
  )

  # Test with mismatched lengths
  expect_error(
    validate_parameter_lengths(
      estimate = c(1, 2),
      ci_low = c(0.5, 1.5, 2.5),
      ci_up = c(1.5, 2.5),
      periods_starts = c(1, 2),
      periods_ends = c(2, 3, 4)
    ),
    "All parameters must have the same length"
  )

  # Test with function name
  expect_error(
    validate_parameter_lengths(
      estimate = c(1, 2),
      ci_low = c(0.5, 1.5, 2.5),
      function_name = "test_function"
    ),
    "In test_function: All parameters must have the same length"
  )
})

test_that("validate_parameter_lengths with required_lengths works correctly", {
  # Test with specific parameter groups
  expect_no_error(
    validate_parameter_lengths(
      estimate = c(1, 2, 3),
      ci_low = c(0.5, 1.5, 2.5),
      periods_starts = c(1, 2, 3),
      periods_ends = c(2, 3, 4),
      required_lengths = list(
        "main_params" = c("estimate", "ci_low"),
        "period_params" = c("periods_starts", "periods_ends")
      )
    )
  )

  # Test with mismatched lengths in a group
  expect_error(
    validate_parameter_lengths(
      estimate = c(1, 2),
      ci_low = c(0.5, 1.5, 2.5),
      periods_starts = c(1, 2, 3),
      periods_ends = c(2, 3, 4),
      required_lengths = list(
        "main_params" = c("estimate", "ci_low"),
        "period_params" = c("periods_starts", "periods_ends")
      )
    ),
    "main_params must have the same length"
  )

  # Test with missing parameters
  expect_error(
    validate_parameter_lengths(
      estimate = c(1, 2, 3),
      ci_low = c(0.5, 1.5, 2.5),
      required_lengths = list(
        "main_params" = c("estimate", "ci_low", "missing_param")
      )
    ),
    "Missing required parameters: missing_param"
  )
})

test_that("validate_parameter_lengths provides detailed error messages", {
  # Test detailed error message
  expect_error(
    validate_parameter_lengths(
      estimate = c(1, 2),
      ci_low = c(0.5, 1.5, 2.5),
      ci_up = c(1.5, 2.5),
      periods_starts = c(1, 2),
      periods_ends = c(2, 3, 4)
    ),
    "Lengths: estimate=2, ci_low=3, ci_up=2, periods_starts=2, periods_ends=3"
  )

  # Test detailed error message with function name
  expect_error(
    validate_parameter_lengths(
      estimate = c(1, 2),
      ci_low = c(0.5, 1.5, 2.5),
      function_name = "my_function"
    ),
    "In my_function: All parameters must have the same length"
  )
})
