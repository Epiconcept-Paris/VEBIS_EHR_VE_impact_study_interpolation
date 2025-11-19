# Compare multiple decay models

Compare multiple decay models

Plot model comparison

Get the best model based on lowest AIC

Get the best model name only

Get the best model instance only

## Usage

``` r
compare_decay_models(
  model_names = c("exponential", "logistic", "logistic_simple"),
  estimate,
  ci_low,
  ci_up,
  periods_starts,
  periods_ends,
  decay_funs = list(exponential = exponential_decay(), logistic = logistic_decay(),
    logistic_simple = logistic_simple_decay()),
  param_configs = list(exponential = list(start = c(VE0 = 0.7, decay_rate = 0.1), lower =
    c(0.01, 0.001), upper = c(0.99, 1), names = c("VE0", "decay_rate")), logistic =
    list(start = c(VE0 = 0.7, decay_rate = 0.1, constant = 0.5), lower = c(0.01, 0.001,
    0.01), upper = c(0.99, 1, 2), names = c("VE0", "decay_rate", "constant")),
    logistic_simple = list(start = c(VE0 = 0.7, decay_rate = 0.1), lower = c(0.01,
    0.001), upper = c(0.99, 1), names = c("VE0", "decay_rate"))),
  seed = 123,
  type_of_fit = "midpoint"
)

plot_model_comparison(
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
)

get_best_model(comparison_result, criterion = "aic")

get_best_model_name(comparison_result, criterion = "aic")

get_best_model_instance(comparison_result, criterion = "aic")
```

## Arguments

- model_names:

  Vector of model names to use (default is c("exponential", "logistic"))

- estimate:

  Vector of estimates (in the HR scale)

- ci_low:

  Vector of lower confidence interval bounds (in the HR scale)

- ci_up:

  Vector of upper confidence interval bounds (in the HR scale)

- periods_starts:

  Vector of periods starts

- periods_ends:

  Vector of periods ends

- decay_funs:

  List of decay functions to use (default is list(exponential =
  exponential_decay(), logistic = logistic_decay(), logistic_simple =
  logistic_simple_decay()))

- param_configs:

  List of parameter configurations to use (default is list(exponential =
  list(start = c(VE0 = 0.7, decay_rate = 0.1), lower = c(0.01, 0.001),
  upper = c(0.99, 1.0), names = c("VE0", "decay_rate")), logistic =
  list(start = c(VE0 = 0.7, decay_rate = 0.1, constant = 0.5), lower =
  c(0.01, 0.001, 0.01), upper = c(0.99, 1.0, 2.0), names = c("VE0",
  "decay_rate", "constant")), logistic_simple = list(start = c(VE0 =
  0.7, decay_rate = 0.1), lower = c(0.01, 0.001), upper = c(0.99, 1.0),
  names = c("VE0", "decay_rate"))))

- seed:

  Seed to use for reproducibility (default is 123)

- type_of_fit:

  Type of fit (default is "midpoint")

- comparison_result:

  List containing models, results, comparison, and summary

- colors:

  Vector of colors to use for the models (default is c("red", "blue",
  "green", "purple", "orange"))

- criterion:

  Criterion to use for best model selection (default is "aic")

- age_group:

  Age group

- site:

  Site

- season:

  Season

- outcome:

  Outcome

## Value

List containing models, results, comparison, and summary

ggplot object

List containing best model name, instance, AIC, and all results

Best model name

Best model instance
