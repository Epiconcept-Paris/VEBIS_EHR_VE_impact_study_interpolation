# DecayModel class

A class to represent a decay model

## Public fields

- `name`:

  Model name

- `decay_function`:

  The decay function

- `param_config`:

  Parameter configuration

- `optimization_results`:

  Results from the last optimization

- `nlrob_samples`:

  List of nlrob samples

- `ve_summary`:

  Summary of VE samples

## Methods

### Public methods

- [`DecayModel$new()`](#method-DecayModel-new)

- [`DecayModel$evaluate()`](#method-DecayModel-evaluate)

- [`DecayModel$create_nls_formula()`](#method-DecayModel-create_nls_formula)

- [`DecayModel$fit_to_determine_bic()`](#method-DecayModel-fit_to_determine_bic)

- [`DecayModel$generate_VE_based_on_resampling()`](#method-DecayModel-generate_VE_based_on_resampling)

- [`DecayModel$generate_VE_summary()`](#method-DecayModel-generate_VE_summary)

- [`DecayModel$plot_VE_with_uncertainty()`](#method-DecayModel-plot_VE_with_uncertainty)

- [`DecayModel$plot_VE_with_uncertainty_plus_raw_data()`](#method-DecayModel-plot_VE_with_uncertainty_plus_raw_data)

- [`DecayModel$plot_bootstrap_samples()`](#method-DecayModel-plot_bootstrap_samples)

- [`DecayModel$get_optimal_params()`](#method-DecayModel-get_optimal_params)

- [`DecayModel$get_aic()`](#method-DecayModel-get_aic)

- [`DecayModel$get_bic()`](#method-DecayModel-get_bic)

- [`DecayModel$predict_optimal()`](#method-DecayModel-predict_optimal)

- [`DecayModel$predict()`](#method-DecayModel-predict)

- [`DecayModel$plot_optimal()`](#method-DecayModel-plot_optimal)

- [`DecayModel$plot()`](#method-DecayModel-plot)

- [`DecayModel$print()`](#method-DecayModel-print)

- [`DecayModel$build_prediction_dataset()`](#method-DecayModel-build_prediction_dataset)

- [`DecayModel$clone()`](#method-DecayModel-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize a specific decay model

#### Usage

    DecayModel$new(name, decay_function, param_config)

#### Arguments

- `name`:

  Model name

- `decay_function`:

  Function that implements the decay model

- `param_config`:

  Parameter configuration list

------------------------------------------------------------------------

### Method `evaluate()`

Evaluate the model with given parameters

#### Usage

    DecayModel$evaluate(time_steps, params)

#### Arguments

- `time_steps`:

  Time steps to evaluate

- `params`:

  Named vector of parameters

#### Returns

Vector of decay values

------------------------------------------------------------------------

### Method `create_nls_formula()`

Create formula for nls fitting

#### Usage

    DecayModel$create_nls_formula()

#### Returns

Formula suitable for nls

------------------------------------------------------------------------

### Method `fit_to_determine_bic()`

Fit the model using 3-step strategy (nlxb -\> nls -\> bootstrap)

#### Usage

    DecayModel$fit_to_determine_bic(
      estimate,
      ci_low,
      ci_up,
      periods_starts,
      periods_ends,
      seed = 123,
      maxiteritations = 50,
      type_of_fit = "midpoint"
    )

#### Arguments

- `estimate`:

  Log estimates (in the HR scale)

- `ci_low`:

  Lower confidence interval bounds (in the HR scale)

- `ci_up`:

  Upper confidence interval bounds (in the HR scale)

- `periods_starts`:

  Vector of periods starts

- `periods_ends`:

  Vector of periods ends

- `seed`:

  Seed for reproducibility (default is 123)

- `maxiteritations`:

  Maximum number of iterations for nls (default is 500)

- `type_of_fit`:

  Type of fit (default is "midpoint")

#### Returns

List with fitting results

------------------------------------------------------------------------

### Method `generate_VE_based_on_resampling()`

Generate VE based on boot fit

#### Usage

    DecayModel$generate_VE_based_on_resampling(
      estimate,
      ci_low,
      ci_up,
      periods_starts,
      periods_ends,
      n_weekly_ve_to_generate = 1000,
      type_of_fit = "midpoint",
      seed = 123
    )

#### Arguments

- `estimate`:

  Log estimates (in the HR scale)

- `ci_low`:

  Lower confidence interval bounds (in the HR scale)

- `ci_up`:

  Upper confidence interval bounds (in the HR scale)

- `periods_starts`:

  Vector of periods starts

- `periods_ends`:

  Vector of periods ends

- `n_weekly_ve_to_generate`:

  Number of weekly VE to generate (default is 1000)

- `type_of_fit`:

  Type of fit (default is "midpoint")

- `seed`:

  Seed for reproducibility (default is 123)

#### Returns

Matrix of VE samples

------------------------------------------------------------------------

### Method `generate_VE_summary()`

Generate VE summary statistics from bootstrap samples

#### Usage

    DecayModel$generate_VE_summary(
      periods_starts,
      periods_ends,
      probs = c(0.025, 0.5, 0.975)
    )

#### Arguments

- `periods_starts`:

  Vector of periods starts

- `periods_ends`:

  Vector of periods ends

- `probs`:

  Vector of quantiles to compute (default is c(0.025, 0.5, 0.975))

#### Returns

Data frame with summary statistics

------------------------------------------------------------------------

### Method `plot_VE_with_uncertainty()`

Plot VE with uncertainty bands from bootstrap

#### Usage

    DecayModel$plot_VE_with_uncertainty(
      periods_starts,
      periods_ends,
      title = NULL,
      color = "blue",
      alpha = 0.3,
      show_points = TRUE,
      probs = c(0.025, 0.975)
    )

#### Arguments

- `periods_starts`:

  Vector of periods starts

- `periods_ends`:

  Vector of periods ends

- `title`:

  Plot title (default is NULL)

- `color`:

  Color for the line and ribbon (default is "blue")

- `alpha`:

  Transparency for the ribbon (default is 0.3)

- `show_points`:

  Whether to show points on the line (default is TRUE)

- `probs`:

  Vector of quantiles for confidence intervals (default is c(0.025,
  0.975))

#### Returns

ggplot object

------------------------------------------------------------------------

### Method `plot_VE_with_uncertainty_plus_raw_data()`

Plot VE with uncertainty bands from bootstrap

#### Usage

    DecayModel$plot_VE_with_uncertainty_plus_raw_data(
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
    )

#### Arguments

- `estimate`:

  Log estimates (in the HR scale)

- `ci_low`:

  Lower confidence interval bounds (in the HR scale)

- `ci_up`:

  Upper confidence interval bounds (in the HR scale)

- `periods_starts`:

  Vector of periods starts

- `periods_ends`:

  Vector of periods ends

- `title`:

  Plot title (default is NULL)

- `color`:

  Color for the line and ribbon (default is "blue")

- `alpha`:

  Transparency for the ribbon (default is 0.3)

- `show_points`:

  Whether to show points on the line (default is TRUE)

- `probs`:

  Vector of quantiles for confidence intervals (default is c(0.025,
  0.975))

#### Returns

ggplot object

------------------------------------------------------------------------

### Method `plot_bootstrap_samples()`

Plot individual bootstrap samples

#### Usage

    DecayModel$plot_bootstrap_samples(
      periods_starts,
      periods_ends,
      n_samples = 5,
      title = NULL,
      alpha = 0.7
    )

#### Arguments

- `periods_starts`:

  Vector of periods starts

- `periods_ends`:

  Vector of periods ends

- `n_samples`:

  Number of bootstrap samples to show (default is 5)

- `title`:

  Plot title (default is NULL)

- `alpha`:

  Transparency for the lines (default is 0.7)

#### Returns

ggplot object

------------------------------------------------------------------------

### Method `get_optimal_params()`

Get optimal parameters from last optimization

#### Usage

    DecayModel$get_optimal_params()

#### Returns

Named vector of optimal parameters

------------------------------------------------------------------------

### Method `get_aic()`

Get AIC value from last optimization

#### Usage

    DecayModel$get_aic()

#### Returns

AIC value

------------------------------------------------------------------------

### Method `get_bic()`

Get BIC value from last optimization

#### Usage

    DecayModel$get_bic()

#### Returns

BIC value

------------------------------------------------------------------------

### Method `predict_optimal()`

Predict using optimal parameters

#### Usage

    DecayModel$predict_optimal(
      estimate,
      periods_starts,
      periods_ends,
      aggregate_by_8_weeks = FALSE,
      criterion = "aic"
    )

#### Arguments

- `estimate`:

  estimates (in the HR scale)

- `periods_starts`:

  Vector of periods starts

- `periods_ends`:

  Vector of periods ends

- `aggregate_by_8_weeks`:

  Whether to aggregate the predictions by 8 weeks (default is FALSE)

- `criterion`:

  Criterion to use for prediction (default is "aic")

#### Returns

Vector of predictions

------------------------------------------------------------------------

### Method [`predict()`](https://rdrr.io/r/stats/predict.html)

Predict decay model for given parameters

#### Usage

    DecayModel$predict(
      params,
      estimate,
      periods_starts,
      periods_ends,
      aggregate_by_8_weeks = FALSE,
      criterion = "aic"
    )

#### Arguments

- `params`:

  Named vector of parameters

- `estimate`:

  estimates (in the HR scale)

- `periods_starts`:

  Vector of periods starts

- `periods_ends`:

  Vector of periods ends

- `aggregate_by_8_weeks`:

  Whether to aggregate the predictions by 8 weeks (default is FALSE)

- `criterion`:

  Criterion to use for prediction (default is "aic")

#### Returns

Vector of predictions

------------------------------------------------------------------------

### Method `plot_optimal()`

Plot using optimal parameters

#### Usage

    DecayModel$plot_optimal(
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
    )

#### Arguments

- `estimate`:

  Vector of estimates (in the HR scale)

- `ci_low`:

  Lower confidence interval (in the HR scale)

- `ci_up`:

  Upper confidence interval (in the HR scale)

- `periods_starts`:

  Vector of periods starts

- `periods_ends`:

  Vector of periods ends

- `title`:

  Plot title

- `color`:

  Color of the line

- `age_group`:

  Age group

- `site`:

  Site

- `season`:

  Season

- `outcome`:

  Outcome

#### Returns

ggplot object

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Plot decay model

#### Usage

    DecayModel$plot(
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
    )

#### Arguments

- `params`:

  Named vector of parameters

- `estimate`:

  Vector of estimates (in the HR scale)

- `ci_low`:

  Lower confidence interval (in the HR scale)

- `ci_up`:

  Upper confidence interval (in the HR scale)

- `periods_starts`:

  Vector of periods starts

- `periods_ends`:

  Vector of periods ends

- `title`:

  Plot title

- `color`:

  Color of the line (default is "red")

- `age_group`:

  Age group

- `site`:

  Site

- `season`:

  Season

- `outcome`:

  Outcome

- `criterion`:

  Criterion to use for prediction (default is "aic")

#### Returns

ggplot object

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print model information

#### Usage

    DecayModel$print()

------------------------------------------------------------------------

### Method `build_prediction_dataset()`

Build final dataset with weekly and 8-weeks predictions

#### Usage

    DecayModel$build_prediction_dataset(
      estimate,
      ci_low,
      ci_up,
      periods_starts,
      periods_ends,
      include_weekly = TRUE,
      include_8_weeks = TRUE,
      criterion = "aic"
    )

#### Arguments

- `estimate`:

  Vector of estimates (in the HR scale)

- `ci_low`:

  Vector of lower confidence interval bounds (in the HR scale)

- `ci_up`:

  Vector of upper confidence interval bounds (in the HR scale)

- `periods_starts`:

  Vector of periods starts

- `periods_ends`:

  Vector of periods ends

- `include_weekly`:

  Whether to include weekly predictions (default is TRUE)

- `include_8_weeks`:

  Whether to include 8-weeks predictions (default is TRUE)

- `criterion`:

  Criterion to use for prediction (default is "aic")

#### Returns

List containing weekly_data and data_8_weeks data frames

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    DecayModel$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
