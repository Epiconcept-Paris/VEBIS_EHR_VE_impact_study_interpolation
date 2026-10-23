# VEinterpolation (development version)

* `generate_VE_based_on_resampling()` now also stores the noise-free VE curves in `nlrob_samples_clean`, and `generate_VE_summary()` adds `VE_clean_median`, `VE_clean_quantile_025` and `VE_clean_quantile_975`.
* `plot_VE_with_uncertainty()` and `plot_VE_with_uncertainty_plus_raw_data()` gain a `show_clean` argument (default `TRUE`): the line is the interpolated VE without measurement error, the dark ribbon its bootstrap interval, and the light ribbon the prediction interval including measurement error. `show_clean = FALSE` gives the previous display.
* `plot_VE_with_uncertainty_plus_raw_data()` no longer fails when `show_points = TRUE`.
* Added unit tests.

# VEinterpolation 0.1.0

Methodology modified according to Belgium comments on Nov 19th

# VEinterpolation 0.0.0.9000

Shared with Belgium for evaluation purposes on Nov 5th

* Initial CRAN submission.
