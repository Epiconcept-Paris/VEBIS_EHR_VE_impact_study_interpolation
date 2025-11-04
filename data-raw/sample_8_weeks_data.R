## code to prepare `sample_8_weeks_data` dataset goes here
set.seed(123)
periods_starts <- structure(
  c(
    19631,
    19662,
    19692,
    19723,
    19754,
    19783,
    19814,
    19844,
    19875,
    19905,
    19936,
    19967
  ),
  class = "Date"
)

periods_ends <- structure(
  c(
    19686,
    19717,
    19747,
    19778,
    19809,
    19838,
    19869,
    19899,
    19930,
    19960,
    19991,
    20022
  ),
  class = "Date"
)

estimates <- c(
  0.355703480676933,
  0.388598976376876,
  0.457007169565882,
  0.517875272577604,
  0.651627831105234,
  0.722080187822221,
  0.749280287133703,
  0.864132324700696,
  0.874473193962537,
  0.893509633032125,
  0.965858168415434,
  1.12608899036491
) +
  rnorm(length(periods_ends), 0, 0.05)

log_estimates <- log(estimates)
se <- runif(length(estimates), 0.05, 0.30)

CIlow <- exp(log_estimates - 1.96 * se)
CIup <- exp(log_estimates + 1.96 * se)

sample_8_weeks_data <- data.frame(
  date_min = periods_starts,
  date_max = periods_ends,
  estimate = estimates,
  CIlow = CIlow,
  CIhigh = CIup
)

usethis::use_data(sample_8_weeks_data, overwrite = TRUE)
checkhelper::use_data_doc("sample_8_weeks_data")
devtools::document()
