library(dplyr)
library(hubData)
source("scripts/cdf_to_quantiles.R")
quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)

bucket_name <- "uscdc-flusight-hub-v1"
hub_bucket <- s3_bucket(bucket_name)
hub_con <- hubData::connect_hub(hub_bucket, file_format = "parquet", skip_checks = TRUE)
model_names <- c("hist-avg", "delphi-epicast")
step_ahead_model_output <- hub_con |>
  dplyr::filter(
    model_id %in% model_names,
    target == "ili perc",
    output_type == "cdf"
  ) |>
  hubData::collect_hub() |>
  dplyr::mutate(
    output_type_id = as.numeric(output_type_id)
  )

# save all the component models in hubverse format
model_dates_df <- distinct(step_ahead_model_output, model_id, origin_date)
model_dates_df |>
  purrr::pwalk(.f = function(model_id, origin_date, ...) {
    model_folder <- file.path("model-output", model_id)
    if (!file.exists(model_folder)) dir.create(model_folder, recursive = TRUE)

    results_path <- file.path(model_folder, paste0(origin_date, "-", model_id, ".csv"))

    model <- model_id
    date <- origin_date
    model_outputs <- step_ahead_model_output |>
      dplyr::filter(model_id == model, origin_date == date, output_type == "cdf") |>
      dplyr::select(-"model_id")
    quantile_outputs <- model_outputs |>
      group_by(horizon, location) |>
      group_split() |>
      purrr::map(.f = function(split_outputs) {
        expand.grid(
          stringsAsFactors = FALSE,
          origin_date = date,
          location = split_outputs$location[1],
          horizon = split_outputs$horizon[1],
          target = split_outputs$target[1],
          output_type = "quantile",
          output_type_id = quantile_levels
        ) |>
          mutate(value = cdf_to_quantiles(cdf_df = split_outputs,
                                          quantile_probs = quantile_levels,
                                          x_col = "output_type_id",
                                          cdf_col = "value"))
      }) |>
      purrr::list_rbind() |>
      mutate(target_end_date = origin_date + 7 * horizon) |>
      select(origin_date, location, target,
             horizon, target_end_date,
             output_type, output_type_id, value)
    ## print(quantile_outputs)

    utils::write.csv(quantile_outputs, file = results_path, row.names = FALSE)
  })
