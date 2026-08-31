# Independently reconstruct and validate the nested-CV computation records.

source(file.path("R", "00_common.R"))
assert_packages(c("data.table", "digest", "R.utils"))
source(file.path("R", "09_nested_cv_helpers.R"))
message_rule("10: Independent nested-validation reconstruction QA")

read_gz <- function(directory, filename) {
  data.table::fread(file.path(directory, filename), showProgress = FALSE)
}

common_predictors <- read_gz(
  paths$analysis, "predictive_common_case_predictors_long.csv.gz"
)
outer_assignments <- read_gz(
  paths$analysis, "nested_cv_outer_fold_assignments.csv.gz"
)
inner_assignments <- read_gz(
  paths$analysis, "nested_cv_inner_fold_assignments.csv.gz"
)
predictions <- read_gz(
  paths$analysis, "nested_cv_fold_level_predictions.csv.gz"
)
averaged_predictions <- read_gz(
  paths$analysis, "nested_cv_averaged_predictions.csv.gz"
)
inner_scalers <- read_gz(
  paths$diagnostics, "nested_cv_inner_training_scalers.csv.gz"
)
outer_scalers <- read_gz(
  paths$diagnostics, "nested_cv_outer_training_scalers.csv.gz"
)
aggregate_losses <- read_gz(
  paths$diagnostics, "nested_cv_inner_aggregate_tuning_losses.csv.gz"
)
alpha_choices <- read_gz(
  paths$diagnostics, "nested_cv_alpha_one_se_choices.csv.gz"
)
tuning_choices <- data.table::fread(
  file.path(paths$tables, "Table_S_elastic_net_tuning_choices.csv")
)

inner_group_columns <- c(
  "dataset", "outcome", "model_type", "cv_repeat", "outer_fold",
  "inner_validation_fold"
)
outer_group_columns <- c(
  "dataset", "outcome", "model_type", "cv_repeat", "outer_fold"
)

inner_reconstruction_rows <- list()
outer_reconstruction_rows <- list()
leakage_rows <- list()

for (dataset_name in names(datasets)) {
  for (outcome_name in names(outcome_specs)) {
    predictor_long <- common_predictors[
      dataset == dataset_name & outcome == outcome_name
    ]
    predictor_wide <- data.table::dcast(
      predictor_long, Subject_ID ~ predictor, value.var = "value"
    )
    subjects <- predictor_wide$Subject_ID
    predictor_wide[, Subject_ID := NULL]
    predictor_matrix <- as.matrix(predictor_wide)
    storage.mode(predictor_matrix) <- "double"
    rownames(predictor_matrix) <- subjects

    outer_local <- outer_assignments[
      dataset == dataset_name & outcome == outcome_name
    ]
    inner_local <- inner_assignments[
      dataset == dataset_name & outcome == outcome_name
    ]
    inner_scaler_local <- inner_scalers[
      dataset == dataset_name & outcome == outcome_name
    ]
    outer_scaler_local <- outer_scalers[
      dataset == dataset_name & outcome == outcome_name
    ]

    for (repeat_index in seq_len(cv_settings$outer_repeats)) {
      for (outer_fold_value in seq_len(cv_settings$outer_folds)) {
        # Use deliberately distinct loop-variable names so data.table columns
        # can never mask the scalar values being checked.
        outer_repeat <- outer_local[outer_local$cv_repeat == repeat_index]
        outer_test_ids <- outer_repeat[
          outer_repeat$outer_fold == outer_fold_value, Subject_ID
        ]
        outer_training_ids <- outer_repeat[
          outer_repeat$outer_fold != outer_fold_value, Subject_ID
        ]
        inner_split <- inner_local[
          inner_local$cv_repeat == repeat_index &
            inner_local$outer_fold == outer_fold_value
        ]
        leakage_rows[[length(leakage_rows) + 1L]] <- data.table::data.table(
          dataset = dataset_name,
          outcome = outcome_name,
          cv_repeat = repeat_index,
          outer_fold = outer_fold_value,
          outer_test_n = length(outer_test_ids),
          outer_training_n = length(outer_training_ids),
          inner_assignment_n = nrow(inner_split),
          outer_test_in_inner_count = length(intersect(
            outer_test_ids, inner_split$Subject_ID
          )),
          inner_ids_equal_outer_training = setequal(
            inner_split$Subject_ID, outer_training_ids
          )
        )

        for (model_type_value in model_types) {
          saved_outer <- outer_scaler_local[
            outer_scaler_local$cv_repeat == repeat_index &
              outer_scaler_local$outer_fold == outer_fold_value &
              outer_scaler_local$model_type == model_type_value
          ]
          predictors <- saved_outer$predictor
          x <- predictor_matrix[outer_training_ids, predictors, drop = FALSE]
          reconstructed_center <- colMeans(x)
          reconstructed_sd <- apply(x, 2L, stats::sd)
          reconstructed_scale <- reconstructed_sd
          reconstructed_scale[
            !is.finite(reconstructed_scale) | reconstructed_scale == 0
          ] <- 1
          outer_reconstruction_rows[[length(outer_reconstruction_rows) + 1L]] <-
            data.table::data.table(
              dataset = dataset_name,
              outcome = outcome_name,
              model_type = model_type_value,
              cv_repeat = repeat_index,
              outer_fold = outer_fold_value,
              training_n_saved = unique(saved_outer$outer_training_n),
              training_n_reconstructed = nrow(x),
              predictor_n = length(predictors),
              maximum_absolute_center_error = max(abs(
                reconstructed_center[predictors] - saved_outer$center
              )),
              maximum_absolute_raw_sd_error = max(abs(
                reconstructed_sd[predictors] - saved_outer$raw_sd
              )),
              maximum_absolute_scale_used_error = max(abs(
                reconstructed_scale[predictors] - saved_outer$scale_used
              ))
            )

          if (model_type_value %in% c(
            "Elastic Net pathways", "Biomarker + Elastic Net pathways"
          )) {
            for (inner_validation_fold_value in seq_len(cv_settings$inner_folds)) {
              inner_training_ids <- inner_split[
                inner_split$inner_fold != inner_validation_fold_value, Subject_ID
              ]
              inner_validation_ids <- inner_split[
                inner_split$inner_fold == inner_validation_fold_value, Subject_ID
              ]
              saved_inner <- inner_scaler_local[
                inner_scaler_local$cv_repeat == repeat_index &
                  inner_scaler_local$outer_fold == outer_fold_value &
                  inner_scaler_local$model_type == model_type_value &
                  inner_scaler_local$inner_validation_fold == inner_validation_fold_value
              ]
              predictors_inner <- saved_inner$predictor
              x_inner <- predictor_matrix[
                inner_training_ids, predictors_inner, drop = FALSE
              ]
              reconstructed_center <- colMeans(x_inner)
              reconstructed_sd <- apply(x_inner, 2L, stats::sd)
              reconstructed_scale <- reconstructed_sd
              reconstructed_scale[
                !is.finite(reconstructed_scale) | reconstructed_scale == 0
              ] <- 1
              inner_reconstruction_rows[[length(inner_reconstruction_rows) + 1L]] <-
                data.table::data.table(
                  dataset = dataset_name,
                  outcome = outcome_name,
                  model_type = model_type_value,
                  cv_repeat = repeat_index,
                  outer_fold = outer_fold_value,
                  inner_validation_fold = inner_validation_fold_value,
                  training_n_saved = unique(saved_inner$inner_training_n),
                  training_n_reconstructed = nrow(x_inner),
                  validation_n_saved = unique(saved_inner$inner_validation_n),
                  validation_n_reconstructed = length(inner_validation_ids),
                  predictor_n = length(predictors_inner),
                  outer_test_in_inner_training_count = length(intersect(
                    outer_test_ids, inner_training_ids
                  )),
                  outer_test_in_inner_validation_count = length(intersect(
                    outer_test_ids, inner_validation_ids
                  )),
                  inner_validation_in_scaler_training_count = length(intersect(
                    inner_validation_ids, inner_training_ids
                  )),
                  maximum_absolute_center_error = max(abs(
                    reconstructed_center[predictors_inner] - saved_inner$center
                  )),
                  maximum_absolute_raw_sd_error = max(abs(
                    reconstructed_sd[predictors_inner] - saved_inner$raw_sd
                  )),
                  maximum_absolute_scale_used_error = max(abs(
                    reconstructed_scale[predictors_inner] - saved_inner$scale_used
                  ))
                )
            }
          }
        }
      }
    }
  }
}

inner_reconstruction <- data.table::rbindlist(inner_reconstruction_rows)
outer_reconstruction <- data.table::rbindlist(outer_reconstruction_rows)
leakage_audit <- data.table::rbindlist(leakage_rows)

# Reconstruct the one-standard-error lambda choice for every alpha.
recalculated_alpha_choices <- aggregate_losses[, {
  minimum_index <- which.min(mean_fold_mse)
  threshold <- mean_fold_mse[minimum_index] + se_fold_mse[minimum_index]
  eligible <- which(mean_fold_mse <= threshold)
  selected_index <- eligible[which.max(lambda_fraction[eligible])]
  .(
    recalculated_minimum_lambda_index = lambda_index[minimum_index],
    recalculated_minimum_lambda_fraction = lambda_fraction[minimum_index],
    recalculated_minimum_mean_fold_mse = mean_fold_mse[minimum_index],
    recalculated_one_se_threshold = threshold,
    recalculated_selected_lambda_index = lambda_index[selected_index],
    recalculated_selected_lambda_fraction = lambda_fraction[selected_index],
    recalculated_selected_mean_fold_mse = mean_fold_mse[selected_index],
    recalculated_selected_pooled_mse = pooled_mse[selected_index]
  )
}, by = .(
  dataset, outcome, model_type, cv_repeat, outer_fold, alpha
)]

alpha_choice_check <- merge(
  alpha_choices,
  recalculated_alpha_choices,
  by = c("dataset", "outcome", "model_type", "cv_repeat", "outer_fold", "alpha"),
  all = TRUE
)
alpha_choice_check[, `:=`(
  lambda_index_match = selected_lambda_index == recalculated_selected_lambda_index,
  lambda_fraction_error = abs(
    selected_lambda_fraction - recalculated_selected_lambda_fraction
  ),
  pooled_mse_error = abs(
    selected_pooled_mse - recalculated_selected_pooled_mse
  ),
  one_se_threshold_error = abs(
    one_se_threshold - recalculated_one_se_threshold
  )
)]

# Reconstruct deterministic alpha selection after the per-alpha one-SE choice.
selected_alpha_recalculated <- recalculated_alpha_choices[, {
  best_pooled_mse <- min(recalculated_selected_pooled_mse)
  tie_tolerance <- 1e-12 * max(1, abs(best_pooled_mse))
  candidates <- .SD[
    abs(recalculated_selected_pooled_mse - best_pooled_mse) <= tie_tolerance
  ]
  data.table::setorder(
    candidates, -alpha, -recalculated_selected_lambda_fraction
  )
  chosen <- candidates[1L]
  chosen[, `:=`(
    recalculated_alpha_tie_tolerance = tie_tolerance,
    recalculated_alpha_tie_candidate_n = nrow(candidates)
  )]
  chosen
}, by = .(dataset, outcome, model_type, cv_repeat, outer_fold)]
tuning_check <- merge(
  tuning_choices,
  selected_alpha_recalculated[, .(
    dataset, outcome, model_type, cv_repeat, outer_fold,
    recalculated_alpha = alpha,
    recalculated_lambda_fraction = recalculated_selected_lambda_fraction,
    recalculated_pooled_mse = recalculated_selected_pooled_mse,
    recalculated_alpha_tie_tolerance,
    recalculated_alpha_tie_candidate_n
  )],
  by = c("dataset", "outcome", "model_type", "cv_repeat", "outer_fold"),
  all = TRUE
)
tuning_check[, `:=`(
  alpha_match = alpha == recalculated_alpha,
  alpha_tie_tolerance_error = abs(
    alpha_tie_tolerance - recalculated_alpha_tie_tolerance
  ),
  alpha_tie_candidate_n_match =
    alpha_tie_candidate_n == recalculated_alpha_tie_candidate_n,
  lambda_fraction_error = abs(lambda_fraction - recalculated_lambda_fraction),
  pooled_mse_error = abs(selected_inner_pooled_mse - recalculated_pooled_mse)
)]

# Fold identity, held-out counts and five-model common-case equality.
prediction_fold_identity <- predictions[, .(
  model_rows = .N,
  model_types = data.table::uniqueN(model_type),
  outer_fold_count = data.table::uniqueN(outer_fold),
  outer_fold = unique(outer_fold)
), by = .(dataset, outcome, cv_repeat, Subject_ID)]

held_out_count_check <- predictions[, .(
  held_out_rows = .N,
  outer_fold_count = data.table::uniqueN(outer_fold)
), by = .(dataset, outcome, model_type, cv_repeat, Subject_ID)]

averaged_count_check <- averaged_predictions[, .(
  rows = .N,
  held_out_prediction_count = unique(held_out_prediction_count)
), by = .(dataset, outcome, model_type, Subject_ID)]

qa_summary <- data.table::data.table(
  check = c(
    "inner_scaler_groups_reconstructed",
    "inner_scaler_max_center_error",
    "inner_scaler_max_sd_error",
    "inner_validation_in_scaler_training",
    "outer_test_in_inner_records",
    "outer_scaler_groups_reconstructed",
    "outer_scaler_max_center_error",
    "outer_scaler_max_sd_error",
    "outer_test_in_any_inner_assignment",
    "inner_ids_equal_outer_training",
    "per_alpha_one_se_choices",
    "per_alpha_choice_mismatches",
    "deployed_alpha_lambda_choices",
    "deployed_choice_mismatches",
    "participant_repeat_five_model_fold_rows",
    "participant_repeat_fold_identity_failures",
    "held_out_participant_repeat_model_rows",
    "held_out_count_failures",
    "averaged_participant_model_rows",
    "ten_prediction_count_failures",
    "finite_negative_repeat_level_r_squared_retained"
  ),
  value = c(
    nrow(inner_reconstruction),
    max(inner_reconstruction$maximum_absolute_center_error),
    max(inner_reconstruction$maximum_absolute_scale_used_error),
    sum(inner_reconstruction$inner_validation_in_scaler_training_count),
    sum(inner_reconstruction$outer_test_in_inner_training_count) +
      sum(inner_reconstruction$outer_test_in_inner_validation_count),
    nrow(outer_reconstruction),
    max(outer_reconstruction$maximum_absolute_center_error),
    max(outer_reconstruction$maximum_absolute_scale_used_error),
    sum(leakage_audit$outer_test_in_inner_count),
    sum(!leakage_audit$inner_ids_equal_outer_training),
    nrow(alpha_choice_check),
    sum(!alpha_choice_check$lambda_index_match |
          alpha_choice_check$lambda_fraction_error > 1e-12 |
          alpha_choice_check$pooled_mse_error > 1e-12),
    nrow(tuning_check),
    sum(!tuning_check$alpha_match |
          !tuning_check$alpha_tie_candidate_n_match |
          tuning_check$alpha_tie_tolerance_error > 1e-12 |
          tuning_check$lambda_fraction_error > 1e-12 |
          tuning_check$pooled_mse_error > 1e-12),
    nrow(prediction_fold_identity),
    sum(prediction_fold_identity$model_rows != 5L |
          prediction_fold_identity$model_types != 5L |
          prediction_fold_identity$outer_fold_count != 1L),
    nrow(held_out_count_check),
    sum(held_out_count_check$held_out_rows != 1L |
          held_out_count_check$outer_fold_count != 1L),
    nrow(averaged_count_check),
    sum(averaged_count_check$rows != 1L |
          averaged_count_check$held_out_prediction_count != 10L),
    sum(data.table::fread(
      file.path(paths$tables, "Table_S_nested_cv_repeat_metrics.csv")
    )$model_scale_r_squared < 0)
  )
)

tolerance <- 1e-12
failures <- c(
  max(inner_reconstruction$maximum_absolute_center_error) > tolerance,
  max(inner_reconstruction$maximum_absolute_scale_used_error) > tolerance,
  sum(inner_reconstruction$inner_validation_in_scaler_training_count) != 0L,
  sum(inner_reconstruction$outer_test_in_inner_training_count) != 0L,
  sum(inner_reconstruction$outer_test_in_inner_validation_count) != 0L,
  max(outer_reconstruction$maximum_absolute_center_error) > tolerance,
  max(outer_reconstruction$maximum_absolute_scale_used_error) > tolerance,
  sum(leakage_audit$outer_test_in_inner_count) != 0L,
  any(!leakage_audit$inner_ids_equal_outer_training),
  any(!alpha_choice_check$lambda_index_match),
  any(alpha_choice_check$lambda_fraction_error > tolerance),
  any(alpha_choice_check$pooled_mse_error > tolerance),
  any(!tuning_check$alpha_match),
  any(!tuning_check$alpha_tie_candidate_n_match),
  any(tuning_check$alpha_tie_tolerance_error > tolerance),
  any(tuning_check$lambda_fraction_error > tolerance),
  any(tuning_check$pooled_mse_error > tolerance),
  any(prediction_fold_identity$model_rows != 5L),
  any(prediction_fold_identity$model_types != 5L),
  any(prediction_fold_identity$outer_fold_count != 1L),
  any(held_out_count_check$held_out_rows != 1L),
  any(held_out_count_check$outer_fold_count != 1L),
  any(averaged_count_check$rows != 1L),
  any(averaged_count_check$held_out_prediction_count != 10L)
)
failure_names <- c(
  "inner_center_error", "inner_scale_error",
  "inner_validation_in_scaler_training",
  "outer_test_in_inner_training", "outer_test_in_inner_validation",
  "outer_center_error", "outer_scale_error",
  "outer_test_in_inner_assignment", "inner_ids_not_outer_training",
  "one_se_lambda_index", "one_se_lambda_fraction", "one_se_pooled_mse",
  "deployed_alpha", "deployed_alpha_tie_candidate_n",
  "deployed_alpha_tie_tolerance", "deployed_lambda_fraction",
  "deployed_pooled_mse",
  "five_model_row_count", "five_model_type_count", "five_model_fold_identity",
  "held_out_row_count", "held_out_fold_identity",
  "averaged_row_count", "ten_held_out_predictions"
)
failure_audit <- data.table::data.table(
  check = failure_names,
  failed = failures
)
data.table::fwrite(
  qa_summary,
  file.path(paths$validation, "NESTED_CV_QA_SUMMARY_PRECHECK.csv")
)
data.table::fwrite(
  failure_audit,
  file.path(paths$validation, "NESTED_CV_QA_FAILURE_PRECHECK.csv")
)
data.table::fwrite(
  tuning_check,
  file.path(paths$validation, "NESTED_CV_DEPLOYED_TUNING_PRECHECK.csv")
)
if (any(failures)) {
  print(qa_summary)
  print(failure_audit[failed == TRUE])
  stop("One or more independent nested-CV QA checks failed.")
}

data.table::fwrite(
  inner_reconstruction,
  file.path(paths$validation, "NESTED_CV_INNER_SCALER_RECONSTRUCTION.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  outer_reconstruction,
  file.path(paths$validation, "NESTED_CV_OUTER_SCALER_RECONSTRUCTION.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  leakage_audit,
  file.path(paths$validation, "NESTED_CV_LEAKAGE_AUDIT.csv")
)
data.table::fwrite(
  alpha_choice_check,
  file.path(paths$validation, "NESTED_CV_ONE_SE_TUNING_RECONSTRUCTION.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  tuning_check,
  file.path(paths$validation, "NESTED_CV_DEPLOYED_TUNING_RECONSTRUCTION.csv")
)
data.table::fwrite(
  qa_summary,
  file.path(paths$validation, "NESTED_CV_QA_SUMMARY.csv")
)

writeLines(
  c(
    "INDEPENDENT NESTED-CV RECONSTRUCTION QA: PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste("Inner scaler groups reconstructed:", nrow(inner_reconstruction)),
    paste("Maximum inner center error:",
          format(max(inner_reconstruction$maximum_absolute_center_error), scientific = TRUE)),
    paste("Maximum inner scale error:",
          format(max(inner_reconstruction$maximum_absolute_scale_used_error), scientific = TRUE)),
    paste("Outer scaler groups reconstructed:", nrow(outer_reconstruction)),
    paste("Maximum outer center error:",
          format(max(outer_reconstruction$maximum_absolute_center_error), scientific = TRUE)),
    paste("Maximum outer scale error:",
          format(max(outer_reconstruction$maximum_absolute_scale_used_error), scientific = TRUE)),
    "No inner-validation participant contributed to its scaler training set.",
    "No outer-test participant appeared in any corresponding inner assignment.",
    paste("Per-alpha one-SE tuning choices reconstructed:", nrow(alpha_choice_check)),
    paste("Deployed alpha/lambda choices reconstructed:", nrow(tuning_check)),
    "All five models used identical participants and outer folds within each comparison.",
    "Every participant had exactly one held-out prediction per repeat/model and ten before averaging.",
    paste("Finite negative repeat-level R-squared values retained:",
          qa_summary[check == "finite_negative_repeat_level_r_squared_retained", value])
  ),
  file.path(paths$validation, "NESTED_CV_QA_STATUS.txt")
)
message("Independent nested-validation QA completed successfully.")
