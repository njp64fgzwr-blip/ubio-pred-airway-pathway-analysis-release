# Five-model repeated nested cross-validation with strictly fold-contained
# preprocessing and identical common cases/outer folds.

source(file.path("R", "00_common.R"))
assert_packages(c("data.table", "digest", "glmnet"))
source(file.path("R", "09_nested_cv_helpers.R"))
message_rule("09: Corrected repeated nested cross-validation")

pathway_names <- names(readRDS(
  file.path(paths$analysis, "airway_39_signatures_harmonised.rds")
))
outer_repeats <- cv_settings$outer_repeats
outer_folds <- cv_settings$outer_folds
inner_folds <- cv_settings$inner_folds

prediction_rows <- list()
outer_assignment_rows <- list()
inner_assignment_rows <- list()
inner_scaler_rows <- list()
outer_scaler_rows <- list()
fold_loss_rows <- list()
aggregate_loss_rows <- list()
alpha_choice_rows <- list()
tuning_choice_rows <- list()
coefficient_rows <- list()
outer_model_rows <- list()
common_case_rows <- list()
common_predictor_rows <- list()

append_metadata <- function(table, dataset_name, outcome_name, model_type,
                            repeat_index, outer_fold) {
  table[, `:=`(
    dataset = dataset_name,
    outcome = outcome_name,
    model_type = model_type,
    cv_repeat = repeat_index,
    outer_fold = outer_fold
  )]
  table
}

for (dataset_name in names(datasets)) {
  merged <- readRDS(file.path(
    paths$analysis, paste0(dataset_name, "_pathway_clinical_analysis.rds")
  ))

  for (outcome_name in names(outcome_specs)) {
    message("Nested CV: ", datasets[[dataset_name]]$label, " / ", outcome_name)
    outcome <- outcome_specs[[outcome_name]]
    common <- prepare_predictive_common_cases(merged, outcome_name, pathway_names)
    model_matrices <- make_model_matrices(common, outcome_name)
    n <- length(common$Subject_ID)

    common_case_rows[[paste(dataset_name, outcome_name, sep = "__")]] <-
      data.table::data.table(
        dataset = dataset_name,
        dataset_label = datasets[[dataset_name]]$label,
        compartment = datasets[[dataset_name]]$compartment,
        platform = datasets[[dataset_name]]$platform,
        omic_type = datasets[[dataset_name]]$omic_type,
        outcome = outcome_name,
        outcome_label = outcome$label,
        Subject_ID = common$Subject_ID,
        observed_original = common$y_original,
        observed_model = common$y_model,
        outcome_transform = outcome$transform
      )

    predictor_values <- cbind(common$x_biomarkers, common$x_pathways)
    predictor_types <- c(
      rep("Clinical biomarker (log1p)", ncol(common$x_biomarkers)),
      rep("Normalized ssGSEA pathway score", ncol(common$x_pathways))
    )
    common_predictor_rows[[paste(dataset_name, outcome_name, sep = "__")]] <-
      data.table::data.table(
        dataset = dataset_name,
        outcome = outcome_name,
        Subject_ID = rep(common$Subject_ID, times = ncol(predictor_values)),
        predictor = rep(colnames(predictor_values), each = n),
        predictor_type = rep(predictor_types, each = n),
        value = as.numeric(predictor_values)
      )

    for (repeat_index in seq_len(outer_repeats)) {
      outer_seed <- seed_from_key(
        cv_settings$seed, "outer", dataset_name, outcome_name, repeat_index
      )
      outer_fold_id <- make_stratified_folds(
        common$y_model, outer_folds, outer_seed
      )
      outer_assignment_rows[[length(outer_assignment_rows) + 1L]] <-
        data.table::data.table(
          dataset = dataset_name,
          outcome = outcome_name,
          cv_repeat = repeat_index,
          Subject_ID = common$Subject_ID,
          outer_fold = outer_fold_id,
          assignment_seed = outer_seed
        )

      for (outer_fold in seq_len(outer_folds)) {
        test_index <- outer_fold_id == outer_fold
        train_index <- !test_index
        if (!any(test_index) || !any(train_index)) stop("Empty outer split.")
        inner_seed <- seed_from_key(
          cv_settings$seed, "inner", dataset_name, outcome_name,
          repeat_index, outer_fold
        )
        inner_fold_id <- make_stratified_folds(
          common$y_model[train_index], inner_folds, inner_seed
        )
        inner_assignment_rows[[length(inner_assignment_rows) + 1L]] <-
          data.table::data.table(
            dataset = dataset_name,
            outcome = outcome_name,
            cv_repeat = repeat_index,
            outer_fold = outer_fold,
            Subject_ID = common$Subject_ID[train_index],
            inner_fold = inner_fold_id,
            assignment_seed = inner_seed
          )

        for (model_type in model_types) {
          model_spec <- model_matrices[[model_type]]
          x <- model_spec$x
          x_train <- x[train_index, , drop = FALSE]
          x_test <- x[test_index, , drop = FALSE]
          y_train <- common$y_model[train_index]

          if (model_spec$engine == "OLS") {
            fitted <- fit_outer_ols(x_train, y_train, x_test)
            selected_alpha <- selected_fraction <- selected_lambda <- NA_real_
            inner_selected_mse <- NA_real_
          } else {
            tuned <- tune_elastic_net_inner(
              x_train = x_train,
              y_train = y_train,
              inner_fold_id = inner_fold_id,
              penalty_factor = model_spec$penalty,
              alpha_grid = cv_settings$alpha_grid,
              lambda_fraction_grid = cv_settings$lambda_fraction_grid
            )

            inner_scalers <- append_metadata(
              tuned$scaler_rows, dataset_name, outcome_name, model_type,
              repeat_index, outer_fold
            )
            inner_scalers[, inner_assignment_seed := inner_seed]
            inner_scaler_rows[[length(inner_scaler_rows) + 1L]] <- inner_scalers

            losses <- append_metadata(
              tuned$fold_losses, dataset_name, outcome_name, model_type,
              repeat_index, outer_fold
            )
            fold_loss_rows[[length(fold_loss_rows) + 1L]] <- losses
            aggregate_losses <- append_metadata(
              tuned$aggregate_losses, dataset_name, outcome_name, model_type,
              repeat_index, outer_fold
            )
            aggregate_loss_rows[[length(aggregate_loss_rows) + 1L]] <- aggregate_losses
            alpha_choices <- append_metadata(
              tuned$alpha_choices, dataset_name, outcome_name, model_type,
              repeat_index, outer_fold
            )
            alpha_choice_rows[[length(alpha_choice_rows) + 1L]] <- alpha_choices

            fitted <- fit_outer_elastic_net(
              x_train = x_train,
              y_train = y_train,
              x_test = x_test,
              alpha = tuned$chosen_alpha,
              lambda_fraction = tuned$chosen_lambda_fraction,
              penalty_factor = model_spec$penalty
            )
            selected_alpha <- tuned$chosen_alpha
            selected_fraction <- tuned$chosen_lambda_fraction
            selected_lambda <- fitted$lambda_actual
            inner_selected_mse <- tuned$chosen_inner_pooled_mse
            tuning_choice_rows[[length(tuning_choice_rows) + 1L]] <-
              data.table::data.table(
                dataset = dataset_name,
                outcome = outcome_name,
                model_type = model_type,
                cv_repeat = repeat_index,
                outer_fold = outer_fold,
                outer_training_n = sum(train_index),
                outer_test_n = sum(test_index),
                inner_folds = inner_folds,
                alpha = selected_alpha,
                lambda_fraction = selected_fraction,
                outer_training_lambda_max = fitted$lambda_max,
                deployed_lambda = selected_lambda,
                selected_inner_pooled_mse = inner_selected_mse,
                alpha_tie_tolerance = tuned$alpha_tie_tolerance,
                alpha_tie_candidate_n = tuned$alpha_tie_candidate_n,
                outer_assignment_seed = outer_seed,
                inner_assignment_seed = inner_seed,
                tuning_rule = cv_settings$tuning_rule
              )
          }

          outer_scalers <- append_metadata(
            scaler_table(fitted$scaler), dataset_name, outcome_name, model_type,
            repeat_index, outer_fold
          )
          outer_scalers[, `:=`(
            outer_training_n = sum(train_index),
            outer_test_n = sum(test_index)
          )]
          outer_scaler_rows[[length(outer_scaler_rows) + 1L]] <- outer_scalers

          coefficient <- as.numeric(fitted$coefficients)
          names(coefficient) <- names(fitted$coefficients)
          feature_type <- ifelse(
            names(coefficient) %in% pathway_names, "Pathway", "Clinical biomarker"
          )
          coefficient_rows[[length(coefficient_rows) + 1L]] <-
            data.table::data.table(
              dataset = dataset_name,
              outcome = outcome_name,
              model_type = model_type,
              cv_repeat = repeat_index,
              outer_fold = outer_fold,
              predictor = names(coefficient),
              predictor_type = feature_type,
              standardized_predictor_coefficient = coefficient,
              selected = abs(coefficient) > 1e-10,
              coefficient_sign = sign(coefficient)
            )

          prediction_model <- fitted$prediction
          prediction_original <- inverse_outcome_transform(
            prediction_model, outcome$transform
          )
          prediction_rows[[length(prediction_rows) + 1L]] <-
            data.table::data.table(
              dataset = dataset_name,
              dataset_label = datasets[[dataset_name]]$label,
              compartment = datasets[[dataset_name]]$compartment,
              platform = datasets[[dataset_name]]$platform,
              omic_type = datasets[[dataset_name]]$omic_type,
              outcome = outcome_name,
              outcome_label = outcome$label,
              outcome_transform = outcome$transform,
              model_type = model_type,
              cv_repeat = repeat_index,
              outer_fold = outer_fold,
              Subject_ID = common$Subject_ID[test_index],
              observed_model = common$y_model[test_index],
              predicted_model = prediction_model,
              observed_original = common$y_original[test_index],
              predicted_original = prediction_original,
              predicted_original_below_zero = prediction_original < 0
            )

          selected_pathways <- sum(
            abs(coefficient[names(coefficient) %in% pathway_names]) > 1e-10
          )
          selected_biomarkers <- sum(
            abs(coefficient[!names(coefficient) %in% pathway_names]) > 1e-10
          )
          outer_model_rows[[length(outer_model_rows) + 1L]] <-
            data.table::data.table(
              dataset = dataset_name,
              outcome = outcome_name,
              model_type = model_type,
              engine = model_spec$engine,
              cv_repeat = repeat_index,
              outer_fold = outer_fold,
              outer_training_n = sum(train_index),
              outer_test_n = sum(test_index),
              candidate_predictors = ncol(x),
              selected_predictors = sum(abs(coefficient) > 1e-10),
              selected_pathways = selected_pathways,
              selected_biomarkers = selected_biomarkers,
              alpha = selected_alpha,
              lambda_fraction = selected_fraction,
              deployed_lambda = selected_lambda,
              inner_selected_pooled_mse = inner_selected_mse,
              outer_assignment_seed = outer_seed,
              inner_assignment_seed = inner_seed
            )
        }
      }
    }
  }
}

message("Combining and writing nested-validation records...")
predictions <- data.table::rbindlist(prediction_rows)
outer_assignments <- data.table::rbindlist(outer_assignment_rows)
inner_assignments <- data.table::rbindlist(inner_assignment_rows)
inner_scalers <- data.table::rbindlist(inner_scaler_rows)
outer_scalers <- data.table::rbindlist(outer_scaler_rows)
fold_losses <- data.table::rbindlist(fold_loss_rows)
aggregate_losses <- data.table::rbindlist(aggregate_loss_rows)
alpha_choices <- data.table::rbindlist(alpha_choice_rows)
tuning_choices <- data.table::rbindlist(tuning_choice_rows)
coefficients <- data.table::rbindlist(coefficient_rows)
outer_models <- data.table::rbindlist(outer_model_rows)
common_cases <- data.table::rbindlist(common_case_rows)
common_predictors <- data.table::rbindlist(common_predictor_rows)

repeat_metrics <- predictions[, {
  model_metrics <- calculate_metrics(observed_model, predicted_model)
  original_metrics <- calculate_metrics(observed_original, predicted_original)
  .(
    n = .N,
    model_scale_r_squared = unname(model_metrics["r_squared"]),
    model_scale_rmse = unname(model_metrics["rmse"]),
    model_scale_mae = unname(model_metrics["mae"]),
    model_scale_pearson_r = suppressWarnings(stats::cor(
      observed_model, predicted_model, method = "pearson"
    )),
    original_scale_r_squared = unname(original_metrics["r_squared"]),
    original_scale_rmse = unname(original_metrics["rmse"]),
    original_scale_mae = unname(original_metrics["mae"]),
    below_zero_predictions = sum(predicted_original_below_zero)
  )
}, by = .(
  dataset, dataset_label, compartment, platform, omic_type,
  outcome, outcome_label, outcome_transform, model_type, cv_repeat
)]

model_complexity <- outer_models[, .(
  candidate_predictors = unique(candidate_predictors),
  mean_selected_predictors = mean(selected_predictors),
  sd_selected_predictors = stats::sd(selected_predictors),
  mean_selected_pathways = mean(selected_pathways),
  mean_selected_biomarkers = mean(selected_biomarkers)
), by = .(dataset, outcome, model_type)]

performance_summary <- repeat_metrics[, .(
  n = unique(n),
  repeats = .N,
  mean_cv_r_squared = mean(model_scale_r_squared),
  sd_cv_r_squared = stats::sd(model_scale_r_squared),
  mean_cv_rmse_model_scale = mean(model_scale_rmse),
  sd_cv_rmse_model_scale = stats::sd(model_scale_rmse),
  mean_cv_mae_model_scale = mean(model_scale_mae),
  sd_cv_mae_model_scale = stats::sd(model_scale_mae),
  mean_cv_pearson_r = mean(model_scale_pearson_r),
  sd_cv_pearson_r = stats::sd(model_scale_pearson_r),
  mean_cv_rmse_original_scale = mean(original_scale_rmse),
  sd_cv_rmse_original_scale = stats::sd(original_scale_rmse),
  mean_cv_mae_original_scale = mean(original_scale_mae),
  sd_cv_mae_original_scale = stats::sd(original_scale_mae),
  finite_negative_r_squared_repeats = sum(
    is.finite(model_scale_r_squared) & model_scale_r_squared < 0
  )
), by = .(
  dataset, dataset_label, compartment, platform, omic_type,
  outcome, outcome_label, outcome_transform, model_type
)]
performance_summary <- merge(
  performance_summary, model_complexity,
  by = c("dataset", "outcome", "model_type"), all.x = TRUE
)

averaged_predictions <- predictions[, .(
  observed_model = unique(observed_model),
  predicted_model = mean(predicted_model),
  predicted_model_sd = stats::sd(predicted_model),
  observed_original = unique(observed_original),
  predicted_original = mean(predicted_original),
  predicted_original_sd = stats::sd(predicted_original),
  held_out_prediction_count = .N,
  below_zero_repeat_predictions = sum(predicted_original_below_zero)
), by = .(
  dataset, dataset_label, compartment, platform, omic_type,
  outcome, outcome_label, outcome_transform, model_type, Subject_ID
)]

selection_stability <- coefficients[
  predictor_type == "Pathway" &
    model_type %in% c("Elastic Net pathways", "Biomarker + Elastic Net pathways"),
  .(
    outer_models = .N,
    selection_count = sum(selected),
    selection_frequency = mean(selected),
    positive_selected_count = sum(selected & coefficient_sign > 0),
    negative_selected_count = sum(selected & coefficient_sign < 0),
    sign_consistency = if (sum(selected) > 0) {
      max(sum(selected & coefficient_sign > 0),
          sum(selected & coefficient_sign < 0)) / sum(selected)
    } else NA_real_,
    predominant_sign = data.table::fcase(
      sum(selected & coefficient_sign > 0) > sum(selected & coefficient_sign < 0), "Positive",
      sum(selected & coefficient_sign < 0) > sum(selected & coefficient_sign > 0), "Negative",
      sum(selected) == 0, "Never selected",
      default = "Mixed/tied"
    ),
    median_coefficient_all_outer_models = stats::median(
      standardized_predictor_coefficient
    ),
    median_nonzero_coefficient = if (sum(selected)) {
      stats::median(standardized_predictor_coefficient[selected])
    } else NA_real_
  ), by = .(dataset, outcome, model_type, pathway = predictor)
]
selection_stability[, pathway_label := unname(pathway_labels[pathway])]
selection_stability[, caution := "Selection frequency is not a p value"]

# Core structural checks before any file is accepted.
expected_common_n <- data.table::data.table(
  dataset = rep(names(datasets), each = 5L),
  outcome = rep(names(outcome_specs), times = 4L),
  expected_n = c(
    rep(113L, 4L), 99L,
    rep(64L, 4L), 40L,
    rep(44L, 4L), 31L,
    rep(54L, 4L), 36L
  )
)
observed_common_n <- common_cases[, .(observed_n = .N), by = .(dataset, outcome)]
n_check <- merge(expected_common_n, observed_common_n, by = c("dataset", "outcome"))
if (any(n_check$expected_n != n_check$observed_n)) {
  stop("Common-case n differs from the prespecified audit counts.")
}
if (nrow(outer_assignments) != 13060L) {
  stop("Expected 13,060 outer assignment rows; found ", nrow(outer_assignments), ".")
}
if (nrow(predictions) != 65300L) {
  stop("Expected 65,300 held-out prediction rows; found ", nrow(predictions), ".")
}
if (nrow(averaged_predictions) != 6530L) {
  stop("Expected 6,530 averaged prediction rows; found ", nrow(averaged_predictions), ".")
}
if (any(averaged_predictions$held_out_prediction_count != outer_repeats)) {
  stop("Not every participant has ten held-out predictions per model.")
}
duplicate_check <- predictions[, .N, by = .(
  dataset, outcome, model_type, cv_repeat, Subject_ID
)]
if (any(duplicate_check$N != 1L)) {
  stop("A participant does not have exactly one held-out prediction per repeat/model.")
}
if (any(!is.finite(repeat_metrics$model_scale_r_squared))) {
  stop("Non-finite repeat-level R-squared encountered.")
}

data.table::fwrite(
  common_cases,
  file.path(paths$analysis, "predictive_common_case_outcomes.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  common_predictors,
  file.path(paths$analysis, "predictive_common_case_predictors_long.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  outer_assignments,
  file.path(paths$analysis, "nested_cv_outer_fold_assignments.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  inner_assignments,
  file.path(paths$analysis, "nested_cv_inner_fold_assignments.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  predictions,
  file.path(paths$analysis, "nested_cv_fold_level_predictions.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  averaged_predictions,
  file.path(paths$analysis, "nested_cv_averaged_predictions.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  inner_scalers,
  file.path(paths$diagnostics, "nested_cv_inner_training_scalers.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  outer_scalers,
  file.path(paths$diagnostics, "nested_cv_outer_training_scalers.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  fold_losses,
  file.path(paths$diagnostics, "nested_cv_inner_fold_tuning_losses.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  aggregate_losses,
  file.path(paths$diagnostics, "nested_cv_inner_aggregate_tuning_losses.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  alpha_choices,
  file.path(paths$diagnostics, "nested_cv_alpha_one_se_choices.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  tuning_choices,
  file.path(paths$tables, "Table_S_elastic_net_tuning_choices.csv")
)
data.table::fwrite(
  coefficients,
  file.path(paths$diagnostics, "nested_cv_outer_model_coefficients.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  outer_models,
  file.path(paths$diagnostics, "nested_cv_outer_model_records.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  repeat_metrics,
  file.path(paths$tables, "Table_S_nested_cv_repeat_metrics.csv")
)
data.table::fwrite(
  performance_summary,
  file.path(paths$tables, "Table_4_repeated_nested_CV_model_performance.csv")
)
data.table::fwrite(
  selection_stability,
  file.path(paths$tables, "Table_5_elastic_net_selection_stability.csv")
)

settings <- data.table::data.table(
  setting = c(
    "model_types", "outer_repeats", "outer_folds", "inner_folds",
    "alpha_grid", "lambda_fraction_grid", "tuning_rule",
    "predictor_preprocessing", "biomarker_penalty_in_combined_model",
    "primary_r_squared_scale", "negative_r_squared_rule", "seed"
  ),
  value = c(
    paste(model_types, collapse = "; "),
    outer_repeats, outer_folds, inner_folds,
    paste(cv_settings$alpha_grid, collapse = ";"),
    paste(cv_settings$lambda_fraction_grid, collapse = ";"),
    cv_settings$tuning_rule,
    paste(
      "Centre/scale fitted on each inner-training fold for tuning;",
      "refitted on complete outer-training data before outer-test prediction"
    ),
    "Clinical biomarkers unpenalised (penalty.factor=0); pathways penalised",
    "Transformed outcome/model scale",
    "Retain all finite negative values; never truncate to zero",
    cv_settings$seed
  )
)
data.table::fwrite(
  settings,
  file.path(paths$tables, "Table_S_elastic_net_model_settings.csv")
)

writeLines(
  c(
    "REPEATED NESTED CROSS-VALIDATION: COMPUTATION PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste("Common-case dataset-outcome cohorts:", nrow(n_check)),
    paste("Outer assignment rows:", nrow(outer_assignments)),
    paste("Inner assignment rows:", nrow(inner_assignments)),
    paste("Held-out participant-repeat-model rows:", nrow(predictions)),
    paste("Averaged participant-model rows:", nrow(averaged_predictions)),
    paste("Inner scaler rows:", nrow(inner_scalers)),
    paste("Inner fold tuning-loss rows:", nrow(fold_losses)),
    paste("Finite negative repeat-level R-squared values:",
          sum(repeat_metrics$model_scale_r_squared < 0)),
    "All five models use identical common cases and outer folds within each dataset-outcome comparison.",
    "Outer-test participants were not used in inner training, validation, scaling or tuning.",
    "Inner-fold scalers were fitted only on the corresponding inner-training participants.",
    "Clinical biomarkers were unpenalised in biomarker-plus-Elastic-Net models.",
    "Negative cross-validated R-squared values were retained.",
    "Selection frequency is an outer-model stability descriptor, not a p value.",
    "This is exploratory internal validation, not external confirmation or clinical utility."
  ),
  file.path(paths$validation, "NESTED_CV_COMPUTATION_STATUS.txt")
)
message("Nested validation computation completed successfully.")
