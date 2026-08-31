# Summarise the verified repeated nested-CV predictions, add the biomarker-only
# comparator tables, and quantify conditional participant-bootstrap uncertainty.
#
# The bootstrap deliberately resamples the fixed participant-level averages of
# ten held-out predictions. It does not refit models, retune hyperparameters or
# regenerate folds, so its intervals are conditional and exploratory.

source(file.path("R", "00_common.R"))
assert_packages(c("data.table", "digest", "R.utils"))
source(file.path("R", "09_nested_cv_helpers.R"))
message_rule("11: Predictive comparison, calibration and sensitivity")

bootstrap_replicates <- predictive_inference_settings$bootstrap_replicates
bootstrap_seed_base <- predictive_inference_settings$bootstrap_seed

averaged <- data.table::fread(
  file.path(paths$analysis, "nested_cv_averaged_predictions.csv.gz"),
  showProgress = FALSE
)
fold_predictions <- data.table::fread(
  file.path(paths$analysis, "nested_cv_fold_level_predictions.csv.gz"),
  showProgress = FALSE
)
performance <- data.table::fread(file.path(
  paths$tables, "Table_4_repeated_nested_CV_model_performance.csv"
))

expected_model_types <- model_types
standalone_models <- c(
  "Biomarker-only", "Curated pathways", "Elastic Net pathways"
)
combined_models <- c(
  "Biomarker + curated pathways", "Biomarker + Elastic Net pathways"
)

calibration_values <- function(observed, predicted) {
  if (length(observed) < 3L || stats::sd(predicted) < 1e-12) {
    return(c(intercept = NA_real_, slope = NA_real_))
  }
  slope <- stats::cov(observed, predicted) / stats::var(predicted)
  intercept <- mean(observed) - slope * mean(predicted)
  c(intercept = unname(intercept), slope = unname(slope))
}

metric_bundle <- function(observed_model, predicted_model,
                          observed_original, predicted_original) {
  model_metrics <- calculate_metrics(observed_model, predicted_model)
  original_metrics <- calculate_metrics(observed_original, predicted_original)
  calibration <- calibration_values(observed_model, predicted_model)
  c(
    model_r_squared = unname(model_metrics["r_squared"]),
    model_rmse = unname(model_metrics["rmse"]),
    model_mae = unname(model_metrics["mae"]),
    original_rmse = unname(original_metrics["rmse"]),
    original_mae = unname(original_metrics["mae"]),
    calibration_intercept = unname(calibration["intercept"]),
    calibration_slope = unname(calibration["slope"])
  )
}

percentile_interval <- function(x) {
  stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

conditional_two_sided_tail_probability <- function(draws) {
  draws <- draws[is.finite(draws)]
  if (!length(draws)) return(NA_real_)
  probability_lower <- (sum(draws <= 0) + 1) / (length(draws) + 1)
  probability_upper <- (sum(draws >= 0) + 1) / (length(draws) + 1)
  min(1, 2 * min(probability_lower, probability_upper))
}

metric_draw_rows <- list()
metric_summary_rows <- list()
delta_draw_rows <- list()
delta_summary_rows <- list()
calibration_bin_rows <- list()

combination_index <- 0L
for (dataset_name in names(datasets)) {
  for (outcome_name in names(outcome_specs)) {
    combination_index <- combination_index + 1L
    local <- averaged[
      averaged$dataset == dataset_name & averaged$outcome == outcome_name
    ]
    present_models <- sort(unique(local$model_type))
    if (!identical(present_models, sort(expected_model_types))) {
      stop("Unexpected model set for ", dataset_name, " / ", outcome_name)
    }

    model_tables <- lapply(expected_model_types, function(model_name) {
      x <- local[local$model_type == model_name]
      data.table::setorder(x, Subject_ID)
      x
    })
    names(model_tables) <- expected_model_types
    reference_subjects <- model_tables[[1L]]$Subject_ID
    if (any(!vapply(
      model_tables,
      function(x) identical(x$Subject_ID, reference_subjects),
      logical(1)
    ))) {
      stop("Model subject sets differ for ", dataset_name, " / ", outcome_name)
    }
    if (any(!vapply(
      model_tables,
      function(x) isTRUE(all.equal(
        x$observed_model, model_tables[[1L]]$observed_model,
        tolerance = 0
      )),
      logical(1)
    ))) {
      stop("Observed outcomes differ by model for ", dataset_name, " / ", outcome_name)
    }

    n <- length(reference_subjects)
    bootstrap_seed <- seed_from_key(
      bootstrap_seed_base, "conditional-participant-bootstrap",
      dataset_name, outcome_name
    )
    set.seed(bootstrap_seed)
    bootstrap_indices <- matrix(
      sample.int(n, n * bootstrap_replicates, replace = TRUE),
      nrow = bootstrap_replicates,
      ncol = n
    )

    model_draw_matrices <- list()
    for (model_name in expected_model_types) {
      x <- model_tables[[model_name]]
      point <- metric_bundle(
        x$observed_model, x$predicted_model,
        x$observed_original, x$predicted_original
      )
      draws <- matrix(
        NA_real_, nrow = bootstrap_replicates, ncol = length(point),
        dimnames = list(NULL, names(point))
      )
      for (bootstrap_index in seq_len(bootstrap_replicates)) {
        indexes <- bootstrap_indices[bootstrap_index, ]
        draws[bootstrap_index, ] <- metric_bundle(
          x$observed_model[indexes], x$predicted_model[indexes],
          x$observed_original[indexes], x$predicted_original[indexes]
        )
      }
      model_draw_matrices[[model_name]] <- draws

      metric_draw_rows[[length(metric_draw_rows) + 1L]] <-
        data.table::data.table(
          dataset = dataset_name,
          outcome = outcome_name,
          model_type = model_name,
          bootstrap_replicate = seq_len(bootstrap_replicates),
          bootstrap_seed = bootstrap_seed,
          as.data.frame(draws, check.names = FALSE)
        )

      summary <- data.table::data.table(
        dataset = dataset_name,
        dataset_label = x$dataset_label[1L],
        compartment = x$compartment[1L],
        platform = x$platform[1L],
        omic_type = x$omic_type[1L],
        outcome = outcome_name,
        outcome_label = x$outcome_label[1L],
        outcome_transform = x$outcome_transform[1L],
        model_type = model_name,
        n = n,
        bootstrap_replicates = bootstrap_replicates,
        bootstrap_seed = bootstrap_seed,
        uncertainty_scope = paste(
          "Participant bootstrap conditional on fixed averages of ten repeated",
          "held-out predictions; excludes model-refitting, tuning and",
          "fold-generation uncertainty"
        )
      )
      for (metric_name in names(point)) {
        interval <- percentile_interval(draws[, metric_name])
        summary[[paste0(metric_name, "_point")]] <- point[metric_name]
        summary[[paste0(metric_name, "_lower_95")]] <- interval[1L]
        summary[[paste0(metric_name, "_upper_95")]] <- interval[2L]
      }
      metric_summary_rows[[length(metric_summary_rows) + 1L]] <- summary

      # Five equal-sized bins are descriptive calibration summaries, not an
      # independent validation sample.
      rank_index <- rank(x$predicted_model, ties.method = "first")
      calibration_bin <- pmin(
        predictive_inference_settings$calibration_bins,
        ceiling(
          rank_index * predictive_inference_settings$calibration_bins / n
        )
      )
      bins <- data.table::data.table(
        calibration_bin = calibration_bin,
        observed_model = x$observed_model,
        predicted_model = x$predicted_model,
        observed_original = x$observed_original,
        predicted_original = x$predicted_original
      )[, .(
        bin_n = .N,
        mean_observed_model = mean(observed_model),
        mean_predicted_model = mean(predicted_model),
        mean_observed_original = mean(observed_original),
        mean_predicted_original = mean(predicted_original),
        se_observed_model = stats::sd(observed_model) / sqrt(.N),
        se_predicted_model = stats::sd(predicted_model) / sqrt(.N)
      ), by = calibration_bin]
      bins[, `:=`(
        dataset = dataset_name,
        outcome = outcome_name,
        model_type = model_name
      )]
      calibration_bin_rows[[length(calibration_bin_rows) + 1L]] <- bins
    }

    baseline_table <- model_tables[["Biomarker-only"]]
    baseline_point <- metric_bundle(
      baseline_table$observed_model, baseline_table$predicted_model,
      baseline_table$observed_original, baseline_table$predicted_original
    )
    baseline_draws <- model_draw_matrices[["Biomarker-only"]]

    for (combined_name in combined_models) {
      comparison_table <- model_tables[[combined_name]]
      comparison_point <- metric_bundle(
        comparison_table$observed_model, comparison_table$predicted_model,
        comparison_table$observed_original, comparison_table$predicted_original
      )
      comparison_draws <- model_draw_matrices[[combined_name]]
      point_delta <- c(
        delta_model_r_squared =
          unname(comparison_point["model_r_squared"] -
                   baseline_point["model_r_squared"]),
        delta_model_rmse_improvement =
          unname(baseline_point["model_rmse"] - comparison_point["model_rmse"]),
        delta_model_mae_improvement =
          unname(baseline_point["model_mae"] - comparison_point["model_mae"]),
        delta_original_rmse_improvement =
          unname(baseline_point["original_rmse"] -
                   comparison_point["original_rmse"]),
        delta_original_mae_improvement =
          unname(baseline_point["original_mae"] -
                   comparison_point["original_mae"])
      )
      delta_draws <- cbind(
        delta_model_r_squared =
          comparison_draws[, "model_r_squared"] -
            baseline_draws[, "model_r_squared"],
        delta_model_rmse_improvement =
          baseline_draws[, "model_rmse"] - comparison_draws[, "model_rmse"],
        delta_model_mae_improvement =
          baseline_draws[, "model_mae"] - comparison_draws[, "model_mae"],
        delta_original_rmse_improvement =
          baseline_draws[, "original_rmse"] -
            comparison_draws[, "original_rmse"],
        delta_original_mae_improvement =
          baseline_draws[, "original_mae"] -
            comparison_draws[, "original_mae"]
      )
      delta_draw_rows[[length(delta_draw_rows) + 1L]] <-
        data.table::data.table(
          dataset = dataset_name,
          outcome = outcome_name,
          comparator_model = "Biomarker-only",
          augmented_model = combined_name,
          bootstrap_replicate = seq_len(bootstrap_replicates),
          bootstrap_seed = bootstrap_seed,
          as.data.frame(delta_draws, check.names = FALSE)
        )
      delta_summary <- data.table::data.table(
        dataset = dataset_name,
        dataset_label = comparison_table$dataset_label[1L],
        compartment = comparison_table$compartment[1L],
        platform = comparison_table$platform[1L],
        omic_type = comparison_table$omic_type[1L],
        outcome = outcome_name,
        outcome_label = comparison_table$outcome_label[1L],
        outcome_transform = comparison_table$outcome_transform[1L],
        comparator_model = "Biomarker-only",
        augmented_model = combined_name,
        n = n,
        bootstrap_replicates = bootstrap_replicates,
        bootstrap_seed = bootstrap_seed
      )
      for (metric_name in names(point_delta)) {
        interval <- percentile_interval(delta_draws[, metric_name])
        delta_summary[[paste0(metric_name, "_point")]] <- point_delta[metric_name]
        delta_summary[[paste0(metric_name, "_lower_95")]] <- interval[1L]
        delta_summary[[paste0(metric_name, "_upper_95")]] <- interval[2L]
      }
      delta_summary$conditional_bootstrap_two_sided_p_delta_r_squared <-
        conditional_two_sided_tail_probability(
          delta_draws[, "delta_model_r_squared"]
        )
      delta_summary_rows[[length(delta_summary_rows) + 1L]] <- delta_summary
    }
  }
}

metric_draws <- data.table::rbindlist(metric_draw_rows, fill = TRUE)
metric_summary <- data.table::rbindlist(metric_summary_rows, fill = TRUE)
delta_draws <- data.table::rbindlist(delta_draw_rows, fill = TRUE)
delta_summary <- data.table::rbindlist(delta_summary_rows, fill = TRUE)
calibration_bins <- data.table::rbindlist(calibration_bin_rows, fill = TRUE)

delta_summary[, exploratory_BH_q_across_40_delta_r_squared_comparisons :=
  stats::p.adjust(
    conditional_bootstrap_two_sided_p_delta_r_squared,
    method = "BH"
  )]
delta_summary[, interpretation := data.table::fcase(
  delta_model_r_squared_point > 0 &
    delta_model_r_squared_lower_95 > 0 &
    exploratory_BH_q_across_40_delta_r_squared_comparisons < 0.05,
  paste(
    "Exploratory internally held-out improvement with a positive conditional",
    "participant-bootstrap interval and BH sensitivity; not external confirmation"
  ),
  delta_model_r_squared_point > 0 & delta_model_r_squared_lower_95 > 0,
  paste(
    "Exploratory internally held-out improvement with a positive conditional",
    "participant-bootstrap interval; not supported after the 40-comparison BH sensitivity"
  ),
  delta_model_r_squared_point > 0,
  "Positive point estimate but conditional interval includes zero",
  default = "No internally held-out improvement over the biomarker baseline"
)]
delta_summary[, inference_caution := paste(
  "Post hoc exploratory internal validation. Bootstrap intervals condition on",
  "fixed averaged held-out predictions and exclude refitting, tuning and",
  "fold-generation uncertainty."
)]
delta_summary[
  outcome == "sputum_eosinophils",
  inference_caution := paste(
    inference_caution,
    "The outcome is a concurrent same-sample phenotype, not a prognostic endpoint."
  )
]

# Requested biomarker-only and three-way standalone comparison tables. These use
# the mean of the ten complete repeated-CV metrics, matching Table 4.
biomarker_performance <- performance[
  performance$model_type == "Biomarker-only"
]
biomarker_performance[, biomarker_predictors := vapply(
  outcome,
  function(x) paste(biomarker_models[[x]], collapse = " + "),
  character(1)
)]
biomarker_performance[, interpretation := data.table::fcase(
  mean_cv_r_squared < 0,
  "Mean repeated-CV R-squared is negative; worse than predicting the held-out sample mean",
  mean_cv_r_squared < 0.10,
  "Weak positive internally held-out performance",
  mean_cv_r_squared < 0.25,
  "Modest positive internally held-out performance",
  default = "Positive internally held-out performance; external validity is untested"
)]

comparison <- performance[performance$model_type %in% standalone_models]
comparison[, best_r_squared_in_three_model_comparison := max(mean_cv_r_squared),
           by = .(dataset, outcome)]
comparison[, difference_from_best_r_squared :=
  mean_cv_r_squared - best_r_squared_in_three_model_comparison]
comparison[, number_of_predictors := data.table::fcase(
  model_type == "Elastic Net pathways",
  sprintf("39 candidates; mean %.1f selected", mean_selected_predictors),
  default = as.character(candidate_predictors)
)]
comparison[, interpretation := data.table::fcase(
  mean_cv_r_squared < 0,
  "Negative mean repeated-CV R-squared; no useful internal predictive evidence",
  abs(difference_from_best_r_squared) < 1e-12,
  "Highest mean repeated-CV R-squared among these three standalone models; exploratory internal comparison",
  difference_from_best_r_squared >= -0.02,
  "Mean repeated-CV R-squared is within 0.02 of the best standalone model",
  default = "Lower mean repeated-CV R-squared than the best standalone model"
)]
comparison[, `:=`(
  mean_cross_validated_r_squared = mean_cv_r_squared,
  RMSE = mean_cv_rmse_original_scale,
  MAE = mean_cv_mae_original_scale
)]
comparison <- comparison[, .(
  dataset, dataset_label, compartment, platform, omic_type,
  outcome, outcome_label, outcome_transform, model_type, n,
  mean_cross_validated_r_squared, RMSE, MAE,
  number_of_predictors, candidate_predictors, mean_selected_predictors,
  interpretation
)]

# Locked clipping-at-zero sensitivity for outcomes whose natural original scale
# is non-negative. Model-scale R-squared is never clipped or recalculated here.
clipping_rows <- list()
for (prediction_level in c("Averaged repeated OOF", "Participant-repeat OOF")) {
  source_table <- if (prediction_level == "Averaged repeated OOF") {
    averaged
  } else {
    fold_predictions
  }
  source_table <- source_table[source_table$outcome != "FEV1"]
  clipping_rows[[prediction_level]] <- source_table[, {
    clipped <- pmax(predicted_original, 0)
    before <- calculate_metrics(observed_original, predicted_original)
    after <- calculate_metrics(observed_original, clipped)
    .(
      prediction_rows = .N,
      below_zero_prediction_rows = sum(predicted_original < 0),
      below_zero_prediction_percent = 100 * mean(predicted_original < 0),
      original_scale_RMSE_unclipped = unname(before["rmse"]),
      original_scale_RMSE_clipped_at_zero = unname(after["rmse"]),
      RMSE_improvement_after_clipping =
        unname(before["rmse"] - after["rmse"]),
      original_scale_MAE_unclipped = unname(before["mae"]),
      original_scale_MAE_clipped_at_zero = unname(after["mae"]),
      MAE_improvement_after_clipping =
        unname(before["mae"] - after["mae"]),
      model_scale_R_squared_was_not_clipped = TRUE
    )
  }, by = .(
    dataset, dataset_label, compartment, platform, omic_type,
    outcome, outcome_label, model_type
  )][, prediction_level := prediction_level]
}
clipping_sensitivity <- data.table::rbindlist(clipping_rows, fill = TRUE)

# Core checks.
if (nrow(metric_summary) != 100L) stop("Expected 100 model metric summaries.")
if (nrow(delta_summary) != 40L) stop("Expected 40 combined-model comparisons.")
if (nrow(biomarker_performance) != 20L) stop("Expected 20 biomarker-only rows.")
if (nrow(comparison) != 60L) stop("Expected 60 three-model comparison rows.")
if (any(metric_summary$n != performance$n[match(
  paste(metric_summary$dataset, metric_summary$outcome, metric_summary$model_type),
  paste(performance$dataset, performance$outcome, performance$model_type)
)])) stop("Metric summary n differs from repeated-CV performance n.")

data.table::fwrite(
  biomarker_performance,
  file.path(paths$tables, "Table_6_biomarker_only_model_performance.csv")
)
data.table::fwrite(
  comparison,
  file.path(paths$tables, "Table_7_pathway_vs_biomarker_model_comparison.csv")
)
data.table::fwrite(
  metric_summary,
  file.path(paths$tables, "Table_S_averaged_OOF_metrics_and_calibration.csv")
)
data.table::fwrite(
  calibration_bins,
  file.path(paths$tables, "Table_S_predictive_calibration_bins.csv")
)
data.table::fwrite(
  delta_summary,
  file.path(paths$tables, "Table_S_combined_vs_biomarker_delta_inference.csv")
)
data.table::fwrite(
  clipping_sensitivity,
  file.path(paths$tables, "Table_S_below_zero_prediction_clipping_sensitivity.csv")
)
data.table::fwrite(
  metric_draws,
  file.path(paths$diagnostics, "predictive_conditional_bootstrap_metric_draws.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  delta_draws,
  file.path(paths$diagnostics, "predictive_conditional_bootstrap_delta_draws.csv.gz"),
  compress = "gzip"
)

writeLines(
  c(
    "PREDICTIVE COMPARATOR AND CONDITIONAL INFERENCE: PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste("Biomarker-only dataset-outcome rows:", nrow(biomarker_performance)),
    paste("Three-model standalone comparison rows:", nrow(comparison)),
    paste("Combined-versus-biomarker delta comparisons:", nrow(delta_summary)),
    paste("Conditional participant-bootstrap replicates per comparison:",
          bootstrap_replicates),
    paste("Averaged below-zero predictions for non-FEV1 outcomes:",
          clipping_sensitivity[
            prediction_level == "Averaged repeated OOF",
            sum(below_zero_prediction_rows)
          ]),
    paste("Participant-repeat below-zero predictions for non-FEV1 outcomes:",
          clipping_sensitivity[
            prediction_level == "Participant-repeat OOF",
            sum(below_zero_prediction_rows)
          ]),
    "All 40 combined-model delta-R-squared comparisons receive a BH multiplicity sensitivity.",
    "All intervals are participant bootstraps conditional on fixed averages of ten held-out predictions.",
    "Intervals do not propagate model-refitting, tuning or fold-generation uncertainty.",
    "Predictive comparisons are exploratory internal validation, not external confirmation or clinical utility.",
    "Sputum-eosinophil estimation is concurrent same-sample phenotype capture, not prognosis.",
    "No model-scale R-squared value was clipped."
  ),
  file.path(paths$validation, "PREDICTIVE_INFERENCE_STATUS.txt")
)
message("Predictive comparison and sensitivity completed successfully.")
