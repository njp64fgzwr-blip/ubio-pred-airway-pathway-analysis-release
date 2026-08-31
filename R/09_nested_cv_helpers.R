# Helper functions for strictly fold-contained repeated nested validation.

fit_scaler <- function(x) {
  x <- as.matrix(x)
  centers <- colMeans(x)
  raw_scales <- apply(x, 2L, stats::sd)
  used_scales <- raw_scales
  used_scales[!is.finite(used_scales) | used_scales == 0] <- 1
  list(center = centers, scale = used_scales, raw_scale = raw_scales)
}

apply_scaler <- function(x, scaler) {
  x <- as.matrix(x)
  sweep(sweep(x, 2L, scaler$center, "-"), 2L, scaler$scale, "/")
}

scaler_table <- function(scaler) {
  data.table::data.table(
    predictor = names(scaler$center),
    center = as.numeric(scaler$center),
    raw_sd = as.numeric(scaler$raw_scale),
    scale_used = as.numeric(scaler$scale),
    zero_variance = !is.finite(scaler$raw_scale) | scaler$raw_scale == 0
  )
}

seed_from_key <- function(...) {
  key <- paste(..., collapse = "||")
  value <- as.double(digest::digest2int(key))
  as.integer(((value %% 2000000000) + 2000000000) %% 2000000000 + 1)
}

transform_biomarker_predictors <- function(data, variables) {
  output <- as.matrix(data[, ..variables])
  storage.mode(output) <- "double"
  if (any(output < 0, na.rm = TRUE)) {
    stop("Negative biomarker predictor value encountered.")
  }
  output <- log1p(output)
  colnames(output) <- variables
  output
}

prepare_predictive_common_cases <- function(merged, outcome_name, pathway_names) {
  outcome <- outcome_specs[[outcome_name]]
  biomarkers <- biomarker_models[[outcome_name]]
  required <- c("Subject_ID", outcome$column, biomarkers, pathway_names)
  source <- data.table::copy(merged[, ..required])
  complete <- stats::complete.cases(source)
  source <- source[complete]
  y_original <- as.numeric(source[[outcome$column]])
  y_model <- apply_outcome_transform(y_original, outcome$transform)
  x_biomarkers <- transform_biomarker_predictors(source, biomarkers)
  x_pathways <- as.matrix(source[, ..pathway_names])
  storage.mode(x_pathways) <- "double"
  list(
    Subject_ID = source$Subject_ID,
    y_original = y_original,
    y_model = y_model,
    x_biomarkers = x_biomarkers,
    x_pathways = x_pathways,
    biomarker_names = biomarkers,
    pathway_names = pathway_names,
    complete_source_rows = which(complete)
  )
}

make_model_matrices <- function(common, outcome_name) {
  curated <- curated_pathway_models[[outcome_name]]
  x_curated <- common$x_pathways[, curated, drop = FALSE]
  list(
    `Biomarker-only` = list(
      x = common$x_biomarkers,
      penalty = NULL,
      engine = "OLS"
    ),
    `Curated pathways` = list(
      x = x_curated,
      penalty = NULL,
      engine = "OLS"
    ),
    `Elastic Net pathways` = list(
      x = common$x_pathways,
      penalty = rep(1, ncol(common$x_pathways)),
      engine = "Elastic Net"
    ),
    `Biomarker + curated pathways` = list(
      x = cbind(common$x_biomarkers, x_curated),
      penalty = NULL,
      engine = "OLS"
    ),
    `Biomarker + Elastic Net pathways` = list(
      x = cbind(common$x_biomarkers, common$x_pathways),
      penalty = c(
        rep(0, ncol(common$x_biomarkers)),
        rep(1, ncol(common$x_pathways))
      ),
      engine = "Elastic Net"
    )
  )
}

fit_outer_ols <- function(x_train, y_train, x_test) {
  scaler <- fit_scaler(x_train)
  train_scaled <- apply_scaler(x_train, scaler)
  test_scaled <- apply_scaler(x_test, scaler)
  training <- data.frame(outcome = y_train, train_scaled, check.names = FALSE)
  testing <- data.frame(test_scaled, check.names = FALSE)
  formula <- stats::reformulate(colnames(train_scaled), response = "outcome")
  fit <- stats::lm(formula, data = training)
  prediction <- suppressWarnings(stats::predict(fit, newdata = testing))
  coefficients <- stats::coef(fit)[colnames(train_scaled)]
  if (any(!is.finite(prediction)) || any(!is.finite(coefficients))) {
    stop("Non-finite OLS prediction or coefficient in an outer fold.")
  }
  list(
    prediction = as.numeric(prediction),
    coefficients = coefficients,
    scaler = scaler,
    rank = fit$rank,
    residual_df = stats::df.residual(fit)
  )
}

get_glmnet_lambda_max <- function(x, y, alpha, penalty_factor) {
  probe <- suppressWarnings(glmnet::glmnet(
    x = x,
    y = y,
    family = "gaussian",
    alpha = alpha,
    penalty.factor = penalty_factor,
    standardize = FALSE,
    intercept = TRUE,
    nlambda = 5L,
    lambda.min.ratio = 0.10,
    thresh = 1e-8,
    maxit = 100000L
  ))
  maximum <- max(probe$lambda)
  if (!is.finite(maximum) || maximum <= 0) stop("Invalid glmnet lambda maximum.")
  maximum
}

tune_elastic_net_inner <- function(x_train, y_train, inner_fold_id,
                                   penalty_factor, alpha_grid,
                                   lambda_fraction_grid) {
  fold_loss_rows <- list()
  scaler_rows <- list()
  fold_levels <- sort(unique(inner_fold_id))

  for (inner_validation_fold in fold_levels) {
    inner_training <- inner_fold_id != inner_validation_fold
    inner_validation <- !inner_training
    scaler <- fit_scaler(x_train[inner_training, , drop = FALSE])
    scaled_training <- apply_scaler(
      x_train[inner_training, , drop = FALSE], scaler
    )
    scaled_validation <- apply_scaler(
      x_train[inner_validation, , drop = FALSE], scaler
    )
    scaler_rows[[length(scaler_rows) + 1L]] <- cbind(
      data.table::data.table(
        inner_validation_fold = inner_validation_fold,
        inner_training_n = sum(inner_training),
        inner_validation_n = sum(inner_validation)
      ),
      scaler_table(scaler)
    )

    for (alpha_value in alpha_grid) {
      lambda_max <- get_glmnet_lambda_max(
        scaled_training, y_train[inner_training], alpha_value, penalty_factor
      )
      lambda_values <- lambda_max * lambda_fraction_grid
      fit <- suppressWarnings(glmnet::glmnet(
        x = scaled_training,
        y = y_train[inner_training],
        family = "gaussian",
        alpha = alpha_value,
        lambda = lambda_values,
        penalty.factor = penalty_factor,
        standardize = FALSE,
        intercept = TRUE,
        thresh = 1e-8,
        maxit = 100000L
      ))
      prediction <- as.matrix(stats::predict(
        fit, newx = scaled_validation, s = lambda_values
      ))
      if (ncol(prediction) != length(lambda_values)) {
        stop("Elastic Net prediction path length differs from lambda grid.")
      }
      error <- sweep(prediction, 1L, y_train[inner_validation], "-")
      sse <- colSums(error^2)
      fold_loss_rows[[length(fold_loss_rows) + 1L]] <- data.table::data.table(
        inner_validation_fold = inner_validation_fold,
        alpha = alpha_value,
        lambda_index = seq_along(lambda_fraction_grid),
        lambda_fraction = lambda_fraction_grid,
        lambda_actual = lambda_values,
        inner_training_n = sum(inner_training),
        inner_validation_n = sum(inner_validation),
        sse = sse,
        mse = sse / sum(inner_validation)
      )
    }
  }

  fold_losses <- data.table::rbindlist(fold_loss_rows)
  aggregate_losses <- fold_losses[, .(
    inner_folds = .N,
    validation_n_total = sum(inner_validation_n),
    sse_total = sum(sse),
    pooled_mse = sum(sse) / sum(inner_validation_n),
    mean_fold_mse = mean(mse),
    sd_fold_mse = stats::sd(mse),
    se_fold_mse = stats::sd(mse) / sqrt(.N),
    lambda_actual_min = min(lambda_actual),
    lambda_actual_median = stats::median(lambda_actual),
    lambda_actual_max = max(lambda_actual)
  ), by = .(alpha, lambda_index, lambda_fraction)]

  alpha_choices <- aggregate_losses[, {
    minimum_index <- which.min(mean_fold_mse)
    threshold <- mean_fold_mse[minimum_index] + se_fold_mse[minimum_index]
    eligible <- which(mean_fold_mse <= threshold)
    selected_index <- eligible[which.max(lambda_fraction[eligible])]
    .(
      minimum_lambda_index = lambda_index[minimum_index],
      minimum_lambda_fraction = lambda_fraction[minimum_index],
      minimum_mean_fold_mse = mean_fold_mse[minimum_index],
      one_se_threshold = threshold,
      selected_lambda_index = lambda_index[selected_index],
      selected_lambda_fraction = lambda_fraction[selected_index],
      selected_mean_fold_mse = mean_fold_mse[selected_index],
      selected_pooled_mse = pooled_mse[selected_index],
      selected_se_fold_mse = se_fold_mse[selected_index]
    )
  }, by = alpha]

  # Treat losses that differ only at floating-point noise as ties. This keeps
  # the declared tie-break reconstructable after loss records are written to
  # CSV, and avoids allowing machine-level rounding to determine alpha.
  best_pooled_mse <- min(alpha_choices$selected_pooled_mse)
  alpha_tie_tolerance <- 1e-12 * max(1, abs(best_pooled_mse))
  alpha_choices[, alpha_within_tie_tolerance :=
    abs(selected_pooled_mse - best_pooled_mse) <= alpha_tie_tolerance]
  eligible_alpha_choices <- alpha_choices[alpha_within_tie_tolerance == TRUE]
  data.table::setorder(
    eligible_alpha_choices,
    -alpha,
    -selected_lambda_fraction
  )
  chosen <- eligible_alpha_choices[1L]

  list(
    chosen_alpha = chosen$alpha,
    chosen_lambda_fraction = chosen$selected_lambda_fraction,
    chosen_inner_pooled_mse = chosen$selected_pooled_mse,
    chosen_inner_mean_fold_mse = chosen$selected_mean_fold_mse,
    alpha_tie_tolerance = alpha_tie_tolerance,
    alpha_tie_candidate_n = nrow(eligible_alpha_choices),
    fold_losses = fold_losses,
    aggregate_losses = aggregate_losses,
    alpha_choices = alpha_choices,
    scaler_rows = data.table::rbindlist(scaler_rows)
  )
}

fit_outer_elastic_net <- function(x_train, y_train, x_test, alpha,
                                  lambda_fraction, penalty_factor) {
  scaler <- fit_scaler(x_train)
  scaled_training <- apply_scaler(x_train, scaler)
  scaled_test <- apply_scaler(x_test, scaler)
  lambda_max <- get_glmnet_lambda_max(
    scaled_training, y_train, alpha, penalty_factor
  )
  lambda_actual <- lambda_max * lambda_fraction
  fit <- suppressWarnings(glmnet::glmnet(
    x = scaled_training,
    y = y_train,
    family = "gaussian",
    alpha = alpha,
    lambda = lambda_actual,
    penalty.factor = penalty_factor,
    standardize = FALSE,
    intercept = TRUE,
    thresh = 1e-8,
    maxit = 100000L
  ))
  prediction <- as.numeric(stats::predict(
    fit, newx = scaled_test, s = lambda_actual
  ))
  coefficient_matrix <- as.matrix(stats::coef(fit, s = lambda_actual))
  coefficients <- coefficient_matrix[colnames(scaled_training), 1L]
  names(coefficients) <- colnames(scaled_training)
  if (any(!is.finite(prediction)) || any(!is.finite(coefficients))) {
    stop("Non-finite Elastic Net prediction or coefficient in an outer fold.")
  }
  list(
    prediction = prediction,
    coefficients = coefficients,
    scaler = scaler,
    lambda_max = lambda_max,
    lambda_actual = lambda_actual
  )
}
