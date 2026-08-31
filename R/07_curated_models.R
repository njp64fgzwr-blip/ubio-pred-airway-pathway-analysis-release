# Fit prespecified multivariable pathway models and sensitivity analyses.

source(file.path("R", "00_common.R"))
assert_packages(c("data.table", "MASS", "car", "broom", "lmtest"))
message_rule("07: Fitting curated pathway models and sensitivities")

coefficient_rows <- list()
performance_rows <- list()
formula_rows <- list()
vif_rows <- list()
diagnostic_rows <- list()
adjusted_rows <- list()
adjusted_performance_rows <- list()
nb_rows <- list()
nb_performance_rows <- list()

standardize_numeric <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) stop("Cannot standardize a zero-variance variable.")
  as.numeric((x - mean(x)) / s)
}

extract_lm_terms <- function(fit, pathways, dataset_name, outcome_name,
                             model_name, n) {
  coefs <- summary(fit)$coefficients
  df_residual <- stats::df.residual(fit)
  critical <- stats::qt(0.975, df = df_residual)
  missing <- setdiff(pathways, rownames(coefs))
  if (length(missing)) {
    stop(dataset_name, " / ", outcome_name, ": pathway terms missing: ",
         paste(missing, collapse = ", "))
  }
  data.table::data.table(
    dataset = dataset_name,
    dataset_label = datasets[[dataset_name]]$label,
    compartment = datasets[[dataset_name]]$compartment,
    platform = datasets[[dataset_name]]$platform,
    omic_type = datasets[[dataset_name]]$omic_type,
    outcome = outcome_name,
    outcome_label = outcome_specs[[outcome_name]]$label,
    model = model_name,
    pathway = pathways,
    pathway_label = unname(pathway_labels[pathways]),
    n = n,
    residual_df = df_residual,
    standardized_beta = coefs[pathways, "Estimate"],
    standard_error = coefs[pathways, "Std. Error"],
    ci_low = coefs[pathways, "Estimate"] - critical * coefs[pathways, "Std. Error"],
    ci_high = coefs[pathways, "Estimate"] + critical * coefs[pathways, "Std. Error"],
    t_statistic = coefs[pathways, "t value"],
    p_value = coefs[pathways, "Pr(>|t|)"],
    confidence_interval_method = paste0(
      "95% residual-df t interval (df=", df_residual, ")"
    )
  )
}

calculate_in_sample_metrics <- function(observed, predicted) {
  base <- calculate_metrics(observed, predicted)
  c(
    r_squared = unname(base["r_squared"]),
    rmse = unname(base["rmse"]),
    mae = unname(base["mae"]),
    pearson_r = suppressWarnings(stats::cor(observed, predicted, method = "pearson"))
  )
}

for (dataset_name in names(datasets)) {
  merged <- readRDS(file.path(
    paths$analysis, paste0(dataset_name, "_pathway_clinical_analysis.rds")
  ))

  for (outcome_name in names(outcome_specs)) {
    outcome <- outcome_specs[[outcome_name]]
    pathways <- curated_pathway_models[[outcome_name]]
    variables <- c(outcome$column, pathways)
    complete <- stats::complete.cases(merged[, ..variables])
    model_source <- data.table::copy(merged[complete, ..variables])
    n <- nrow(model_source)
    if (n <= length(pathways) + 5L) {
      stop(dataset_name, " / ", outcome_name, ": insufficient primary complete cases.")
    }

    y_model <- apply_outcome_transform(model_source[[outcome$column]], outcome$transform)
    primary_data <- data.frame(outcome_z = standardize_numeric(y_model))
    for (pathway in pathways) {
      primary_data[[pathway]] <- standardize_numeric(model_source[[pathway]])
    }
    primary_formula <- stats::reformulate(pathways, response = "outcome_z")
    primary_fit <- stats::lm(primary_formula, data = primary_data)
    primary_terms <- extract_lm_terms(
      primary_fit, pathways, dataset_name, outcome_name,
      "Primary curated pathways", n
    )
    coefficient_rows[[paste(dataset_name, outcome_name, sep = "__")]] <- primary_terms

    fitted <- stats::fitted(primary_fit)
    metrics <- calculate_in_sample_metrics(primary_data$outcome_z, fitted)
    fit_summary <- summary(primary_fit)
    performance_rows[[paste(dataset_name, outcome_name, sep = "__")]] <-
      data.table::data.table(
        dataset = dataset_name,
        dataset_label = datasets[[dataset_name]]$label,
        compartment = datasets[[dataset_name]]$compartment,
        platform = datasets[[dataset_name]]$platform,
        omic_type = datasets[[dataset_name]]$omic_type,
        outcome = outcome_name,
        outcome_label = outcome$label,
        model = "Primary curated pathways",
        n = n,
        predictors = length(pathways),
        outcome_transform = outcome$transform,
        coefficient_scale = "Model-standardized outcome and pathway predictors",
        r_squared = fit_summary$r.squared,
        adjusted_r_squared = fit_summary$adj.r.squared,
        rmse_standardized_outcome = unname(metrics["rmse"]),
        mae_standardized_outcome = unname(metrics["mae"]),
        pearson_r = unname(metrics["pearson_r"]),
        model_f_statistic = unname(fit_summary$fstatistic["value"]),
        model_df1 = unname(fit_summary$fstatistic["numdf"]),
        model_df2 = unname(fit_summary$fstatistic["dendf"]),
        model_p_value = stats::pf(
          fit_summary$fstatistic["value"],
          fit_summary$fstatistic["numdf"],
          fit_summary$fstatistic["dendf"],
          lower.tail = FALSE
        )
      )

    formula_rows[[paste(dataset_name, outcome_name, "primary", sep = "__")]] <-
      data.table::data.table(
        dataset = dataset_name,
        dataset_label = datasets[[dataset_name]]$label,
        compartment = datasets[[dataset_name]]$compartment,
        platform = datasets[[dataset_name]]$platform,
        omic_type = datasets[[dataset_name]]$omic_type,
        outcome = outcome_name,
        outcome_label = outcome$label,
        analysis = "Primary curated OLS",
        n = n,
        outcome_transform = outcome$transform,
        formula = paste(deparse(primary_formula), collapse = " "),
        pathway_subset = paste(pathways, collapse = ";"),
        covariates = "None"
      )

    # VIF is calculated for the primary pathway terms only.
    vif_values <- car::vif(primary_fit)
    if (is.matrix(vif_values)) {
      vif_table <- data.table::as.data.table(vif_values, keep.rownames = "term")
      data.table::setnames(vif_table, names(vif_table),
                           c("term", "gvif", "df", "gvif_adjusted"))
    } else {
      vif_table <- data.table::data.table(
        term = names(vif_values),
        gvif = as.numeric(vif_values),
        df = 1,
        gvif_adjusted = sqrt(as.numeric(vif_values))
      )
    }
    vif_table[, `:=`(
      dataset = dataset_name,
      outcome = outcome_name,
      model = "Primary curated pathways"
    )]
    vif_rows[[paste(dataset_name, outcome_name, sep = "__")]] <- vif_table

    residuals <- stats::residuals(primary_fit)
    cooks <- stats::cooks.distance(primary_fit)
    leverage <- stats::hatvalues(primary_fit)
    bp <- tryCatch(lmtest::bptest(primary_fit), error = function(e) NULL)
    diagnostic_rows[[paste(dataset_name, outcome_name, sep = "__")]] <-
      data.table::data.table(
        dataset = dataset_name,
        outcome = outcome_name,
        model = "Primary curated pathways",
        n = n,
        residual_mean = mean(residuals),
        residual_sd = stats::sd(residuals),
        residual_minimum = min(residuals),
        residual_maximum = max(residuals),
        max_cooks_distance = max(cooks),
        influential_cooks_over_4_n = sum(cooks > 4 / n),
        max_leverage = max(leverage),
        breusch_pagan_statistic = if (is.null(bp)) NA_real_ else unname(bp$statistic),
        breusch_pagan_p_value = if (is.null(bp)) NA_real_ else bp$p.value
      )

    saveRDS(
      list(
        fit = primary_fit,
        source_subject_ids = merged$Subject_ID[complete],
        outcome_transform = outcome$transform,
        pathways = pathways
      ),
      file.path(paths$models, paste0(dataset_name, "_", outcome_name,
                                     "_curated_primary_model.rds")),
      compress = "gzip"
    )

    # Covariate-adjusted sensitivity analysis.
    covariates <- c("age", "sex", "smoking_status", "bmi", "ocs_current")
    adjusted_variables <- c(outcome$column, pathways, covariates)
    adjusted_complete <- stats::complete.cases(merged[, ..adjusted_variables])
    adjusted_source <- data.table::copy(merged[adjusted_complete, ..adjusted_variables])
    adjusted_n <- nrow(adjusted_source)
    adjusted_y <- apply_outcome_transform(
      adjusted_source[[outcome$column]], outcome$transform
    )
    adjusted_data <- data.frame(
      outcome_z = standardize_numeric(adjusted_y),
      age_z = standardize_numeric(adjusted_source$age),
      sex = stats::relevel(factor(adjusted_source$sex), ref = "male"),
      smoking_status = stats::relevel(
        factor(adjusted_source$smoking_status), ref = "non_smoker"
      ),
      bmi_z = standardize_numeric(adjusted_source$bmi),
      ocs_current = stats::relevel(factor(adjusted_source$ocs_current), ref = "No")
    )
    for (pathway in pathways) {
      adjusted_data[[pathway]] <- standardize_numeric(adjusted_source[[pathway]])
    }
    adjusted_formula <- stats::reformulate(
      c(pathways, "age_z", "sex", "smoking_status", "bmi_z", "ocs_current"),
      response = "outcome_z"
    )
    adjusted_fit <- stats::lm(adjusted_formula, data = adjusted_data)
    adjusted_terms <- extract_lm_terms(
      adjusted_fit, pathways, dataset_name, outcome_name,
      "Covariate-adjusted sensitivity", adjusted_n
    )
    adjusted_rows[[paste(dataset_name, outcome_name, sep = "__")]] <- adjusted_terms
    adjusted_summary <- summary(adjusted_fit)
    adjusted_metrics <- calculate_in_sample_metrics(
      adjusted_data$outcome_z, stats::fitted(adjusted_fit)
    )
    adjusted_performance_rows[[paste(dataset_name, outcome_name, sep = "__")]] <-
      data.table::data.table(
        dataset = dataset_name,
        dataset_label = datasets[[dataset_name]]$label,
        compartment = datasets[[dataset_name]]$compartment,
        platform = datasets[[dataset_name]]$platform,
        omic_type = datasets[[dataset_name]]$omic_type,
        outcome = outcome_name,
        outcome_label = outcome$label,
        model = "Covariate-adjusted sensitivity",
        n = adjusted_n,
        pathway_predictors = length(pathways),
        covariates = paste(covariates, collapse = ";"),
        r_squared = adjusted_summary$r.squared,
        adjusted_r_squared = adjusted_summary$adj.r.squared,
        rmse_standardized_outcome = unname(adjusted_metrics["rmse"]),
        mae_standardized_outcome = unname(adjusted_metrics["mae"]),
        pearson_r = unname(adjusted_metrics["pearson_r"])
      )
    formula_rows[[paste(dataset_name, outcome_name, "adjusted", sep = "__")]] <-
      data.table::data.table(
        dataset = dataset_name,
        dataset_label = datasets[[dataset_name]]$label,
        compartment = datasets[[dataset_name]]$compartment,
        platform = datasets[[dataset_name]]$platform,
        omic_type = datasets[[dataset_name]]$omic_type,
        outcome = outcome_name,
        outcome_label = outcome$label,
        analysis = "Covariate-adjusted OLS sensitivity",
        n = adjusted_n,
        outcome_transform = outcome$transform,
        formula = paste(deparse(adjusted_formula), collapse = " "),
        pathway_subset = paste(pathways, collapse = ";"),
        covariates = paste(covariates, collapse = ";")
      )
    saveRDS(
      list(
        fit = adjusted_fit,
        source_subject_ids = merged$Subject_ID[adjusted_complete],
        outcome_transform = outcome$transform,
        pathways = pathways,
        covariates = covariates
      ),
      file.path(paths$models, paste0(dataset_name, "_", outcome_name,
                                     "_curated_adjusted_sensitivity.rds")),
      compress = "gzip"
    )

    # Count-model sensitivity for prior-year exacerbation frequency.
    if (outcome_name == "exacerbation_frequency") {
      nb_source <- model_source
      nb_data <- data.frame(outcome_count = nb_source[[outcome$column]])
      if (any(abs(nb_data$outcome_count - round(nb_data$outcome_count)) > 1e-8)) {
        stop(dataset_name, ": exacerbation outcome is not integer-valued.")
      }
      for (pathway in pathways) {
        nb_data[[pathway]] <- standardize_numeric(nb_source[[pathway]])
      }
      nb_formula <- stats::reformulate(pathways, response = "outcome_count")
      nb_fit <- MASS::glm.nb(nb_formula, data = nb_data, link = log)
      nb_coef <- summary(nb_fit)$coefficients
      critical <- stats::qnorm(0.975)
      nb_terms <- data.table::data.table(
        dataset = dataset_name,
        dataset_label = datasets[[dataset_name]]$label,
        compartment = datasets[[dataset_name]]$compartment,
        platform = datasets[[dataset_name]]$platform,
        omic_type = datasets[[dataset_name]]$omic_type,
        outcome = outcome_name,
        outcome_label = outcome$label,
        model = "Negative-binomial count sensitivity",
        pathway = pathways,
        pathway_label = unname(pathway_labels[pathways]),
        n = nrow(nb_data),
        log_rate_coefficient = nb_coef[pathways, "Estimate"],
        standard_error = nb_coef[pathways, "Std. Error"],
        z_statistic = nb_coef[pathways, "z value"],
        p_value = nb_coef[pathways, "Pr(>|z|)"],
        incidence_rate_ratio = exp(nb_coef[pathways, "Estimate"]),
        irr_ci_low = exp(nb_coef[pathways, "Estimate"] -
                           critical * nb_coef[pathways, "Std. Error"]),
        irr_ci_high = exp(nb_coef[pathways, "Estimate"] +
                            critical * nb_coef[pathways, "Std. Error"]),
        interval_method = "95% asymptotic normal interval"
      )
      nb_rows[[dataset_name]] <- nb_terms
      nb_performance_rows[[dataset_name]] <- data.table::data.table(
        dataset = dataset_name,
        outcome = outcome_name,
        n = nrow(nb_data),
        theta = nb_fit$theta,
        theta_standard_error = nb_fit$SE.theta,
        log_likelihood = as.numeric(stats::logLik(nb_fit)),
        aic = stats::AIC(nb_fit),
        null_deviance = nb_fit$null.deviance,
        residual_deviance = nb_fit$deviance
      )
      formula_rows[[paste(dataset_name, outcome_name, "negative_binomial", sep = "__")]] <-
        data.table::data.table(
          dataset = dataset_name,
          dataset_label = datasets[[dataset_name]]$label,
          compartment = datasets[[dataset_name]]$compartment,
          platform = datasets[[dataset_name]]$platform,
          omic_type = datasets[[dataset_name]]$omic_type,
          outcome = outcome_name,
          outcome_label = outcome$label,
          analysis = "Negative-binomial count sensitivity",
          n = nrow(nb_data),
          outcome_transform = "Untransformed count with log link",
          formula = paste(deparse(nb_formula), collapse = " "),
          pathway_subset = paste(pathways, collapse = ";"),
          covariates = "None"
        )
      saveRDS(
        list(
          fit = nb_fit,
          source_subject_ids = merged$Subject_ID[complete],
          pathways = pathways
        ),
        file.path(paths$models, paste0(dataset_name,
                                       "_exacerbation_negative_binomial_sensitivity.rds")),
        compress = "gzip"
      )
    }
  }
}

coefficients <- data.table::rbindlist(coefficient_rows)
coefficients[, fdr_bh := stats::p.adjust(p_value, method = "BH"),
             by = .(dataset, outcome)]
coefficients[, absolute_standardized_beta := abs(standardized_beta)]
coefficients[, relative_contribution_percent :=
               100 * absolute_standardized_beta / sum(absolute_standardized_beta),
             by = .(dataset, outcome)]
coefficients[, relative_contribution_caution :=
  "Absolute standardized coefficient share; not variance decomposition"]

adjusted_coefficients <- data.table::rbindlist(adjusted_rows)
adjusted_coefficients[, fdr_bh := stats::p.adjust(p_value, method = "BH"),
                      by = .(dataset, outcome)]

negative_binomial <- data.table::rbindlist(nb_rows)
negative_binomial[, fdr_bh := stats::p.adjust(p_value, method = "BH"),
                  by = .(dataset, outcome)]

data.table::fwrite(
  coefficients,
  file.path(paths$tables, "Table_3_curated_model_coefficients_95CI.csv")
)
data.table::fwrite(
  data.table::rbindlist(performance_rows),
  file.path(paths$tables, "Table_S_curated_model_in_sample_performance.csv")
)
data.table::fwrite(
  adjusted_coefficients,
  file.path(paths$tables, "Table_S_covariate_adjusted_curated_coefficients.csv")
)
data.table::fwrite(
  data.table::rbindlist(adjusted_performance_rows),
  file.path(paths$tables, "Table_S_covariate_adjusted_model_performance.csv")
)
data.table::fwrite(
  data.table::rbindlist(formula_rows, fill = TRUE),
  file.path(paths$tables, "Table_S_curated_model_formulas.csv")
)
data.table::fwrite(
  data.table::rbindlist(vif_rows, fill = TRUE),
  file.path(paths$tables, "Table_S_curated_model_VIF.csv")
)
data.table::fwrite(
  data.table::rbindlist(diagnostic_rows, fill = TRUE),
  file.path(paths$tables, "Table_S_curated_model_diagnostics.csv")
)
data.table::fwrite(
  negative_binomial,
  file.path(paths$tables, "Table_S_exacerbation_negative_binomial_coefficients.csv")
)
data.table::fwrite(
  data.table::rbindlist(nb_performance_rows),
  file.path(paths$tables, "Table_S_exacerbation_negative_binomial_performance.csv")
)

if (nrow(coefficients) != 104L) {
  stop("Expected 104 primary curated pathway coefficients; found ", nrow(coefficients), ".")
}
if (nrow(adjusted_coefficients) != 104L) {
  stop("Expected 104 adjusted pathway coefficients; found ", nrow(adjusted_coefficients), ".")
}

writeLines(
  c(
    "CURATED MODELS: PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste("Primary prespecified pathway coefficients:", nrow(coefficients)),
    paste("Primary coefficients with BH q<0.05:", sum(coefficients$fdr_bh < 0.05)),
    paste("Covariate-adjusted pathway coefficients:", nrow(adjusted_coefficients)),
    paste("Adjusted coefficients with BH q<0.05:",
          sum(adjusted_coefficients$fdr_bh < 0.05)),
    "Primary OLS coefficients use standardized transformed outcomes and standardized pathway predictors.",
    "Confidence intervals use the residual-df t distribution.",
    "Only prespecified pathway subsets were fitted; these are not selected top pathways.",
    "Adjusted sensitivity models add age, sex, smoking status, BMI and current/ongoing OCS exposure.",
    "Negative-binomial sensitivity models use the untransformed prior-year event count.",
    "A negative adjusted IL-5 term can reflect collinearity/suppression and is not biological protection.",
    "All estimates are associations, not causal effects."
  ),
  file.path(paths$validation, "CURATED_MODELS_STATUS.txt")
)
message("Curated models and sensitivities completed successfully.")
