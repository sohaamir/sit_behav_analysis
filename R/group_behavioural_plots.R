# group_behavioural_plots.R

################## LOAD REQUIRED PACKAGES ###################
required_packages <- c(
  "ggplot2", 
  "tidyverse", 
  "ggpubr", 
  "rstatix", 
  "ez", 
  "lme4", 
  "lmerTest", 
  "car", 
  "emmeans", 
  "here", 
  "MuMIn", 
  "conflicted"
)

# Install if missing
install_if_missing <- required_packages[!required_packages %in% installed.packages()]
if (length(install_if_missing) > 0) {
  install.packages(install_if_missing, quietly = TRUE)
}

# Load packages
for (pkg in required_packages) {
  library(pkg, character.only = TRUE)
}

# Set conflicts preferences
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflicts_prefer(lmerTest::lmer)

################## MODEL DIAGNOSTICS ###################
check_model_diagnostics <- function(model) {
  # Calculate basic diagnostics
  residuals <- residuals(model)
  fitted_vals <- fitted(model)
  
  # Handle Shapiro test sample size
  if(length(residuals) > 5000) {
    set.seed(123)
    shapiro_residuals <- sample(residuals, 5000)
  } else if(length(residuals) < 3) {
    shapiro_test <- list(p.value = NA)
  } else {
    shapiro_residuals <- residuals
  }
  
  shapiro_test <- shapiro.test(shapiro_residuals)
  outliers <- boxplot.stats(residuals)$out
  vif_values <- car::vif(model)
  
  # Calculate R² values
  r2_values <- tryCatch({
    MuMIn::r.squaredGLMM(model)
  }, error = function(e) {
    c(NA, NA)
  })
  
  # Check convergence
  convergence <- !is.null(model@optinfo$conv$lme4$messages)
  
  # Check homoscedasticity
  homo_test <- cor.test(abs(residuals), fitted_vals)
  
  list(
    residuals = list(
      mean = mean(residuals),
      sd = sd(residuals),
      normality_p = shapiro_test$p.value,
      n_residuals = length(residuals)
    ),
    model_fit = list(
      r2_marginal = r2_values[1],
      r2_conditional = r2_values[2],
      aic = AIC(model),
      bic = BIC(model)
    ),
    assumptions = list(
      vif = vif_values,
      homoscedasticity_cor = homo_test$estimate,
      homoscedasticity_p = homo_test$p.value,
      convergence = convergence
    ),
    outliers = list(
      n_outliers = length(outliers),
      outlier_values = outliers
    )
  )
}

################## SHARED PLOTTING THEME ###################
theme_custom <- theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

################## HELPER FUNCTIONS ###################
# Function to get shared y-axis limits for consistent plotting
get_shared_y_limits <- function(data_low, data_high, value_col, buffer = 0.05) {
  # Combine both datasets for min/max calculation
  all_values <- c(
    data_low[[value_col]] - data_low$se,
    data_low[[value_col]] + data_low$se,
    data_high[[value_col]] - data_high$se,
    data_high[[value_col]] + data_high$se
  )
  
  min_val <- min(all_values, na.rm = TRUE)
  max_val <- max(all_values, na.rm = TRUE)
  
  # Add buffer
  range <- max_val - min_val
  min_val <- min_val - (range * buffer)
  max_val <- max_val + (range * buffer)
  
  return(list(min = min_val, max = max_val))
}

################## EFFECT SIZE CALCULATIONS ###################
calculate_partial_eta_squared <- function(anova_results) {
  # Create a data frame to store the results
  eta_sq <- data.frame(
    Effect = rownames(anova_results),
    partial_eta_sq = NA
  )
  
  # Calculate partial eta-squared for each effect
  for (i in 1:nrow(anova_results)) {
    # Check if F value exists
    if ("F" %in% colnames(anova_results)) {
      F_val <- anova_results[i, "F"]
      df_num <- anova_results[i, "Df"][1]  # numerator df
      df_denom <- anova_results[i, "Df"][2]  # denominator df
      
      # Calculate partial eta-squared
      eta_sq$partial_eta_sq[i] <- (F_val * df_num) / (F_val * df_num + df_denom)
    } else {
      # For Chi-square tests, use an alternative calculation
      chisq_val <- anova_results[i, "Chisq"]
      df_num <- anova_results[i, "Df"]
      
      # Calculate partial eta-squared (approximation)
      eta_sq$partial_eta_sq[i] <- chisq_val / (chisq_val + nrow(data))
    }
  }
  
  return(eta_sq)
}

# Function to calculate Cohen's d from t-statistic
calculate_cohens_d <- function(t_value, df) {
  return(2 * t_value / sqrt(df))
}

################## RESULTS FORMATTING ###################
format_results <- function(model_results, anova_results, posthoc_results, 
                           plot_data, analysis_name, diagnostics) {
  
  # Calculate partial eta-squared
  eta_sq <- calculate_partial_eta_squared(anova_results)
  
  output_text <- sprintf("%s - STATISTICAL ANALYSIS RESULTS\n\n", 
                         toupper(analysis_name))
  output_text <- paste0(output_text, "Analysis run on: ", 
                        format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # Model comparison results
  output_text <- paste0(output_text, "MODEL COMPARISON RESULTS:\n",
                        "========================\n")
  model_comp_df <- data.frame(
    Model = paste0("m", 1:length(model_results$aic_values)),
    AIC = round(model_results$aic_values, 3),
    BIC = round(model_results$bic_values, 3),
    Is_Winner = ifelse(1:length(model_results$aic_values) == 
                         model_results$winning_idx, "Yes", "No")
  )
  output_text <- paste0(output_text, 
                        paste(capture.output(print.data.frame(model_comp_df)), 
                              collapse = "\n"), "\n\n")
  
  # Diagnostics section
  output_text <- paste0(output_text, "MODEL DIAGNOSTICS:\n",
                        "=================\n",
                        "Residuals:\n",
                        "  Mean: ", round(diagnostics$residuals$mean, 3), "\n",
                        "  SD: ", round(diagnostics$residuals$sd, 3), "\n",
                        "  Normality test p-value: ", 
                        format.pval(diagnostics$residuals$normality_p, 
                                    digits = 3), "\n",
                        "  Number of residuals: ", 
                        diagnostics$residuals$n_residuals, "\n\n",
                        "Model Fit:\n",
                        "  R² marginal: ", 
                        round(diagnostics$model_fit$r2_marginal, 3), "\n",
                        "  R² conditional: ", 
                        round(diagnostics$model_fit$r2_conditional, 3), "\n",
                        "  AIC: ", round(diagnostics$model_fit$aic, 2), "\n",
                        "  BIC: ", round(diagnostics$model_fit$bic, 2), "\n\n",
                        "Assumptions:\n",
                        "  VIF range: ", round(min(diagnostics$assumptions$vif), 2), 
                        " to ", round(max(diagnostics$assumptions$vif), 2), "\n",
                        "  Homoscedasticity (residual-fitted correlation): ", 
                        round(diagnostics$assumptions$homoscedasticity_cor, 3),
                        " (p = ", format.pval(diagnostics$assumptions$homoscedasticity_p, 
                                              digits = 3), ")\n",
                        "  Model convergence achieved: ", 
                        !diagnostics$assumptions$convergence, "\n\n",
                        "Outliers:\n",
                        "  Number of outliers: ", 
                        diagnostics$outliers$n_outliers, "\n\n")
  
  # Effect sizes section
  output_text <- paste0(output_text, "EFFECT SIZES:\n",
                        "===========\n",
                        "R² marginal: ", round(diagnostics$model_fit$r2_marginal, 3), 
                        " (variance explained by fixed effects)\n",
                        "R² conditional: ", round(diagnostics$model_fit$r2_conditional, 3), 
                        " (variance explained by both fixed and random effects)\n\n")
  
  # Add partial eta-squared values for complex interactions
  output_text <- paste0(output_text, "PARTIAL ETA-SQUARED VALUES:\n",
                        "==========================\n")
  
  # Format partial eta-squared table
  eta_sq_table <- data.frame(
    Effect = eta_sq$Effect,
    Partial_Eta_Sq = round(eta_sq$partial_eta_sq, 3)
  )
  
  output_text <- paste0(output_text,
                        paste(capture.output(print.data.frame(eta_sq_table)),
                              collapse = "\n"), "\n\n")
  
  # ANOVA results
  output_text <- paste0(output_text, "ANOVA RESULTS:\n",
                        "==============\n",
                        paste(capture.output(anova_results), 
                              collapse = "\n"), "\n\n")
  
  # Post-hoc tests
  output_text <- paste0(output_text, "POST-HOC TESTS:\n",
                        "==============\n",
                        paste(capture.output(posthoc_results), 
                              collapse = "\n"), "\n\n")
  
  return(output_text)
}


################## CORE PIPELINE ###################
run_mixed_model_pipeline <- function(
    data,
    preprocessing_function,
    model_formulas,
    plot_function,
    analysis_name,
    output_path,
    additional_params = list()
) {
  # Process data
  processed_data <- preprocessing_function(data)
  
  # Run all models
  models <- lapply(model_formulas, function(formula) {
    lmer(formula, data = processed_data, 
         control = lmerControl(optimizer = "bobyqa"))
  })
  
  # Compare models and select winner
  aic_values <- sapply(models, AIC)
  bic_values <- sapply(models, BIC)
  winning_idx <- which.min(aic_values)
  winning_model <- models[[winning_idx]]
  
  # Compile model results
  model_results <- list(
    models = models,
    winning_model = winning_model,
    aic_values = aic_values,
    bic_values = bic_values,
    winning_idx = winning_idx
  )
  
  # Run diagnostics
  diagnostics <- check_model_diagnostics(winning_model)
  
  # Get statistical results
  anova_results <- car::Anova(winning_model, type = 2)
  
  # Calculate partial eta-squared values
  eta_squared_values <- calculate_partial_eta_squared(anova_results)
  
  # Run post-hoc tests
  emmeans_results <- emmeans(winning_model, 
                             specs = additional_params$emmeans_specs)
  posthoc_results <- pairs(emmeans_results, adjust = "tukey")
  
  # Create a properly merged list of parameters
  plot_params <- additional_params
  plot_params$eta_squared <- eta_squared_values
  plot_params$diagnostics <- diagnostics
  
  # Create plot with all parameters
  plot <- plot_function(processed_data, model_results, 
                        anova_results, plot_params)
  
  # Format and save results
  results_text <- format_results(
    model_results,
    anova_results,
    posthoc_results,
    processed_data,
    analysis_name,
    diagnostics
  )
  
  # Create output directory if needed
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(results_text, output_path)
  
  # Print plot
  print(plot)
  
  # Return all results
  return(list(
    model_results = model_results,
    anova_results = anova_results,
    posthoc_results = posthoc_results,
    plot = plot,
    processed_data = processed_data,
    diagnostics = diagnostics,
    eta_squared = eta_squared_values
  ))
}