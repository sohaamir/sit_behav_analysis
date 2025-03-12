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

create_consensus_plot <- function(data, model_results, anova_results, params) {
  interaction_p <- anova_results["direction:consensus_numeric", "Pr(>Chisq)"]
  
  # Get the effect size from the eta_squared values
  interaction_effect <- params$eta_squared %>%
    filter(Effect == "direction:consensus_numeric") %>%
    pull(partial_eta_sq)
  
  # Aggregate data for plotting purposes only (after model fitting)
  plot_data <- data %>%
    mutate(
      consensus_level = case_when(
        consensus_numeric == 0 ~ "2:2",
        consensus_numeric == 1 ~ "3:1",
        consensus_numeric == 2 ~ "4:0"
      ),
      consensus_level = factor(consensus_level, levels = c("2:2", "3:1", "4:0"))
    ) %>%
    group_by(consensus_level, direction) %>%
    summarise(
      mean_outcome = mean(outcome_value),
      se = sd(outcome_value) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Get the 2:2 values from "With group"
  with_2_2 <- plot_data %>%
    filter(direction == "With group", consensus_level == "2:2") %>%
    mutate(direction = "Against group")  # Create a copy with "Against group" direction
  
  # Add this to plot_data
  plot_data <- bind_rows(plot_data, with_2_2)
  
  # Calculate data range
  data_min <- min(plot_data$mean_outcome - plot_data$se)
  data_max <- max(plot_data$mean_outcome + plot_data$se)
  data_range <- data_max - data_min
  
  # Determine the appropriate axis settings based on data characteristics
  if(data_range < 1 && abs(data_max) < 1 && abs(data_min) < 1) {
    # Small value range (likely bet difference data)
    # Round to nearest 0.1 with some padding
    y_min <- floor(data_min * 10) / 10 - 0.1
    y_max <- ceiling(data_max * 10) / 10 + 0.1
    
    # Add more padding if range is very narrow
    if(y_max - y_min < 0.5) {
      y_min <- max(-0.2, y_min - 0.1)
      y_max <- min(0.5, y_max + 0.1)
    }
    
    y_breaks <- seq(y_min, y_max, by = 0.1)
    
  } else if(max(plot_data$mean_outcome) > 20) {
    # Likely accuracy data (values in 0-100 range)
    y_min <- floor(data_min / 10) * 10
    y_max <- ceiling(data_max / 10) * 10
    y_breaks <- seq(y_min, y_max, by = 10)
  } else {
    # For other data types, use increments of 1 or 5
    if(data_range < 10) {
      y_min <- floor(data_min)
      y_max <- ceiling(data_max)
      y_breaks <- seq(y_min, y_max, by = 1)
    } else {
      y_min <- floor(data_min / 5) * 5
      y_max <- ceiling(data_max / 5) * 5
      y_breaks <- seq(y_min, y_max, by = 5)
    }
  }
  
  # Create the base plot
  p <- ggplot(plot_data, 
              aes(x = consensus_level, 
                  y = mean_outcome, 
                  color = direction,
                  group = direction)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_outcome - se, 
                      ymax = mean_outcome + se), 
                  width = 0.2,
                  linewidth = 1) +
    geom_ribbon(aes(ymin = mean_outcome - se, 
                    ymax = mean_outcome + se,
                    fill = direction), 
                alpha = 0.2,
                color = NA) +
    scale_color_manual(values = c("Against group" = "red", "With group" = "blue"),
                       name = NULL) +
    scale_fill_manual(values = c("Against group" = "red", "With group" = "blue"),
                      name = NULL) +
    scale_y_continuous(limits = c(y_min, y_max), breaks = y_breaks) +
    labs(x = "Group consensus",
         y = params$y_label) +
    theme_custom +
    # Position legend in top-left with no box
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.position = c(0.15, 0.95),    # Position in the top-left
      legend.justification = c(0, 1),     # Anchor point at top-left of legend box
      legend.background = element_blank(), # Remove background box
      legend.key = element_blank(),       # Remove box around keys in legend
      legend.text = element_text(size = 14),
      legend.key.size = unit(0.8, "cm"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
    )
  
  # Now add the interaction p-value and effect size on two lines with no box
  # Position depends on the plot type
  text_x_pos <- 1  # Just after the "2:2" position
  
  if(y_max > 20) {
    # For the switch probability plot (larger y range)
    text_y_pos <- y_min + (y_max - y_min) * 0.85  # 85% up from bottom
    
    # Add lines of text
    p <- p + 
      # First line: "Interaction:"
      annotate(
        "text",
        x = text_x_pos,
        y = text_y_pos,
        label = "Interaction:",
        hjust = 0,
        fontface = "bold",
        size = 5
      ) +
      # Second line: "P = value, η²ₚ = value"
      annotate(
        "text",
        x = text_x_pos,
        y = text_y_pos - (y_max - y_min) * 0.03,  # Slightly below first line
        label = sprintf("P = %.3e, η²ₚ = %.3f", interaction_p, interaction_effect),
        hjust = 0,
        fontface = "bold",
        size = 5
      )
  } else {
    # For the bet difference plot (smaller y range)
    text_y_pos <- y_min + (y_max - y_min) * 0.75  # 75% up from bottom
    
    # Add lines of text
    p <- p + 
      # First line: "Interaction:"
      annotate(
        "text",
        x = text_x_pos,
        y = text_y_pos,
        label = "Interaction:",
        hjust = 0,
        fontface = "bold",
        size = 5
      ) +
      # Second line: "P = value, η²ₚ = value"
      annotate(
        "text",
        x = text_x_pos,
        y = text_y_pos - (y_max - y_min) * 0.03,  # Slightly below first line
        label = sprintf("P = %.3f, η²ₚ = %.3f", interaction_p, interaction_effect),
        hjust = 0,
        fontface = "bold",
        size = 5
      )
  }
  
  return(p)
}

create_switch_accuracy_plot <- function(data, model_results, anova_results, params) {
  interaction_p <- as.numeric(anova_results[
    "direction:consensus_numeric:choice_switch_across_trials", 
    "Pr(>Chisq)"])
  
  # Get the effect size for the three-way interaction
  interaction_effect <- params$eta_squared %>%
    filter(Effect == "direction:consensus_numeric:choice_switch_across_trials") %>%
    pull(partial_eta_sq)
  
  # Aggregate data for plotting purposes only
  plot_data <- data %>%
    group_by(consensus_numeric, direction, choice_switch_across_trials) %>%
    summarise(
      value = mean(!!sym(params$outcome_var)),
      se = sd(!!sym(params$outcome_var)) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    mutate(
      consensus_category = case_when(
        consensus_numeric == 2 ~ "4:0",
        consensus_numeric == 1 ~ "3:1",
        consensus_numeric == 0 ~ "2:2"
      ),
      consensus_category = factor(consensus_category, 
                                  levels = c("2:2", "3:1", "4:0")),
      # Create a new variable for trial display with the desired order
      trial_display = ifelse(choice_switch_across_trials == 1, 
                             "Switch trials", "Stay trials"),
      # Set the factor levels to control the order (Switch first, then Stay)
      trial_display = factor(trial_display, 
                             levels = c("Switch trials", "Stay trials"))
    )
  
  if(params$is_accuracy) {
    plot_data <- plot_data %>%
      mutate(
        value = value * 100,
        se = se * 100
      )
  }
  
  # Get the 2:2 values from "Against group" for both switch and stay trials
  against_2_2 <- plot_data %>%
    filter(direction == "Against group", consensus_category == "2:2") %>%
    mutate(direction = "With group")
  
  # Add this to plot_data
  plot_data <- bind_rows(plot_data, against_2_2)
  
  p <- ggplot(plot_data, 
              aes(x = consensus_category, 
                  y = value, 
                  color = direction, 
                  group = direction)) +
    # Use the trial_display variable to control facet order
    facet_wrap(~trial_display) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = value - se,
                      ymax = value + se),
                  width = 0.2) +
    geom_ribbon(aes(ymin = value - se,
                    ymax = value + se,
                    fill = direction),
                alpha = 0.2,
                color = NA) +
    scale_color_manual(values = c("Against group" = "red", 
                                  "With group" = "blue")) +
    scale_fill_manual(values = c("Against group" = "red", 
                                 "With group" = "blue")) +
    labs(
      x = "Group consensus",
      y = params$y_label,
      color = "Direction",
      fill = "Direction",
      title = "Three-way interaction:",
      subtitle = sprintf("P < %.2f, η²ₚ = %.3f", interaction_p, interaction_effect)
    ) +
    theme_custom +
    # Make axes labels and titles larger
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14)
    )
  
  # Set y-axis specifically for accuracy plots to 30-80 in increments of 10
  if(params$is_accuracy) {
    p <- p + scale_y_continuous(limits = c(30, 75), 
                                breaks = seq(30, 75, by = 10))
  }
  
  return(p)
}

create_bet_difference_plot <- function(data, model_results, anova_results, params) {
  interaction_p <- anova_results[
    "direction:consensus_numeric:switch_vs_stay", 
    "Pr(>Chisq)"]
  
  # Get effect size for the interaction
  interaction_effect <- params$eta_squared %>%
    filter(Effect == "direction:consensus_numeric:switch_vs_stay") %>%
    pull(partial_eta_sq)
  
  # Aggregate data for plotting purposes only
  plot_data <- data %>%
    group_by(consensus_numeric, direction, switch_vs_stay) %>%
    summarise(
      mean_diff = mean(bet_difference),
      se = sd(bet_difference) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    mutate(
      consensus_category = case_when(
        consensus_numeric == 2 ~ "4:0",
        consensus_numeric == 1 ~ "3:1",
        consensus_numeric == 0 ~ "2:2"
      ),
      consensus_category = factor(consensus_category, 
                                  levels = c("2:2", "3:1", "4:0")),
      # Create a new variable for display ordering
      trial_display = ifelse(switch_vs_stay == 1, 
                             '"Switch" trials:\nC2(t) ≠ C1(t)', 
                             '"Stay" trials:\nC2(t) = C1(t)'),
      # Set factor levels to control the order (Switch first, then Stay)
      trial_display = factor(trial_display, 
                             levels = c('"Switch" trials:\nC2(t) ≠ C1(t)', 
                                        '"Stay" trials:\nC2(t) = C1(t)'))
    )
  
  # Get the 2:2 values from "Against group" for both switch and stay trials
  against_2_2 <- plot_data %>%
    filter(direction == "Against group", consensus_category == "2:2") %>%
    mutate(direction = "With group")
  
  # Add this to plot_data to extend With group line to include 2:2
  plot_data <- bind_rows(plot_data, against_2_2)
  
  # Calculate appropriate y-axis bounds based on data
  y_min <- floor(min(plot_data$mean_diff - plot_data$se) * 10) / 10
  y_max <- ceiling(max(plot_data$mean_diff + plot_data$se) * 10) / 10
  
  # Add some padding to ensure data points aren't at the edges
  y_min <- max(y_min - 0.1, -0.3)
  y_max <- min(y_max + 0.1, 0.6)
  
  ggplot(plot_data, 
         aes(x = consensus_category, 
             y = mean_diff, 
             color = direction, 
             group = direction)) +
    facet_wrap(~trial_display) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_diff - se,
                      ymax = mean_diff + se),
                  width = 0.2) +
    geom_ribbon(aes(ymin = mean_diff - se,
                    ymax = mean_diff + se,
                    fill = direction),
                alpha = 0.2,
                color = NA) +
    scale_color_manual(values = c("Against group" = "red", 
                                  "With group" = "blue")) +
    scale_fill_manual(values = c("Against group" = "red", 
                                 "With group" = "blue")) +
    scale_y_continuous(limits = c(y_min, y_max), 
                       breaks = seq(y_min, y_max, by = 0.1),
                       labels = function(x) sprintf("%.1f", x)) +
    labs(
      x = "Group consensus",
      y = "Bet difference (Bet 2 - Bet 1)",
      color = "Direction",
      fill = "Direction",
      title = "Direction by Choice:",
      subtitle = sprintf("P = %.2e, η²ₚ = %.3f", interaction_p, interaction_effect)
    ) +
    theme_custom +
    # Make axes text and labels larger
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 14),
      plot.subtitle = element_text(size = 12)
    )
}

create_reversal_plot <- function(data, model_results, anova_results, params) {
  trial_p <- anova_results["trial_to_reversal", "Pr(>Chisq)"]
  type_col <- params$type_col
  
  # Get R² values from the diagnostics
  r2_marginal <- params$diagnostics$model_fit$r2_marginal
  r2_conditional <- params$diagnostics$model_fit$r2_conditional
  
  # Aggregate data for plotting purposes only
  plot_data <- data %>%
    group_by(trial_to_reversal, !!sym(type_col)) %>%
    summarise(
      value = mean(!!sym(params$value_col)) * params$scale_factor,
      se = sd(!!sym(params$value_col)) / sqrt(n()) * params$scale_factor,
      .groups = 'drop'
    )
  
  # Calculate data range
  data_min <- min(plot_data$value - plot_data$se)
  data_max <- max(plot_data$value + plot_data$se)
  data_range <- data_max - data_min
  
  # Set appropriate y-axis bounds based on data type
  if(params$is_accuracy) {
    # For accuracy data
    y_min <- max(20, floor(data_min / 10) * 10)
    y_max <- min(100, ceiling(data_max / 10) * 10)
    y_breaks <- seq(y_min, y_max, by = 10)
  } else if(data_min >= 2 && data_max <= 3 && data_range < 0.5) {
    # For bet magnitude data (narrow range between 2-3)
    y_min <- floor(data_min * 10) / 10
    y_max <- ceiling(data_max * 10) / 10
    
    # Ensure we have a reasonable range even with small variations
    if(y_max - y_min < 0.3) {
      y_min <- floor((data_min - 0.1) * 10) / 10
      y_max <- ceiling((data_max + 0.1) * 10) / 10
    }
    
    y_breaks <- seq(y_min, y_max, by = 0.1)
  } else {
    # For other data types
    y_min <- floor(data_min)
    y_max <- ceiling(data_max)
    y_breaks <- seq(y_min, y_max, by = 1)
  }
  
  # Set dodge width for jittering
  dodge_width <- 0.3
  
  ggplot(plot_data, 
         aes(x = trial_to_reversal, 
             y = value, 
             color = !!sym(type_col), 
             group = !!sym(type_col))) +
    # Add jittering with position_dodge
    geom_line(linewidth = 1, position = position_dodge(width = dodge_width)) +
    geom_point(size = 3, position = position_dodge(width = dodge_width)) +
    geom_errorbar(aes(ymin = value - se,
                      ymax = value + se),
                  width = 0.2,
                  position = position_dodge(width = dodge_width)) +
    scale_color_manual(values = params$colors) +
    scale_x_discrete(labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
    # Use coord_cartesian instead of scale_y_continuous limits
    coord_cartesian(ylim = c(y_min, y_max), expand = FALSE) +
    scale_y_continuous(breaks = y_breaks) +
    labs(
      x = "Trial (relative to reversal)",
      y = params$y_label,
      color = params$legend_title,
      title = "",
      subtitle = sprintf("", 
                         trial_p, r2_marginal, r2_conditional)
    ) +
    theme_custom +
    # Make axes text and labels larger
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 14),
      plot.subtitle = element_text(size = 12)
    )
}