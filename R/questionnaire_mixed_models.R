# questionnaire_mixed_models.R

# Add this at the start of your script
options(immediate.write = TRUE)

################## SETUP AND CONSTANTS ###################
required_packages <- c(
  "ggplot2", "tidyverse", "ggpubr", "rstatix", "ez", 
  "lme4", "lmerTest", "car", "emmeans", "MuMIn", 
  "psych", "interactions", "effects", "here", 
  "lm.beta", "effectsize", "corrplot", "reshape2"
)

# Install if missing
install_if_missing <- required_packages[!required_packages %in% installed.packages()]
if (length(install_if_missing) > 0) {
  install.packages(install_if_missing, quietly = TRUE)
}

# Load all packages
for (pkg in required_packages) {
  library(pkg, character.only = TRUE)
}

# Set conflicts preferences
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflicts_prefer(lmerTest::lmer)
conflicts_prefer(effectsize::eta_squared)

# Analysis types for directory structure
# Valid analysis types for validation
VALID_ANALYSIS_TYPES <- c(
  "choice_consensus",    # Original choice/bet by consensus analysis
  "switch_difference"    # Switch difference analysis
)

# Plot types for directory structure
PLOT_TYPES <- c(
  "continuous",      # For continuous relationship plots
  "median_split",    # For median split plots
  "moderation",      # For moderation effect plots
  "simple_slopes"    # For simple slopes analysis plots
)

# Add this function to handle plot saving
save_analysis_plots <- function(results, analysis_type, base_dir) {
  for(scale in names(results)) {
    scale_results <- results[[scale]]
    
    # Map plot names to directories
    plot_dir_mapping <- list(
      continuous = "continuous",
      median_split = "median_split",
      moderation = "moderation",
      simple_slopes = "simple_slopes"
    )
    
    # Save main scale plots
    for(plot_name in names(scale_results$plots)) {
      if(!is.null(scale_results$plots[[plot_name]])) {
        # Determine which subfolder to use
        plot_type <- NA
        for(type in names(plot_dir_mapping)) {
          if(grepl(type, plot_name)) {
            plot_type <- plot_dir_mapping[[type]]
            break
          }
        }
        
        if(!is.na(plot_type)) {
          # Construct the full path using the new directory structure
          output_dir <- file.path(base_dir, "plots", analysis_type, plot_type)
          file_name <- paste0(scale, "_", plot_name, ".png")
          
          # Create directory if it doesn't exist
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
          
          ggsave(
            filename = file.path(output_dir, file_name),
            plot = scale_results$plots[[plot_name]],
            width = 10,
            height = 7
          )
        }
      }
    }
  }
}

################## CORE MODEL PIPELINE ###################
run_model_pipeline <- function(data, params, output_path) {
  # Print trial counts once at the start if this is switch_difference analysis
  if(params$analysis_type == "switch_difference") {
    cat("\n===== Trial Type Summary =====\n")
    stay_count <- sum(data$choice_switch_across_trials == 0, na.rm = TRUE)
    switch_count <- sum(data$choice_switch_across_trials == 1, na.rm = TRUE)
    cat(sprintf("Stay trials: %d\n", stay_count))
    cat(sprintf("Switch trials: %d\n", switch_count))
    cat("============================\n")
  }
  
  # First: Run pooled model selection
  cat("\n===== Running Model Selection =====\n")
  pooled_results <- run_pooled_model_selection(data, params)
  
  cat("\n===== Participant Summary =====\n")
  cat(sprintf("Total participants: %d\n", n_distinct(data$participant.id_in_session)))
  cat("============================\n")
  
  all_results <- list()
  interaction_p_values <- list()
  
  # Apply winning model structure to each questionnaire
  message("\nAnalyzing questionnaires...")
  for(var in params$questionnaire_vars) {
    message(sprintf("\nProcessing: %s", var))
    
    # Process data
    processed_data <- params$preprocessing_fn(data, var)
    
    # Run analysis with winning model structure
    analysis_results <- run_individual_analysis(
      processed_data,
      pooled_results$winning_model,
      var,
      params
    )
    
    # Store raw p-values for three-way interactions
    interaction_p_values[[paste0(var, "_main")]] <- analysis_results$main_results$raw_p
    
    # Store subscale p-values if they exist
    if(!is.null(analysis_results$subscale_results)) {
      for(subscale in names(analysis_results$subscale_results)) {
        interaction_p_values[[paste0(var, "_", subscale)]] <- 
          analysis_results$subscale_results[[subscale]]$raw_p
      }
    }
    
    # Store results without adjusted p-values for now
    all_results[[var]] <- list(
      model_results = pooled_results,
      analysis_results = analysis_results$main_results,
      moderation_results = analysis_results,
      processed_data = processed_data
    )
  }
  
  # Perform FDR correction across all interaction p-values
  message("\nPerforming FDR correction...")
  all_p_values <- unlist(interaction_p_values)
  all_adj_p_values <- p.adjust(all_p_values, method = "fdr")
  names(all_adj_p_values) <- names(all_p_values)
  
  # Update results with corrected p-values and generate plots
  for(var in params$questionnaire_vars) {
    # Update main scale adjusted p-value
    all_results[[var]]$moderation_results$main_results$adj_p <- 
      all_adj_p_values[paste0(var, "_main")]
    
    # Update subscale adjusted p-values if they exist
    if(!is.null(all_results[[var]]$moderation_results$subscale_results)) {
      for(subscale in names(all_results[[var]]$moderation_results$subscale_results)) {
        all_results[[var]]$moderation_results$subscale_results[[subscale]]$adj_p <- 
          all_adj_p_values[paste0(var, "_", subscale)]
      }
    }
    
    # Generate plots
    all_results[[var]]$plots <- generate_all_plots(
      all_results[[var]]$processed_data,
      all_results[[var]]$moderation_results,
      var,
      params
    )
  }
  
  # Save results
  save_results(all_results, params$analysis_name, output_path, params)
  
  return(all_results)
}

################## POOLED MODEL SELECTION ###################
run_pooled_model_selection <- function(data, params) {
  # Create pooled dataset for initial model selection
  pooled_data <- map_dfr(params$questionnaire_vars, function(var) {
    processed <- params$preprocessing_fn(data, var)
    
    # Base columns to select
    base_cols <- c("participant.id_in_session", "consensus_level", 
                   "direction", "outcome_value", "age", "gender", 
                   "scale_name")
    
    # Add switch_difference if it's a switch_difference analysis
    if(params$analysis_type == "switch_difference") {
      base_cols <- c(base_cols, "switch_difference")
    }
    
    processed %>%
      select(all_of(base_cols)) %>%
      distinct()
  })
  
  # Run initial model comparison on pooled data
  message("Running model comparison on pooled data...")
  model_results <- run_model_comparison(params$model_formulas, pooled_data, params)
  
  return(model_results)
}

################## MODEL FUNCTIONS ###################
run_model_comparison <- function(model_formulas, data, params) {
  message("\nComparing models:")
  models <- vector("list", length(model_formulas))
  names(models) <- names(model_formulas)
  
  for(model_name in names(model_formulas)) {
    models[[model_name]] <- tryCatch({
      lmer(
        formula = model_formulas[[model_name]]$formula, 
        data = data,
        control = lmerControl(optimizer = "bobyqa",
                              optCtrl = list(maxfun = 100000))
      )
    }, error = function(e) {
      message(sprintf("Failed to fit %s: %s", model_name, e$message))
      NULL
    })
  }
  
  # Remove NULL models 
  models <- models[!sapply(models, is.null)]
  if(length(models) == 0) stop("No models successfully fitted")
  
  # Compare models
  aic_values <- sapply(models, AIC)
  bic_values <- sapply(models, BIC)
  winning_idx <- which.min(aic_values)
  
  # Print concise comparison
  comparison_df <- data.frame(
    Model = names(models),
    Description = sapply(model_formulas[names(models)], function(x) x$description),
    AIC = round(aic_values, 2),
    BIC = round(bic_values, 2),
    Is_Winner = seq_along(models) == winning_idx,
    stringsAsFactors = FALSE
  )
  
  message("\nModel comparison results:")
  print(comparison_df)
  
  message("\nWinning model formula:")
  print(formula(models[[winning_idx]]))
  
  return(list(
    models = models,
    winning_model = models[[winning_idx]],
    aic_values = aic_values,
    bic_values = bic_values,
    winning_idx = winning_idx,
    model_names = names(models),
    descriptions = sapply(model_formulas[names(models)], function(x) x$description)
  ))
}

################## INDIVIDUAL ANALYSIS ###################
run_individual_analysis <- function(processed_data, winning_model, var, params) {
  # Initialize results lists
  main_results <- list()
  subscale_results <- list()
  
  # Fit main scale model
  main_model <- tryCatch({
    model <- lmer(formula(winning_model), data = processed_data,
                  control = lmerControl(optimizer = "bobyqa"))
    model
  }, error = function(e) {
    message(sprintf("Error fitting model for %s: %s", var, e$message))
    return(NULL)
  })
  
  if(is.null(main_model)) {
    stop(paste("Failed to fit main model for variable:", var))
  }
  
  # Get ANOVA results for main model
  main_anova <- car::Anova(main_model, type = 2)
  
  # Extract p-value for three-way interaction
  interaction_term <- "consensus_level:direction:scale_name"
  raw_p <- main_anova[interaction_term, "Pr(>Chisq)"]
  
  # Only print if significant
  if(raw_p < 0.05) {
    message(sprintf("  Significant interaction found (p = %.3f)", raw_p))
  }
  
  # Create main results list
  main_results <- list(
    model = main_model,
    anova_results = main_anova,
    raw_p = raw_p,
    diagnostics = check_model_diagnostics(main_model, params)
  )
  
  # Run subscale analyses if applicable
  if(!is.null(params$subscale_mapping[[var]])) {
    message(sprintf("  Processing subscales for %s", var))
    
    for(subscale in params$subscale_mapping[[var]]) {
      if(subscale %in% colnames(processed_data)) {
        # Create subscale data
        subscale_data <- processed_data
        subscale_data$scale_name <- as.numeric(scale(processed_data[[subscale]]))
        
        # Fit subscale model
        subscale_model <- tryCatch({
          model <- lmer(formula(winning_model), data = subscale_data,
                        control = lmerControl(optimizer = "bobyqa"))
          model
        }, error = function(e) {
          message(sprintf("    Error fitting model for subscale %s: %s", 
                          subscale, e$message))
          return(NULL)
        })
        
        if(!is.null(subscale_model)) {
          # Get ANOVA results for subscale
          sub_anova <- car::Anova(subscale_model, type = 2)
          sub_raw_p <- sub_anova[interaction_term, "Pr(>Chisq)"]
          
          if(sub_raw_p < 0.05) {
            message(sprintf("    Significant interaction for %s (p = %.3f)", 
                            subscale, sub_raw_p))
          }
          
          # Store subscale results
          subscale_results[[subscale]] <- list(
            model = subscale_model,
            anova_results = sub_anova,
            raw_p = sub_raw_p,
            diagnostics = check_model_diagnostics(subscale_model, params)
          )
        }
      }
    }
  }
  
  return(list(
    main_results = main_results,
    subscale_results = subscale_results,
    processed_data = processed_data
  ))
}

################## ANALYZE WINNING MODEL ###################
analyze_winning_model <- function(model, data, params) {
  # Run diagnostics
  diagnostics <- check_model_diagnostics(model, params)  # Add params here
  
  # Run ANOVA
  anova_results <- car::Anova(model, type = 2)
  
  # Get p-values and adjust
  p_values <- anova_results[,"Pr(>Chisq)"]
  
  # Get interaction term
  interaction_term <- grep("consensus_level:direction:scale_name", 
                           rownames(anova_results), 
                           value = TRUE)
  
  # Extract raw p-value for three-way interaction
  raw_p <- p_values[interaction_term]
  
  list(
    diagnostics = diagnostics,
    anova_results = anova_results,
    raw_p = raw_p,
    model = model
  )
}

################## PLOTTING THEME ###################
theme_custom <- theme_classic() +  # Changed from theme_minimal()
  theme(
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

################## PLOTTING FUNCTIONS ###################
generate_all_plots <- function(processed_data, moderation_results, var, params) {
  current_plots <- list()
  
  # Main scale plots
  current_plots[["continuous"]] <- tryCatch({
    plot_continuous_relationship(
      data = processed_data,
      var = "scale_name",
      outcome_var = "outcome_value",
      results = moderation_results$main_results,
      params = params,
      display_name = var
    )
  }, error = function(e) {
    print(paste("Error in continuous plot:", e$message))
    return(NULL)
  })
  
  current_plots[["median_split"]] <- tryCatch({
    plot_median_split(
      data = processed_data,
      var = var,
      results = moderation_results$main_results,
      median_var = "scale_name",
      params = params,
      display_name = var
    )
  }, error = function(e) {
    print(paste("Error in median split plot:", e$message))
    return(NULL)
  })
  
  current_plots[["moderation"]] <- tryCatch({
    plot_moderation_effects(
      model = moderation_results$main_results$model,
      data = processed_data,
      var = "scale_name",
      raw_p = moderation_results$main_results$raw_p,
      adj_p = moderation_results$main_results$adj_p,
      params = params,
      display_name = var
    )$plot
  }, error = function(e) {
    print(paste("Error in moderation plot:", e$message))
    return(NULL)
  })
  
  current_plots[["simple_slopes"]] <- tryCatch({
    plot_simple_slopes(
      model = moderation_results$main_results$model,
      var = "scale_name",
      data = processed_data,
      raw_p = moderation_results$main_results$raw_p,
      adj_p = moderation_results$main_results$adj_p,
      params = params,
      display_name = var
    )
  }, error = function(e) {
    print(paste("Error in simple slopes plot:", e$message))
    return(NULL)
  })
  
  # Add subscale plots
  if(!is.null(moderation_results$subscale_results)) {
    for(subscale in names(moderation_results$subscale_results)) {
      sub_results <- moderation_results$subscale_results[[subscale]]
      
      # Create temporary data with subscale as scale_name
      temp_data <- processed_data
      temp_data$scale_name <- temp_data[[subscale]]
      
      current_plots[[paste0(subscale, "_continuous")]] <- tryCatch({
        plot_continuous_relationship(
          data = temp_data,
          var = "scale_name",
          outcome_var = "outcome_value",
          results = sub_results,
          params = params,
          display_name = subscale
        )
      }, error = function(e) {
        print(paste("Error in subscale continuous plot:", e$message))
        return(NULL)
      })
      
      current_plots[[paste0(subscale, "_moderation")]] <- tryCatch({
        plot_moderation_effects(
          model = sub_results$model,
          data = temp_data,
          var = "scale_name",
          raw_p = sub_results$raw_p,
          adj_p = sub_results$adj_p,
          params = params,
          display_name = subscale
        )$plot
      }, error = function(e) {
        print(paste("Error in subscale moderation plot:", e$message))
        return(NULL)
      })
      
      current_plots[[paste0(subscale, "_simple_slopes")]] <- tryCatch({
        plot_simple_slopes(
          model = sub_results$model,
          var = "scale_name",
          data = temp_data,
          raw_p = sub_results$raw_p,
          adj_p = sub_results$adj_p,
          params = params,
          display_name = subscale
        )
      }, error = function(e) {
        print(paste("Error in subscale simple slopes plot:", e$message))
        return(NULL)
      })
    }
  }
  
  return(current_plots)
}

plot_consensus_results <- function(all_results, params) {
  all_plots <- list()
  
  for(var in names(all_results)) {
    results <- all_results[[var]]
    current_plots <- list()
    
    # Main scale plots with updated results structure
    current_plots[["continuous"]] <- plot_continuous_relationship(
      data = results$processed_data,
      var = "scale_name",
      outcome_var = "outcome_value",
      results = results$moderation_results$main_results,  # Updated
      params = params
    )
    
    current_plots[["median_split"]] <- plot_median_split(
      data = results$processed_data,
      var = var,
      results = results$moderation_results$main_results,  # Updated
      median_var = "scale_name",
      params = params
    )
    
    current_plots[["simple_slopes"]] <- plot_simple_slopes(
      model = results$moderation_results$main_results$model,  # Updated
      var = "scale_name",
      data = results$processed_data,
      params = params
    )
    
    # Handle subscale plots with updated results structure
    if(!is.null(results$moderation_results$subscale_results)) {
      for(subscale in names(results$moderation_results$subscale_results)) {
        sub_results <- results$moderation_results$subscale_results[[subscale]]
        
        current_plots[[paste0(subscale, "_continuous")]] <- plot_continuous_relationship(
          data = results$processed_data,
          var = subscale,
          outcome_var = "outcome_value",
          results = sub_results,
          params = params
        )
        
        current_plots[[paste0(subscale, "_simple_slopes")]] <- plot_simple_slopes(
          model = sub_results$model,
          var = subscale,
          data = results$processed_data,
          params = params
        )
      }
    }
    
    all_plots[[var]] <- current_plots
  }
  
  return(all_plots)
}

plot_continuous_relationship <- function(data, var, outcome_var, results, params, display_name) {
  # Print counts per condition
  participant_counts <- data %>% 
    group_by(consensus_level, direction) %>% 
    summarise(
      n = n_distinct(participant.id_in_session),
      .groups = 'drop'
    )
  
  # Prepare data for plotting
  plot_data <- data
  if(params$analysis_type == "switch_difference") {
    plot_data <- plot_data %>%
      mutate(
        # Create trial_type factor based on switch_difference
        trial_type = factor(if_else(switch_difference == 0, "Stay", "Switch"), 
                            levels = c("Stay", "Switch")),
        # Ensure consensus_level is properly ordered
        consensus_level = factor(consensus_level, levels = c("2:2", "3:1", "4:0")),
        direction = factor(direction, levels = c("Against group", "With group"))
      )
  }
  
  # Create subtitle text
  subtitle_text <- if(!is.null(results$raw_p) && !is.na(results$raw_p)) {
    paste("Three-way interaction:\nRaw p =",
          format.pval(results$raw_p, digits = 3),
          "\nFDR-adjusted p =",
          format.pval(results$adj_p, digits = 3))
  } else {
    "P-values not available"
  }
  
  # Calculate correlations
  effects <- plot_data %>%
    group_by(direction, consensus_level) %>%
    summarise(
      correlation = cor(.data[[var]], .data[["outcome_value"]]),
      .groups = 'drop'
    )
  
  # Update y-axis label based on analysis type and plot type
  y_label <- if(params$analysis_type == "switch_difference") {
    if(grepl("accuracy", params$y_label, ignore.case = TRUE)) {
      "Choice 1 accuracy on trial t + 1"
    } else if(grepl("bet", params$y_label, ignore.case = TRUE)) {
      "Bet magnitude on trial t + 1"
    } else {
      params$y_label
    }
  } else {
    params$y_label
  }
  
  # Base plot setup
  p <- ggplot(plot_data, 
              aes(x = .data[[var]], 
                  y = .data[["outcome_value"]], 
                  color = direction)) +
    geom_smooth(method = "lm", 
                formula = y ~ x,
                se = TRUE) +
    labs(x = display_name,
         y = y_label,
         title = paste("Relationship between", display_name, "and", y_label),
         subtitle = subtitle_text,
         caption = paste("Effect sizes (r) range:", 
                         round(min(effects$correlation), 3), "to",
                         round(max(effects$correlation), 3))) +
    scale_color_manual(values = c("Against group" = "red", "With group" = "blue")) +
    theme_custom
  
  # Add specific faceting based on analysis type
  if(params$analysis_type == "switch_difference") {
    p <- p + facet_grid(trial_type ~ consensus_level) +
      theme(panel.spacing = unit(1, "lines"))
  } else if(params$analysis_type == "choice_consensus") {
    p <- p + facet_wrap(~consensus_level)
  }
  
  return(p)
}

plot_median_split <- function(data, var, results, median_var, params, display_name) {
  subtitle_text <- if(!is.null(results$raw_p) && !is.na(results$raw_p)) {
    paste("Three-way interaction:\nRaw p =", format.pval(results$raw_p, digits = 3),
          "\nFDR-adjusted p =", format.pval(results$adj_p, digits = 3))
  } else {
    "P-values not available"
  }
  
  # Different plotting logic based on analysis type
  if(params$analysis_type == "switch_difference") {
    summary_data <- data %>%
      mutate(quest_group = ifelse(.data[[median_var]] > median(.data[[median_var]]), "High", "Low")) %>%
      group_by(consensus_level, direction, quest_group, trial_type) %>%
      summarise(
        mean_outcome = mean(outcome_value),
        se = sd(outcome_value) / sqrt(n()),
        .groups = 'drop'
      )
    
    p <- ggplot(summary_data, 
                aes(x = consensus_level, 
                    y = mean_outcome, 
                    color = direction,
                    linetype = quest_group,
                    group = interaction(direction, quest_group))) +
      geom_line() +
      geom_point() +
      geom_errorbar(aes(ymin = mean_outcome - se, 
                        ymax = mean_outcome + se), 
                    width = 0.2) +
      scale_color_manual(values = c("Against group" = "red", "With group" = "blue")) +
      facet_wrap(~trial_type, ncol = 2) +  # Side by side plots
      labs(x = "Consensus Level",
           y = params$y_label,
           title = paste("Effect of", display_name, "(Median-split)"),
           subtitle = subtitle_text,
           linetype = paste(display_name, "Group")) +
      theme_custom
    
  } else {
    # Original plotting code for non-switch_difference analysis
    summary_data <- data %>%
      mutate(quest_group = ifelse(.data[[median_var]] > median(.data[[median_var]]), "High", "Low")) %>%
      group_by(consensus_level, direction, quest_group) %>%
      summarise(
        mean_outcome = mean(outcome_value),
        se = sd(outcome_value) / sqrt(n()),
        .groups = 'drop'
      )
    
    p <- ggplot(summary_data, 
                aes(x = consensus_level, 
                    y = mean_outcome, 
                    color = direction,
                    linetype = quest_group,
                    group = interaction(direction, quest_group))) +
      geom_line() +
      geom_point() +
      geom_errorbar(aes(ymin = mean_outcome - se, 
                        ymax = mean_outcome + se), 
                    width = 0.2) +
      scale_color_manual(values = c("Against group" = "red", "With group" = "blue")) +
      labs(x = "Consensus Level",
           y = params$y_label,
           title = paste("Effect of", display_name, "(Median-split)"),
           subtitle = subtitle_text,
           linetype = paste(display_name, "Group")) +
      theme_custom
  }
  
  return(p)
}

plot_moderation_effects <- function(model, data, var, raw_p = NULL, adj_p = NULL, params, display_name) {
  if(!params$analysis_type %in% VALID_ANALYSIS_TYPES) {
    stop(paste("Invalid analysis type. Must be one of:", 
               paste(VALID_ANALYSIS_TYPES, collapse = ", ")))
  }
  
  subtitle_text <- if(!is.null(raw_p) && !is.na(raw_p)) {
    paste("Interaction test:",
          "\nRaw p =", format.pval(raw_p, digits = 3),
          "\nFDR-adjusted p =", format.pval(adj_p, digits = 3))
  } else {
    "P-values not available"
  }
  
  common_theme <- theme_custom
  color_scheme <- scale_color_manual(values = c("Against group" = "red", "With group" = "blue"))
  
  if(params$analysis_type == "switch_difference") {
    # Create prediction grid with fewer switch_difference values
    pred_data <- expand.grid(
      consensus_level = levels(data$consensus_level),
      direction = levels(data$direction),
      switch_difference = c(0, 1),  # Just binary switch/stay values
      scale_name = seq(from = -2, to = 2, length.out = 100),
      age = 0
    ) %>%
      mutate(trial_type = factor(ifelse(switch_difference == 0, "Stay", "Switch"), 
                                 levels = c("Stay", "Switch")))
    
    # Generate predictions
    pred_data$predicted <- predict(model, newdata = pred_data, re.form = NA)
    if(params$is_percentage) pred_data$predicted <- pred_data$predicted * 100
    
    # Create plot with trial_type faceting
    p <- ggplot(pred_data, 
                aes(x = scale_name, 
                    y = predicted, 
                    color = direction)) +
      geom_line() +
      facet_grid(trial_type ~ consensus_level) +  # Removed ncol argument
      labs(title = paste("Moderation effect of", display_name),
           subtitle = subtitle_text,
           x = paste(display_name, "score (standardized)"),
           y = params$y_label) +
      common_theme +
      color_scheme
    
  } else if(params$analysis_type == "choice_consensus") {
    # Original code for choice_consensus
    pred_data <- expand.grid(
      consensus_level = levels(data$consensus_level),
      direction = levels(data$direction),
      scale_name = seq(from = -2, to = 2, length.out = 100),
      age = 0
    )
    
    pred_data$predicted <- predict(model, newdata = pred_data, re.form = NA)
    if(params$is_percentage) pred_data$predicted <- pred_data$predicted * 100
    
    p <- ggplot(pred_data, 
                aes(x = scale_name, 
                    y = predicted, 
                    color = direction)) +
      geom_line() +
      facet_wrap(~consensus_level) +
      labs(title = paste("Moderation effect of", display_name),
           subtitle = subtitle_text,
           x = paste(display_name, "score (standardized)"),
           y = params$y_label) +
      common_theme +
      color_scheme
  }
  
  return(list(
    plot = p,
    predictions = pred_data,
    raw_p = raw_p,
    adj_p = adj_p
  ))
}

plot_simple_slopes <- function(model, var, data, raw_p = NULL, adj_p = NULL, params, display_name) {
  if(!params$analysis_type %in% VALID_ANALYSIS_TYPES) {
    stop(paste("Invalid analysis type. Must be one of:", 
               paste(VALID_ANALYSIS_TYPES, collapse = ", ")))
  }
  
  subtitle_text <- if(!is.null(raw_p) && !is.na(raw_p)) {
    paste("Interaction test:",
          "\nRaw p =", format.pval(raw_p, digits = 3),
          "\nFDR-adjusted p =", format.pval(adj_p, digits = 3))
  } else {
    "P-values not available"
  }
  
  common_theme <- theme_custom
  color_scheme <- scale_color_manual(values = c("Against group" = "red", "With group" = "blue"))
  
  if(params$analysis_type == "switch_difference") {
    # Create prediction grid including trial_type
    pred_data <- expand.grid(
      consensus_level = levels(data$consensus_level),
      direction = levels(data$direction),
      switch_difference = c(0, 1),  # Binary values for Stay (0) and Switch (1)
      scale_name = c(-1, 0, 1),
      age = 0
    ) %>%
      mutate(trial_type = factor(ifelse(switch_difference == 0, "Stay", "Switch"), 
                                 levels = c("Stay", "Switch")))
    
    # Generate predictions
    pred_data$predicted <- predict(model, newdata = pred_data, re.form = NA)
    if(params$is_percentage) pred_data$predicted <- pred_data$predicted * 100
    
    # Create plot with trial_type faceting
    p <- ggplot(pred_data, 
                aes(x = consensus_level, 
                    y = predicted, 
                    color = direction)) +
      geom_line(aes(linetype = factor(scale_name,
                                      labels = c("-1 SD", "Mean", "+1 SD")),
                    group = interaction(direction, scale_name))) +
      geom_point() +
      facet_wrap(~trial_type, ncol = 2) +
      labs(title = paste("Simple slopes analysis for", display_name),
           subtitle = subtitle_text,
           x = "Consensus Level",
           y = params$y_label,
           color = "Direction",
           linetype = paste(display_name, "Score")) +
      common_theme +
      color_scheme
    
  } else if(params$analysis_type == "choice_consensus") {
    # Original code for choice_consensus remains unchanged
    pred_data <- expand.grid(
      consensus_level = levels(data$consensus_level),
      direction = levels(data$direction),
      scale_name = c(-1, 0, 1),
      age = 0
    )
    
    pred_data$predicted <- predict(model, newdata = pred_data, re.form = NA)
    if(params$is_percentage) pred_data$predicted <- pred_data$predicted * 100
    
    p <- ggplot(pred_data, 
                aes(x = consensus_level, 
                    y = predicted, 
                    color = direction,
                    linetype = factor(scale_name,
                                      labels = c("-1 SD", "Mean", "+1 SD")),
                    group = interaction(direction, scale_name))) +
      geom_line() +
      geom_point() +
      labs(title = paste("Simple slopes analysis for", display_name),
           subtitle = subtitle_text,
           x = "Consensus Level",
           y = params$y_label,
           color = "Direction",
           linetype = paste(display_name, "Score")) +
      common_theme +
      color_scheme
  }
  
  return(p)
}

################## CHECK MODEL DIAGNOSTICS ###################
check_model_diagnostics <- function(model, params) {
  # Get model frame and terms
  model_frame <- model@frame
  model_terms <- terms(model)
  
  # Extract residuals with appropriate handling for large samples
  residuals <- residuals(model)
  if(length(residuals) > 5000) {
    set.seed(123)
    shapiro_residuals <- sample(residuals, 5000)
  } else {
    shapiro_residuals <- residuals
  }
  
  # Calculate basic diagnostics
  shapiro_test <- shapiro.test(shapiro_residuals)
  outliers <- boxplot.stats(residuals)$out
  fitted_vals <- fitted(model)
  
  # Calculate VIF for fixed effects
  fixed_terms <- attr(terms(model), "term.labels")
  fixed_formula <- reformulate(fixed_terms)
  vif_values <- tryCatch({
    fixed_model <- lm(fixed_formula, data = model_frame)
    car::vif(fixed_model)
  }, error = function(e) {
    warning("Could not calculate VIF values")
    NA
  })
  
  # Calculate R² values
  r2_values <- tryCatch(
    MuMIn::r.squaredGLMM(model),
    error = function(e) c(NA, NA)
  )
  
  # Check convergence
  convergence <- !is.null(model@optinfo$conv$lme4$messages)
  
  # Test homoscedasticity
  homo_test <- cor.test(abs(residuals), fitted_vals)
  
  # Get random effects structure
  ranef_summary <- summary(ranef(model))
  
  # Analysis-specific checks
  if(params$analysis_type == "switch_difference") {
    # Additional checks specific to switch difference analysis
    switch_diff_diagnostics <- tryCatch({
      list(
        switch_diff_range = range(model_frame$switch_difference),
        switch_diff_distribution = shapiro.test(model_frame$switch_difference)
      )
    }, error = function(e) NULL)
  }
  
  # Compile results
  diagnostics <- list(
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
    ),
    random_effects = ranef_summary
  )
  
  # Add analysis-specific diagnostics
  if(params$analysis_type == "switch_difference" && !is.null(switch_diff_diagnostics)) {
    diagnostics$switch_difference <- switch_diff_diagnostics
  }
  
  return(diagnostics)
}

################## MODERATION ANALYSIS ###################
run_moderation_analysis <- function(data, analysis_results, params) {
  print("Starting moderation analysis...")
  print(paste("Analyzing scale:", params$current_var))
  
  # Main scale models
  print("Fitting base model...")
  base_model <- tryCatch({
    lmer(
      outcome_value ~ consensus_level + direction + scale_name + 
        age + (1|participant.id_in_session),
      data = data,
      control = lmerControl(optimizer = "bobyqa")
    )
  }, error = function(e) {
    print(paste("Error fitting base model:", e$message))
    return(NULL)
  })
  
  print("Fitting moderation model...")
  mod_model <- tryCatch({
    lmer(
      outcome_value ~ consensus_level * direction * scale_name + 
        age + (1|participant.id_in_session),
      data = data,
      control = lmerControl(optimizer = "bobyqa")
    )
  }, error = function(e) {
    print(paste("Error fitting moderation model:", e$message))
    return(NULL)
  })
  
  if(is.null(base_model) || is.null(mod_model)) {
    stop("Failed to fit one or both models")
  }
  
  # Model comparisons and diagnostics
  print("Running model diagnostics...")
  main_anova <- car::Anova(mod_model, type = 2)
  model_comparison <- anova(base_model, mod_model)
  model_diagnostics <- check_model_diagnostics(mod_model, params)  # Add params here
  
  print("Model comparison results:")
  print(model_comparison)
  
  print("ANOVA results:")
  print(main_anova)
  
  # Calculate effect sizes
  print("Calculating effect sizes...")
  effect_sizes <- tryCatch({
    effectsize::eta_squared(mod_model, partial = TRUE)
  }, error = function(e) {
    print(paste("Error calculating effect sizes:", e$message))
    return(NULL)
  })
  
  print("Effect sizes:")
  print(effect_sizes)
  
  # Get interaction terms and p-values
  interaction_term <- "consensus_level:direction:scale_name"
  raw_p <- main_anova[interaction_term, "Pr(>Chisq)"]
  all_p_values <- main_anova[, "Pr(>Chisq)"]
  
  print(paste("Three-way interaction p-value:", raw_p))
  
  # Create moderation plot
  print("Creating moderation plot...")
  main_moderation_plot <- plot_moderation_effects(
    model = mod_model,
    data = data,
    var = "scale_name",
    params = params,
    raw_p = raw_p
  )
  
  # Store main results
  main_results <- list(
    model = mod_model,
    base_model = base_model,
    anova_results = main_anova,
    comparison = model_comparison,
    raw_p = raw_p,
    all_p_values = all_p_values,
    interaction_term = interaction_term,
    moderation_plot = main_moderation_plot,
    diagnostics = model_diagnostics,
    effect_sizes = effect_sizes
  )
  
  # Store subscale results
  subscale_results <- list()
  
  # Run subscale analyses if applicable
  current_scale <- params$current_var
  if(!is.null(params$subscale_mapping[[current_scale]])) {
    print(paste("Processing subscales for", current_scale))
    for(subscale in params$subscale_mapping[[current_scale]]) {
      print(paste("Analyzing subscale:", subscale))
      
      if(subscale %in% colnames(data)) {
        # Fit subscale models
        print(paste("Fitting models for subscale:", subscale))
        sub_base_model <- tryCatch({
          lmer(
            as.formula(paste0(
              "outcome_value ~ consensus_level + direction + ", subscale,
              " + age + (1|participant.id_in_session)"
            )),
            data = data,
            control = lmerControl(optimizer = "bobyqa")
          )
        }, error = function(e) {
          print(paste("Error fitting subscale base model:", e$message))
          return(NULL)
        })
        
        sub_mod_model <- tryCatch({
          lmer(
            as.formula(paste0(
              "outcome_value ~ consensus_level * direction * ", subscale,
              " + age + (1|participant.id_in_session)"
            )),
            data = data,
            control = lmerControl(optimizer = "bobyqa")
          )
        }, error = function(e) {
          print(paste("Error fitting subscale moderation model:", e$message))
          return(NULL)
        })
        
        if(!is.null(sub_base_model) && !is.null(sub_mod_model)) {
          # Calculate subscale diagnostics and statistics
          sub_anova <- car::Anova(sub_mod_model, type = 2)
          sub_comparison <- anova(sub_base_model, sub_mod_model)
          sub_diagnostics <- check_model_diagnostics(sub_mod_model, params)  # Add params here
          
          print(paste("Calculating effect sizes for subscale:", subscale))
          sub_effect_sizes <- tryCatch({
            effectsize::eta_squared(sub_mod_model, partial = TRUE)
          }, error = function(e) {
            print(paste("Error calculating subscale effect sizes:", e$message))
            return(NULL)
          })
          
          # Get subscale interaction terms and p-values
          sub_interaction_term <- paste0("consensus_level:direction:", subscale)
          sub_raw_p <- sub_anova[sub_interaction_term, "Pr(>Chisq)"]
          sub_all_p_values <- sub_anova[, "Pr(>Chisq)"]
          
          print(paste("Subscale three-way interaction p-value:", sub_raw_p))
          
          # Create subscale moderation plot
          sub_moderation_plot <- plot_moderation_effects(
            model = sub_mod_model,
            data = data,
            var = subscale,
            params = params,
            raw_p = sub_raw_p
          )
          
          subscale_results[[subscale]] <- list(
            model = sub_mod_model,
            base_model = sub_base_model,
            anova_results = sub_anova,
            comparison = sub_comparison,
            raw_p = sub_raw_p,
            all_p_values = sub_all_p_values,
            interaction_term = sub_interaction_term,
            diagnostics = sub_diagnostics,
            effect_sizes = sub_effect_sizes,
            moderation_plot = sub_moderation_plot
          )
        }
      }
    }
  }
  
  print("Moderation analysis complete.")
  return(list(
    main_results = main_results,
    subscale_results = subscale_results
  ))
}

################## HELPER FUNCTIONS FOR RESULTS ###################

# Calculate effect sizes using partial eta-squared
calculate_effect_sizes <- function(model, params) {
  tryCatch({
    # Get ANOVA results
    anova_results <- anova(model)
    
    # Calculate partial eta-squared
    eta_sq <- effectsize::F_to_eta2(
      f = anova_results$`F value`,
      df = anova_results$NumDF,
      df_error = anova_results$DenDF,
      ci = 0.95
    )
    
    # Get relevant terms based on analysis type
    if(params$analysis_type == "switch_difference") {
      main_terms <- c("switch_difference", "consensus_level", "direction", "scale_name")
    } else if(params$analysis_type == "choice_consensus") {
      main_terms <- c("consensus_level", "direction", "scale_name")
    }
    
    # Filter and format results
    formatted_eta <- data.frame(
      Term = rownames(anova_results),
      Eta2_partial = sprintf("%.3f", eta_sq$Eta2_partial),
      CI_low = sprintf("%.3f", eta_sq$CI_low),
      CI_high = sprintf("%.3f", eta_sq$CI_high),
      stringsAsFactors = FALSE
    )
    
    # Add indicator for main effects vs interactions
    formatted_eta$Effect_Type <- sapply(formatted_eta$Term, function(x) {
      if(x %in% main_terms) "Main Effect" else "Interaction"
    })
    
    return(formatted_eta)
    
  }, error = function(e) {
    warning(paste("Error calculating effect sizes:", e$message))
    return(NULL)
  })
}

# Calculate simple slopes
compute_simple_slopes <- function(model, data) {
  tryCatch({
    # Create list to store results
    slopes_results <- list()
    
    # For each consensus level
    for(cons in levels(data$consensus_level)) {
      # For each direction
      for(dir in levels(data$direction)) {
        # Create new data for this combination
        test_data <- data.frame(
          consensus_level = cons,
          direction = dir,
          scale_name = c(-1, 1),  # Test at ±1 SD
          age = 0  # Set to mean
        )
        
        # Get predictions
        preds <- predict(model, newdata = test_data, re.form = NA)
        
        # Calculate slope
        slope <- (preds[2] - preds[1]) / 2  # Divide by difference in scale_name
        
        # Store results
        slopes_results[[paste(cons, dir)]] <- slope
      }
    }
    
    return(slopes_results)
  }, error = function(e) {
    return(NULL)
  })
}

################## RESULTS FORMATTING ###################
format_results <- function(all_results, analysis_name, params) {
  output_text <- sprintf("%s - STATISTICAL ANALYSIS RESULTS\n\n", toupper(analysis_name))
  output_text <- paste0(output_text, "Analysis run on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  output_text <- paste0(output_text, "Analysis type: ", params$analysis_type, "\n\n")
  
  # Get first result for pooled comparison
  first_result <- all_results[[1]]
  
  # Add pooled model comparison section
  output_text <- paste0(output_text, "POOLED MODEL COMPARISON:\n=======================\n")
  
  # Create model comparison table
  model_comp_df <- data.frame(
    Model = sprintf("m%d", seq_along(first_result$model_results$aic_values)),
    AIC = round(first_result$model_results$aic_values, 3),
    BIC = round(first_result$model_results$bic_values, 3),
    Is_Winner = ifelse(seq_along(first_result$model_results$aic_values) == 
                         first_result$model_results$winning_idx, "Yes", "No")
  )
  rownames(model_comp_df) <- NULL
  
  # Model descriptions
  model_explanations <- sprintf("m%d: %s", 
                                seq_along(first_result$model_results$descriptions),
                                first_result$model_results$descriptions)
  
  # Add comparison section
  output_text <- paste0(output_text, 
                        "MODEL COMPARISON:\n================\n",
                        paste(capture.output(print.data.frame(model_comp_df, row.names = FALSE)), 
                              collapse = "\n"),
                        "\n\nModel Descriptions:\n=================\n",
                        paste(model_explanations, collapse = "\n"),
                        "\n\nWinning model formula:\n",
                        paste(deparse(formula(first_result$model_results$winning_model)), 
                              collapse = "\n"),
                        "\n\n")
  
  # Add winning model diagnostics with analysis type specific checks
  diag <- first_result$analysis_results$diagnostics
  output_text <- paste0(output_text,
                        "WINNING MODEL DIAGNOSTICS (POOLED):\n================================\n",
                        "R² marginal/conditional: ", 
                        round(diag$model_fit$r2_marginal, 3), "/",
                        round(diag$model_fit$r2_conditional, 3), "\n",
                        "Residual normality p: ", 
                        format.pval(diag$residuals$normality_p, digits = 3), "\n",
                        "VIF range: ", format_vif_range(diag$assumptions$vif), "\n",
                        "Homoscedasticity: p = ", 
                        format.pval(diag$assumptions$homoscedasticity_p, digits = 3), "\n",
                        "Outliers: ", diag$outliers$n_outliers, "\n")
  
  # Add switch difference diagnostics if applicable
  if(params$analysis_type == "switch_difference" && !is.null(diag$switch_difference)) {
    output_text <- paste0(output_text,
                          "Switch difference range: ", 
                          paste(round(diag$switch_difference$switch_diff_range, 3), collapse = " to "), "\n",
                          "Switch difference normality p: ",
                          format.pval(diag$switch_difference$switch_diff_distribution$p.value, digits = 3), "\n")
  }
  output_text <- paste0(output_text, "\n")
  
  # Individual questionnaire results
  for(var in names(all_results)) {
    results <- all_results[[var]]
    
    output_text <- paste0(output_text, 
                          "\n", var, " ANALYSIS\n",
                          "=================\n")
    
    # Add ANOVA results
    output_text <- paste0(output_text,
                          "\nANOVA RESULTS:\n==============\n",
                          paste(capture.output(results$analysis_results$anova_results), 
                                collapse = "\n"),
                          "\n")
    
    # Add effect sizes
    if(!is.null(results$moderation_results$main_results$model)) {
      effect_sizes <- calculate_effect_sizes(results$moderation_results$main_results$model, params)
      if(!is.null(effect_sizes) && nrow(effect_sizes) > 0) {
        effect_size_table <- format_effect_size_table(effect_sizes)
        output_text <- paste0(output_text, effect_size_table, "\n")
      }
    }
    
    # Add interaction results based on analysis type
    if(params$analysis_type == "switch_difference") {
      output_text <- paste0(output_text,
                            "\nINTERACTION RESULTS:\n===================\n")
      # Extract p-values from results structure
      raw_p <- results$moderation_results$main_results$raw_p
      adj_p <- results$moderation_results$main_results$adj_p
      
      output_text <- paste0(output_text,
                            "Switch difference x Consensus x Direction interaction:\n",
                            "Raw p = ", format.pval(raw_p, digits = 3), "\n",
                            "FDR-adjusted p = ", format.pval(adj_p, digits = 3),
                            "\n")
    } else {
      output_text <- paste0(output_text,
                            "\nTHREE-WAY INTERACTION:\n=====================\n")
      # Extract p-values from results structure
      raw_p <- results$moderation_results$main_results$raw_p
      adj_p <- results$moderation_results$main_results$adj_p
      
      output_text <- paste0(output_text,
                            "Raw p = ", format.pval(raw_p, digits = 3), "\n",
                            "FDR-adjusted p = ", format.pval(adj_p, digits = 3),
                            "\n")
    }
    
    # Add model diagnostics
    output_text <- paste0(output_text, 
                          "\n", format_diagnostics(results$moderation_results$main_results$diagnostics))
    
    # Add subscale results
    if(!is.null(results$moderation_results$subscale_results)) {
      output_text <- paste0(output_text, format_subscale_results(
        results$moderation_results$subscale_results, params))
    }
  }
  
  return(output_text)
}

format_diagnostics <- function(diagnostics) {
  output <- paste0(
    "MODEL DIAGNOSTICS:\n",
    "=================\n",
    "Residuals:\n",
    "  Mean: ", round(diagnostics$residuals$mean, 3), "\n",
    "  SD: ", round(diagnostics$residuals$sd, 3), "\n",
    "  Normality test p-value: ", 
    format.pval(diagnostics$residuals$normality_p, digits = 3), "\n",
    "  N: ", diagnostics$residuals$n_residuals, "\n\n",
    
    "Model Fit:\n",
    "  R² marginal: ", round(diagnostics$model_fit$r2_marginal, 3), "\n",
    "  R² conditional: ", round(diagnostics$model_fit$r2_conditional, 3), "\n",
    "  AIC: ", round(diagnostics$model_fit$aic, 2), "\n",
    "  BIC: ", round(diagnostics$model_fit$bic, 2), "\n\n",
    
    "Assumptions:\n",
    "  VIF range: ", format_vif_range(diagnostics$assumptions$vif), "\n",
    "  Homoscedasticity correlation: ", 
    round(diagnostics$assumptions$homoscedasticity_cor, 3),
    " (p = ", format.pval(diagnostics$assumptions$homoscedasticity_p, digits = 3), ")\n",
    "  Model converged: ", !diagnostics$assumptions$convergence, "\n",
    "  Random effects: ", paste(diagnostics$assumptions$random_effects, collapse = ", "), "\n\n",
    
    "Outliers:\n",
    "  Count: ", diagnostics$outliers$n_outliers, "\n\n"
  )
  
  # Add switch difference diagnostics if they exist
  if(!is.null(diagnostics$switch_difference)) {
    output <- paste0(output,
                     "Switch Difference:\n",
                     "  Range: ", paste(round(diagnostics$switch_difference$switch_diff_range, 3), 
                                        collapse = " to "), "\n",
                     "  Normality p: ", 
                     format.pval(diagnostics$switch_difference$switch_diff_distribution$p.value,
                                 digits = 3), "\n\n")
  }
  
  return(output)
}

format_vif_range <- function(vif_values) {
  if(is.null(vif_values) || all(is.na(vif_values))) return("NA")
  paste(round(range(vif_values, na.rm = TRUE), 2), collapse = " - ")
}

# Helper function for effect size table formatting
format_effect_size_table <- function(effect_sizes) {
  table <- paste0(
    "\nEFFECT SIZES:\n============\n",
    sprintf("%-40s %10s %10s %10s\n", "Term", "Eta2", "CI_low", "CI_high")
  )
  
  for(i in 1:nrow(effect_sizes)) {
    table <- paste0(
      table,
      sprintf("%-40s %10s %10s %10s\n",
              effect_sizes$Term[i],
              effect_sizes$Eta2_partial[i],
              effect_sizes$CI_low[i],
              effect_sizes$CI_high[i])
    )
  }
  
  return(table)
}

# Helper function for subscale results formatting
format_subscale_results <- function(subscale_results, params) {
  output <- paste0("\nSUBSCALE RESULTS:\n================\n")
  
  for(subscale_name in names(subscale_results)) {
    sub_result <- subscale_results[[subscale_name]]
    
    output <- paste0(output,
                     "\nSubscale: ", subscale_name, "\n",
                     "----------------------\n")
    
    if(!is.null(sub_result$anova_results)) {
      output <- paste0(output,
                       "ANOVA Results:\n-------------\n",
                       paste(capture.output(sub_result$anova_results), 
                             collapse = "\n"),
                       "\n")
    }
    
    if(!is.null(sub_result$effect_sizes)) {
      output <- paste0(output, 
                       format_effect_size_table(sub_result$effect_sizes))
    }
    
    if(!is.null(sub_result$raw_p) && !is.null(sub_result$adj_p)) {
      if(params$analysis_type == "switch_difference") {
        output <- paste0(output,
                         "\nSwitch Difference Interaction:\n",
                         "Raw p = ", format.pval(sub_result$raw_p, digits = 3), "\n",
                         "FDR-adjusted p = ", format.pval(sub_result$adj_p, digits = 3),
                         "\n")
      } else {
        output <- paste0(output,
                         "\nThree-way Interaction:\n",
                         "Raw p = ", format.pval(sub_result$raw_p, digits = 3), "\n",
                         "FDR-adjusted p = ", format.pval(sub_result$adj_p, digits = 3),
                         "\n")
      }
    }
    
    if(!is.null(sub_result$diagnostics)) {
      output <- paste0(output, "\n", format_diagnostics(sub_result$diagnostics))
    }
  }
  
  return(output)
}


format_moderation_results <- function(moderation_results, params) {
  output <- "MODERATION ANALYSIS RESULTS:\n"
  output <- paste0(output, "==========================\n\n")
  output <- paste0(output, "Analysis type: ", params$analysis_type, "\n\n")
  
  # Add main moderation results
  output <- paste0(output,
                   "Main Model Comparison:\n",
                   paste(capture.output(moderation_results$comparison), 
                         collapse = "\n"),
                   "\n\n")
  
  # Add analysis-specific results
  if(params$analysis_type == "switch_difference") {
    output <- paste0(output,
                     "Switch Difference Effects:\n",
                     "------------------------\n",
                     paste(capture.output(moderation_results$switch_effects), 
                           collapse = "\n"),
                     "\n\n")
  }
  
  # Add subscale results if present
  if(!is.null(moderation_results$subscales)) {
    output <- paste0(output, "SUBSCALE RESULTS:\n",
                     "=================\n")
    
    for(sub_result in moderation_results$subscales) {
      if(!is.null(sub_result)) {
        output <- paste0(output,
                         "\nSubscale: ", sub_result$subscale, "\n",
                         paste(capture.output(sub_result$comparison), 
                               collapse = "\n"),
                         "\n")
      }
    }
  }
  
  return(output)
}

################## SAVE RESULTS ###################
save_results <- function(all_results, analysis_name, output_path, params) {
  # Generate formatted results
  output_text <- format_results(all_results, analysis_name, params)
  
  # Create output directory if it doesn't exist
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  # Save to file
  writeLines(output_text, output_path)
  
  # Print confirmation
  print(paste("Results saved to:", output_path))
}

# Helper function to create appropriate file names based on analysis type
create_output_filename <- function(params, prefix) {
  paste0(params$analysis_type, "_", prefix, ".txt")
}