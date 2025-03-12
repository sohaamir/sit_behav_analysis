# questionnaire_behavioural_plots.R

# Add this at the start of your script
options(immediate.write = TRUE)

################## SETUP AND CONSTANTS ###################
required_packages <- c(
  "ggplot2", "tidyverse", "ggpubr", "rstatix", "ez", 
  "lme4", "lmerTest", "car", "emmeans", "MuMIn", 
  "psych", "interactions", "effects", "here", 
  "lm.beta", "effectsize", "corrplot", "reshape2", "gridExtra", "grid", "conflicted"
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
conflicts_prefer(effectsize::cohens_d)

# Analysis types for directory structure
# Valid analysis types for validation
VALID_ANALYSIS_TYPES <- c(
  "choice_consensus",    # Original choice/bet by consensus analysis
  "switch_difference",   # Switch difference analysis
  "within_trial_switch"  # Within-trial switch analysis
)

# Plot types for directory structure
PLOT_TYPES <- c(
  "continuous",      # For continuous relationship plots
  "median_split",    # For median split plots
  "moderation",      # For moderation effect plots
  "simple_slopes"    # For simple slopes analysis plots
)

# Subscale mapping for questionnaires
SUBSCALE_MAPPING <- list(
  lsas = c("lsas_p", "lsas_s"),
  dass = c("dass_a", "dass_d", "dass_s"),
  ssms = c("ssms_cd", "ssms_ia"),
  srp_sf = c("srp_sf_ipm", "srp_sf_ca", "srp_sf_els", "srp_sf_ct"),
  ami = c("ami_es", "ami_sm", "ami_ba"),
  aq_10 = NULL
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
          
          # Extract plot object and handle different plot types
          plot_obj <- if(is.list(scale_results$plots[[plot_name]]) && 
                         "plot" %in% names(scale_results$plots[[plot_name]])) {
            scale_results$plots[[plot_name]]$plot
          } else if(inherits(scale_results$plots[[plot_name]], "gtable")) {
            # Handle gridExtra arranged plots
            scale_results$plots[[plot_name]]
          } else {
            scale_results$plots[[plot_name]]
          }
          
          # Save plot with appropriate handling
          if(inherits(plot_obj, "gtable")) {
            # For gridExtra arranged plots
            png(file.path(output_dir, file_name), width = 12, height = 6, units = "in", res = 300)
            grid::grid.draw(plot_obj)
            dev.off()
          } else if(inherits(plot_obj, "ggplot")) {
            # For regular ggplot objects
            ggsave(
              filename = file.path(output_dir, file_name),
              plot = plot_obj,
              width = 10,
              height = 7
            )
          }
        }
      }
    }
  }
}

################## CORE MODEL PIPELINE ###################
run_model_pipeline <- function(data, params, output_path) {
  # Print trial counts based on analysis type
  if(params$analysis_type == "switch_difference") {
    cat("\n===== Trial Type Summary =====\n")
    stay_count <- sum(data$choice_switch_across_trials == 0, na.rm = TRUE)
    switch_count <- sum(data$choice_switch_across_trials == 1, na.rm = TRUE)
    cat(sprintf("Stay trials: %d\n", stay_count))
    cat(sprintf("Switch trials: %d\n", switch_count))
    cat("============================\n")
  } else if(params$analysis_type == "within_trial_switch") {
    cat("\n===== Trial Type Summary =====\n")
    stay_count <- sum(data$player.switch_vs_stay == 0, na.rm = TRUE)
    switch_count <- sum(data$player.switch_vs_stay == 1, na.rm = TRUE)
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
  median_split_p_values <- list()
  
  # Apply winning model structure to each MAIN questionnaire only
  message("\nAnalyzing main questionnaires...")
  for(var in params$questionnaire_vars) {
    message(sprintf("\nProcessing: %s", var))
    
    # Process data
    processed_data <- params$preprocessing_fn(data, var, params)
    
    # Run analysis with winning model structure
    analysis_results <- run_individual_analysis(
      processed_data,
      pooled_results$winning_model,
      var,
      params
    )
    
    # Store raw p-values for main questionnaires only (simplified naming)
    interaction_p_values[[var]] <- analysis_results$main_results$raw_p
    
    # Store median split p-values for main questionnaires only
    if(!is.null(analysis_results$main_results$median_split) && 
       !is.na(analysis_results$main_results$median_split$p_value)) {
      median_split_p_values[[var]] <- analysis_results$main_results$median_split$p_value
    }
    
    # Store results
    all_results[[var]] <- list(
      model_results = pooled_results,
      analysis_results = analysis_results$main_results,
      moderation_results = analysis_results,
      processed_data = processed_data
    )
  }
  
  # Perform FDR correction - now only across main questionnaires
  message("\nPerforming FDR correction (main questionnaires only)...")
  
  # Correct continuous model p-values
  all_p_values <- unlist(interaction_p_values)
  all_adj_p_values <- p.adjust(all_p_values, method = "fdr")
  names(all_adj_p_values) <- names(all_p_values)
  
  # Correct median split p-values if any exist
  if(length(median_split_p_values) > 0) {
    all_median_p_values <- unlist(median_split_p_values)
    all_median_adj_p_values <- p.adjust(all_median_p_values, method = "fdr")
    names(all_median_adj_p_values) <- names(all_median_p_values)
  } else {
    all_median_adj_p_values <- NULL
  }
  
  # Update results with corrected p-values and generate plots
  for(var in params$questionnaire_vars) {
    # Update main scale adjusted p-value
    all_results[[var]]$moderation_results$main_results$adj_p <- all_adj_p_values[var]
    
    # Update median split adjusted p-value
    if(!is.null(all_results[[var]]$moderation_results$main_results$median_split) &&
       !is.null(all_median_adj_p_values) && var %in% names(all_median_adj_p_values)) {
      all_results[[var]]$moderation_results$main_results$median_split$adj_p <- 
        all_median_adj_p_values[var]
    }
    
    # Generate plots for main questionnaires only
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
    processed <- params$preprocessing_fn(data, var, params)
    
    # Base columns to select
    base_cols <- c("participant.id_in_session", "consensus_level", 
                   "direction", "age", "gender", 
                   "scale_name", "outcome_value")
    
    # Add appropriate switch variable based on analysis type
    if(params$analysis_type == "switch_difference") {
      base_cols <- c(base_cols, "switch_difference")
    } else if(params$analysis_type == "within_trial_switch") {
      base_cols <- c(base_cols, "switch_vs_stay")
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
  
  # Extract p-value for interaction based on analysis type
  if(params$analysis_type == "within_trial_switch") {
    interaction_term <- "consensus_level:direction:switch_vs_stay:scale_name"
  } else if(params$analysis_type == "switch_difference") {
    interaction_term <- "consensus_level:direction:switch_difference:scale_name"
  } else {
    interaction_term <- "consensus_level:direction:scale_name"
  }
  
  raw_p <- tryCatch({
    if(interaction_term %in% rownames(main_anova)) {
      main_anova[interaction_term, "Pr(>Chisq)"]
    } else {
      message(sprintf("Interaction term '%s' not found in ANOVA results", interaction_term))
      NA
    }
  }, error = function(e) {
    message(sprintf("Error extracting p-value: %s", e$message))
    NA
  })
  
  if(!is.na(raw_p) && raw_p < 0.05) {
    message(sprintf("  Significant interaction found (p = %.3f)", raw_p))
  }
  
  # Run median split analysis
  message("  Running median split analysis with winning model...")
  median_analysis <- run_median_split_analysis(processed_data, winning_model, params)
  
  # Store main results including median split analysis
  main_results <- list(
    model = main_model,
    anova_results = main_anova,
    raw_p = raw_p,
    diagnostics = check_model_diagnostics(main_model, params),
    median_split = median_analysis
  )
  
  # No subscale results
  return(list(
    main_results = main_results,
    processed_data = processed_data
  ))
}
run_median_split_analysis <- function(data, winning_model, params) {
  # Create median-split data
  prepped_data <- data %>%
    mutate(quest_group = factor(
      ifelse(!is.na(scale_name), 
             ifelse(scale_name > median(scale_name, na.rm = TRUE), "High", "Low"),
             NA),
      levels = c("Low", "High")
    )) %>%
    filter(!is.na(quest_group))
  
  # Extract the formula from the winning model and adapt it
  model_formula <- formula(winning_model)
  
  # Convert formula to string
  formula_str <- as.character(model_formula)
  formula_rhs <- formula_str[3] # Right-hand side of the formula
  
  # Replace scale_name with quest_group in the formula
  new_formula_str <- gsub("scale_name", "quest_group", formula_rhs)
  
  # Construct new formula
  new_formula <- as.formula(paste("outcome_value ~", new_formula_str))
  
  # Fit the model with adapted formula
  median_model <- tryCatch({
    lmer(new_formula, data = prepped_data)
  }, error = function(e) {
    message(paste("Error fitting median split model:", e$message))
    return(NULL)
  })
  
  # If model fitting failed, return NA results
  if(is.null(median_model)) {
    return(list(p_value = NA, model = NULL, anova = NULL, data = prepped_data))
  }
  
  # Get ANOVA results
  median_anova <- car::Anova(median_model, type = 2)
  
  # Get appropriate interaction term based on analysis type
  if(params$analysis_type == "within_trial_switch") {
    interaction_term <- "consensus_level:direction:switch_vs_stay:quest_group"
  } else if(params$analysis_type == "switch_difference") {
    interaction_term <- "consensus_level:direction:switch_difference:quest_group"
  } else {
    interaction_term <- "consensus_level:direction:quest_group"
  }
  
  # Extract p-value
  p_value <- if(interaction_term %in% rownames(median_anova)) {
    median_anova[interaction_term, "Pr(>Chisq)"]
  } else {
    NA
  }
  
  return(list(
    p_value = p_value,
    model = median_model,
    anova = median_anova,
    data = prepped_data
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
generate_all_plots <- function(processed_data, moderation_results, var, params, display_name) {
  print("Levels of factors:")
  print(levels(processed_data$switch_vs_stay))
  print(levels(processed_data$consensus_level))
  print(levels(processed_data$direction))
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
  
  # Add subscale plots if they exist
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
    
    # Main scale plots only
    current_plots[["continuous"]] <- plot_continuous_relationship(
      data = results$processed_data,
      var = "scale_name",
      outcome_var = "outcome_value",
      results = results$moderation_results$main_results,
      params = params
    )
    
    current_plots[["median_split"]] <- plot_median_split(
      data = results$processed_data,
      var = var,
      results = results$moderation_results$main_results,
      median_var = "scale_name",
      params = params
    )
    
    current_plots[["simple_slopes"]] <- plot_simple_slopes(
      model = results$moderation_results$main_results$model,
      var = "scale_name",
      data = results$processed_data,
      params = params
    )
    
    # No subscale plots
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
        trial_type = factor(if_else(switch_difference == 0, "Stay", "Switch"), 
                            levels = c("Stay", "Switch")),
        consensus_level = factor(consensus_level, levels = c("2:2", "3:1", "4:0")),
        direction = factor(direction, levels = c("Against group", "With group"))
      )
  } else if(params$analysis_type == "within_trial_switch") {
    plot_data <- plot_data %>%
      mutate(
        switch_vs_stay = factor(switch_vs_stay, 
                                levels = c(0, 1),
                                labels = c("Stay trials", "Switch trials")),
        consensus_level = factor(consensus_level, levels = c("2:2", "3:1", "4:0")),
        direction = factor(direction, levels = c("Against group", "With group"))
      )
  } else {
    plot_data <- plot_data %>%
      mutate(
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
  
  # Base plot setup
  p <- ggplot(plot_data, 
              aes(x = .data[[var]], 
                  y = .data[["outcome_value"]], 
                  color = direction)) +
    geom_smooth(method = "lm", 
                formula = y ~ x,
                se = TRUE) +
    labs(x = display_name,
         y = params$y_label,
         title = paste("Relationship between", display_name, "and", params$y_label),
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
  } else if(params$analysis_type == "within_trial_switch") {
    p <- p + facet_grid(switch_vs_stay ~ consensus_level) +
      theme(panel.spacing = unit(1, "lines"))
  } else if(params$analysis_type == "choice_consensus") {
    p <- p + facet_wrap(~consensus_level)
  }
  
  return(p)
}

plot_median_split <- function(data, var, results, median_var, params, display_name) {
  # Check if we have dedicated median split results
  median_results <- results$median_split
  
  # Use the data from median split analysis if available
  if (!is.null(median_results) && !is.null(median_results$data)) {
    data_for_plot <- median_results$data
  } else {
    # Fallback to creating the median split data directly
    data_for_plot <- data %>%
      mutate(quest_group = factor(
        ifelse(!is.na(.data[[median_var]]), 
               ifelse(.data[[median_var]] > median(.data[[median_var]], na.rm = TRUE), "High", "Low"),
               NA),
        levels = c("Low", "High")
      )) %>%
      filter(!is.na(quest_group))
  }
  
  # Get p-value from median split analysis
  median_p <- if(!is.null(median_results) && !is.na(median_results$p_value)) {
    median_results$p_value
  } else {
    NA
  }
  
  # Get adjusted p-value if available
  median_adj_p <- if(!is.null(median_results) && !is.null(median_results$adj_p)) {
    median_results$adj_p
  } else {
    NA
  }
  
  # Create subtitle text using median split p-value if available
  subtitle_text <- if(!is.na(median_p)) {
    if(!is.na(median_adj_p)) {
      paste("Median-split interaction:\nRaw p =", format.pval(median_p, digits = 3),
            "\nFDR-adjusted p =", format.pval(median_adj_p, digits = 3))
    } else {
      paste("Median-split interaction:\np =", format.pval(median_p, digits = 3))
    }
  } else {
    if(!is.null(results$raw_p) && !is.na(results$raw_p)) {
      paste("Three-way interaction:\nRaw p =", format.pval(results$raw_p, digits = 3),
            "\nFDR-adjusted p =", format.pval(results$adj_p, digits = 3))
    } else {
      "P-values not available"
    }
  }
  
  # Create plots using appropriate grouping variables based on analysis type
  if(params$analysis_type == "switch_difference") {
    summary_data <- data_for_plot %>%
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
                    group = direction)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = mean_outcome - se, 
                        ymax = mean_outcome + se), 
                    width = 0.1,
                    linewidth = 0.8) +
      scale_color_manual(values = c("Against group" = "red", "With group" = "blue")) +
      facet_grid(trial_type ~ quest_group) +  
      labs(x = "Consensus Level",
           y = params$y_label,
           title = paste("Effect of", display_name, "(Median-split)"),
           subtitle = subtitle_text) +
      theme_custom
    
  } else if(params$analysis_type == "within_trial_switch") {
    summary_data <- data_for_plot %>%
      group_by(consensus_level, direction, quest_group, switch_vs_stay) %>%
      summarise(
        mean_outcome = mean(outcome_value),
        se = sd(outcome_value) / sqrt(n()),
        .groups = 'drop'
      )
    
    # Create labels for switch/stay
    switch_labels <- c("0" = "Stay trials", "1" = "Switch trials")
    
    p <- ggplot(summary_data, 
                aes(x = consensus_level, 
                    y = mean_outcome, 
                    color = direction,
                    group = direction)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = mean_outcome - se, 
                        ymax = mean_outcome + se), 
                    width = 0.1,
                    linewidth = 0.8) +
      scale_color_manual(values = c("Against group" = "red", "With group" = "blue")) +
      facet_grid(switch_vs_stay ~ quest_group, 
                 labeller = labeller(switch_vs_stay = switch_labels)) +  
      labs(x = "Consensus Level",
           y = params$y_label,
           title = paste("Effect of", display_name, "(Median-split)"),
           subtitle = subtitle_text) +
      theme_custom
    
  } else {  # choice_consensus
    summary_data <- data_for_plot %>%
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
                    group = direction)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = mean_outcome - se, 
                        ymax = mean_outcome + se), 
                    width = 0.1,
                    linewidth = 0.8) +
      scale_color_manual(values = c("Against group" = "red", "With group" = "blue")) +
      facet_wrap(~ quest_group, ncol = 2) +  # Side by side plots for Low/High
      labs(x = "Consensus Level",
           y = params$y_label,
           title = paste("Effect of", display_name, "(Median-split)"),
           subtitle = subtitle_text) +
      theme_custom
  }
  
  return(p)
}

plot_moderation_effects <- function(model, data, var, raw_p = NULL, adj_p = NULL, params, display_name) {
  # Create diagnostic storage list
  diagnostics <- list(
    input_data = data,
    model = model,
    pred_data = NULL,
    predictions = NULL,
    final_plot_data = NULL
  )
  
  subtitle_text <- if(!is.null(raw_p) && !is.na(raw_p)) {
    paste("Interaction test:",
          "\nRaw p =", format.pval(raw_p, digits = 3),
          "\nFDR-adjusted p =", format.pval(adj_p, digits = 3))
  } else {
    "P-values not available"
  }
  
  common_theme <- theme_custom
  color_scheme <- scale_color_manual(values = c("Against group" = "red", "With group" = "blue"))
  
  # Create prediction grid based on analysis type
  if(params$analysis_type == "switch_difference") {
    pred_data <- expand.grid(
      consensus_level = levels(data$consensus_level),
      direction = levels(data$direction),
      switch_difference = c(0, 1),
      scale_name = seq(from = -2, to = 2, length.out = 100),
      age = 0
    ) %>%
      mutate(trial_type = factor(ifelse(switch_difference == 0, "Stay", "Switch"), 
                                 levels = c("Stay", "Switch"))) %>%
      filter(!(consensus_level == "2:2" & direction == "With group"))
    
  } else if(params$analysis_type == "within_trial_switch") {
    pred_data <- expand.grid(
      consensus_level = levels(data$consensus_level),
      direction = levels(data$direction),
      switch_vs_stay = levels(data$switch_vs_stay),
      scale_name = seq(from = -2, to = 2, length.out = 100),
      age = 0
    ) %>%
      filter(!(consensus_level == "2:2" & direction == "With group"))
    
  } else {
    pred_data <- expand.grid(
      consensus_level = levels(data$consensus_level),
      direction = levels(data$direction),
      scale_name = seq(from = -2, to = 2, length.out = 100),
      age = 0
    ) %>%
      filter(!(consensus_level == "2:2" & direction == "With group"))
  }
  
  diagnostics$pred_data <- pred_data
  
  # Generate predictions
  pred_data$predicted <- predict(model, newdata = pred_data, re.form = NA)
  
  diagnostics$predictions <- pred_data$predicted
  diagnostics$final_plot_data <- pred_data
  
  # Create appropriate plot based on analysis type
  if(params$analysis_type == "switch_difference") {
    p <- ggplot(pred_data, 
                aes(x = scale_name, 
                    y = predicted, 
                    color = direction)) +
      geom_line() +
      facet_grid(trial_type ~ consensus_level)
    
  } else if(params$analysis_type == "within_trial_switch") {
    p <- ggplot(pred_data, 
                aes(x = scale_name, 
                    y = predicted, 
                    color = direction)) +
      geom_line() +
      facet_grid(switch_vs_stay ~ consensus_level,
                 labeller = labeller(switch_vs_stay = c("0" = "Stay trials", 
                                                        "1" = "Switch trials")))
    
  } else {
    p <- ggplot(pred_data, 
                aes(x = scale_name, 
                    y = predicted, 
                    color = direction)) +
      geom_line() +
      facet_wrap(~consensus_level)
  }
  
  p <- p +
    labs(title = paste("Moderation effect of", display_name),
         subtitle = subtitle_text,
         x = paste(display_name, "score (standardized)"),
         y = params$y_label) +
    common_theme +
    color_scheme
  
  return(list(
    plot = p,
    predictions = pred_data,
    raw_p = raw_p,
    adj_p = adj_p,
    diagnostics = diagnostics
  ))
}

plot_simple_slopes <- function(model, var, data, raw_p = NULL, adj_p = NULL, params, display_name) {
  # Create diagnostic storage list
  diagnostics <- list(
    input_data = data,
    model = model,
    pred_data = NULL,
    predictions = NULL,
    ribbon_data = NULL
  )
  
  subtitle_text <- if(!is.null(raw_p) && !is.na(raw_p)) {
    paste("Interaction test:",
          "\nRaw p =", format.pval(raw_p, digits = 3),
          "\nFDR-adjusted p =", format.pval(adj_p, digits = 3))
  } else {
    "P-values not available"
  }
  
  common_theme <- theme_custom
  color_scheme <- scale_color_manual(values = c("Against group" = "red", "With group" = "blue"))
  fill_scheme <- scale_fill_manual(values = c("Against group" = "#f36a7b", "With group" = "#3f9fef"))
  
  if(params$analysis_type == "choice_consensus") {
    # Create prediction grid for choice_consensus
    pred_data <- expand.grid(
      scale_name = seq(from = -2, to = 2, length.out = 100),
      consensus_level = levels(data$consensus_level),
      direction = levels(data$direction),
      age = 0
    ) %>%
      filter(!(consensus_level == "2:2" & direction == "With group"))
    
    # Generate predictions
    pred_data$predicted <- predict(model, newdata = pred_data, re.form = NA)
    
    # Create plot
    p <- ggplot(pred_data, 
                aes(x = scale_name, 
                    y = predicted, 
                    color = consensus_level,
                    group = consensus_level)) +
      geom_line(linewidth = 1) +
      facet_wrap(~direction) +
      labs(title = paste("Simple slopes for", display_name),
           subtitle = subtitle_text,
           x = paste(display_name, "score (standardized)"),
           y = params$y_label,
           color = "Consensus Level") +
      common_theme +
      scale_color_manual(values = c("2:2" = "#E69F00", 
                                    "3:1" = "#56B4E9", 
                                    "4:0" = "#009E73"))
    
    diagnostics$pred_data <- pred_data
    diagnostics$predictions <- pred_data$predicted
    
    return(list(
      plot = p,
      diagnostics = diagnostics,
      individual_plots = NULL
    ))
    
  } else if(params$analysis_type == "switch_difference") {
    # Create prediction grid with questionnaire scores as x-axis
    all_pred_data <- list()
    
    for(this_type in unique(data$trial_type)) {
      pred_data <- expand.grid(
        scale_name = seq(from = -2, to = 2, length.out = 100),
        consensus_level = levels(data$consensus_level),
        direction = levels(data$direction),
        trial_type = this_type,
        switch_difference = ifelse(this_type == "Stay", 0, 1),
        age = 0
      ) %>%
        filter(!(consensus_level == "2:2" & direction == "With group"))
      
      # Generate predictions
      pred_data$predicted <- predict(model, newdata = pred_data, re.form = NA)
      
      all_pred_data[[this_type]] <- pred_data
    }
    
    # Combine prediction data
    combined_pred_data <- bind_rows(all_pred_data)
    
    # Create plot
    p <- ggplot(combined_pred_data, 
                aes(x = scale_name, 
                    y = predicted, 
                    color = consensus_level,
                    group = consensus_level)) +
      geom_line(linewidth = 1) +
      facet_grid(trial_type ~ direction) +
      labs(title = paste("Simple slopes for", display_name),
           subtitle = subtitle_text,
           x = paste(display_name, "score (standardized)"),
           y = params$y_label,
           color = "Consensus Level") +
      common_theme +
      scale_color_manual(values = c("2:2" = "#E69F00", 
                                    "3:1" = "#56B4E9", 
                                    "4:0" = "#009E73"))
    
    diagnostics$pred_data <- combined_pred_data
    diagnostics$predictions <- combined_pred_data$predicted
    
    return(list(
      plot = p,
      diagnostics = diagnostics,
      individual_plots = NULL
    ))
    
  } else if(params$analysis_type == "within_trial_switch") {
    # Create prediction grid with questionnaire scores as x-axis
    all_pred_data <- list()
    
    for(trial_type in c("0", "1")) {
      pred_data <- expand.grid(
        scale_name = seq(from = -2, to = 2, length.out = 100),
        consensus_level = levels(data$consensus_level),
        direction = levels(data$direction),
        switch_vs_stay = trial_type,
        age = 0
      ) %>%
        filter(!(consensus_level == "2:2" & direction == "With group"))
      
      # Generate predictions
      pred_data$predicted <- predict(model, newdata = pred_data, re.form = NA)
      
      # Add trial type labels
      pred_data$trial_type <- ifelse(trial_type == "0", "Stay trials", "Switch trials")
      
      all_pred_data[[trial_type]] <- pred_data
    }
    
    # Combine prediction data
    combined_pred_data <- bind_rows(all_pred_data)
    
    # Create plot
    p <- ggplot(combined_pred_data, 
                aes(x = scale_name, 
                    y = predicted, 
                    color = consensus_level,
                    group = consensus_level)) +
      geom_line(linewidth = 1) +
      facet_grid(trial_type ~ direction) +
      labs(title = paste("Simple slopes for", display_name),
           subtitle = subtitle_text,
           x = paste(display_name, "score (standardized)"),
           y = params$y_label,
           color = "Consensus Level") +
      common_theme +
      scale_color_manual(values = c("2:2" = "#E69F00", 
                                    "3:1" = "#56B4E9", 
                                    "4:0" = "#009E73"))
    
    diagnostics$pred_data <- combined_pred_data
    diagnostics$predictions <- combined_pred_data$predicted
  }
  
  return(list(
    plot = p,
    diagnostics = diagnostics,
    individual_plots = NULL
  ))
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
  } else if(params$analysis_type == "within_trial_switch") {
    # Additional checks specific to within-trial switch analysis
    within_trial_diagnostics <- tryCatch({
      list(
        switch_stay_counts = table(model_frame$switch_vs_stay),
        switch_stay_distribution = chisq.test(table(model_frame$switch_vs_stay))
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
  } else if(params$analysis_type == "within_trial_switch" && !is.null(within_trial_diagnostics)) {
    diagnostics$within_trial <- within_trial_diagnostics
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
  model_diagnostics <- check_model_diagnostics(mod_model, params)
  
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
  
  # Empty subscale results list
  subscale_results <- list()
  
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
    main_terms <- if(params$analysis_type == "switch_difference") {
      c("switch_difference", "consensus_level", "direction", "scale_name")
    } else if(params$analysis_type == "within_trial_switch") {
      c("switch_vs_stay", "consensus_level", "direction", "scale_name")
    } else {
      c("consensus_level", "direction", "scale_name")
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
  
  # Add analysis-specific diagnostics
  if(params$analysis_type == "switch_difference" && !is.null(diag$switch_difference)) {
    output_text <- paste0(output_text,
                          "Switch difference range: ", 
                          paste(round(diag$switch_difference$switch_diff_range, 3), collapse = " to "), "\n",
                          "Switch difference normality p: ",
                          format.pval(diag$switch_difference$switch_diff_distribution$p.value, digits = 3), "\n")
  } else if(params$analysis_type == "within_trial_switch" && !is.null(diag$within_trial)) {
    output_text <- paste0(output_text,
                          "Switch vs Stay counts:\n",
                          paste(capture.output(diag$within_trial$switch_stay_counts), collapse = "\n"), "\n",
                          "Switch vs Stay distribution chi-square p: ",
                          format.pval(diag$within_trial$switch_stay_distribution$p.value, digits = 3), "\n")
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
      
    } else if(params$analysis_type == "within_trial_switch") {
      output_text <- paste0(output_text,
                            "\nINTERACTION RESULTS:\n===================\n")
      # Extract p-values from results structure
      raw_p <- results$moderation_results$main_results$raw_p
      adj_p <- results$moderation_results$main_results$adj_p
      
      output_text <- paste0(output_text,
                            "Within-trial switch x Consensus x Direction interaction:\n",
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
    
    # Add section for continuous model interaction
    output_text <- paste0(output_text,
                          "\nCONTINUOUS MODEL INTERACTION:\n",
                          "=============================\n",
                          "Raw p = ", format.pval(results$moderation_results$main_results$raw_p, digits = 3), "\n",
                          "FDR-adjusted p = ", format.pval(results$moderation_results$main_results$adj_p, digits = 3),
                          "\n\n")
    
    # Add section for median split results
    median_results <- results$moderation_results$main_results$median_split
    if(!is.null(median_results) && !is.na(median_results$p_value)) {
      output_text <- paste0(output_text,
                            "MEDIAN SPLIT INTERACTION:\n",
                            "==========================\n",
                            "Raw p = ", format.pval(median_results$p_value, digits = 3), "\n")
      
      if(!is.null(median_results$adj_p)) {
        output_text <- paste0(output_text,
                              "FDR-adjusted p = ", format.pval(median_results$adj_p, digits = 3), "\n")
      }
      output_text <- paste0(output_text, "\n")
    }
    
    # Add model diagnostics
    output_text <- paste0(output_text, 
                          "\n", format_diagnostics(results$moderation_results$main_results$diagnostics))
    
    # No subscale results section
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
  # Since we're not analyzing subscales anymore, return empty string
  return("")
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

################## REVERSAL ANALYSIS FUNCTIONS ###################
plot_median_split_reversal <- function(data, var, results, params, display_name) {
  
  # Calculate median directly
  median_value <- median(data[[var]], na.rm = TRUE)
  
  # Determine if this is choice or bet analysis based on the data
  analysis_type <- if(all(c("player.bet1", "player.bet2") %in% names(data))) "BET ANALYSIS" else "CHOICE ANALYSIS"
  
  # Create variable-specific group column name
  group_col_name <- paste0("quest_group_", var)
  
  # Create groups with variable-specific column using !!sym() 
  data <- data %>%
    # First create a temporary column with current variable's values
    mutate(
      temp_values = !!sym(var),
      !!group_col_name := factor(
        if_else(
          is.na(temp_values), 
          NA_character_,
          if_else(temp_values > median_value, "High", "Low")
        ),
        levels = c("Low", "High")
      )
    ) %>%
    # Remove temporary column
    select(-temp_values)
  
  # Get and print participant IDs for each group
  grouping_info <- data %>%
    select(participant.id_in_session, !!sym(group_col_name)) %>%
    distinct()
  
  # Modified summary calculation using variable-specific group column
  summary_data <- data %>%
    filter(!is.na(.data[[group_col_name]])) %>%
    group_by(trial_to_reversal, !!sym(group_col_name)) %>%
    {
      # Debug: Print grouped data before summarizing
      print(head(.))
      .
    } %>%
    summarise(
      # Debug: Print group sizes
      n_obs = n(),
      
      # Choice 1
      choice1_mean = mean(player.choice1_accuracy, na.rm = TRUE) * 100,
      choice1_sd = sd(player.choice1_accuracy, na.rm = TRUE) * 100,
      choice1_n = sum(!is.na(player.choice1_accuracy)),
      choice1_se = (choice1_sd / sqrt(choice1_n)),
      
      # Choice 2 
      choice2_mean = mean(player.choice2_accuracy, na.rm = TRUE) * 100,
      choice2_sd = sd(player.choice2_accuracy, na.rm = TRUE) * 100,
      choice2_n = sum(!is.na(player.choice2_accuracy)),
      choice2_se = (choice2_sd / sqrt(choice2_n)),
      
      # Bet 1
      bet1_mean = mean(player.bet1, na.rm = TRUE),
      bet1_sd = sd(player.bet1, na.rm = TRUE),
      bet1_n = sum(!is.na(player.bet1)),
      bet1_se = (bet1_sd / sqrt(bet1_n)),
      
      # Bet 2
      bet2_mean = mean(player.bet2, na.rm = TRUE),
      bet2_sd = sd(player.bet2, na.rm = TRUE),
      bet2_n = sum(!is.na(player.bet2)),
      bet2_se = (bet2_sd / sqrt(bet2_n)),
      
      .groups = 'drop'
    )
  
  # Run t-tests with effect sizes
  test_results <- data %>%
    group_by(trial_to_reversal) %>%
    summarise(
      # Choice 1
      choice1_t = tryCatch({
        t.test(player.choice1_accuracy ~ .data[[group_col_name]])$statistic
      }, error = function(e) NA),
      choice1_p = tryCatch({
        t.test(player.choice1_accuracy ~ .data[[group_col_name]])$p.value
      }, error = function(e) NA),
      choice1_effect = tryCatch({
        fisher_result <- fisher.test(table(.data[[group_col_name]], player.choice1_accuracy))
        sqrt(log(fisher_result$estimate)^2 / sum(table(.data[[group_col_name]], player.choice1_accuracy)))
      }, error = function(e) NA),
      
      # Choice 2
      choice2_t = tryCatch({
        t.test(player.choice2_accuracy ~ .data[[group_col_name]])$statistic
      }, error = function(e) NA),
      choice2_p = tryCatch({
        t.test(player.choice2_accuracy ~ .data[[group_col_name]])$p.value
      }, error = function(e) NA),
      choice2_effect = tryCatch({
        fisher_result <- fisher.test(table(.data[[group_col_name]], player.choice2_accuracy))
        sqrt(log(fisher_result$estimate)^2 / sum(table(.data[[group_col_name]], player.choice2_accuracy)))
      }, error = function(e) NA),
      
      # Bet 1
      bet1_t = tryCatch({
        t.test(player.bet1 ~ .data[[group_col_name]])$statistic
      }, error = function(e) NA),
      bet1_p = tryCatch({
        t.test(player.bet1 ~ .data[[group_col_name]])$p.value
      }, error = function(e) NA),
      bet1_d = tryCatch({
        effectsize::cohens_d(player.bet1 ~ .data[[group_col_name]])$Cohens_d
      }, error = function(e) NA),
      
      # Bet 2
      bet2_t = tryCatch({
        t.test(player.bet2 ~ .data[[group_col_name]])$statistic
      }, error = function(e) NA),
      bet2_p = tryCatch({
        t.test(player.bet2 ~ .data[[group_col_name]])$p.value
      }, error = function(e) NA),
      bet2_d = tryCatch({
        effectsize::cohens_d(player.bet2 ~ .data[[group_col_name]])$Cohens_d
      }, error = function(e) NA),
      .groups = 'drop'
    ) %>%
    mutate(across(ends_with("_p"), 
                  list(adj = ~p.adjust(., method = "fdr")), 
                  .names = "{.col}_adj"))
  
  # Function to get significance stars (unchanged)
  get_stars <- function(p) {
    if(p < 0.0001) return("****")
    if(p < 0.001) return("***")
    if(p < 0.01) return("**")
    if(p < 0.05) return("*")
    return("")
  }
  
  # Base plot theme (unchanged)
  theme_custom_no_legend <- theme_custom + theme(legend.position = "none")
  
  # Choice plots
  p_choice1 <- ggplot(summary_data, aes(x = trial_to_reversal, group = !!sym(group_col_name))) +
    geom_line(aes(y = choice1_mean, color = !!sym(group_col_name)), linewidth = 1) +
    geom_errorbar(aes(ymin = choice1_mean - choice1_se, 
                      ymax = choice1_mean + choice1_se, 
                      color = !!sym(group_col_name)), 
                  width = 0.2, alpha = 0.5) +
    scale_color_manual(values = c("Low" = "blue", "High" = "red")) +
    labs(title = "Choice 1 Accuracy",
         x = "Trial relative to reversal",
         y = "Accuracy (%)",
         color = "Group") +
    theme_custom_no_legend
  
  p_choice2 <- ggplot(summary_data, aes(x = trial_to_reversal, group = !!sym(group_col_name))) +
    geom_line(aes(y = choice2_mean, color = !!sym(group_col_name)), linewidth = 1) +
    geom_errorbar(aes(ymin = choice2_mean - choice2_se, 
                      ymax = choice2_mean + choice2_se, 
                      color = !!sym(group_col_name)), 
                  width = 0.2, alpha = 0.5) +
    scale_color_manual(values = c("Low" = "blue", "High" = "red")) +
    labs(title = "Choice 2 Accuracy",
         x = "Trial relative to reversal",
         y = "Accuracy (%)",
         color = "Group") +
    theme_custom_no_legend
  
  # Bet plots
  p_bet1 <- ggplot(summary_data, aes(x = trial_to_reversal, group = !!sym(group_col_name))) +
    geom_line(aes(y = bet1_mean, color = !!sym(group_col_name)), linewidth = 1) +
    geom_errorbar(aes(ymin = bet1_mean - bet1_se, 
                      ymax = bet1_mean + bet1_se, 
                      color = !!sym(group_col_name)), 
                  width = 0.2, alpha = 0.5) +
    scale_color_manual(values = c("Low" = "blue", "High" = "red")) +
    labs(title = "Bet 1 Magnitude",
         x = "Trial relative to reversal",
         y = "Bet magnitude",
         color = "Group") +
    theme_custom_no_legend
  
  p_bet2 <- ggplot(summary_data, aes(x = trial_to_reversal, group = !!sym(group_col_name))) +
    geom_line(aes(y = bet2_mean, color = !!sym(group_col_name)), linewidth = 1) +
    geom_errorbar(aes(ymin = bet2_mean - bet2_se, 
                      ymax = bet2_mean + bet2_se, 
                      color = !!sym(group_col_name)), 
                  width = 0.2, alpha = 0.5) +
    scale_color_manual(values = c("Low" = "blue", "High" = "red")) +
    labs(title = "Bet 2 Magnitude",
         x = "Trial relative to reversal",
         y = "Bet magnitude",
         color = "Group") +
    theme_custom_no_legend
  
  # Add significance stars for all plots
  for(i in 1:nrow(test_results)) {
    # Choice 1
    if(!is.na(test_results$choice1_p_adj[i]) && test_results$choice1_p_adj[i] < 0.05) {
      stars <- get_stars(test_results$choice1_p_adj[i])
      p_choice1 <- p_choice1 + 
        annotate("text", x = test_results$trial_to_reversal[i], 
                 y = max(summary_data$choice1_mean) + 5, 
                 label = stars, size = 5)
    }
    
    # Choice 2
    if(!is.na(test_results$choice2_p_adj[i]) && test_results$choice2_p_adj[i] < 0.05) {
      stars <- get_stars(test_results$choice2_p_adj[i])
      p_choice2 <- p_choice2 + 
        annotate("text", x = test_results$trial_to_reversal[i], 
                 y = max(summary_data$choice2_mean) + 5, 
                 label = stars, size = 5)
    }
    
    # Bet 1
    if(!is.na(test_results$bet1_p_adj[i]) && test_results$bet1_p_adj[i] < 0.05) {
      stars <- get_stars(test_results$bet1_p_adj[i])
      p_bet1 <- p_bet1 + 
        annotate("text", x = test_results$trial_to_reversal[i], 
                 y = max(summary_data$bet1_mean) + 0.5, 
                 label = stars, size = 5)
    }
    
    # Bet 2
    if(!is.na(test_results$bet2_p_adj[i]) && test_results$bet2_p_adj[i] < 0.05) {
      stars <- get_stars(test_results$bet2_p_adj[i])
      p_bet2 <- p_bet2 + 
        annotate("text", x = test_results$trial_to_reversal[i], 
                 y = max(summary_data$bet2_mean) + 0.5, 
                 label = stars, size = 5)
    }
  }
  
  # Create legend plot
  legend_plot <- ggplot(summary_data, aes(x = trial_to_reversal, group = !!sym(group_col_name))) +
    geom_line(aes(y = choice1_mean, color = !!sym(group_col_name))) +
    scale_color_manual(values = c("Low" = "blue", "High" = "red")) +
    labs(color = "Group") +
    theme_custom
  legend <- cowplot::get_legend(legend_plot)
  
  # Combine plots with side legend
  p_choice <- gridExtra::grid.arrange(
    gridExtra::arrangeGrob(p_choice1, p_choice2, ncol = 1),
    legend,
    ncol = 2,
    widths = c(4, 1),
    top = grid::textGrob(paste("Choice Accuracy by", display_name, "(Median-split)"),
                         gp = grid::gpar(fontsize = 12, fontface = "bold"))
  )
  
  p_bet <- gridExtra::grid.arrange(
    gridExtra::arrangeGrob(p_bet1, p_bet2, ncol = 1),
    legend,
    ncol = 2,
    widths = c(4, 1),
    top = grid::textGrob(paste("Bet Magnitude by", display_name, "(Median-split)"),
                         gp = grid::gpar(fontsize = 12, fontface = "bold"))
  )
  
  # Debug: Print final return object structure
  cat("\nReturning results for", var, "\n")
  
  return(list(
    choice_plot = p_choice,
    bet_plot = p_bet,
    test_results = test_results,
    summary_data = summary_data
  ))
}

# Reversal results have a different format
format_results_reversal <- function(test_results, summary_data, var_name) {
  # Create variable-specific group column name
  group_col_name <- paste0("quest_group_", var_name)
  
  choice_text <- sprintf("CHOICE ACCURACY ANALYSIS FOR %s\n\n", toupper(var_name))
  choice_text <- paste0(choice_text, "Analysis run on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # Add group sizes using the variable-specific group column
  choice_text <- paste0(choice_text, "GROUP SIZES:\n",
                        "Low group n = ", sum(summary_data[[group_col_name]] == "Low")/7, "\n",
                        "High group n = ", sum(summary_data[[group_col_name]] == "High")/7, "\n\n",
                        "SUMMARY OF SIGNIFICANT DIFFERENCES:\n",
                        "================================\n")
  
  # Add significant results summary for choices
  choice1_sig <- test_results %>% 
    filter(choice1_p_adj < 0.05) %>%
    arrange(trial_to_reversal)
  
  choice2_sig <- test_results %>% 
    filter(choice2_p_adj < 0.05) %>%
    arrange(trial_to_reversal)
  
  if(nrow(choice1_sig) > 0) {
    choice_text <- paste0(choice_text, "\nChoice 1 significant differences found at trials:\n")
    for(i in 1:nrow(choice1_sig)) {
      # Get detailed statistics for this trial
      trial_stats <- summary_data %>%
        filter(trial_to_reversal == choice1_sig$trial_to_reversal[i]) %>%
        select(!!sym(group_col_name), choice1_mean, choice1_sd, choice1_n) %>%
        rename(
          mean_acc = choice1_mean,
          sd_acc = choice1_sd,
          n = choice1_n
        )
      
      # Calculate confidence interval
      df <- sum(trial_stats$n) - 2
      t_value <- qt(0.975, df = df)
      pooled_sd <- sqrt(((trial_stats$n[1] - 1) * trial_stats$sd_acc[1]^2 + 
                           (trial_stats$n[2] - 1) * trial_stats$sd_acc[2]^2) / 
                          (df))
      se_diff <- pooled_sd * sqrt(1/trial_stats$n[1] + 1/trial_stats$n[2])
      mean_diff <- trial_stats$mean_acc[1] - trial_stats$mean_acc[2]
      ci_lower <- mean_diff - t_value * se_diff
      ci_upper <- mean_diff + t_value * se_diff
      
      choice_text <- paste0(choice_text,
                            "Trial ", choice1_sig$trial_to_reversal[i], ":\n",
                            "  High group: ", round(trial_stats$mean_acc[1], 1), 
                            "% (SD = ", round(trial_stats$sd_acc[1], 1), ")\n",
                            "  Low group: ", round(trial_stats$mean_acc[2], 1), 
                            "% (SD = ", round(trial_stats$sd_acc[2], 1), ")\n",
                            "  Mean difference: ", round(mean_diff, 1),
                            "% [", round(ci_lower, 1), ", ", round(ci_upper, 1), "]\n",
                            "  Effect size (phi) = ", round(choice1_sig$choice1_effect[i], 3), "\n",
                            "  t(", df, ") = ", round(choice1_sig$choice1_t[i], 3), 
                            ", adj-p = ", format.pval(choice1_sig$choice1_p_adj[i], digits = 3), "\n\n")
    }
  }
  
  if(nrow(choice2_sig) > 0) {
    choice_text <- paste0(choice_text, "\nChoice 2 significant differences found at trials:\n")
    for(i in 1:nrow(choice2_sig)) {
      trial_stats <- summary_data %>%
        filter(trial_to_reversal == choice2_sig$trial_to_reversal[i]) %>%
        select(!!sym(group_col_name), choice2_mean, choice2_sd, choice2_n) %>%
        rename(
          mean_acc = choice2_mean,
          sd_acc = choice2_sd,
          n = choice2_n
        )
      
      # Calculate confidence interval
      df <- sum(trial_stats$n) - 2
      t_value <- qt(0.975, df = df)
      pooled_sd <- sqrt(((trial_stats$n[1] - 1) * trial_stats$sd_acc[1]^2 + 
                           (trial_stats$n[2] - 1) * trial_stats$sd_acc[2]^2) / 
                          (df))
      se_diff <- pooled_sd * sqrt(1/trial_stats$n[1] + 1/trial_stats$n[2])
      mean_diff <- trial_stats$mean_acc[1] - trial_stats$mean_acc[2]
      ci_lower <- mean_diff - t_value * se_diff
      ci_upper <- mean_diff + t_value * se_diff
      
      choice_text <- paste0(choice_text,
                            "Trial ", choice2_sig$trial_to_reversal[i], ":\n",
                            "  High group: ", round(trial_stats$mean_acc[1], 1), 
                            "% (SD = ", round(trial_stats$sd_acc[1], 1), ")\n",
                            "  Low group: ", round(trial_stats$mean_acc[2], 1), 
                            "% (SD = ", round(trial_stats$sd_acc[2], 1), ")\n",
                            "  Mean difference: ", round(mean_diff, 1),
                            "% [", round(ci_lower, 1), ", ", round(ci_upper, 1), "]\n",
                            "  Effect size (phi) = ", round(choice2_sig$choice2_effect[i], 3), "\n",
                            "  t(", df, ") = ", round(choice2_sig$choice2_t[i], 3), 
                            ", adj-p = ", format.pval(choice2_sig$choice2_p_adj[i], digits = 3), "\n\n")
    }
  }
  
  bet_text <- sprintf("BET MAGNITUDE ANALYSIS FOR %s\n\n", toupper(var_name))
  bet_text <- paste0(bet_text, "Analysis run on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  bet_text <- paste0(bet_text, "GROUP SIZES:\n",
                     "Low group n = ", sum(summary_data[[group_col_name]] == "Low")/7, "\n",
                     "High group n = ", sum(summary_data[[group_col_name]] == "High")/7, "\n\n",
                     "SUMMARY OF SIGNIFICANT DIFFERENCES:\n",
                     "================================\n")
  
  bet1_sig <- test_results %>% 
    filter(bet1_p_adj < 0.05) %>%
    arrange(trial_to_reversal)
  
  bet2_sig <- test_results %>% 
    filter(bet2_p_adj < 0.05) %>%
    arrange(trial_to_reversal)
  
  if(nrow(bet1_sig) > 0) {
    bet_text <- paste0(bet_text, "\nBet 1 significant differences found at trials:\n")
    for(i in 1:nrow(bet1_sig)) {
      trial_stats <- summary_data %>%
        filter(trial_to_reversal == bet1_sig$trial_to_reversal[i]) %>%
        select(!!sym(group_col_name), bet1_mean, bet1_sd, bet1_n) %>%
        rename(
          mean_bet = bet1_mean,
          sd_bet = bet1_sd,
          n = bet1_n
        )
      
      # Calculate confidence interval
      df <- sum(trial_stats$n) - 2
      t_value <- qt(0.975, df = df)
      pooled_sd <- sqrt(((trial_stats$n[1] - 1) * trial_stats$sd_bet[1]^2 + 
                           (trial_stats$n[2] - 1) * trial_stats$sd_bet[2]^2) / 
                          (df))
      se_diff <- pooled_sd * sqrt(1/trial_stats$n[1] + 1/trial_stats$n[2])
      mean_diff <- trial_stats$mean_bet[1] - trial_stats$mean_bet[2]
      ci_lower <- mean_diff - t_value * se_diff
      ci_upper <- mean_diff + t_value * se_diff
      
      bet_text <- paste0(bet_text,
                         "Trial ", bet1_sig$trial_to_reversal[i], ":\n",
                         "  High group: ", round(trial_stats$mean_bet[1], 2), 
                         " (SD = ", round(trial_stats$sd_bet[1], 2), ")\n",
                         "  Low group: ", round(trial_stats$mean_bet[2], 2), 
                         " (SD = ", round(trial_stats$sd_bet[2], 2), ")\n",
                         "  Mean difference: ", round(mean_diff, 2),
                         " [", round(ci_lower, 2), ", ", round(ci_upper, 2), "]\n",
                         "  Effect size (Cohen's d) = ", round(bet1_sig$bet1_d[i], 3), "\n",
                         "  t(", df, ") = ", round(bet1_sig$bet1_t[i], 3), 
                         ", adj-p = ", format.pval(bet1_sig$bet1_p_adj[i], digits = 3), "\n\n")
    }
  }
  
  if(nrow(bet2_sig) > 0) {
    bet_text <- paste0(bet_text, "\nBet 2 significant differences found at trials:\n")
    for(i in 1:nrow(bet2_sig)) {
      trial_stats <- summary_data %>%
        filter(trial_to_reversal == bet2_sig$trial_to_reversal[i]) %>%
        select(!!sym(group_col_name), bet2_mean, bet2_sd, bet2_n) %>%
        rename(
          mean_bet = bet2_mean,
          sd_bet = bet2_sd,
          n = bet2_n
        )
      
      # Calculate confidence interval
      df <- sum(trial_stats$n) - 2
      t_value <- qt(0.975, df = df)
      pooled_sd <- sqrt(((trial_stats$n[1] - 1) * trial_stats$sd_bet[1]^2 + 
                           (trial_stats$n[2] - 1) * trial_stats$sd_bet[2]^2) / 
                          (df))
      se_diff <- pooled_sd * sqrt(1/trial_stats$n[1] + 1/trial_stats$n[2])
      mean_diff <- trial_stats$mean_bet[1] - trial_stats$mean_bet[2]
      ci_lower <- mean_diff - t_value * se_diff
      ci_upper <- mean_diff + t_value * se_diff
      
      bet_text <- paste0(bet_text,
                         "Trial ", bet2_sig$trial_to_reversal[i], ":\n",
                         "  High group: ", round(trial_stats$mean_bet[1], 2), 
                         " (SD = ", round(trial_stats$sd_bet[1], 2), ")\n",
                         "  Low group: ", round(trial_stats$mean_bet[2], 2), 
                         " (SD = ", round(trial_stats$sd_bet[2], 2), ")\n",
                         "  Mean difference: ", round(mean_diff, 2),
                         " [", round(ci_lower, 2), ", ", round(ci_upper, 2), "]\n",
                         "  Effect size (Cohen's d) = ", round(bet2_sig$bet2_d[i], 3), "\n",
                         "  t(", df, ") = ", round(bet2_sig$bet2_t[i], 3), 
                         ", adj-p = ", format.pval(bet2_sig$bet2_p_adj[i], digits = 3), "\n\n")
    }
  }
  
  return(list(choice_text = choice_text, bet_text = bet_text))
}

# Specific formatting function for bet and choice accuracy differences
format_results_overall <- function(main_results, subscale_results = NULL, params) {
  output_text <- sprintf("%s - STATISTICAL ANALYSIS RESULTS\n\n", 
                         toupper(params$analysis_name))
  output_text <- paste0(output_text, 
                        "Analysis run on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # Format main results
  output_text <- paste0(output_text, "PRIMARY ANALYSIS RESULTS:\n",
                        "================================\n")
  
  for(i in 1:nrow(main_results)) {
    result <- main_results[i,]
    output_text <- paste0(output_text,
                          sprintf("\nQuestionnaire: %s\n", result$questionnaire),
                          "\nChoice Analysis:\n",
                          sprintf("p-value: %s\n", format.pval(result$choice_p, digits = 3)),
                          sprintf("FDR-adjusted p-value: %s\n", format.pval(result$choice_p_adj, digits = 3)),
                          sprintf("Correlation (r): %.3f [%.3f, %.3f]\n", 
                                  result$choice_r, result$choice_r_ci_low, result$choice_r_ci_high),
                          sprintf("Standardized beta: %.3f\n", result$choice_beta),
                          sprintf("R²: %.3f\n", result$choice_r2),
                          sprintf("F-statistic: F(1,%d) = %.2f\n", result$choice_df - 2, result$choice_f),
                          sprintf("Residual normality: %s\n", format.pval(result$choice_resid_norm, digits = 3)),
                          "\nBet Analysis:\n",
                          sprintf("p-value: %s\n", format.pval(result$bet_p, digits = 3)),
                          sprintf("FDR-adjusted p-value: %s\n", format.pval(result$bet_p_adj, digits = 3)),
                          sprintf("Correlation (r): %.3f [%.3f, %.3f]\n", 
                                  result$bet_r, result$bet_r_ci_low, result$bet_r_ci_high),
                          sprintf("Standardized beta: %.3f\n", result$bet_beta),
                          sprintf("R²: %.3f\n", result$bet_r2),
                          sprintf("F-statistic: F(1,%d) = %.2f\n", result$bet_df - 2, result$bet_f),
                          sprintf("Residual normality: %s\n\n", format.pval(result$bet_resid_norm, digits = 3)))
  }
  
  # No subscale results section
  
  return(output_text)
}

################## POST-HOC TESTING FUNCTIONS ###################
# Function to run post-hoc tests for significant relationships
run_posthoc_tests <- function(model, data, var_name, analysis_type, output_dir) {
  # Create more descriptive analysis name
  analysis_name <- if(analysis_type == "choice_consensus") {
    "choice_switch_by_consensus"
  } else if(analysis_type == "within_trial_switch") {
    "bet_difference_by_switch"
  } else {
    analysis_type
  }
  
  # Create output directory with clearer naming
  posthoc_dir <- file.path(output_dir, "posthoc", analysis_name, var_name, "txt")
  dir.create(posthoc_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get total participant count
  total_participants <- length(unique(data$participant.id_in_session))
  
  # Initialize results storage
  results_text <- sprintf("POST-HOC ANALYSIS: %s - %s\n\n", toupper(var_name), toupper(analysis_name))
  results_text <- paste0(results_text, "Analysis run on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  results_text <- paste0(results_text, "Total participants in dataset: ", total_participants, "\n\n")
  
  # Initialize storage for statistics
  p_values <- numeric()
  test_details <- character()
  test_results <- list()
  participants_per_condition <- list()
  
  # Run tests based on analysis type
  if(analysis_type == "choice_consensus") {
    results_text <- paste0(results_text, "SIMPLE SLOPES ANALYSIS BY CONSENSUS LEVEL:\n========================================\n\n")
    
    # For each consensus level
    for(cons in c("2:2", "3:1", "4:0")) {
      # For each direction
      for(dir in c("Against group", "With group")) {
        # Skip invalid combination
        if(cons == "2:2" && dir == "With group") next
        
        # Create model for this specific combination
        subset_data <- data %>% 
          filter(consensus_level == cons, direction == dir)
        
        # Count unique participants in this condition
        participant_ids <- unique(subset_data$participant.id_in_session)
        n_participants <- length(participant_ids)
        n_trials <- nrow(subset_data)
        participants_per_condition[[paste(cons, dir)]] <- participant_ids
        
        # Only proceed if we have enough data
        if(nrow(subset_data) < 30) {
          results_text <- paste0(results_text, 
                                 sprintf("Insufficient data for %s, %s (n=%d trials from %d participants)\n\n", 
                                         cons, dir, n_trials, n_participants))
          next
        }
        
        # Run the simple model
        simple_model <- tryCatch({
          lm(outcome_value ~ scale_name + age, data = subset_data)
        }, error = function(e) {
          return(NULL)
        })
        
        if(is.null(simple_model)) {
          results_text <- paste0(results_text, 
                                 sprintf("Failed to fit model for %s, %s\n\n", cons, dir))
          next
        }
        
        # Get model results
        model_summary <- summary(simple_model)
        
        # Get correlation 
        cor_test <- cor.test(subset_data$scale_name, subset_data$outcome_value)
        
        # Get key statistics
        slope <- model_summary$coefficients["scale_name", "Estimate"]
        slope_se <- model_summary$coefficients["scale_name", "Std. Error"]
        slope_t <- model_summary$coefficients["scale_name", "t value"]
        slope_p <- model_summary$coefficients["scale_name", "Pr(>|t|)"]
        slope_df <- model_summary$df[2]
        r_value <- cor_test$estimate
        r_ci_low <- cor_test$conf.int[1]
        r_ci_high <- cor_test$conf.int[2]
        
        # Get standardized beta
        beta <- lm.beta(simple_model)$standardized.coefficients["scale_name"]
        
        # Store p-value and test details for correction
        p_values <- c(p_values, slope_p)
        test_details <- c(test_details, sprintf("Consensus level: %s, Direction: %s", cons, dir))
        
        # Store test results in a named list element
        test_key <- paste(cons, dir, sep="_")
        test_results[[test_key]] <- list(
          p_value = slope_p,
          r_value = r_value,
          r_ci_low = r_ci_low,
          r_ci_high = r_ci_high,
          beta = beta,
          slope = slope,
          t_value = slope_t,
          df = slope_df,
          n = nrow(subset_data),
          n_participants = n_participants
        )
        
        # Add to results
        results_text <- paste0(results_text,
                               sprintf("Consensus level: %s, Direction: %s\n", cons, dir),
                               sprintf("Sample size: %d trials from %d unique participants\n", 
                                       n_trials, n_participants),
                               sprintf("Correlation: r = %.3f [%.3f, %.3f]\n", r_value, r_ci_low, r_ci_high),
                               sprintf("Slope: β = %.3f (SE = %.3f)\n", slope, slope_se),
                               sprintf("Standardized beta: %.3f\n", beta),
                               sprintf("t(%d) = %.3f, p = %.4f\n\n", 
                                       slope_df, slope_t, slope_p))
      }
    }
  } else if(analysis_type == "within_trial_switch") {
    results_text <- paste0(results_text, "SIMPLE SLOPES ANALYSIS BY CONSENSUS LEVEL AND SWITCH STATUS:\n=========================================================\n\n")
    
    # For each consensus level
    for(cons in c("2:2", "3:1", "4:0")) {
      # For each direction
      for(dir in c("Against group", "With group")) {
        # Skip invalid combination
        if(cons == "2:2" && dir == "With group") next
        
        # For each switch status
        for(switch_status in c("0", "1")) {
          # Create model for this specific combination
          subset_data <- data %>% 
            filter(consensus_level == cons, 
                   direction == dir,
                   switch_vs_stay == switch_status)
          
          # Count unique participants in this condition
          participant_ids <- unique(subset_data$participant.id_in_session)
          n_participants <- length(participant_ids)
          n_trials <- nrow(subset_data)
          condition_key <- paste(cons, dir, switch_status, sep="_")
          participants_per_condition[[condition_key]] <- participant_ids
          
          # Only proceed if we have enough data
          if(nrow(subset_data) < 30) {
            results_text <- paste0(results_text, 
                                   sprintf("Insufficient data for %s, %s, switch=%s (n=%d trials from %d participants)\n\n", 
                                           cons, dir, switch_status, n_trials, n_participants))
            next
          }
          
          # Run the simple model
          simple_model <- tryCatch({
            lm(outcome_value ~ scale_name + age, data = subset_data)
          }, error = function(e) {
            return(NULL)
          })
          
          if(is.null(simple_model)) {
            results_text <- paste0(results_text, 
                                   sprintf("Failed to fit model for %s, %s, switch=%s\n\n", 
                                           cons, dir, switch_status))
            next
          }
          
          # Get model results
          model_summary <- summary(simple_model)
          
          # Get correlation 
          cor_test <- cor.test(subset_data$scale_name, subset_data$outcome_value)
          
          # Get key statistics
          slope <- model_summary$coefficients["scale_name", "Estimate"]
          slope_se <- model_summary$coefficients["scale_name", "Std. Error"]
          slope_t <- model_summary$coefficients["scale_name", "t value"]
          slope_p <- model_summary$coefficients["scale_name", "Pr(>|t|)"]
          slope_df <- model_summary$df[2]
          r_value <- cor_test$estimate
          r_ci_low <- cor_test$conf.int[1]
          r_ci_high <- cor_test$conf.int[2]
          
          # Get standardized beta
          beta <- lm.beta(simple_model)$standardized.coefficients["scale_name"]
          
          # Store p-value and test details for correction
          p_values <- c(p_values, slope_p)
          switch_label <- ifelse(switch_status == "1", "Switch", "Stay")
          test_details <- c(test_details, 
                            sprintf("Consensus: %s, Direction: %s, Trial: %s", 
                                    cons, dir, switch_label))
          
          # Store test results in a named list element
          test_key <- paste(cons, dir, switch_status, sep="_")
          test_results[[test_key]] <- list(
            p_value = slope_p,
            r_value = r_value,
            r_ci_low = r_ci_low,
            r_ci_high = r_ci_high,
            beta = beta,
            slope = slope,
            t_value = slope_t,
            df = slope_df,
            n = nrow(subset_data),
            n_participants = n_participants
          )
          
          # Add to results
          results_text <- paste0(results_text,
                                 sprintf("Consensus level: %s, Direction: %s, Trial type: %s\n", 
                                         cons, dir, switch_label),
                                 sprintf("Sample size: %d trials from %d unique participants\n", 
                                         n_trials, n_participants),
                                 sprintf("Correlation: r = %.3f [%.3f, %.3f]\n", r_value, r_ci_low, r_ci_high),
                                 sprintf("Slope: β = %.3f (SE = %.3f)\n", slope, slope_se),
                                 sprintf("Standardized beta: %.3f\n", beta),
                                 sprintf("t(%d) = %.3f, p = %.4f\n\n", 
                                         slope_df, slope_t, slope_p))
        }
      }
    }
  }
  
  # Apply FDR correction and add to results
  adj_p_values <- p.adjust(p_values, method = "fdr")
  
  results_text <- paste0(results_text, 
                         "FDR CORRECTED P-VALUES:\n========================\n")
  
  for(i in seq_along(p_values)) {
    results_text <- paste0(results_text,
                           sprintf("%s: Original p = %.4f, FDR-adjusted p = %.4f\n", 
                                   test_details[i], p_values[i], adj_p_values[i]))
    
    # Also update the test_results with adjusted p-values
    if(i <= length(test_details)) {
      # Extract condition details
      cond_parts <- strsplit(test_details[i], ", ")[[1]]
      cons_part <- sub("Consensus level: ", "", cond_parts[1])
      cons_part <- sub("Consensus: ", "", cons_part) # Handle both formats
      dir_part <- sub("Direction: ", "", cond_parts[2])
      
      if(length(cond_parts) > 2) {
        # Handle within-trial case with switch status
        switch_part <- sub("Trial: ", "", cond_parts[3])
        switch_val <- ifelse(switch_part == "Switch", "1", "0")
        key <- paste(cons_part, dir_part, switch_val, sep="_")
      } else {
        # Handle consensus case
        key <- paste(cons_part, dir_part, sep="_")
      }
      
      # Update the adj_p_value if the key exists
      if(key %in% names(test_results)) {
        test_results[[key]]$adj_p_value <- adj_p_values[i]
      }
    }
  }
  
  # Add summary of participant inclusion
  results_text <- paste0(results_text, 
                         "\n\nPARTICIPANT INCLUSION SUMMARY:\n============================\n")
  
  # Count participants included in at least one analysis
  all_included_participants <- unique(unlist(participants_per_condition))
  
  results_text <- paste0(results_text,
                         sprintf("Total participants in dataset: %d\n", total_participants),
                         sprintf("Participants included in at least one condition: %d (%.1f%%)\n", 
                                 length(all_included_participants),
                                 100 * length(all_included_participants) / total_participants))
  
  # Add condition-specific counts
  results_text <- paste0(results_text, "\nParticipants per condition:\n")
  
  # Sort condition names for better readability
  sorted_conditions <- sort(names(participants_per_condition))
  for(cond in sorted_conditions) {
    if(is.null(participants_per_condition[[cond]])) next
    
    # Format the condition name for display
    display_name <- if(grepl("_", cond)) {
      parts <- strsplit(cond, "_")[[1]]
      if(length(parts) == 3) {
        # Handle within-trial case
        switch_label <- ifelse(parts[3] == "1", "Switch", "Stay")
        sprintf("%s, %s, %s", parts[1], parts[2], switch_label)
      } else {
        # Handle consensus case
        sprintf("%s, %s", parts[1], parts[2])
      }
    } else {
      cond
    }
    
    results_text <- paste0(results_text,
                           sprintf("  %s: %d participants\n", 
                                   display_name, length(participants_per_condition[[cond]])))
  }
  
  # Save results with descriptive filename
  output_file <- file.path(posthoc_dir, paste0(analysis_name, "_", var_name, "_posthoc.txt"))
  writeLines(results_text, output_file)
  
  return(list(
    text = results_text, 
    file = output_file, 
    p_values = p_values, 
    adj_p_values = adj_p_values,
    test_results = test_results,
    analysis_name = analysis_name,
    participant_counts = participants_per_condition
  ))
}

# Function to create improved visualizations of the post-hoc tests using participant averages
plot_posthoc_simple_slopes <- function(data, var_name, analysis_type, output_dir, test_results) {
  # Create more descriptive analysis name
  analysis_name <- if(analysis_type == "choice_consensus") {
    "choice_switch_by_consensus"
  } else if(analysis_type == "within_trial_switch") {
    "bet_difference_by_switch"
  } else {
    analysis_type
  }
  
  # Create output directory with clearer naming
  posthoc_plot_dir <- file.path(output_dir, "posthoc", analysis_name, var_name, "plots")
  dir.create(posthoc_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Define y-axis label based on analysis type
  y_label <- if(analysis_type == "choice_consensus") {
    "Choice switch probability (%)"
  } else if(analysis_type == "within_trial_switch") {
    "Bet difference (Bet 2 - Bet 1)"
  } else {
    "Outcome value"
  }
  
  if(analysis_type == "choice_consensus") {
    # For each consensus level
    for(cons in c("2:2", "3:1", "4:0")) {
      # Get test results for this consensus level
      against_key <- paste(cons, "Against group", sep="_")
      with_key <- paste(cons, "With group", sep="_")
      
      # Create subtitle with statistics
      against_stats <- if(against_key %in% names(test_results)) {
        result <- test_results[[against_key]]
        sprintf("Against group: r = %.3f, β = %.3f, p = %.4f, adj-p = %.4f", 
                result$r_value, result$beta, result$p_value, result$adj_p_value)
      } else { "Against group: No data" }
      
      with_stats <- if(with_key %in% names(test_results) && cons != "2:2") {
        result <- test_results[[with_key]]
        sprintf("With group: r = %.3f, β = %.3f, p = %.4f, adj-p = %.4f", 
                result$r_value, result$beta, result$p_value, result$adj_p_value)
      } else { "With group: N/A" }
      
      subtitle <- paste(against_stats, with_stats, sep="\n")
      
      # Filter data for this consensus level
      subset_data <- data %>% 
        filter(consensus_level == cons)
      
      # CHANGE: Aggregate by participant for each direction
      against_data <- subset_data %>%
        filter(direction == "Against group") %>%
        group_by(participant.id_in_session) %>%
        summarise(
          outcome_value = mean(outcome_value, na.rm = TRUE),
          scale_name = first(scale_name),  # Should be identical for all trials
          n_trials = n(),
          direction = "Against group",
          .groups = 'drop'
        )
      
      with_data <- subset_data %>%
        filter(direction == "With group") %>%
        group_by(participant.id_in_session) %>%
        summarise(
          outcome_value = mean(outcome_value, na.rm = TRUE),
          scale_name = first(scale_name),  # Should be identical for all trials
          n_trials = n(),
          direction = "With group",
          .groups = 'drop'
        )
      
      # Combine the aggregated data
      aggregated_data <- bind_rows(against_data, with_data)
      
      # DEBUG: Print point counts
      message(sprintf("\nPlotting %s at %s consensus:", var_name, cons))
      message(sprintf("Against group points: %d", nrow(against_data)))
      message(sprintf("With group points: %d", nrow(with_data)))
      message(sprintf("Total points: %d", nrow(aggregated_data)))
      
      # Add additional info for title
      point_count_info <- sprintf("Counts: Against=%d, With=%d", 
                                  nrow(against_data), 
                                  nrow(with_data))
      
      # Create plot with participant averages and jittering
      p <- ggplot(aggregated_data, aes(x = scale_name, y = outcome_value, color = direction)) +
        # Add jittering to reveal overlapping points
        geom_point(aes(size = n_trials), alpha = 0.6, 
                   position = position_jitter(width = 0.1, height = 0.05, seed = 42), 
                   shape = 16) +
        scale_size_continuous(name = "Number of trials", range = c(2, 6)) +
        geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
        scale_color_manual(values = c("Against group" = "red", "With group" = "blue")) +
        # Ensure axes include all points with some padding
        coord_cartesian(
          ylim = c(min(aggregated_data$outcome_value, na.rm = TRUE) - 0.5, 
                   max(aggregated_data$outcome_value, na.rm = TRUE) + 0.5)
        ) +
        labs(
          title = paste("Simple slope for", var_name, "at", cons, "consensus"),
          subtitle = paste(subtitle, point_count_info, sep="\n"),
          x = paste(var_name, "score (standardized)"),
          y = y_label,
          color = "Direction",
          caption = "Note: Points are jittered to reduce overlap. Each point is one participant's average."
        ) +
        theme_custom
      
      # Save plot with descriptive filename
      ggsave(
        filename = file.path(posthoc_plot_dir, paste0(analysis_name, "_", var_name, "_", cons, ".png")),
        plot = p,
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  } else if(analysis_type == "within_trial_switch") {
    # For each consensus level
    for(cons in c("2:2", "3:1", "4:0")) {
      # For each switch status
      for(switch_status in c("0", "1")) {
        # Get test results for this consensus level and switch status
        switch_label <- ifelse(switch_status == "0", "Stay", "Switch")
        
        against_key <- paste(cons, "Against group", switch_status, sep="_")
        with_key <- paste(cons, "With group", switch_status, sep="_")
        
        # Create subtitle with statistics
        against_stats <- if(against_key %in% names(test_results)) {
          result <- test_results[[against_key]]
          sprintf("Against group: r = %.3f, β = %.3f, p = %.4f, adj-p = %.4f", 
                  result$r_value, result$beta, result$p_value, result$adj_p_value)
        } else { "Against group: No data" }
        
        with_stats <- if(with_key %in% names(test_results) && cons != "2:2") {
          result <- test_results[[with_key]]
          sprintf("With group: r = %.3f, β = %.3f, p = %.4f, adj-p = %.4f", 
                  result$r_value, result$beta, result$p_value, result$adj_p_value)
        } else { "With group: N/A" }
        
        subtitle <- paste(against_stats, with_stats, sep="\n")
        
        # Filter data for this consensus level and switch status
        subset_data <- data %>% 
          filter(consensus_level == cons, switch_vs_stay == switch_status)
        
        # Skip if not enough data
        if(nrow(subset_data) < 5) {
          message(sprintf("Skipping %s plot for consensus %s, switch=%s (n=%d)", 
                          var_name, cons, switch_status, nrow(subset_data)))
          next
        }
        
        # CHANGE: Aggregate by participant for each direction
        against_data <- subset_data %>%
          filter(direction == "Against group") %>%
          group_by(participant.id_in_session) %>%
          summarise(
            outcome_value = mean(outcome_value, na.rm = TRUE),
            scale_name = first(scale_name),  # Should be identical for all trials
            n_trials = n(),
            direction = "Against group",
            .groups = 'drop'
          )
        
        with_data <- subset_data %>%
          filter(direction == "With group") %>%
          group_by(participant.id_in_session) %>%
          summarise(
            outcome_value = mean(outcome_value, na.rm = TRUE),
            scale_name = first(scale_name),  # Should be identical for all trials
            n_trials = n(),
            direction = "With group",
            .groups = 'drop'
          )
        
        # Combine the aggregated data
        aggregated_data <- bind_rows(against_data, with_data)
        
        # DEBUG: Print point counts
        message(sprintf("\nPlotting %s at %s consensus, %s trials:", 
                        var_name, cons, switch_label))
        message(sprintf("Against group points: %d", nrow(against_data)))
        message(sprintf("With group points: %d", nrow(with_data)))
        message(sprintf("Total points: %d", nrow(aggregated_data)))
        
        # Add additional info for title
        point_count_info <- sprintf("Counts: Against=%d, With=%d", 
                                    nrow(against_data), 
                                    nrow(with_data))
        
        # Create plot with participant averages and jittering
        p <- ggplot(aggregated_data, aes(x = scale_name, y = outcome_value, color = direction)) +
          # Add jittering to reveal overlapping points
          geom_point(aes(size = n_trials), alpha = 0.6, 
                     position = position_jitter(width = 0.1, height = 0.05, seed = 42), 
                     shape = 16) +
          scale_size_continuous(name = "Number of trials", range = c(2, 6)) +
          geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
          scale_color_manual(values = c("Against group" = "red", "With group" = "blue")) +
          # Ensure axes include all points with some padding
          coord_cartesian(
            ylim = c(min(aggregated_data$outcome_value, na.rm = TRUE) - 0.5, 
                     max(aggregated_data$outcome_value, na.rm = TRUE) + 0.5)
          ) +
          labs(
            title = paste("Simple slope for", var_name, "at", cons, "consensus,", switch_label, "trials"),
            subtitle = paste(subtitle, point_count_info, sep="\n"),
            x = paste(var_name, "score (standardized)"),
            y = y_label,
            color = "Direction",
            caption = "Note: Points are jittered to reduce overlap. Each point is one participant's average."
          ) +
          theme_custom
        
        # Save plot with descriptive filename
        ggsave(
          filename = file.path(posthoc_plot_dir, 
                               paste0(analysis_name, "_", var_name, "_", cons, "_", switch_status, ".png")),
          plot = p,
          width = 8,
          height = 6,
          dpi = 300
        )
      }
    }
  }
  
  # Create combined plots for quick comparison with participant averages
  if(analysis_type == "choice_consensus") {
    # Get data for Against group only
    subset_data <- data %>% 
      filter(direction == "Against group")
    
    # CHANGE: Aggregate by participant for each consensus level
    aggregated_data <- subset_data %>%
      group_by(participant.id_in_session, consensus_level) %>%
      summarise(
        outcome_value = mean(outcome_value, na.rm = TRUE),
        scale_name = first(scale_name),
        n_trials = n(),
        .groups = 'drop'
      )
    
    # DEBUG: Print point counts for combined plot
    message("\nPlotting combined Against group comparison:")
    for(cons in c("2:2", "3:1", "4:0")) {
      count <- sum(aggregated_data$consensus_level == cons)
      message(sprintf("%s consensus: %d points", cons, count))
    }
    
    combined_p <- ggplot(aggregated_data, aes(x = scale_name, y = outcome_value, color = consensus_level)) +
      geom_point(aes(size = n_trials), alpha = 0.6, 
                 position = position_jitter(width = 0.1, height = 0.05, seed = 42), 
                 shape = 16) +
      scale_size_continuous(name = "Number of trials", range = c(2, 6)) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
      scale_color_manual(values = c("2:2" = "#E69F00", "3:1" = "#56B4E9", "4:0" = "#009E73")) +
      coord_cartesian(
        ylim = c(min(aggregated_data$outcome_value, na.rm = TRUE) - 0.5, 
                 max(aggregated_data$outcome_value, na.rm = TRUE) + 0.5)
      ) +
      labs(
        title = paste("Comparison of", var_name, "effect across consensus levels"),
        subtitle = "Against group direction only",
        x = paste(var_name, "score (standardized)"),
        y = y_label,
        color = "Consensus Level",
        caption = "Note: Points are jittered to reduce overlap. Each point is one participant's average."
      ) +
      theme_custom
    
    ggsave(
      filename = file.path(posthoc_plot_dir, paste0(analysis_name, "_", var_name, "_combined.png")),
      plot = combined_p,
      width = 10,
      height = 7,
      dpi = 300
    )
  } else if(analysis_type == "within_trial_switch") {
    # Create combined plots by switch status
    for(switch_status in c("0", "1")) {
      switch_label <- ifelse(switch_status == "0", "Stay", "Switch")
      
      # Get data for Against group and current switch status
      subset_data <- data %>% 
        filter(direction == "Against group", switch_vs_stay == switch_status)
      
      # Skip if not enough data
      if(nrow(subset_data) < 10) {
        message(sprintf("Skipping combined plot for switch=%s (insufficient data)", switch_status))
        next
      }
      
      # CHANGE: Aggregate by participant for each consensus level
      aggregated_data <- subset_data %>%
        group_by(participant.id_in_session, consensus_level) %>%
        summarise(
          outcome_value = mean(outcome_value, na.rm = TRUE),
          scale_name = first(scale_name),
          n_trials = n(),
          .groups = 'drop'
        )
      
      # DEBUG: Print point counts for combined plot
      message(sprintf("\nPlotting combined Against group comparison for %s trials:", switch_label))
      for(cons in c("2:2", "3:1", "4:0")) {
        count <- sum(aggregated_data$consensus_level == cons)
        message(sprintf("%s consensus: %d points", cons, count))
      }
      
      combined_p <- ggplot(aggregated_data, aes(x = scale_name, y = outcome_value, color = consensus_level)) +
        geom_point(aes(size = n_trials), alpha = 0.6, 
                   position = position_jitter(width = 0.1, height = 0.05, seed = 42), 
                   shape = 16) +
        scale_size_continuous(name = "Number of trials", range = c(2, 6)) +
        geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
        scale_color_manual(values = c("2:2" = "#E69F00", "3:1" = "#56B4E9", "4:0" = "#009E73")) +
        coord_cartesian(
          ylim = c(min(aggregated_data$outcome_value, na.rm = TRUE) - 0.5, 
                   max(aggregated_data$outcome_value, na.rm = TRUE) + 0.5)
        ) +
        labs(
          title = paste("Comparison of", var_name, "effect across consensus levels"),
          subtitle = paste("Against group direction,", switch_label, "trials"),
          x = paste(var_name, "score (standardized)"),
          y = y_label,
          color = "Consensus Level",
          caption = "Note: Points are jittered to reduce overlap. Each point is one participant's average."
        ) +
        theme_custom
      
      ggsave(
        filename = file.path(posthoc_plot_dir, 
                             paste0(analysis_name, "_", var_name, "_combined_", switch_status, ".png")),
        plot = combined_p,
        width = 10,
        height = 7,
        dpi = 300
      )
    }
  }
  
  return(posthoc_plot_dir)
}