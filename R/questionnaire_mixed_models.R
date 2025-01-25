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
  "switch_difference",   # Switch difference analysis
  "within_trial_switch", # Within-trial switch analysis
  "reversal_learning"    # New reversal learning analysis
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
          
          # Handle different plot types based on analysis type
          plot_obj <- if(results[[scale]]$moderation_results$main_results$model %>%
                         inherits("list") && analysis_type == "reversal_learning") {
            # For reversal learning, combine choice1 and choice2 plots
            cowplot::plot_grid(
              scale_results$plots[[plot_name]]$choice1,
              scale_results$plots[[plot_name]]$choice2,
              ncol = 2,
              labels = c("Choice 1", "Choice 2")
            )
          } else if(is.list(scale_results$plots[[plot_name]]) && 
                    "plot" %in% names(scale_results$plots[[plot_name]])) {
            scale_results$plots[[plot_name]]$plot
          } else if(inherits(scale_results$plots[[plot_name]], "gtable")) {
            scale_results$plots[[plot_name]]
          } else {
            scale_results$plots[[plot_name]]
          }
          
          # Save plot with appropriate handling
          if(inherits(plot_obj, "gtable") || inherits(plot_obj, "arrange")) {
            png(file.path(output_dir, file_name), 
                width = 12, height = 6, units = "in", res = 300)
            grid::grid.draw(plot_obj)
            dev.off()
          } else if(inherits(plot_obj, "ggplot")) {
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
  if(params$analysis_type == "reversal_learning") {
    cat("\n===== Trial Position Summary =====\n")
    trial_counts <- table(data$trial_to_reversal)
    print(trial_counts)
    cat("================================\n")
  } else if(params$analysis_type == "switch_difference") {
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
  
  # Apply winning model structure to each questionnaire
  message("\nAnalyzing questionnaires...")
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

    # Store raw p-values for interactions
    if(params$analysis_type == "reversal_learning") {
      interaction_p_values[[paste0(var, "_choice1")]] <- analysis_results$main_results$raw_p$choice1
      interaction_p_values[[paste0(var, "_choice2")]] <- analysis_results$main_results$raw_p$choice2
    } else {
      interaction_p_values[[paste0(var, "_main")]] <- analysis_results$main_results$raw_p
    }
    
    print(str(analysis_results$main_results$raw_p))
    
    # Store subscale p-values if they exist
    if(!is.null(analysis_results$subscale_results)) {
      for(subscale in names(analysis_results$subscale_results)) {
        if(params$analysis_type == "reversal_learning") {
          interaction_p_values[[paste0(var, "_", subscale, "_choice1")]] <- 
            analysis_results$subscale_results[[subscale]]$raw_p$choice1
          interaction_p_values[[paste0(var, "_", subscale, "_choice2")]] <- 
            analysis_results$subscale_results[[subscale]]$raw_p$choice2
        } else {
          interaction_p_values[[paste0(var, "_", subscale)]] <- 
            analysis_results$subscale_results[[subscale]]$raw_p
        }
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
    if(params$analysis_type == "reversal_learning") {
      # Update main scale adjusted p-values for both choices
      all_results[[var]]$moderation_results$main_results$adj_p <- list(
        choice1 = all_adj_p_values[paste0(var, "_choice1")],
        choice2 = all_adj_p_values[paste0(var, "_choice2")]
      )
    } else {
      # Update main scale adjusted p-value
      all_results[[var]]$moderation_results$main_results$adj_p <- 
        all_adj_p_values[paste0(var, "_main")]
    }
    
    # Update subscale adjusted p-values if they exist
    if(!is.null(all_results[[var]]$moderation_results$subscale_results)) {
      for(subscale in names(all_results[[var]]$moderation_results$subscale_results)) {
        if(params$analysis_type == "reversal_learning") {
          all_results[[var]]$moderation_results$subscale_results[[subscale]]$adj_p <- list(
            choice1 = all_adj_p_values[paste0(var, "_", subscale, "_choice1")],
            choice2 = all_adj_p_values[paste0(var, "_", subscale, "_choice2")]
          )
        } else {
          all_results[[var]]$moderation_results$subscale_results[[subscale]]$adj_p <- 
            all_adj_p_values[paste0(var, "_", subscale)]
        }
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
  pooled_data <- if(params$analysis_type == "reversal_learning") {
    map_dfr(params$questionnaire_vars, function(var) {
      processed <- params$preprocessing_fn(data, var, params)
      
      # For reversal learning, we need the reshaped data structure
      base_cols <- c("participant.id_in_session", "trial_to_reversal",
                     "choice_type", "outcome_value",
                     "age", "gender", "scale_name")
      
      processed %>%
        select(all_of(base_cols)) %>%
        distinct()
    })
  } else {
    map_dfr(params$questionnaire_vars, function(var) {
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
  }
  
  # Add diagnostic prints
  message("\nDiagnostic information for pooled data:")
  message("Structure of pooled data:")
  print(str(pooled_data))
  message("\nFirst few rows:")
  print(head(pooled_data))
  message("\nSummary of key variables:")
  print(summary(pooled_data))
  message("\nCounts of trial positions and choice types:")
  if(params$analysis_type == "reversal_learning") {
    print(table(pooled_data$trial_to_reversal, pooled_data$choice_type))
  }
  
  # Continue with model comparison...
  message("\nRunning model comparison on pooled data...")
  model_results <- run_model_comparison(params$model_formulas, pooled_data, params)
  
  return(model_results)
}

################## MODEL FUNCTIONS ###################
run_model_comparison <- function(model_formulas, data, params) {
  message("\nComparing models:")
  models <- vector("list", length(model_formulas))
  names(models) <- names(model_formulas)
  
  for(model_name in names(model_formulas)) {
    message(sprintf("\nTrying to fit %s...", model_name))
    message("Formula:", deparse(model_formulas[[model_name]]$formula))
    
    if(params$analysis_type == "reversal_learning") {
      # For reversal learning, fit a single model with choice_type as predictor
      models[[model_name]] <- tryCatch({
        model <- lmer(
          formula = model_formulas[[model_name]]$formula,
          data = data,
          control = lmerControl(optimizer = "bobyqa",
                                optCtrl = list(maxfun = 100000))
        )
        message("Model successfully fitted")
        message("\nModel summary:")
        print(summary(model)$coefficients)
        model
      }, error = function(e) {
        message(sprintf("Failed to fit model: %s", e$message))
        NULL
      }, warning = function(w) {
        message(sprintf("Warning in model fitting: %s", w$message))
      })
    } else {
      # Original code for other analysis types
      models[[model_name]] <- tryCatch({
        model <- lmer(
          formula = model_formulas[[model_name]]$formula,
          data = data,
          control = lmerControl(optimizer = "bobyqa",
                                optCtrl = list(maxfun = 100000))
        )
        message("Model successfully fitted")
        message("\nModel summary:")
        print(summary(model)$coefficients)
        model
      }, error = function(e) {
        message(sprintf("Failed to fit %s: %s", model_name, e$message))
        NULL
      }, warning = function(w) {
        message(sprintf("Warning in model fitting: %s", w$message))
      })
    }
  }
  
  # Remove NULL models
  models <- models[!sapply(models, is.null)]
  
  if(length(models) == 0) {
    message("\nDetailed model fitting diagnostics:")
    message("Number of observations:", nrow(data))
    message("Number of participants:", length(unique(data$participant.id_in_session)))
    if(params$analysis_type == "reversal_learning") {
      message("\nChoice type distribution:")
      print(table(data$choice_type))
      message("\nTrial position distribution:")
      print(table(data$trial_to_reversal))
    }
    stop("No models successfully fitted")
  }
  
  # Compare models using AIC/BIC
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
  
  # Extract p-value for interaction based on analysis type
  if(params$analysis_type == "reversal_learning") {
    interaction_term <- "trial_to_reversal:choice_type:scale_name"
    raw_p <- list(
      choice1 = tryCatch({
        choice1_interaction <- "trial_to_reversal:scale_name"
        main_anova[choice1_interaction, "Pr(>Chisq)"]
      }, error = function(e) NA),
      choice2 = tryCatch({
        choice2_interaction <- "choice_typechoice2_accuracy:scale_name"
        main_anova[choice2_interaction, "Pr(>Chisq)"]
      }, error = function(e) NA)
    )
  } else if(params$analysis_type == "within_trial_switch") {
    interaction_term <- "consensus_level:direction:switch_vs_stay:scale_name"
    raw_p <- main_anova[interaction_term, "Pr(>Chisq)"]
  } else if(params$analysis_type == "switch_difference") {
    interaction_term <- "consensus_level:direction:switch_difference:scale_name"
    raw_p <- main_anova[interaction_term, "Pr(>Chisq)"]
  } else {
    interaction_term <- "consensus_level:direction:scale_name"
    raw_p <- main_anova[interaction_term, "Pr(>Chisq)"]
  }
  
  message("\nDebugging p-values:")
  message("Analysis type: ", params$analysis_type)
  message("raw_p structure:")
  print(str(raw_p))
  message("raw_p content:")
  print(raw_p)
  
  # Print significant results, handling both single p-values and lists
  if(params$analysis_type == "reversal_learning") {
    # For reversal learning, check both choice types
    for(choice in c("choice1", "choice2")) {
      if(!is.na(raw_p[[choice]]) && raw_p[[choice]] < 0.05) {
        message(sprintf("  Significant interaction found for %s (p = %.3f)", 
                        choice, raw_p[[choice]]))
      }
    }
  } else {
    # For other analysis types, check single p-value
    if(!is.na(raw_p) && raw_p < 0.05) {
      message(sprintf("  Significant interaction found (p = %.3f)", raw_p))
    }
  }
  
  # Create main results list
  main_results <- list(
    model = main_model,
    anova_results = main_anova,
    raw_p = raw_p,
    diagnostics = check_model_diagnostics(main_model, params)
  )
  
  # Store subscale results
  subscale_results <- list()
  
  # Run subscale analyses if applicable
  current_scale <- var
  if(!is.null(params$subscale_mapping[[current_scale]])) {
    message(sprintf("  Processing subscales for %s", current_scale))
    
    for(subscale in params$subscale_mapping[[current_scale]]) {
      if(subscale %in% colnames(processed_data)) {
        message(sprintf("    Analyzing subscale: %s", subscale))
        
        # Create subscale data
        subscale_data <- processed_data
        subscale_data$scale_name <- as.numeric(scale(processed_data[[subscale]]))
        
        # Fit subscale model
        sub_model <- tryCatch({
          lmer(formula(winning_model), data = subscale_data,
               control = lmerControl(optimizer = "bobyqa"))
        }, error = function(e) NULL)
        
        if(!is.null(sub_model)) {
          # Get ANOVA results for subscale
          sub_anova <- car::Anova(sub_model, type = 2)
          
          # Extract p-value with same logic as main model
          if(params$analysis_type == "reversal_learning") {
            sub_raw_p <- list(
              choice1 = tryCatch({
                sub_anova[interaction_term, "Pr(>Chisq)"]
              }, error = function(e) NA),
              choice2 = tryCatch({
                sub_anova[interaction_term, "Pr(>Chisq)"]
              }, error = function(e) NA)
            )
            
            # Print significant results for subscales
            for(choice in c("choice1", "choice2")) {
              if(!is.na(sub_raw_p[[choice]]) && sub_raw_p[[choice]] < 0.05) {
                message(sprintf("    Significant interaction for %s in %s (p = %.3f)", 
                                subscale, choice, sub_raw_p[[choice]]))
              }
            }
          } else {
            sub_raw_p <- sub_anova[interaction_term, "Pr(>Chisq)"]
            if(!is.na(sub_raw_p) && sub_raw_p < 0.05) {
              message(sprintf("    Significant interaction for %s (p = %.3f)", 
                              subscale, sub_raw_p))
            }
          }
          
          # Store subscale results
          subscale_results[[subscale]] <- list(
            model = sub_model,
            anova_results = sub_anova,
            raw_p = sub_raw_p,
            diagnostics = check_model_diagnostics(sub_model, params)
          )
        }
      }
    }
  }
  
  return(list(
    main_results = main_results,
    subscale_results = if(length(subscale_results) > 0) subscale_results else NULL,
    processed_data = processed_data
  ))
}

################## ANALYZE WINNING MODEL ###################
analyze_winning_model <- function(model, data, params) {
  if(params$analysis_type == "reversal_learning") {
    # Run diagnostics
    diagnostics <- check_model_diagnostics(model, params)
    
    # Run ANOVA
    anova_results <- car::Anova(model, type = 2)
    
    # Get interaction terms
    three_way_term <- "trial_to_reversal:choice_type:scale_name"
    two_way_terms <- c(
      "trial_to_reversal:choice_type",
      "trial_to_reversal:scale_name",
      "choice_type:scale_name"
    )
    
    # Extract raw p-values
    raw_p <- tryCatch({
      if(three_way_term %in% rownames(anova_results)) {
        anova_results[three_way_term, "Pr(>Chisq)"]
      } else {
        # Get minimum p-value from two-way interactions
        min(sapply(two_way_terms, function(term) {
          if(term %in% rownames(anova_results)) {
            anova_results[term, "Pr(>Chisq)"]
          } else {
            1.0
          }
        }))
      }
    }, error = function(e) NA)
    
    return(list(
      diagnostics = diagnostics,
      anova_results = anova_results,
      raw_p = raw_p,  # Now a single value
      model = model
    ))
  } else {
    # Run diagnostics
    diagnostics <- check_model_diagnostics(model, params)
    
    # Run ANOVA
    anova_results <- car::Anova(model, type = 2)
    
    # Get p-values and adjust
    p_values <- anova_results[,"Pr(>Chisq)"]
    
    # Get interaction term based on analysis type
    if(params$analysis_type == "within_trial_switch") {
      interaction_term <- "consensus_level:direction:switch_vs_stay:scale_name"
    } else if(params$analysis_type == "switch_difference") {
      interaction_term <- "consensus_level:direction:switch_difference:scale_name"
    } else {
      interaction_term <- "consensus_level:direction:scale_name"
    }
    
    # Extract raw p-value for interaction
    raw_p <- p_values[interaction_term]
    
    list(
      diagnostics = diagnostics,
      anova_results = anova_results,
      raw_p = raw_p,
      model = model
    )
  }
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
  print("Data structure for plotting:")
  print(str(processed_data))
  
  current_plots <- list()
  
  if(params$analysis_type == "reversal_learning") {
    print("Generating reversal learning plots...")
    
    # For reversal learning, we need separate plots for each choice/bet
    current_plots[["continuous"]] <- list(
      choice1 = tryCatch({
        plot_continuous_relationship(
          data = processed_data,
          var = "choice1_accuracy",
          outcome_var = "choice1_accuracy",
          results = moderation_results$main_results,
          params = params,
          display_name = paste(var, "- Choice 1")
        )
      }, error = function(e) NULL),
      
      choice2 = tryCatch({
        plot_continuous_relationship(
          data = processed_data,
          var = "choice2_accuracy",
          outcome_var = "choice2_accuracy",
          results = moderation_results$main_results,
          params = params,
          display_name = paste(var, "- Choice 2")
        )
      }, error = function(e) NULL)
    )
    
    current_plots[["median_split"]] <- tryCatch({
      plot_median_split(
        data = processed_data,
        var = var,
        results = moderation_results$main_results,
        median_var = "scale_name",
        params = params,
        display_name = var
      )
    }, error = function(e) NULL)
    
    current_plots[["moderation"]] <- list(
      choice1 = tryCatch({
        plot_moderation_effects(
          model = moderation_results$main_results$models$choice1,
          data = processed_data,
          var = "scale_name",
          raw_p = moderation_results$main_results$raw_p$choice1,
          adj_p = moderation_results$main_results$adj_p$choice1,
          params = params,
          display_name = paste(var, "- Choice 1")
        )$plot
      }, error = function(e) NULL),
      
      choice2 = tryCatch({
        plot_moderation_effects(
          model = moderation_results$main_results$models$choice2,
          data = processed_data,
          var = "scale_name",
          raw_p = moderation_results$main_results$raw_p$choice2,
          adj_p = moderation_results$main_results$adj_p$choice2,
          params = params,
          display_name = paste(var, "- Choice 2")
        )$plot
      }, error = function(e) NULL)
    )
    
    current_plots[["simple_slopes"]] <- list(
      choice1 = tryCatch({
        plot_simple_slopes(
          model = moderation_results$main_results$models$choice1,
          var = "scale_name",
          data = processed_data,
          raw_p = moderation_results$main_results$raw_p$choice1,
          adj_p = moderation_results$main_results$adj_p$choice1,
          params = params,
          display_name = paste(var, "- Choice 1")
        )
      }, error = function(e) NULL),
      
      choice2 = tryCatch({
        plot_simple_slopes(
          model = moderation_results$main_results$models$choice2,
          var = "scale_name",
          data = processed_data,
          raw_p = moderation_results$main_results$raw_p$choice2,
          adj_p = moderation_results$main_results$adj_p$choice2,
          params = params,
          display_name = paste(var, "- Choice 2")
        )
      }, error = function(e) NULL)
    )
    
  } else {
    # Original plotting code for other analysis types
    current_plots[["continuous"]] <- tryCatch({
      plot_continuous_relationship(
        data = processed_data,
        var = "scale_name",
        outcome_var = "outcome_value",
        results = moderation_results$main_results,
        params = params,
        display_name = var
      )
    }, error = function(e) NULL)
    
    current_plots[["median_split"]] <- tryCatch({
      plot_median_split(
        data = processed_data,
        var = var,
        results = moderation_results$main_results,
        median_var = "scale_name",
        params = params,
        display_name = var
      )
    }, error = function(e) NULL)
    
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
    }, error = function(e) NULL)
    
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
    }, error = function(e) NULL)
  }
  
  # Add subscale plots if they exist and are significant
  if(!is.null(moderation_results$subscale_results)) {
    for(subscale in names(moderation_results$subscale_results)) {
      sub_result <- moderation_results$subscale_results[[subscale]]
      
      if(params$analysis_type == "reversal_learning") {
        current_plots[[paste0(subscale, "_moderation")]] <- list(
          choice1 = tryCatch({
            plot_moderation_effects(
              model = sub_result$models$choice1,
              data = processed_data,
              var = subscale,
              raw_p = sub_result$raw_p$choice1,
              adj_p = sub_result$adj_p$choice1,
              params = params,
              display_name = subscale
            )$plot
          }, error = function(e) NULL),
          
          choice2 = tryCatch({
            plot_moderation_effects(
              model = sub_result$models$choice2,
              data = processed_data,
              var = subscale,
              raw_p = sub_result$raw_p$choice2,
              adj_p = sub_result$adj_p$choice2,
              params = params,
              display_name = subscale
            )$plot
          }, error = function(e) NULL)
        )
        
      } else {
        current_plots[[paste0(subscale, "_moderation")]] <- tryCatch({
          plot_moderation_effects(
            model = sub_result$model,
            data = processed_data,
            var = subscale,
            raw_p = sub_result$raw_p,
            adj_p = sub_result$adj_p,
            params = params,
            display_name = subscale
          )$plot
        }, error = function(e) NULL)
      }
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
  if(params$analysis_type == "reversal_learning") {
    # Print counts per condition
    participant_counts <- data %>% 
      group_by(trial_to_reversal, choice_type) %>% 
      summarise(
        n = n_distinct(participant.id_in_session),
        .groups = 'drop'
      )
    
    # Calculate correlations
    effects <- data %>%
      group_by(trial_to_reversal, choice_type) %>%
      summarise(
        correlation = cor(.data[[var]], outcome_value),
        .groups = 'drop'
      )
    
    # Extract p-values for subtitle
    if(!is.null(results$raw_p)) {
      subtitle_text <- sprintf("Choice 1 interaction p = %.3f (adj: %.3f)\nChoice 2 interaction p = %.3f (adj: %.3f)",
                               results$raw_p$choice1, results$adj_p$choice1,
                               results$raw_p$choice2, results$adj_p$choice2)
    } else {
      subtitle_text <- "P-values not available"
    }
    
    # Create plot
    p <- ggplot(data, 
                aes(x = .data[[var]], 
                    y = outcome_value, 
                    color = choice_type)) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", 
                  formula = y ~ x,
                  se = TRUE) +
      facet_wrap(~trial_to_reversal) +
      scale_color_manual(values = c("Choice 1" = "#56B4E9", 
                                    "Choice 2" = "#E69F00")) +
      labs(x = display_name,
           y = params$y_label,
           title = paste("Relationship between", display_name, "and", params$y_label),
           subtitle = subtitle_text,
           caption = paste("Effect sizes (r) range:", 
                           round(min(effects$correlation, na.rm = TRUE), 3), "to",
                           round(max(effects$correlation, na.rm = TRUE), 3))) +
      theme_custom
    
  } else if(params$analysis_type == "switch_difference") {
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
  
  # Extract p-values for subtitle
  if(!is.null(results$raw_p)) {
    if(params$analysis_type == "reversal_learning") {
      subtitle_text <- paste("Interaction p-values:\n",
                             sprintf("Choice 1: %.3f (adj: %.3f)\n", 
                                     results$raw_p$choice1, results$adj_p$choice1),
                             sprintf("Choice 2: %.3f (adj: %.3f)", 
                                     results$raw_p$choice2, results$adj_p$choice2))
    } else {
      subtitle_text <- paste("Interaction:\nRaw p =",
                             format.pval(results$raw_p, digits = 3),
                             "\nFDR-adjusted p =",
                             format.pval(results$adj_p, digits = 3))
    }
  } else {
    subtitle_text <- "P-values not available"
  }
  
  # Calculate correlations
  if(params$analysis_type == "reversal_learning") {
    effects <- plot_data %>%
      group_by(trial_to_reversal) %>%
      summarise(
        correlation = cor(.data[[var]], .data[[outcome_var]]),
        .groups = 'drop'
      )
  } else {
    effects <- plot_data %>%
      group_by(direction, consensus_level) %>%
      summarise(
        correlation = cor(.data[[var]], .data[[outcome_var]]),
        .groups = 'drop'
      )
  }
  
  # Base plot setup
  if(params$analysis_type == "reversal_learning") {
    p <- ggplot(plot_data, 
                aes(x = .data[[var]], 
                    y = .data[[outcome_var]])) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", 
                  formula = y ~ x,
                  se = TRUE) +
      facet_wrap(~trial_to_reversal) +
      labs(x = display_name,
           y = params$y_label,
           title = paste("Relationship between", display_name, "and", params$y_label),
           subtitle = subtitle_text,
           caption = paste("Effect sizes (r) range:", 
                           round(min(effects$correlation, na.rm = TRUE), 3), "to",
                           round(max(effects$correlation, na.rm = TRUE), 3))) +
      theme_custom
    
  } else {
    p <- ggplot(plot_data, 
                aes(x = .data[[var]], 
                    y = .data[[outcome_var]], 
                    color = direction)) +
      geom_smooth(method = "lm", 
                  formula = y ~ x,
                  se = TRUE) +
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
    
    p <- p + labs(x = display_name,
                  y = params$y_label,
                  title = paste("Relationship between", display_name, "and", params$y_label),
                  subtitle = subtitle_text,
                  caption = paste("Effect sizes (r) range:", 
                                  round(min(effects$correlation, na.rm = TRUE), 3), "to",
                                  round(max(effects$correlation, na.rm = TRUE), 3)))
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
  
  if(params$analysis_type == "reversal_learning") {
    # Work with the long format data
    summary_data <- data %>%
      mutate(quest_group = factor(
        ifelse(.data[[median_var]] > median(.data[[median_var]]), "High", "Low"),
        levels = c("Low", "High")
      )) %>%
      group_by(trial_to_reversal, quest_group, choice_type) %>%
      summarise(
        mean = mean(outcome_value),  # Use outcome_value instead of choice1_accuracy
        se = sd(outcome_value) / sqrt(n()),
        .groups = 'drop'
      )
    
    # Create plot
    p <- ggplot(summary_data, 
                aes(x = trial_to_reversal, 
                    y = mean, 
                    color = choice_type,
                    group = interaction(choice_type, quest_group))) +
      geom_line(aes(linetype = quest_group), linewidth = 1) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = mean - se, 
                        ymax = mean + se), 
                    width = 0.2) +
      scale_color_manual(values = c("Choice 1" = "#56B4E9", 
                                    "Choice 2" = "#E69F00")) +
      scale_linetype_manual(values = c("Low" = "dashed", "High" = "solid")) +
      labs(x = "Trial relative to reversal",
           y = params$y_label,
           title = paste("Effect of", display_name, "(Median-split)"),
           subtitle = subtitle_text,
           color = "Choice Type",
           linetype = "Scale Group") +
      theme_custom
    
  } else {  # Non-reversal learning analysis types
    # Create quest_group factor with ordered levels
    data <- data %>%
      mutate(quest_group = factor(
        ifelse(.data[[median_var]] > median(.data[[median_var]]), "High", "Low"),
        levels = c("Low", "High")
      ))
    
    summary_data <- data %>%
      group_by(consensus_level, direction, quest_group) %>%
      summarise(
        mean_outcome = mean(outcome_value),
        se = sd(outcome_value) / sqrt(n()),
        .groups = 'drop'
      )
    
    if(params$is_percentage) {
      summary_data$mean_outcome <- summary_data$mean_outcome * 100
      summary_data$se <- summary_data$se * 100
    }
    
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
      facet_wrap(~ quest_group, ncol = 2) +
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
  
  # Handle subtitle text differently for reversal learning
  subtitle_text <- if(params$analysis_type == "reversal_learning") {
    if(!is.null(raw_p) && !is.null(adj_p)) {
      paste("Interaction test:",
            "\nChoice 1: Raw p =", format.pval(raw_p$choice1, digits = 3),
            ", FDR-adjusted p =", format.pval(adj_p$choice1, digits = 3),
            "\nChoice 2: Raw p =", format.pval(raw_p$choice2, digits = 3),
            ", FDR-adjusted p =", format.pval(adj_p$choice2, digits = 3))
    } else {
      "P-values not available"
    }
  } else {
    if(!is.null(raw_p) && !is.na(raw_p)) {
      paste("Interaction test:",
            "\nRaw p =", format.pval(raw_p, digits = 3),
            "\nFDR-adjusted p =", format.pval(adj_p, digits = 3))
    } else {
      "P-values not available"
    }
  }
  
  common_theme <- theme_custom
  color_scheme <- scale_color_manual(values = c("Against group" = "red", "With group" = "blue"))
  
  # Rest of your prediction grid code remains the same
  if(params$analysis_type == "reversal_learning") {
    pred_data <- expand.grid(
      trial_to_reversal = levels(data$trial_to_reversal),
      scale_name = seq(from = -2, to = 2, length.out = 100),
      choice_type = c("Choice 1", "Choice 2"),
      age = 0
    )
  } else if(params$analysis_type == "switch_difference") {
    pred_data <- expand.grid(
      consensus_level = levels(data$consensus_level),
      direction = levels(data$direction),
      switch_difference = c(0, 1),
      scale_name = seq(from = -2, to = 2, length.out = 100),
      age = 0
    ) %>%
      mutate(trial_type = factor(ifelse(switch_difference == 0, "Stay", "Switch")))
  } else if(params$analysis_type == "within_trial_switch") {
    pred_data <- expand.grid(
      consensus_level = levels(data$consensus_level),
      direction = levels(data$direction),
      switch_vs_stay = levels(data$switch_vs_stay),
      scale_name = seq(from = -2, to = 2, length.out = 100),
      age = 0
    )
  } else {
    pred_data <- expand.grid(
      consensus_level = levels(data$consensus_level),
      direction = levels(data$direction),
      scale_name = seq(from = -2, to = 2, length.out = 100),
      age = 0
    )
  }
  
  diagnostics$pred_data <- pred_data
  
  # Generate predictions
  pred_data$predicted <- predict(model, newdata = pred_data, re.form = NA)
  
  diagnostics$predictions <- pred_data$predicted
  diagnostics$final_plot_data <- pred_data
  
  # Create appropriate plot based on analysis type
  if(params$analysis_type == "reversal_learning") {
    p <- ggplot(pred_data, 
                aes(x = scale_name, 
                    y = predicted, 
                    color = choice_type,
                    group = choice_type)) +
      geom_line(linewidth = 1) +
      facet_wrap(~trial_to_reversal) +
      scale_color_manual(values = c("Choice 1" = "#56B4E9", "Choice 2" = "#E69F00"))
  } else if(params$analysis_type == "switch_difference") {
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
      facet_grid(switch_vs_stay ~ consensus_level)
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
         x = ifelse(params$analysis_type == "reversal_learning",
                    paste(display_name, "score (standardized)"),
                    paste(display_name, "score (standardized)")),
         y = params$y_label) +
    common_theme
  
  if(params$analysis_type != "reversal_learning") {
    p <- p + color_scheme
  }
  
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
  
  # Handle subtitle text differently for reversal learning
  subtitle_text <- if(params$analysis_type == "reversal_learning") {
    if(!is.null(raw_p) && !is.null(adj_p)) {
      paste("Interaction test:",
            "\nChoice 1: Raw p =", format.pval(raw_p$choice1, digits = 3),
            ", FDR-adjusted p =", format.pval(adj_p$choice1, digits = 3),
            "\nChoice 2: Raw p =", format.pval(raw_p$choice2, digits = 3),
            ", FDR-adjusted p =", format.pval(adj_p$choice2, digits = 3))
    } else {
      "P-values not available"
    }
  } else {
    if(!is.null(raw_p) && !is.na(raw_p)) {
      paste("Interaction test:",
            "\nRaw p =", format.pval(raw_p, digits = 3),
            "\nFDR-adjusted p =", format.pval(adj_p, digits = 3))
    } else {
      "P-values not available"
    }
  }
  
  common_theme <- theme_custom
  color_scheme <- scale_color_manual(values = c("Against group" = "red", "With group" = "blue"))
  fill_scheme <- scale_fill_manual(values = c("Against group" = "#f36a7b", "With group" = "#3f9fef"))
  
  # Rest of your existing code remains the same until the plot creation
  if(params$analysis_type == "reversal_learning") {
    pred_data <- expand.grid(
      trial_to_reversal = levels(data$trial_to_reversal),
      scale_name = c(-1, 0, 1),  # -1 SD, Mean, +1 SD
      choice_type = c("Choice 1", "Choice 2"),
      age = 0
    )
    
    # Generate predictions
    pred_data$predicted <- predict(model, newdata = pred_data, re.form = NA)
    
    if(params$is_percentage) {
      pred_data$predicted <- pred_data$predicted * 100
    }
    
    # Create plot
    p <- ggplot(pred_data, 
                aes(x = scale_name, 
                    y = predicted, 
                    color = choice_type,
                    group = interaction(choice_type, trial_to_reversal))) +
      geom_line(aes(linetype = trial_to_reversal), 
                linewidth = 1) +
      scale_color_manual(values = c("Choice 1" = "#56B4E9", "Choice 2" = "#E69F00"))
    
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
      scale_color_manual(values = c("2:2" = "#E69F00", 
                                    "3:1" = "#56B4E9", 
                                    "4:0" = "#009E73"))
    
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
      scale_color_manual(values = c("2:2" = "#E69F00", 
                                    "3:1" = "#56B4E9", 
                                    "4:0" = "#009E73"))
    
  } else {  # choice_consensus
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
      scale_color_manual(values = c("2:2" = "#E69F00", 
                                    "3:1" = "#56B4E9", 
                                    "4:0" = "#009E73"))
  }
  
  diagnostics$pred_data <- pred_data
  diagnostics$predictions <- pred_data$predicted
  
  p <- p +
    labs(title = paste("Simple slopes for", display_name),
         subtitle = subtitle_text,
         x = ifelse(params$analysis_type == "reversal_learning",
                    paste(display_name, "score (standardized)"),
                    paste(display_name, "score (standardized)")),
         y = params$y_label,
         color = "Choice Type",
         linetype = "Trial Position") +
    common_theme
  
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
  if(params$analysis_type == "reversal_learning") {
    # Additional checks specific to reversal learning analysis
    reversal_diagnostics <- tryCatch({
      list(
        trial_distribution = table(model_frame$trial_to_reversal),
        complete_cases = sum(complete.cases(model_frame))
      )
    }, error = function(e) NULL)
  } else if(params$analysis_type == "switch_difference") {
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
  if(params$analysis_type == "reversal_learning" && !is.null(reversal_diagnostics)) {
    diagnostics$reversal <- reversal_diagnostics
  } else if(params$analysis_type == "switch_difference" && !is.null(switch_diff_diagnostics)) {
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
  
  if(params$analysis_type == "reversal_learning") {
    # Main scale models for both choices/bets
    print("Fitting base models...")
    base_models <- list(
      choice1 = lmer(
        choice1_accuracy ~ trial_to_reversal + scale_name + 
          age + (1|participant.id_in_session),
        data = data,
        control = lmerControl(optimizer = "bobyqa")
      ),
      choice2 = lmer(
        choice2_accuracy ~ trial_to_reversal + scale_name + 
          age + (1|participant.id_in_session),
        data = data,
        control = lmerControl(optimizer = "bobyqa")
      )
    )
    
    print("Fitting moderation models...")
    mod_models <- list(
      choice1 = lmer(
        choice1_accuracy ~ trial_to_reversal * scale_name + 
          age + (1|participant.id_in_session),
        data = data,
        control = lmerControl(optimizer = "bobyqa")
      ),
      choice2 = lmer(
        choice2_accuracy ~ trial_to_reversal * scale_name + 
          age + (1|participant.id_in_session),
        data = data,
        control = lmerControl(optimizer = "bobyqa")
      )
    )
    
    # Model comparisons and diagnostics
    main_anova <- list(
      choice1 = car::Anova(mod_models$choice1, type = 2),
      choice2 = car::Anova(mod_models$choice2, type = 2)
    )
    
    model_comparison <- list(
      choice1 = anova(base_models$choice1, mod_models$choice1),
      choice2 = anova(base_models$choice2, mod_models$choice2)
    )
    
    model_diagnostics <- list(
      choice1 = check_model_diagnostics(mod_models$choice1, params),
      choice2 = check_model_diagnostics(mod_models$choice2, params)
    )
    
  } else {
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
    main_anova <- car::Anova(mod_model, type = 2)
    model_comparison <- anova(base_model, mod_model)
    model_diagnostics <- check_model_diagnostics(mod_model, params)
  }
  
  # Calculate effect sizes
  print("Calculating effect sizes...")
  if(params$analysis_type == "reversal_learning") {
    effect_sizes <- list(
      choice1 = tryCatch({
        effectsize::eta_squared(mod_models$choice1, partial = TRUE)
      }, error = function(e) NULL),
      choice2 = tryCatch({
        effectsize::eta_squared(mod_models$choice2, partial = TRUE)
      }, error = function(e) NULL)
    )
  } else {
    effect_sizes <- tryCatch({
      effectsize::eta_squared(mod_model, partial = TRUE)
    }, error = function(e) NULL)
  }
  
  # Store main results
  if(params$analysis_type == "reversal_learning") {
    main_results <- list(
      models = mod_models,
      base_models = base_models,
      anova_results = main_anova,
      comparison = model_comparison,
      diagnostics = model_diagnostics,
      effect_sizes = effect_sizes
    )
  } else {
    main_results <- list(
      model = mod_model,
      base_model = base_model,
      anova_results = main_anova,
      comparison = model_comparison,
      diagnostics = model_diagnostics,
      effect_sizes = effect_sizes
    )
  }
  
  # Store subscale results
  subscale_results <- list()
  
  # Run subscale analyses if applicable
  current_scale <- params$current_var
  if(!is.null(params$subscale_mapping[[current_scale]])) {
    print(paste("Processing subscales for", current_scale))
    for(subscale in params$subscale_mapping[[current_scale]]) {
      print(paste("Analyzing subscale:", subscale))
      
      if(subscale %in% colnames(data)) {
        if(params$analysis_type == "reversal_learning") {
          # Fit subscale models for both choices/bets
          sub_models <- list(
            choice1 = tryCatch({
              lmer(
                choice1_accuracy ~ trial_to_reversal * scale_name * subscale + 
                  age + (1|participant.id_in_session),
                data = data,
                control = lmerControl(optimizer = "bobyqa")
              )
            }, error = function(e) NULL),
            choice2 = tryCatch({
              lmer(
                choice2_accuracy ~ trial_to_reversal * scale_name * subscale + 
                  age + (1|participant.id_in_session),
                data = data,
                control = lmerControl(optimizer = "bobyqa")
              )
            }, error = function(e) NULL)
          )
          
          if(!is.null(sub_models$choice1) && !is.null(sub_models$choice2)) {
            sub_anova <- list(
              choice1 = car::Anova(sub_models$choice1, type = 2),
              choice2 = car::Anova(sub_models$choice2, type = 2)
            )
            
            sub_diagnostics <- list(
              choice1 = check_model_diagnostics(sub_models$choice1, params),
              choice2 = check_model_diagnostics(sub_models$choice2, params)
            )
            
            subscale_results[[subscale]] <- list(
              models = sub_models,
              anova_results = sub_anova,
              diagnostics = sub_diagnostics
            )
          }
        } else {
          # Fit subscale model
          sub_model <- tryCatch({
            lmer(
              outcome_value ~ consensus_level * direction * scale_name * subscale + 
                age + (1|participant.id_in_session),
              data = data,
              control = lmerControl(optimizer = "bobyqa")
            )
          }, error = function(e) NULL)
          
          if(!is.null(sub_model)) {
            sub_anova <- car::Anova(sub_model, type = 2)
            sub_diagnostics <- check_model_diagnostics(sub_model, params)
            
            subscale_results[[subscale]] <- list(
              model = sub_model,
              anova_results = sub_anova,
              diagnostics = sub_diagnostics
            )
          }
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
    if(params$analysis_type == "reversal_learning" && is.list(model)) {
      # Calculate effect sizes for both choice models
      choice1_effects <- effectsize::F_to_eta2(
        f = car::Anova(model$choice1, type = 2)$`F value`,
        df = car::Anova(model$choice1, type = 2)$NumDF,
        df_error = car::Anova(model$choice1, type = 2)$DenDF,
        ci = 0.95
      )
      
      choice2_effects <- effectsize::F_to_eta2(
        f = car::Anova(model$choice2, type = 2)$`F value`,
        df = car::Anova(model$choice2, type = 2)$NumDF,
        df_error = car::Anova(model$choice2, type = 2)$DenDF,
        ci = 0.95
      )
      
      # Format results for both choices
      choice1_formatted <- data.frame(
        Term = rownames(car::Anova(model$choice1, type = 2)),
        Eta2_partial = sprintf("%.3f", choice1_effects$Eta2_partial),
        CI_low = sprintf("%.3f", choice1_effects$CI_low),
        CI_high = sprintf("%.3f", choice1_effects$CI_high),
        Choice = "Choice 1",
        stringsAsFactors = FALSE
      )
      
      choice2_formatted <- data.frame(
        Term = rownames(car::Anova(model$choice2, type = 2)),
        Eta2_partial = sprintf("%.3f", choice2_effects$Eta2_partial),
        CI_low = sprintf("%.3f", choice2_effects$CI_low),
        CI_high = sprintf("%.3f", choice2_effects$CI_high),
        Choice = "Choice 2",
        stringsAsFactors = FALSE
      )
      
      # Combine results
      formatted_eta <- rbind(choice1_formatted, choice2_formatted)
      
    } else {
      # Get ANOVA results
      anova_results <- car::Anova(model, type = 2)
      
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
      } else if(params$analysis_type == "reversal_learning") {
        c("trial_to_reversal", "scale_name")
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
    }
    
    # Add indicator for main effects vs interactions
    formatted_eta$Effect_Type <- sapply(formatted_eta$Term, function(x) {
      if(params$analysis_type == "reversal_learning") {
        if(x %in% c("trial_to_reversal", "scale_name")) "Main Effect" else "Interaction"
      } else {
        if(x %in% main_terms) "Main Effect" else "Interaction"
      }
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
    
    if("trial_to_reversal" %in% names(data)) {
      # For reversal learning analysis
      for(trial in unique(data$trial_to_reversal)) {
        # Create new data for this trial
        test_data <- data.frame(
          trial_to_reversal = trial,
          scale_name = c(-1, 1),  # Test at ±1 SD
          age = 0  # Set to mean
        )
        
        # Get predictions
        preds <- predict(model, newdata = test_data, re.form = NA)
        
        # Calculate slope
        slope <- (preds[2] - preds[1]) / 2
        
        # Store results
        slopes_results[[as.character(trial)]] <- slope
      }
    } else {
      # Original analysis for consensus levels and directions
      for(cons in levels(data$consensus_level)) {
        for(dir in levels(data$direction)) {
          test_data <- data.frame(
            consensus_level = cons,
            direction = dir,
            scale_name = c(-1, 1),  # Test at ±1 SD
            age = 0  # Set to mean
          )
          
          # Get predictions
          preds <- predict(model, newdata = test_data, re.form = NA)
          
          # Calculate slope
          slope <- (preds[2] - preds[1]) / 2
          
          # Store results
          slopes_results[[paste(cons, dir)]] <- slope
        }
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
  
  if(params$analysis_type == "reversal_learning") {
    # Add trial-by-trial statistics
    output_text <- paste0(output_text, "TRIAL-BY-TRIAL STATISTICS:\n",
                          "==========================\n\n")
    
    for(trial in sort(unique(all_results[[1]]$processed_data$trial_to_reversal))) {
      trial_data <- all_results[[1]]$processed_data %>%
        filter(trial_to_reversal == trial)
      
      if(params$is_percentage) {
        stats <- trial_data %>%
          summarise(
            choice1_mean = mean(choice1_accuracy) * 100,
            choice1_sd = sd(choice1_accuracy) * 100,
            choice2_mean = mean(choice2_accuracy) * 100,
            choice2_sd = sd(choice2_accuracy) * 100
          )
        
        output_text <- paste0(output_text,
                              "Trial ", trial, ":\n",
                              sprintf("Choice 1: %.2f%% (SD = %.2f)\n", 
                                      stats$choice1_mean, stats$choice1_sd),
                              sprintf("Choice 2: %.2f%% (SD = %.2f)\n", 
                                      stats$choice2_mean, stats$choice2_sd),
                              "\n")
      } else {
        stats <- trial_data %>%
          summarise(
            bet1_mean = mean(bet1_magnitude),
            bet1_sd = sd(bet1_magnitude),
            bet2_mean = mean(bet2_magnitude),
            bet2_sd = sd(bet2_magnitude)
          )
        
        output_text <- paste0(output_text,
                              "Trial ", trial, ":\n",
                              sprintf("Bet 1: %.2f (SD = %.2f)\n", 
                                      stats$bet1_mean, stats$bet1_sd),
                              sprintf("Bet 2: %.2f (SD = %.2f)\n", 
                                      stats$bet2_mean, stats$bet2_sd),
                              "\n")
      }
    }
  }
  
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
  
  # Add analysis-specific diagnostics
  if(!is.null(diagnostics$reversal)) {
    output <- paste0(output,
                     "Reversal Learning:\n",
                     "  Trial distribution:\n",
                     paste(capture.output(diagnostics$reversal$trial_distribution), collapse = "\n"), "\n",
                     "  Complete cases: ", diagnostics$reversal$complete_cases, "\n\n")
  } else if(!is.null(diagnostics$switch_difference)) {
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
  if(any(grepl("Choice", names(effect_sizes)))) {
    # Format for reversal learning (separate tables for each choice)
    table <- paste0(
      "\nEFFECT SIZES:\n============\n"
    )
    
    for(choice in unique(effect_sizes$Choice)) {
      choice_data <- effect_sizes[effect_sizes$Choice == choice, ]
      
      table <- paste0(
        table,
        "\n", choice, ":\n",
        sprintf("%-40s %10s %10s %10s\n", "Term", "Eta2", "CI_low", "CI_high")
      )
      
      for(i in 1:nrow(choice_data)) {
        table <- paste0(
          table,
          sprintf("%-40s %10s %10s %10s\n",
                  choice_data$Term[i],
                  choice_data$Eta2_partial[i],
                  choice_data$CI_low[i],
                  choice_data$CI_high[i])
        )
      }
      table <- paste0(table, "\n")
    }
  } else {
    # Original format for other analysis types
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
    
    if(params$analysis_type == "reversal_learning") {
      # Format results for both choices
      for(choice in c("choice1", "choice2")) {
        if(!is.null(sub_result$anova_results[[choice]])) {
          output <- paste0(output,
                           paste0("\n", toupper(choice), " Results:\n"),
                           "ANOVA Results:\n-------------\n",
                           paste(capture.output(sub_result$anova_results[[choice]]), 
                                 collapse = "\n"),
                           "\n")
        }
      }
    } else {
      if(!is.null(sub_result$anova_results)) {
        output <- paste0(output,
                         "ANOVA Results:\n-------------\n",
                         paste(capture.output(sub_result$anova_results), 
                               collapse = "\n"),
                         "\n")
      }
    }
    
    # Add diagnostics for each subscale
    if(!is.null(sub_result$diagnostics)) {
      if(params$analysis_type == "reversal_learning") {
        for(choice in c("choice1", "choice2")) {
          output <- paste0(output, 
                           "\nDiagnostics for ", toupper(choice), ":\n",
                           format_diagnostics(sub_result$diagnostics[[choice]]))
        }
      } else {
        output <- paste0(output, "\n", format_diagnostics(sub_result$diagnostics))
      }
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
# Helper function to create appropriate file names based on analysis type
create_output_filename <- function(params, prefix) {
  paste0(params$analysis_type, "_", prefix, ".txt")
}

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