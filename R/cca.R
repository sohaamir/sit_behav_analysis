# cca_analysis.R

# This script contains all the functions required to perform a Canonical Correlation Analysis (CCA) on two sets of variables.

# Immediate write to console for all print statements and warnings
options(immediate.write = TRUE)

################## SETUP AND CONSTANTS ###################
required_packages <- c(
  "CCA", "ggplot2", "tidyverse", "yacca", "CCP", 
  "PMA", "gridExtra", "dplyr", "knitr", "here", 
  "MVN", "corrplot", "mvtnorm", "boot", "caret", 
  "moments", "psych", "grid"
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

# Import required libraries
library(CCA)
library(ggplot2)
library(tidyverse)
library(yacca)
library(CCP)
library(PMA)
library(gridExtra)
library(dplyr)
library(knitr)
library(here)
library(MVN)
library(corrplot)
library(mvtnorm)
library(boot)
library(caret)
library(moments)
library(psych)

# Set the threshold for multicollinearity - only variables with correlation BELOW this threshold will be retained (i.e., higher value is less stringent)
COLINEAR_THRESHOLD = 0.9

# Set the threshold for loadings - only variables with loadings above/below this threshold will be retained (absolute value)
LOADINGS_THRESHOLD = 0

###################
# TESTING FUNCTIONS
###################

# Function to check normality
check_normality <- function(X, Y) {
  # Univariate tests for X
  X_sw <- data.frame(
    Variable = colnames(X),
    do.call(rbind, lapply(1:ncol(X), function(i) {
      test <- shapiro.test(X[,i])
      data.frame(W = test$statistic, p.value = test$p.value)
    }))
  )
  
  # Univariate tests for Y
  Y_sw <- data.frame(
    Variable = colnames(Y),
    do.call(rbind, lapply(1:ncol(Y), function(i) {
      test <- shapiro.test(Y[,i])
      data.frame(W = test$statistic, p.value = test$p.value)
    }))
  )
  
  # Descriptive statistics
  desc_stats <- function(data) {
    data.frame(
      Variable = colnames(data),
      Mean = colMeans(data),
      SD = apply(data, 2, sd),
      Skewness = apply(data, 2, function(x) moments::skewness(x)),
      Kurtosis = apply(data, 2, function(x) moments::kurtosis(x))
    )
  }
  
  X_desc <- desc_stats(X)
  Y_desc <- desc_stats(Y)
  
  # Mardia's test
  combined_data <- cbind(X, Y)
  mardia_result <- mvn(data = combined_data, mvnTest = "mardia")
  
  return(list(
    univariateNormality = list(
      Behavioral = X_sw,
      Questionnaire = Y_sw
    ),
    Descriptives = list(
      Behavioral = X_desc,
      Questionnaire = Y_desc
    ),
    multivariateNormality = mardia_result$multivariateNormality
  ))
}

# Function for Box's M Test
box_m_test <- function(X, Y) {
  # Get dimensions
  n <- nrow(X)  # sample size
  p <- ncol(X)  # number of X variables
  q <- ncol(Y)  # number of Y variables
  
  # Calculate separate covariance matrices
  S1 <- cov(X)
  S2 <- cov(Y)
  
  # Calculate pooled covariance matrices
  n1 <- n - 1
  n2 <- n - 1
  Sp1 <- S1
  Sp2 <- S2
  
  # Calculate determinants
  det_S1 <- det(S1)
  det_S2 <- det(S2)
  
  # Calculate Box's M statistic
  M <- (n1 * log(det(Sp1)) + n2 * log(det(Sp2))) - 
    ((n1 * log(det_S1)) + (n2 * log(det_S2)))
  
  # Calculate chi-square approximation
  p1 <- p * (p + 1) / 2  
  p2 <- q * (q + 1) / 2  
  df <- p1 + p2
  
  c <- (((2 * p1^2 + 2 * p2^2 - 3) / (6 * (p1 + p2))) * 
          (1/(n1) + 1/(n2) - 1/(n1 + n2)))
  
  chi_sq <- M * (1 - c)
  p_value <- pchisq(chi_sq, df, lower.tail = FALSE)
  
  return(list(
    M_statistic = M,
    chi_square = chi_sq,
    df = df,
    p_value = p_value
  ))
}

# Function to calculate VIF
calculate_vif <- function(data) {
  vif_df <- data.frame(
    Variable = character(),
    VIF = numeric(),
    Tolerance = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(i in 1:ncol(data)) {
    dependent_var <- colnames(data)[i]
    predictor_vars <- colnames(data)[-i]
    
    if(length(predictor_vars) > 0) {
      formula_str <- paste(dependent_var, "~", 
                           paste(predictor_vars, collapse = " + "))
      
      model <- lm(as.formula(formula_str), data = as.data.frame(data))
      
      r2 <- summary(model)$r.squared
      vif <- 1/(1-r2)
      tolerance <- 1-r2
      
      vif_df <- rbind(vif_df, data.frame(
        Variable = dependent_var,
        VIF = vif,
        Tolerance = tolerance,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(vif_df)
}

# Function to handle multicollinearity
handle_multicollinearity <- function(data_matrix, threshold = COLINEAR_THRESHOLD) {
  cat("Starting multicollinearity check with threshold:", threshold, "\n")
  
  # Remove any constant columns first
  var_zeros <- apply(data_matrix, 2, var, na.rm = TRUE) == 0
  if(any(var_zeros)) {
    data_matrix <- data_matrix[, !var_zeros, drop = FALSE]
    cat("Removed constant columns:", 
        paste(colnames(data_matrix)[var_zeros], collapse=", "), "\n")
  }
  
  cor_matrix <- cor(data_matrix, use = "pairwise.complete.obs")
  condition_number <- kappa(cor_matrix)
  
  severity <- if(condition_number > 1e15) {
    "Severe"
  } else if(condition_number > 1e10) {
    "Moderate"
  } else {
    "Low"
  }
  
  vars_to_remove <- c()
  removed_due_to <- list()
  problematic_pairs <- data.frame(
    Var1 = character(), 
    Var2 = character(), 
    Correlation = numeric()
  )
  
  for(i in 1:(ncol(data_matrix)-1)) {
    for(j in (i+1):ncol(data_matrix)) {
      if(abs(cor_matrix[i,j]) > threshold) {
        problematic_pairs <- rbind(
          problematic_pairs,
          data.frame(
            Var1 = colnames(data_matrix)[i],
            Var2 = colnames(data_matrix)[j],
            Correlation = cor_matrix[i,j]
          )
        )
        
        correlations_i <- sum(abs(cor_matrix[i,]) > threshold) - 1
        correlations_j <- sum(abs(cor_matrix[j,]) > threshold) - 1
        
        to_remove <- if(correlations_i > correlations_j) i else j
        vars_to_remove <- c(vars_to_remove, to_remove)
        
        removed_due_to[[length(removed_due_to) + 1]] <- 
          sprintf("%s removed due to r=%.2f correlation with %s",
                  colnames(data_matrix)[to_remove],
                  cor_matrix[i,j],
                  colnames(data_matrix)[if(to_remove == i) j else i])
      }
    }
  }
  
  vars_to_remove <- unique(vars_to_remove)
  vars_to_keep <- setdiff(1:ncol(data_matrix), vars_to_remove)
  
  return(list(
    reduced_matrix = data_matrix[, vars_to_keep, drop = FALSE],
    removed_vars = colnames(data_matrix)[vars_to_remove],
    kept_vars = colnames(data_matrix)[vars_to_keep],
    removal_reasons = removed_due_to,
    correlation_matrix = cor_matrix,
    condition_number = condition_number,
    multicollinearity_severity = severity,
    problematic_pairs = problematic_pairs,
    threshold_used = threshold
  ))
}

###################
# SIGNIFICANCE TESTING
###################

# Wilks' Lambda test
wilks_test <- function(cca_fit, n, p, q) {
  k <- length(cca_fit$cor)
  wilks <- rep(NA, k)
  F_stats <- rep(NA, k)
  p_values <- rep(NA, k)
  df1 <- rep(NA, k)
  df2 <- rep(NA, k)
  
  for(i in 1:k) {
    wilks[i] <- prod((1 - cca_fit$cor[i:k]^2))
    m <- n - 3/2 - (p + q)/2
    p1 <- p - i + 1
    q1 <- q - i + 1
    d <- max(p1, q1)
    df1[i] <- p1 * q1
    df2[i] <- (n - p - q + 1) * d/2 + 1
    
    if(df2[i] > 0) {
      F_stats[i] <- ((1 - wilks[i]^(1/d))/wilks[i]^(1/d)) * df2[i]/df1[i]
      p_values[i] <- pf(F_stats[i], df1[i], df2[i], lower.tail = FALSE)
    }
  }
  
  data.frame(
    Dimension = 1:k,
    Canonical_Correlation = cca_fit$cor,
    Wilks_Lambda = wilks,
    F_statistic = F_stats,
    df1 = df1,
    df2 = df2,
    p_value = p_values
  )
}

# Bartlett's test
bartlett_test <- function(cca_fit, n) {
  k <- length(cca_fit$cor)
  bartlett_stats <- numeric(k)
  p_values <- numeric(k)
  
  for(i in 1:k) {
    remaining_cors <- cca_fit$cor[i:k]
    m <- n - 1/2 * (ncol(X) + ncol(Y) + 3)
    bartlett_stats[i] <- -(m * sum(log(1 - remaining_cors^2)))
    df <- (ncol(X) - i + 1) * (ncol(Y) - i + 1)
    p_values[i] <- pchisq(bartlett_stats[i], df, lower.tail = FALSE)
  }
  
  data.frame(
    Dimension = 1:k,
    Bartlett_Statistic = bartlett_stats,
    P_Value = p_values
  )
}


###################
# BOOTSTRAP FUNCTIONS
###################

# Bootstrap confidence intervals
bootstrap_ci <- function(X, Y, cca_fit, n_boot = 1000, alpha = 0.05) {
  n_dims <- length(cca_fit$cor)
  boot_cors <- matrix(NA, nrow = n_boot, ncol = n_dims)
  
  for(i in 1:n_boot) {
    boot_idx <- sample(nrow(X), replace = TRUE)
    tryCatch({
      boot_cca <- cancor(X[boot_idx,], Y[boot_idx, 1:nrow(cca_fit$ycoef)])
      if(!is.null(boot_cca$cor)) {
        boot_cors[i,] <- boot_cca$cor[1:n_dims]
      }
    }, error = function(e) {
      boot_cors[i,] <- NA
    })
  }
  
  boot_cors <- na.omit(boot_cors)
  ci_lower <- apply(boot_cors, 2, quantile, probs = alpha/2, na.rm = TRUE)
  ci_upper <- apply(boot_cors, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)
  se <- apply(boot_cors, 2, sd, na.rm = TRUE)
  
  data.frame(
    correlation = cca_fit$cor[1:n_dims],
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    se = se
  )
}

###################
# ANALYSIS FUNCTIONS
###################

# Compute redundancy indices
compute_redundancy <- function(X, Y, cca_fit) {
  # Ensure proper dimensions
  n_dims <- min(ncol(X), ncol(Y))
  n_dims <- min(n_dims, ncol(cca_fit$ycoef))
  
  # Use only first dimension for simplicity
  xcoef <- matrix(cca_fit$xcoef[, 1], ncol = 1)
  ycoef <- matrix(cca_fit$ycoef[, 1], ncol = 1)
  
  X_loadings <- cor(X, X %*% xcoef)
  Y_loadings <- cor(Y, Y %*% ycoef)
  
  var_extracted_X <- colMeans(X_loadings^2)
  var_extracted_Y <- colMeans(Y_loadings^2)
  
  redundancy_Y_given_X <- var_extracted_Y * cca_fit$cor[1]^2
  redundancy_X_given_Y <- var_extracted_X * cca_fit$cor[1]^2
  
  list(
    redundancy_Y_given_X = redundancy_Y_given_X,
    redundancy_X_given_Y = redundancy_X_given_Y,
    var_extracted_X = var_extracted_X,
    var_extracted_Y = var_extracted_Y
  )
}

# Enhanced commonality analysis
enhanced_commonality <- function(X, Y, cca_fit, loadings_X, loadings_Y, threshold = 0.3) {
  # Get minimum dimension
  n_dims <- min(ncol(X), ncol(Y))
  n_dims <- min(n_dims, ncol(cca_fit$ycoef))
  
  # Calculate structure coefficients using first dimension only
  struct_coef_X <- cor(X, X %*% cca_fit$xcoef[, 1, drop = FALSE])
  struct_coef_Y <- cor(Y, Y %*% cca_fit$ycoef[, 1, drop = FALSE])
  
  # Calculate function coefficients
  func_coef_X <- cca_fit$xcoef[, 1, drop = FALSE] * apply(X, 2, sd)
  func_coef_Y <- cca_fit$ycoef[, 1, drop = FALSE] * apply(Y, 2, sd)
  
  compare_coefficients <- function(struct, func, names) {
    data.frame(
      Variable = names,
      Loading = round(struct[,1], 3),
      Weight = round(func[,1], 3),
      Discrepancy = round(abs(struct[,1] - func[,1]), 3),
      Suppression = struct[,1] * func[,1] < 0 |
        (abs(struct[,1]) < threshold & abs(func[,1]) > threshold) |
        (abs(struct[,1]) > threshold & abs(func[,1]) < threshold)
    )
  }
  
  X_structure <- compare_coefficients(struct_coef_X, func_coef_X, colnames(X))
  Y_structure <- compare_coefficients(struct_coef_Y, func_coef_Y, colnames(Y))
  
  calc_commonality <- function(data, coef, loadings, vars) {
    unique_var <- sapply(1:ncol(data), function(i) {
      mod <- lm(data[,i] ~ data %*% coef)
      summary(mod)$r.squared
    })
    
    total_var <- loadings[,1]^2
    common_var <- total_var - unique_var
    
    data.frame(
      Variable = vars,
      Unique = round(unique_var, 3),
      Common = round(common_var, 3),
      Total = round(total_var, 3)
    )
  }
  
  X_commonality <- calc_commonality(X, func_coef_X, loadings_X, colnames(X))
  Y_commonality <- calc_commonality(Y, func_coef_Y, loadings_Y, colnames(Y))
  
  list(
    X_structure_function = X_structure,
    Y_structure_function = Y_structure,
    X_commonality = X_commonality,
    Y_commonality = Y_commonality
  )
}

###################
# VISUALIZATION FUNCTIONS
###################

# Create correlation heatmap
create_correlation_heatmap <- function(X, Y, cca_fit) {
  # Get dimensions of canonical vectors
  n_vars_x <- ncol(X)
  n_vars_y <- ncol(Y)
  n_dims <- min(n_vars_x, n_vars_y)  # number of canonical dimensions
  
  # Calculate correlations with properly sized matrices
  X_cors <- cor(X, X %*% cca_fit$xcoef[, 1:n_dims, drop = FALSE])
  Y_cors <- cor(Y, Y %*% cca_fit$ycoef[, 1:n_dims, drop = FALSE])
  
  # Create combined matrix with proper dimensions
  combined_matrix <- matrix(0, 
                            nrow = n_vars_x + n_vars_y,
                            ncol = n_dims)
  
  # Add row names
  rownames(combined_matrix) <- c(colnames(X), colnames(Y))
  colnames(combined_matrix) <- paste0("CV", 1:n_dims)
  
  # Fill in the correlations
  combined_matrix[1:n_vars_x, ] <- X_cors
  combined_matrix[(n_vars_x + 1):(n_vars_x + n_vars_y), ] <- Y_cors
  
  # Check for problematic values
  if(any(is.na(combined_matrix)) || any(!is.finite(combined_matrix))) {
    warning("Problematic values in correlation heatmap matrix. Replacing with zeros.")
    combined_matrix[is.na(combined_matrix)] <- 0
    combined_matrix[!is.finite(combined_matrix)] <- 0
  }
  
  return(combined_matrix)
}

# Create canonical biplot
create_canonical_biplot <- function(X, Y, cca_fit, conf.level = 0.95) {
  scores_X <- scale(X %*% cca_fit$xcoef)
  scores_Y <- scale(Y %*% cca_fit$ycoef)
  
  plot_data <- data.frame(
    CV1_X = scores_X[,1],
    CV1_Y = scores_Y[,1],
    CV2_X = scores_X[,2],
    CV2_Y = scores_Y[,2]
  )
  
  ggplot(plot_data) +
    geom_point(aes(x = CV1_X, y = CV1_Y), color = "blue", alpha = 0.5) +
    stat_ellipse(aes(x = CV1_X, y = CV1_Y), 
                 type = "norm", level = conf.level, color = "red") +
    labs(title = "Canonical Variate Biplot",
         x = "First Canonical Variate",
         y = "Second Canonical Variate") +
    theme_minimal()
}

# Create scree plot
create_scree_plot <- function(cca_fit) {
  variance_explained <- data.frame(
    Dimension = 1:length(cca_fit$cor),
    R_squared = cca_fit$cor^2
  )
  
  ggplot(variance_explained, aes(x = Dimension)) +
    geom_line(aes(y = R_squared), linewidth = 1) +
    geom_point(aes(y = R_squared), size = 3) +
    scale_x_continuous(breaks = seq_along(cca_fit$cor)) +
    labs(title = "Scree Plot of Canonical Correlations",
         x = "Dimension",
         y = "Squared Canonical Correlation") +
    theme_minimal()
}

# Plot loadings structure
plot_loadings_structure <- function(data_frame, title) {
  ggplot(data_frame, aes(x = Loading, y = Weight)) +
    geom_point() +
    geom_text(aes(label = Variable), hjust = -0.1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggtitle(title) +
    theme_minimal()
}

# Plot loadings, cross-loadings, weights, and total contribution
create_loadings_matrix <- function(X, Y, cca_fit, loadings_X, loadings_Y, threshold = LOADINGS_THRESHOLD) {
  # Get correct dimensions
  n_y_vars <- ncol(Y)
  n_dims <- min(ncol(X), n_y_vars)
  
  # Create properly sized coefficient matrix
  ycoef_full <- matrix(0, nrow = n_y_vars, ncol = n_dims)
  ycoef_full[1:nrow(cca_fit$ycoef), 1:ncol(cca_fit$ycoef)] <- cca_fit$ycoef
  
  # Calculate values using full-sized coefficient matrix
  struct_coef_Y <- cor(Y, Y %*% ycoef_full[, 1, drop = FALSE])
  loadings <- struct_coef_Y[,1]
  
  xcoef_full <- matrix(0, nrow = ncol(X), ncol = n_dims)
  xcoef_full[1:nrow(cca_fit$xcoef), 1:ncol(cca_fit$xcoef)] <- cca_fit$xcoef
  
  cross_loadings <- cor(Y, scale(X %*% xcoef_full[, 1, drop = FALSE]))[,1]
  weights <- ycoef_full[,1]
  total_contribution <- abs(loadings * weights)
  
  # Create data frame
  plot_data <- data.frame(
    variable = colnames(Y),
    loadings = loadings,
    cross_loadings = cross_loadings,
    weights = weights,
    total_contribution = total_contribution
  )
  
  # Filter based on significance threshold
  plot_data <- plot_data[abs(plot_data$loadings) >= threshold | 
                           abs(plot_data$cross_loadings) >= threshold, ]
  
  # Check if we have any data to plot
  if(nrow(plot_data) == 0) {
    warning("No variables meet the significance threshold for plotting")
    return(NULL)
  }
  
  # Create the four panels
  p1 <- ggplot(plot_data, aes(x = reorder(variable, loadings), y = loadings)) +
    geom_bar(stat = "identity", fill = "grey") +
    coord_flip() +
    labs(title = "A) Loadings", x = "", y = "Loadings") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  p2 <- ggplot(plot_data, aes(x = reorder(variable, cross_loadings), y = cross_loadings)) +
    geom_bar(stat = "identity", fill = "grey") +
    coord_flip() +
    labs(title = "B) Cross-Loadings", x = "", y = "Cross-Loadings") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  p3 <- ggplot(plot_data, aes(x = reorder(variable, weights), y = weights)) +
    geom_bar(stat = "identity", fill = "grey") +
    coord_flip() +
    labs(title = "C) Weights", x = "", y = "Weights") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  plot_data_sorted <- plot_data %>%
    arrange(desc(total_contribution))
  
  p4 <- ggplot(plot_data_sorted, 
               aes(x = reorder(variable, total_contribution), 
                   y = total_contribution)) +
    geom_bar(stat = "identity", fill = "grey") +
    coord_flip() +
    labs(title = "D) Total Contribution", 
         x = "", 
         y = "Total Contribution on predicting prosocial canonical variate") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  # Arrange all plots in a grid
  grid.arrange(p1, p2, p3, p4, ncol = 2)
}

###################
# OUTPUT FUNCTIONS
###################

# Write detailed results to file
write_detailed_results <- function(results, filename) {
  sink(filename)
  
  # 1. Diagnostic Checks
  cat("\nDIAGNOSTIC CHECKS\n")
  cat("================\n")
  
  # Variables Analysis
  cat("\nVariables in Analysis:\n")
  cat("Original variables set 1: ", 
      paste(colnames(results$data_used$X), collapse=", "), "\n")
  cat("Original variables set 2: ", 
      paste(colnames(results$data_used$Y_original), collapse=", "), "\n")
  cat("Retained variables set 2: ", 
      paste(colnames(results$data_used$Y), collapse=", "), "\n")
  
  # Multicollinearity section
  cat("\nMULTICOLLINEARITY DIAGNOSTICS\n")
  cat("=========================\n")
  cat("Overall Severity:", results$multicollinearity$Y$multicollinearity_severity)
  cat("\nCondition Numbers:")
  cat(sprintf("\nSet 1 Variables: %.2e", 
              kappa(cor(results$data_used$X))))
  cat(sprintf("\nSet 2 Variables: %.2e\n", 
              results$multicollinearity$Y$condition_number))
  
  # Correlation Analysis
  cat("\nCorrelation Analysis:\n")
  cat("\nSet 1 Variables Correlation Matrix:\n")
  print(cor(results$data_used$X), digits = 3)
  cat("\nSet 2 Variables Correlation Matrix:\n")
  print(cor(results$data_used$Y), digits = 3)
  
  # Problematic pairs
  if(!is.null(results$multicollinearity$Y$problematic_pairs) && 
     nrow(results$multicollinearity$Y$problematic_pairs) > 0) {
    cat("\nProblematic Variable Pairs (r >", 
        results$multicollinearity$Y$threshold_used, "):\n")
    problematic_df <- data.frame(
      Variable_1 = results$multicollinearity$Y$problematic_pairs$Var1,
      Variable_2 = results$multicollinearity$Y$problematic_pairs$Var2,
      Correlation = sprintf("%.3f", 
                            results$multicollinearity$Y$problematic_pairs$Correlation)
    )
    print.data.frame(problematic_df, row.names=FALSE)
  }
  
  # Variable Removal Summary
  if(length(results$multicollinearity$Y$removed_vars) > 0) {
    cat("\nVariables Removed Due to Multicollinearity:\n")
    cat(paste("- ", results$multicollinearity$Y$removal_reasons, collapse = "\n"), "\n")
  }
  
  # VIF Results
  cat("\nVARIANCE INFLATION FACTORS:\n")
  cat("\nSet 1 Variables:\n")
  print.data.frame(results$vif_results$X, row.names=FALSE, digits=3)
  cat("\nSet 2 Variables:\n")
  print.data.frame(results$vif_results$Y, row.names=FALSE, digits=3)
  
  # Box's M Test Results
  cat("\nBOX'S M TEST OF HOMOGENEITY:\n")
  cat(sprintf("M statistic: %.3f\n", results$box_m_results$M_statistic))
  cat(sprintf("Chi-square: %.3f\n", results$box_m_results$chi_square))
  cat(sprintf("df: %d\n", results$box_m_results$df))
  cat(sprintf("p-value: %.4f\n", results$box_m_results$p_value))
  
  # Normality Tests
  cat("\nNORMALITY TESTS\n")
  cat("==============\n")
  if(!is.null(results$normality_results$multivariateNormality)) {
    cat("\nMultivariate Normality:\n")
    print.data.frame(results$normality_results$multivariateNormality, row.names=FALSE)
  }
  
  if(!is.null(results$normality_results$univariateNormality)) {
    cat("\nUnivariate Normality:\n")
    cat("\nSet 1 Variables:\n")
    print.data.frame(results$normality_results$univariateNormality$Behavioral, row.names=FALSE)
    cat("\nSet 2 Variables:\n")
    print.data.frame(results$normality_results$univariateNormality$Questionnaire, row.names=FALSE)
  }
  
  # CCA Results
  cat("\nCANONICAL CORRELATION ANALYSIS RESULTS\n")
  cat("=====================================\n")
  
  # Bartlett's Test
  cat("\nBARTLETT'S TEST OF SIGNIFICANCE:\n")
  if(!is.null(results$bartlett_results)) {
    bartlett_df <- data.frame(
      Dimension = results$bartlett_results$Dimension,
      Statistic = sprintf("%.3f", results$bartlett_results$Bartlett_Statistic),
      P_Value = sprintf("%.4f", results$bartlett_results$P_Value)
    )
    print.data.frame(bartlett_df, row.names=FALSE)
  }
  
  # Canonical Correlations
  cat("\nCANONICAL CORRELATIONS AND SIGNIFICANCE:\n")
  for(i in 1:nrow(results$significance_tests)) {
    cat(sprintf("\nDimension %d:\n", i))
    cat(sprintf("Canonical correlation = %.2f (CI: %.2f, %.2f)\n",
                results$significance_tests$Canonical_Correlation[i],
                results$bootstrap_results$ci_lower[i],
                results$bootstrap_results$ci_upper[i]))
    cat(sprintf("Wilks' λ = %.3f, F(%f,%f) = %.2f, p = %.3f\n",
                results$significance_tests$Wilks_Lambda[i],
                results$significance_tests$df1[i],
                results$significance_tests$df2[i],
                results$significance_tests$F_statistic[i],
                results$significance_tests$p_value[i]))
  }
  
  # Add Structure-Function and Commonality Analysis section
  cat("\nSTRUCTURE-FUNCTION AND COMMONALITY ANALYSIS:\n")
  
  # Behavioral Variables
  cat("\nBehavioral Variables:\n")
  cat("\na) Structure-Function Results:\n")
  print.data.frame(results$commonality_results$X_structure_function, 
                   row.names=FALSE, right=FALSE)
  cat("\nb) Commonality Analysis:\n")
  print.data.frame(results$commonality_results$X_commonality, 
                   row.names=FALSE, right=FALSE)
  
  # Questionnaire Variables
  cat("\nQuestionnaire Variables:\n")
  cat("\na) Structure-Function Results:\n")
  print.data.frame(results$commonality_results$Y_structure_function, 
                   row.names=FALSE, right=FALSE)
  cat("\nb) Commonality Analysis:\n")
  print.data.frame(results$commonality_results$Y_commonality, 
                   row.names=FALSE, right=FALSE)
  
  # Add Cross-Loadings section
  cat("\nCROSS-LOADINGS:\n")
  cat("\nBehavioral Variables with Questionnaire Variates:\n")
  print(cor(results$data_used$X, 
            scale(results$data_used$Y %*% results$cca_fit$ycoef)), 
        digits = 3)
  
  cat("\nQuestionnaire Variables with Behavioral Variates:\n")
  print(cor(results$data_used$Y, 
            scale(results$data_used$X %*% results$cca_fit$xcoef)), 
        digits = 3)
  
  # Enhanced Redundancy Analysis
  cat("\nREDUNDANCY ANALYSIS:\n")
  if(!is.null(results$redundancy_results)) {
    # Create redundancy dataframe
    n_dims <- length(results$redundancy_results$redundancy_Y_given_X)
    redundancy_df <- data.frame(
      Dimension = seq_len(n_dims),
      Y_given_X = results$redundancy_results$redundancy_Y_given_X,
      X_given_Y = results$redundancy_results$redundancy_X_given_Y
    )
    names(redundancy_df) <- c("Dimension", 
                              "Questionnaire from Behavioral", 
                              "Behavioral from Questionnaire")
    
    cat("\na) Redundancy Indices:\n")
    redundancy_df[,2:3] <- lapply(redundancy_df[,2:3], 
                                  function(x) sprintf("%.4f", x))
    print.data.frame(redundancy_df, row.names=FALSE)
    
    cat("\nb) Variance Extracted:\n")
    variance_df <- data.frame(
      Set = c("Behavioral", "Questionnaire"),
      Variance_Extracted = c(
        mean(results$redundancy_results$var_extracted_X),
        mean(results$redundancy_results$var_extracted_Y)
      )
    )
    variance_df$Variance_Extracted <- sprintf("%.4f", 
                                              variance_df$Variance_Extracted)
    print.data.frame(variance_df, row.names=FALSE)
  }
  
  # Add Variance and Effect Sizes section
  cat("\nVARIANCE AND EFFECT SIZES:\n")
  
  # Variance Explained table
  cat("\nVariance Explained:\n")
  if(!is.null(results$variance_explained)) {
    variance_df <- results$variance_explained
    variance_df$R_squared <- sprintf("%.3f", variance_df$R_squared)
    variance_df$Proportion <- sprintf("%.4f", variance_df$Proportion)
    variance_df$Cumulative <- sprintf("%.3f", variance_df$Cumulative)
    print.data.frame(variance_df, row.names=FALSE, right=FALSE)
  }
  
  # Effect Sizes table
  if(!is.null(results$effect_sizes)) {
    cat("\nEffect Sizes:\n")
    effect_df <- results$effect_sizes
    effect_df$Correlation <- sprintf("%.3f", effect_df$Correlation)
    print.data.frame(effect_df, row.names=FALSE, right=FALSE)
  }
  
  sink()
}

# Save all plots
save_all_plots <- function(results, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Correlation heatmap
  if (!is.null(results$plots$correlation_heatmap)) {
    png(file.path(output_dir, "correlation_heatmap.png"), 
        width = 12, height = 10, units = "in", res = 300)
    tryCatch({
      corrplot(results$plots$correlation_heatmap,
               method = "color",
               type = "full",
               tl.col = "black",
               tl.cex = 0.7)
    }, error = function(e) {
      warning("Error in correlation heatmap plotting: ", e$message, 
              "\nTrying alternative visualization...")
      # Simple alternative visualization
      image(t(results$plots$correlation_heatmap), 
            axes = FALSE, 
            main = "Correlation Heatmap")
      axis(1, at = seq(0, 1, length.out = ncol(results$plots$correlation_heatmap)),
           labels = colnames(results$plots$correlation_heatmap), las = 2)
      axis(2, at = seq(0, 1, length.out = nrow(results$plots$correlation_heatmap)),
           labels = rownames(results$plots$correlation_heatmap), las = 2)
    })
    dev.off()
  }
  
  # Canonical biplot
  if (!is.null(results$plots$canonical_biplot)) {
    ggsave(file.path(output_dir, "canonical_biplot.png"), 
           results$plots$canonical_biplot, 
           width = 10, height = 8)
  }
  
  # Scree plot
  if (!is.null(results$plots$scree_plot)) {
    ggsave(file.path(output_dir, "scree_plot.png"), 
           results$plots$scree_plot, 
           width = 8, height = 6)
  }
  
  # Structure-function plots
  if (!is.null(results$plots$structure_x)) {
    ggsave(file.path(output_dir, "structure_function_set1.png"), 
           results$plots$structure_x, 
           width = 10, height = 8)
  }
  
  if (!is.null(results$plots$structure_y)) {
    ggsave(file.path(output_dir, "structure_function_set2.png"), 
           results$plots$structure_y, 
           width = 10, height = 8)
  }
  # Save loadings matrix
  if (!is.null(results$plots$loadings_matrix)) {
    png(file.path(output_dir, "loadings_matrix.png"), 
        width = 12, height = 10, units = "in", res = 300)
    grid.arrange(results$plots$loadings_matrix)
    dev.off()
  }
}

###################
# MAIN CCA FUNCTION
###################

# Main CCA analysis function
run_cca_analysis <- function(X, Y, output_dir = NULL) {
  # Add X variable names if missing
  if(is.null(colnames(X))) {
    colnames(X) <- paste0("X", 1:ncol(X))
  }
  
  # Add Y variable names if missing
  if(is.null(colnames(Y))) {
    colnames(Y) <- paste0("Y", 1:ncol(Y))
  }
  
  # Handle multicollinearity
  Y_multi <- handle_multicollinearity(Y)
  Y_reduced <- Y_multi$reduced_matrix
  
  # Calculate VIF
  vif_X <- calculate_vif(X)
  vif_Y <- calculate_vif(Y_reduced)
  
  # Box's M test
  box_m_results <- box_m_test(X, Y_reduced)
  
  # Add these diagnostic prints after Y_reduced is created and before CCA
  cat("Dimensions of X:", dim(X)[1], "x", dim(X)[2], "\n")
  cat("Dimensions of Y_reduced:", dim(Y_reduced)[1], "x", dim(Y_reduced)[2], "\n")
  
  # Run CCA with explicit dimension handling
  cca_fit <- tryCatch({
    # Get minimum dimension
    n_dims <- min(ncol(X), ncol(Y_reduced))
    
    # Run cancor
    cca_result <- cancor(X, Y_reduced)
    
    # Get actual dimensions from cancor result
    n_ycoef_rows <- min(nrow(cca_result$ycoef), ncol(Y_reduced))
    n_ycoef_cols <- min(ncol(cca_result$ycoef), n_dims)
    n_xcoef_rows <- min(nrow(cca_result$xcoef), ncol(X))
    n_xcoef_cols <- min(ncol(cca_result$xcoef), n_dims)
    
    # Ensure coefficient matrices match Y_reduced dimensions
    new_ycoef <- matrix(0, nrow = ncol(Y_reduced), ncol = n_dims)
    if (n_ycoef_rows > 0 && n_ycoef_cols > 0) {
      new_ycoef[1:n_ycoef_rows, 1:n_ycoef_cols] <- cca_result$ycoef[1:n_ycoef_rows, 1:n_ycoef_cols]
    }
    cca_result$ycoef <- new_ycoef
    
    # Also adjust xcoef if needed
    new_xcoef <- matrix(0, nrow = ncol(X), ncol = n_dims)
    if (n_xcoef_rows > 0 && n_xcoef_cols > 0) {
      new_xcoef[1:n_xcoef_rows, 1:n_xcoef_cols] <- cca_result$xcoef[1:n_xcoef_rows, 1:n_xcoef_cols]
    }
    cca_result$xcoef <- new_xcoef
    
    # Truncate correlations if needed
    n_cors <- min(length(cca_result$cor), n_dims)
    if (n_cors > 0) {
      cca_result$cor <- cca_result$cor[1:n_cors]
    }
    
    cca_result
  }, error = function(e) {
    cat("Error in CCA:", e$message, "\n")
    stop(e)
  })
  
  # Add these after CCA
  cat("Dimensions of cca_fit$xcoef:", dim(cca_fit$xcoef)[1], "x", dim(cca_fit$xcoef)[2], "\n")
  cat("Dimensions of cca_fit$ycoef:", dim(cca_fit$ycoef)[1], "x", dim(cca_fit$ycoef)[2], "\n")
  
  # Calculate dimensions for subsequent analyses
  n_dims <- min(ncol(X), ncol(Y_reduced))
  
  # Ensure coefficient matrices have correct dimensions
  cca_fit$xcoef <- cca_fit$xcoef[, 1:n_dims, drop = FALSE]
  cca_fit$ycoef <- cca_fit$ycoef[, 1:n_dims, drop = FALSE]
  cca_fit$cor <- cca_fit$cor[1:n_dims]
  
  # Calculate various test results
  bartlett_results <- bartlett_test(cca_fit, nrow(X))
  normality_results <- check_normality(X, Y_reduced)
  significance_tests <- wilks_test(cca_fit, nrow(X), ncol(X), ncol(Y_reduced))
  bootstrap_results <- bootstrap_ci(X, Y_reduced, cca_fit)
  
  # Calculate loadings with proper scaling
  loadings_X <- cor(X, scale(X %*% cca_fit$xcoef))
  loadings_Y <- cor(Y_reduced, scale(Y_reduced %*% cca_fit$ycoef))
  
  # Additional analyses
  commonality_results <- enhanced_commonality(X, Y_reduced, cca_fit, loadings_X, loadings_Y)
  redundancy_results <- compute_redundancy(X, Y_reduced, cca_fit)
  
  # Calculate variance explained
  variance_explained <- data.frame(
    Dimension = 1:length(cca_fit$cor),
    R_squared = cca_fit$cor^2,
    Proportion = cca_fit$cor^2/sum(cca_fit$cor^2),
    Cumulative = cumsum(cca_fit$cor^2/sum(cca_fit$cor^2))
  )
  
  # Calculate effect sizes
  effect_sizes <- data.frame(
    Dimension = 1:length(cca_fit$cor),
    Correlation = cca_fit$cor,
    Effect_Size = ifelse(cca_fit$cor >= 0.5, "Large",
                         ifelse(cca_fit$cor >= 0.3, "Medium", "Small"))
  )
  
  # Generate plots
  plots <- list(
    correlation_heatmap = create_correlation_heatmap(X, Y_reduced, cca_fit),
    canonical_biplot = create_canonical_biplot(X, Y_reduced, cca_fit),
    scree_plot = create_scree_plot(cca_fit),
    structure_x = plot_loadings_structure(
      commonality_results$X_structure_function, 
      "Structure-Function Plot for Set 1 Variables"
    ),
    structure_y = plot_loadings_structure(
      commonality_results$Y_structure_function,
      "Structure-Function Plot for Set 2 Variables"
    ),
    
    loadings_correlation = corrplot(rbind(loadings_X, loadings_Y)),
    
    loadings_matrix = create_loadings_matrix(X, Y_reduced, cca_fit, loadings_X, loadings_Y)
  )
  
  # Compile results
  results <- list(
    cca_fit = cca_fit,
    normality_results = normality_results,
    multicollinearity = list(Y = Y_multi),
    significance_tests = significance_tests,
    bootstrap_results = bootstrap_results,
    commonality_results = commonality_results,
    loadings = list(X = loadings_X, Y = loadings_Y),
    data_used = list(X = X, Y = Y_reduced, Y_original = Y),
    redundancy_results = redundancy_results,
    vif_results = list(X = vif_X, Y = vif_Y),
    box_m_results = box_m_results,
    bartlett_results = bartlett_results,
    variance_explained = variance_explained,
    effect_sizes = effect_sizes,
    plots = plots
  )
  
  # Save results if output directory provided
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write_detailed_results(results, file.path(output_dir, "cca_detailed_results.txt"))
    save_all_plots(results, file.path(output_dir, "plots"))
  }
  
  return(results)
}






