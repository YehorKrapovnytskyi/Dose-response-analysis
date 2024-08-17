# Load necessary libraries
#library(pracma) # For the detrend function
library(dplyr)
library(car)
library(tidyr)

calculate_LoB <- function(data, dose_col = "dose", response_col = "transformed_response", ks_test_alpha = 0.05) {
  
  # Filter for zero-dose data
  zero_dose_data <- data %>%
    dplyr::filter(!!sym(dose_col) == 0) %>%
    dplyr::select(!!sym(response_col)) %>%
    pull() # Extract the response values as a vector
  
  # Check if zero-dose data is empty
  if (length(zero_dose_data) == 0) {
    stop("No zero-dose data found. Check if the dose column is correct.")
  }
  
  # Calculate degrees of freedom
  df <- length(zero_dose_data) - 1
  
  # Find the critical value from the t-distribution for 95% one-sided confidence level
  t_critical <- qt(1 - ks_test_alpha, df)
  
  # Calculate the LoB
  LoB <- mean(zero_dose_data) + t_critical * sd(zero_dose_data)
  
  # Print the calculated LoB
  cat(sprintf("The Level of Blank (LoB) is: %.3f\n", LoB))
  return(LoB)
}


check_normality_and_homoscedasticity <- function(data, dose_col = "dose", response_col = "response", doses_to_include = NULL) 
{
  
  data <- data %>%
    dplyr::select(all_of(dose_col), all_of(response_col)) %>%
    mutate(dose = as.factor(data[[dose_col]]), response = as.numeric(data[[response_col]]))
  
  if (!is.null(doses_to_include)) {
    data <- data %>%
      filter(dose %in% doses_to_include)
  }

  # Check normality for each dose
  normality_results <- data %>%
    group_by(dose) %>%
    summarize(p_value = ks.test(response, "pnorm", mean = mean(response), sd = sd(response))$p.value) %>%
    ungroup()
  
  # Print KS test results for normality
  print("Normality Test (KS Test) Results:")
  print(normality_results)
  
  # Check for homoscedasticity using Levene's Test
  levene_test_result <- leveneTest(response ~ factor(dose), data = data)
  
  # Print Levene's Test result
  print("Homoscedasticity Test (Levene's Test) Result:")
  print(levene_test_result)
}

calculate_LoD <- function(data, dose_col = "dose", response_col = "transformed_response", low_doses, LoB) {
  
  # Filter for low-dose data
  low_dose_data <- data %>%
    filter(!!sym(dose_col) %in% low_doses)
  
  # Calculate variance for each low dose
  variance_data <- low_dose_data %>%
    group_by(!!sym(dose_col)) %>%
    summarize(
      variance = var(!!sym(response_col)),
      n = n()
    ) %>%
    ungroup()
  
  # Calculate pooled variance
  weighted_variance_sum <- sum(variance_data$n * variance_data$variance)
  total_n <- sum(variance_data$n)
  pooled_variance <- weighted_variance_sum / total_n
  root_pooled_variance <- sqrt(pooled_variance)
  
  # Calculate degrees of freedom (f)
  f <- total_n - length(low_doses)
  
  # Calculate C_beta
  c_beta <- 1.645 / (1 - 1 / (4 * f))
  
  # Calculate the LoD
  LoD <- LoB + c_beta * root_pooled_variance
  
  # Print and return the calculated LoD
  cat(sprintf("The Limit of Detection (LoD) is: %.3f\n", LoD))
  return(LoD)
}


min_max_scale <- function(x, new_min = 1, new_max = 2) 
{
  old_min <- min(x)
  old_max <- max(x)
  scaled_x <- ((x - old_min) / (old_max - old_min)) * (new_max - new_min) + new_min
  return(scaled_x)
}



box_cox_transform <- function(x, lambda) {
  if (lambda == 0) {
    return(log(x))
  } else {
    return((x^lambda - 1) / lambda)
  }
}
