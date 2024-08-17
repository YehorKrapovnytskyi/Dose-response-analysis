# Load necessary libraries
library(readxl)
library(dplyr)
library(drc)
library(MASS)

source("utils.R")

# Switch to a the current working directory a load the data file 
current_directory <- getwd()
print(current_directory)
data <- read_excel("Data_before_norm_all.xlsx", col_names=FALSE, sheet="BaseLine")
colnames(data) <- c("dose", paste0("response", 1:(ncol(data) - 1)))

# Reshape data to a long format and remove NA elements
data_long <- data %>%
  pivot_longer(cols = starts_with("response"), names_to = "replicate", values_to = "response") %>%
  drop_na()

# Check the response for normally and homoscedasity for every dose for accurate response estimation
check_normality_and_homoscedasticity(data_long, dose_col = "dose", response_col = "response")


# Calculation of dose response #
#                              #

bc <- boxcox(response ~ dose, data = data_long, plotit = FALSE)
optimal_lambda <- bc$x[which.max(bc$y)]
cat("Optimal lambda:", optimal_lambda, "\n")

# Transform the data using the optimal lambda
data_long <- data_long %>%
  mutate(transformed_response = box_cox_transform(x=data_long$response, lambda=0))

# Check the response for normally and homoscedasity again after transform
check_normality_and_homoscedasticity(data_long, dose_col = "dose", response_col = "transformed_response", doses_to_include = c(0.05, 0.1, 0.2, 1))

# Fit the dose-response model using the transformed data
model <- drm(transformed_response ~ dose, data = data_long, fct = LL.4(fixed = c(NA, NA, NA, NA)))

# # Summarize the model fit
summary(model)

# Plot the dose-response curve with transformed
par(mfrow = c(1, 1), pin = c(7, 4))
plot(model, type = "all", main = "Dose-Response Curve", xlab = "Dose", ylab = "Scaled Transformed Response")

# Calculate the Limit of Blank in response and dose units
LoB <- calculate_LoB(data_long)
dose_for_LoB <- ED(model, LoB, type="absolute")

LoD <- calculate_LoD(data_long, dose_col = "dose", response_col = "transformed_response", low_doses = c(0.05, 0.1, 0.2, 1), LoB)

dose_for_LoB <- ED(model, LoB, type="absolute")
dose_for_LoD <- ED(model, LoD, type="absolute")






































