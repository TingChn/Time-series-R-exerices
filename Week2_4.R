#4.1
# Load necessary libraries
library(MASS)       # For multivariate normal distribution
library(vars)       # For fitting VAR models
library(tseries)


# Define parameters
A1 <- matrix(c(0.2, 0.0, -0.3, 0.4), nrow=2, byrow=TRUE)
A2 <- matrix(c(-0.1, 0.1, 0.2, -0.3), nrow=2, byrow=TRUE)
num_iter <- 1000
sample_sizes <- c(100, 200, 500)
criteria_results <- list()

# Function to simulate the VAR(2) process
simulate_VAR2 <- function(T) {
  y <- matrix(0, ncol=2, nrow=T)
  errors <- mvrnorm(T, mu=rep(0, 2), Sigma=diag(2))  # iid errors from N(0, I2)
  
  for (t in 3:T) {
    y[t, ] <- A1 %*% y[t - 1, ] + A2 %*% y[t - 2, ] + errors[t, ]
  }
  
  return(y)
}

# Checking stability
set.seed(123) 
for (T in sample_sizes) {

    data <- simulate_VAR2(T)
    var_data <- data.frame(y = data[, 1], x = data[, 2])
    
      var_model1 <- VAR(var_data, p = 1, type = "const")
      var_model2 <- VAR(var_data, p = 2, type = "const")
      var_model3 <- VAR(var_data, p = 3, type = "const")
      var_model4 <- VAR(var_data, p = 4, type = "const")

# Check stability for each VAR model
stability_var_model1 <- roots(var_model1)
stability_var_model2 <- roots(var_model2)
stability_var_model3 <- roots(var_model3)
stability_var_model4 <- roots(var_model4)  

# Check if each model is stable (all roots must be < 1 for stability)
is_stable_model1 <- all(Mod(stability_var_model1) < 1)
is_stable_model2 <- all(Mod(stability_var_model2) < 1)
is_stable_model3 <- all(Mod(stability_var_model3) < 1)
is_stable_model4 <- all(Mod(stability_var_model4) < 1)

# Display stability status
cat("For sample size", T)
cat("\nIs VAR Model 1 stable? ", is_stable_model1, "\n")
cat("Is VAR Model 2 stable? ", is_stable_model2, "\n")
cat("Is VAR Model 3 stable? ", is_stable_model3, "\n")
cat("Is VAR Model 4 stable? ", is_stable_model4, "\n")
}

#Check stationarity
lag.length = 25
pvalues_matrix <- matrix(nrow = 2, ncol = length(sample_sizes))
rownames(pvalues_matrix) <- c("y_bt_pvalues", "x_bt_pvalues")
colnames(pvalues_matrix) <- sample_sizes
for (i in seq_along(sample_sizes)) {
  T <- sample_sizes[i]
  data <- simulate_VAR2(T)
  var_data <- data.frame(y = data[, 1], x = data[, 2])  
  data <- simulate_VAR2(T)
  var_data <- data.frame(y = data[, 1], x = data[, 2])
  y <- var_data[,1]
  x <- var_data[,2]
  
  #plot
  par(mfrow = c(2, 2))
  plot(1:T, var_data$y, type='l', col='red', xlab="time (t)", ylab="Y(t)",
       main="Signal for Y(t)")
  acf(var_data$y, lag.max = T, xlab = "lag #", ylab = 'ACF', main="ACF of Y(t)", col = "red")
  plot(1:T, var_data$x, type='l', col='blue', xlab="time (t)", ylab="X(t)",
       main="Signal for X(t)")
  acf(var_data$x, lag.max = T, xlab = "lag #", ylab = 'ACF', main="ACF of X(t)", col = "blue")
  par(mfrow = c(1, 1))

  # Perform Ljung-Box test and store p-values
  y_bt <- Box.test(y, lag = lag.length, type = "Ljung-Box")
  x_bt <- Box.test(x, lag = lag.length, type = "Ljung-Box")
  pvalues_matrix["y_bt_pvalues", i] <- y_bt$p.value
  pvalues_matrix["x_bt_pvalues", i] <- x_bt$p.value
}

print(pvalues_matrix)


#4.2
results_tables <- list()
num_iterations <- 1000
sample_sizes <- c(100, 200, 500, 1000)
# Iterate over each sample size
for (T in sample_sizes) {  
  # Initialize vectors to store the best model orders for each iteration
  best_model_AIC <- numeric(num_iterations)
  best_model_BIC <- numeric(num_iterations)
  best_model_HQ <- numeric(num_iterations)
  
  # Iterate for each simulation
  for (i in 1:num_iterations) {
    # Step 1: Generate data based on the VAR(2) process
    Y <- matrix(0, nrow = T, ncol = 2)
    eps <- mvrnorm(n = T, mu = c(0, 0), Sigma = diag(2))  # N(0, I2)
    
    # Initialize the first two rows for the VAR(2) process
    Y[1, ] <- eps[1, ]
    Y[2, ] <- A1 %*% Y[1, ] + eps[2, ]
    
    # Generate the remaining time series values
    for (t in 3:T) {
      Y[t, ] <- A1 %*% Y[t-1, ] + A2 %*% Y[t-2, ] + eps[t, ]
    }
    
    # Step 2: Select the best VAR model based on different criteria
    var_model <- VARselect(Y, lag.max = 4, type="const")
    
    # Store the best lag order for each criterion
    best_model_AIC[i] <- var_model[["selection"]][["AIC(n)"]]
    best_model_BIC[i] <- var_model[["selection"]][["SC(n)"]]
    best_model_HQ[i] <- var_model[["selection"]][["HQ(n)"]]
  }
  
  # Create frequency tables for each criterion
  table_AIC <- table(best_model_AIC)
  table_BIC <- table(best_model_BIC)
  table_HQ <- table(best_model_HQ)
  
  # Combine all tables for the current sample size into a single list
  results_tables[[paste("Sample size", T)]] <- list(
    AIC = table_AIC,
    BIC = table_BIC,
    HQ = table_HQ
  )
}

for (sample_size in names(results_tables)) {
  cat("Results for", sample_size, ":\n")
  print(results_tables[[sample_size]])
  cat("\n")
}
