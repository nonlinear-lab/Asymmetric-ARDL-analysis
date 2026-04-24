#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: Set up simulation parameters
#++++++++++++++++++++++++++++++++++++++++++++++++

set.seed(2026) # Set seed for reproducibility

# 1. Simulation Parameters based on Pesaran et al. (2001)
T_obs <- 1000      # Sample size T = 1000
reps  <- 1000      # Number of replications (Change to 40000 for exact replication)
k <- 1             # Number of regressors

# Vectors to store the F-statistics for both bounds
F_stat_I0 <- numeric(reps)
F_stat_I1 <- numeric(reps)

# Degrees of freedom for the F-test
# 2 restrictions (coefficients on y_t-1 and x_t-1 are 0)
df_num <- k + 1
# T_obs minus parameters (intercept, y_lag, x_lag)
df_den <- T_obs - (k + 2)

cat("Starting simulation. This may take a few moments...\n")

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: ARDL Bounds test function
#++++++++++++++++++++++++++++++++++++++++++++++++
ardl_bounds <- function(y,x){
  T_obs   <- length(y)
  k       <- 1    # k is number of regressor
  y_lag  <- y[-T_obs]
  dy     <- diff(y)
  x_lag  <- x[-T_obs]
  dx     <- diff(x)
  # restricted model
  X_r <- cbind(1, dx)
  beta_r <- solve(t(X_r) %*% X_r) %*% t(X_r) %*% dy
  res_r <- dy - (X_r %*% beta_r)
  SSR_R <- sum(res_r^2)
  # unrestricted model
  X_ur <- cbind(1, dx, y_lag, x_lag)
  beta_ur <- solve(t(X_ur) %*% X_ur) %*% t(X_ur) %*% dy
  res_ur <- dy - (X_ur %*% beta_ur)
  SSR_ur <- sum(res_ur^2)
  # degree of freedom
  df_num <- k + 1
  df_den <- length(dy) - ncol(X_ur)

  bounds_stat <- ((SSR_R - SSR_ur) / df_num) / (SSR_ur / df_den)
  return(bounds_stat)
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: Run the Fourier ARDL Bounds function
#++++++++++++++++++++++++++++++++++++++++++++++++
for (i in 1:reps) {

  # A. Generate independent standard normal errors
  e1 <- rnorm(T_obs)
  e2 <- rnorm(T_obs)

  # B. Generate the dependent variable y (always a random walk under H0)
  y <- cumsum(e1)

  # C. Generate regressor x for the I(0) lower bound (P = 0)
  x_I0 <- e2

  # D. Generate regressor x for the I(1) upper bound (P = 1)
  x_I1 <- cumsum(e2)

  # E. Run the function
  F_stat_I0[i] <- ardl_bounds(y,x_I0)
  F_stat_I1[i] <- ardl_bounds(y,x_I1)
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 4: Report the findings
#++++++++++++++++++++++++++++++++++++++++++++++++
# Extract and Format the Asymptotic Critical Values
alpha_levels <- c(0.01, 0.025, 0.05, 0.10,0.50, 0.90, 0.95, 0.975, 0.99)

crit_vals_I0 <- quantile(F_stat_I0, probs = alpha_levels)
crit_vals_I1 <- quantile(F_stat_I1, probs = alpha_levels)

# Create a clean summary table
results_table <- data.frame(
  Significance = c("1%", "2.5%", "5%", "10%","50%", "10%", "5%", "2.5%", "1%"),
  `I(0) Lower Bound` = round(crit_vals_I0, 2),
  `I(1) Upper Bound` = round(crit_vals_I1, 2),
  row.names = NULL,
  check.names = FALSE
)

print("Monte Carlo Simulation Complete. Simulated Critical Values (k=1, Case III):")
print(results_table)
