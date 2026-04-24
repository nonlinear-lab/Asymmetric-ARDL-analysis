#+++++++++++++++++++++++++++++++++
# Step 0: set a dorking directory
#+++++++++++++++++++++++++++++++
setwd("C:/R/ardl")

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: load data
#++++++++++++++++++++++++++++++++++++++++++++++++
my_data <- read.csv("data(2).csv")
my_data <- as.matrix(my_data)
T_obs   <- nrow(my_data)
mout    <- matrix(NA,nrow=5,ncol=1)

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: Asymmetric bounds test function
#++++++++++++++++++++++++++++++++++++++++++++++++
asymmetric_ardl_bounds <- function(y,x){
  y <- as.numeric(y)
  x <- as.numeric(x)
  T_obs <- length(y)
  k       <- 2    # k is number of regressor
  dy     <- diff(y)
  dx     <- c(0,diff(x))
  dx_pos <- ifelse(dx > 0, dx, 0)
  dx_neg <- ifelse(dx < 0, dx, 0)
  x_pos  <- cumsum(dx_pos)
  x_neg  <- cumsum(dx_neg)

  dy_t       <- dy
  dx_pos_t   <- dx_pos[-1]
  dx_neg_t   <- dx_neg[-1]
  y_lag      <- y[-T_obs]
  x_pos_lag  <- x_pos[-T_obs]
  x_neg_lag  <- x_neg[-T_obs]

  # restricted model
  X_r <- cbind(1, dx_pos_t, dx_neg_t)
  beta_r <- solve(t(X_r) %*% X_r) %*% t(X_r) %*% dy_t
  res_r <- dy_t - (X_r %*% beta_r)
  SSR_R <- sum(res_r^2)
  # unrestricted model
  X_ur <- cbind(1, dx_pos_t, dx_neg_t, y_lag, x_pos_lag, x_neg_lag)
  beta_ur <- solve(t(X_ur) %*% X_ur) %*% t(X_ur) %*% dy_t
  res_ur <- dy_t - (X_ur %*% beta_ur)
  SSR_ur <- sum(res_ur^2)
  # degree of freedom
  df_num <- k + 1
  df_den <- length(dy_t) - ncol(X_ur)

  bounds_stat <- ((SSR_R - SSR_ur) / df_num) / (SSR_ur / df_den)
  return(bounds_stat)
}
#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: Run the analysis
#++++++++++++++++++++++++++++++++++++++++++++++++

for (j in 1:5) {
  if (j == 1) {
    y <- as.matrix(my_data[, 3])
    x <- as.matrix(my_data[, 2])
  } else if (j == 2) {
    y <- as.matrix(my_data[, 4])
    x <- as.matrix(my_data[, 2])
  } else if (j == 3) {
    y <- as.matrix(my_data[, 5])
    x <- as.matrix(my_data[, 2])
  } else if (j == 4) {
    y <- as.matrix(my_data[, 6])
    x <- as.matrix(my_data[, 2])
  } else {
    y <- as.matrix(my_data[, 7])
    x <- as.matrix(my_data[, 2])
  }
  results <- asymmetric_ardl_bounds (y,x)
  mout[j,1]=results
}

print(" bounds test (k=1, Case III):")
print(mout)
