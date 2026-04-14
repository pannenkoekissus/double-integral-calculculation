# 1. Define your parameters (example numbers)
r_P <- 0.1
m_P <- 0.01
P_t <- 100

mu_P <- 5
sigma_P <- 1
mu_A <- 6
sigma_A <- 1.5

# New exact parameters from the paper text
c_cost <- 0.1 # proportionality constant for cost 'c'
A_t <- 50     # Animal abundance
M_P0 <- 10    # Half-saturation constant for Plant benefits
b_P <- 1.5    # Max benefit for Plant
sigma_c <- 1  # Attenuation of benefit

# Heaviside function (Eq 47)
H <- function(val) {
  ifelse(val >= 0, 1, 0)
}

# The f function (Eq 98)
f_func <- function(diff) {
  ifelse(diff < 0, exp(-(diff^2) / (2 * sigma_c^2)), 1)
}

# 2. Define the Benefit and Cost functions exactly as per text (1-plant, 1-animal case)
B_P <- function(x_P, y_A) {
  # M_P{i,j} (Eq 92)
  # The denominator represents the total plant population available to the animal y_A, 
  # evaluated at the population mean mu_P (adding 1e-9 to avoid div by zero if unaccessible)
  denom_f <- P_t * f_func(mu_P - y_A) + 1e-9
  M_P_val <- (f_func(x_P - y_A) / denom_f) * A_t
  
  # B_P (Eq 86, Single species case)
  b_P * (M_P_val / (M_P_val + M_P0))
}

C_P <- function(x_P, y_A) {
  # M^C_A (Eq 41)
  # Denominator is total plant population available to animal y_A based on population mean mu_P
  denom_H <- P_t * H(y_A - mu_P) + 1e-9
  M_C_A <- (H(y_A - x_P) / denom_H) * A_t
  
  # Costs (Eq 35)
  c_cost * M_C_A
}

# 3. Define the full integrand function
integrand <- function(x_P, y_A) {
  # The exponential growth factor from your Ricker model
  growth_factor <- exp(r_P + B_P(x_P, y_A) - C_P(x_P, y_A) - m_P * P_t)
  
  # The Normal probability density distributions N(x_P) and N(y_A)
  N_xP <- dnorm(x_P, mean = mu_P, sd = sigma_P)
  N_yA <- dnorm(y_A, mean = mu_A, sd = sigma_A)
  
  return(growth_factor * N_xP * N_yA)
}

# 4. Create a nested integral wrapper
# We use sapply to ensure it is vectorized over y_A for the outer integrate()
inner_integral <- function(y_A) {
  sapply(y_A, function(y) {
    integrate(function(x_P) integrand(x_P, y), lower = -Inf, upper = Inf)$value
  })
}

# 5. Compute the final double integral!
result <- integrate(inner_integral, lower = -Inf, upper = Inf)
print(paste("Numerical Integration Result:", result$value))
