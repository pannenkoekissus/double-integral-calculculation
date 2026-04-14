library(numDeriv)

# New exact parameters from the paper text
c_cost <- 0.1 # proportionality constant for cost 'c'
A_t <- 50     # Animal abundance
M_P0 <- 10    # Half-saturation constant for Plant benefits
b_P <- 1.5    # Max benefit for Plant
sigma_c <- 1  # Attenuation of benefit

r_P <- 0.1
m_P <- 0.01
P_t <- 100

mu_P <- 5
sigma_P <- 1
mu_A <- 6
sigma_A <- 1.5

# Heaviside function (Eq 47)
H <- function(val) {
  ifelse(val >= 0, 1, 0)
}

# The f function (Eq 98)
f_func <- function(diff) {
  ifelse(diff < 0, exp(-(diff^2) / (2 * sigma_c^2)), 1)
}

# 1. Define the Benefit and Cost functions exactly as per text
B_P <- function(x_P, y_A) {
  denom_f <- P_t * f_func(mu_P - y_A) + 1e-9
  M_P_val <- (f_func(x_P - y_A) / denom_f) * A_t
  b_P * (M_P_val / (M_P_val + M_P0))
}

C_P <- function(x_P, y_A) {
  denom_H <- P_t * H(y_A - mu_P) + 1e-9
  M_C_A <- (H(y_A - x_P) / denom_H) * A_t
  c_cost * M_C_A
}

# 2. Define the exponentiated growth function g(x_P, y_A)
g <- function(v) {
  x_P <- v[1]
  y_A <- v[2]
  exp(r_P + B_P(x_P, y_A) - C_P(x_P, y_A) - m_P * P_t)
}

means <- c(mu_P, mu_A)

# --- 0th Order Term ---
term0 <- g(means)

# --- 2nd Order Term ---
H_matrix <- hessian(g, x = means)
term2 <- 0.5 * (H_matrix[1, 1] * sigma_P^2 + H_matrix[2, 2] * sigma_A^2)

# --- 4th Order Term ---
# The 3rd order terms mathematically vanish (evaluate to 0) because normal distributions are symmetric.
# To compute computationally complex 4th derivatives, we wrap the hessian!
g_xx <- function(v) hessian(g, v)[1,1]
g_yy <- function(v) hessian(g, v)[2,2]

H_xx <- hessian(g_xx, means) # Computes g_xxxx and g_xxyy
H_yy <- hessian(g_yy, means) # Computes g_yxyy and g_yyyy

g_xxxx <- H_xx[1,1]
g_xxyy <- H_xx[2,2] 
g_yyyy <- H_yy[2,2]

# The mathematical expectation for independent 4th moments of Normal distributions:
term4 <- (1/8) * g_xxxx * (sigma_P^4) + 
         (1/4) * g_xxyy * (sigma_P^2 * sigma_A^2) + 
         (1/8) * g_yyyy * (sigma_A^4)

approx_value_4th <- term0 + term2 + term4

print(paste("Zero-th order (Means only):", term0))
print(paste("2nd-order Taylor Approx:", term0 + term2))
print(paste("4th-order Taylor Approx:", approx_value_4th))
