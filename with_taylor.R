library(numDeriv)
r_P <- 0.1
m_P <- 0.01
P_t <- 100

mu_P <- 5
sigma_P <- 1
mu_A <- 6
sigma_A <- 1.5
# New exact parameters from the paper text
c_cost <- 0.1 # proportionality constant for cost 'c'
A_t <- 50 # Animal abundance
M_P0 <- 10 # Half-saturation constant for Plant benefits
b_P <- 1.5 # Max benefit for Plant
sigma_c <- 1 # Attenuation of benefit

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

# 3. The mean traits corresponding to the zero-th order point
means <- c(mu_P, mu_A)

# 4. Calculate Zero-th Order Term (Replacing x_P and y_A with means)
term0 <- g(means)

# 5. Calculate Second-Order Term
# We calculate the Hessian matrix (second partial derivatives) at the means
H <- hessian(g, x = means)

# Because x_P and y_A are independent normals, the covariance is 0.
# The expected value of the second order term is: 1/2 * (g_xx * var_x + g_yy * var_y)
term2 <- 0.5 * (H[1, 1] * sigma_P^2 + H[2, 2] * sigma_A^2)

# 6. Final approximated value
approx_value <- term0 + term2

print(paste("Zero-th order (Means only):", term0))
print(paste("Second-order Taylor Approximation:", approx_value))
