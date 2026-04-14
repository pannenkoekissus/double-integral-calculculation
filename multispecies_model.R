# =========================================================
# Multi-Species Analytical Mutualistic Model
# Translating Eq 20 to 148 fully vectorized over S_P and S_A
# =========================================================

# --- 1. Define Example Vectors for Multiple Species ---

# --- 1. Define Example Vectors for Multiple Species ---

mu_P <- c(1:10) # Corollas of varied lengths
S_P <- length(mu_P) # Number of Plant Species

# We use rep() so all vectors automatically match the amount of plant species
P_t <- rep(100, S_P)       
sigma_P <- rep(1, S_P)     
r_P <- rep(0.1, S_P)       
m_P <- rep(0.01, S_P)      
b_P <- rep(1.5, S_P)       
M_P0 <- rep(10, S_P)       
c_cost <- 0.1             

mu_A <- c(1:10) # Proboscises of varied lengths
S_A <- length(mu_A) # Number of Animal Species

# Automatically scale vectors to match the amount of animal species
A_t <- rep(50, S_A)
sigma_A <- rep(1.2, S_A)
b_A <- rep(2.0, S_A)
M_A0 <- rep(15, S_A)
C_A0 <- 0.05
mu_A0 <- 10

sigma_c <- 1 # Benefit attenuation standard deviation

# --- 2. Calculate Continuous Interaction Matrices (H_int & f_int) ---
# H_int[j, i] represents animal j's continuous access to plant i
d_AP <- outer(mu_A, mu_P, "-") # S_A x S_P matrix
var_sum_AP <- outer(sigma_A^2, sigma_P^2, "+")
sd_sum_AP <- sqrt(var_sum_AP)
H_int <- pnorm(d_AP / sd_sum_AP)

# f_int[i, j] represents trait-overlap benefits plant i gets from animal j
d_PA <- outer(mu_P, mu_A, "-") # S_P x S_A matrix
var_sum_PA <- outer(sigma_P^2, sigma_A^2, "+")
var_total <- sigma_c^2 + var_sum_PA
sd_sum_PA <- sqrt(var_sum_PA)
sd_total <- sqrt(var_total)

part1 <- pnorm(d_PA / sd_sum_PA)
part2 <- (sigma_c / sd_total) * exp(-(d_PA^2) / (2 * var_total)) * pnorm((d_PA * sigma_c) / (sd_sum_PA * sd_total))
f_int <- part1 + part2


# --- 3. Compute Plant Fitness ---
# Plant Costs (Eq 35, 41)
# Denom_H calculates the competitive sum over all plants available to animal j
denom_H <- as.vector(H_int %*% P_t) + 1e-9 # Length S_A

# M^C_A Matrix (S_A x S_P): sweeping along rows (MARGIN=1) for animal-specific denominators
M_C_A <- sweep(H_int, 1, denom_H, "/")
M_C_A <- sweep(M_C_A, 1, A_t, "*")
C_P_val <- c_cost * colSums(M_C_A) # Summing all animal j interactions for each plant i

# Plant Benefits (Eq 86, 92)
# Denom_f calculates competitive sum over all plants accessible to animal j
denom_f <- as.vector(t(f_int) %*% P_t) + 1e-9 # Length S_A

# M_P Matrix (S_P x S_A): sweeping along cols (MARGIN=2) for animal-specific denominators
M_P_matrix <- sweep(f_int, 2, denom_f, "/")
M_P_matrix <- sweep(M_P_matrix, 2, A_t, "*")
M_P_total <- rowSums(M_P_matrix) # Summing all animal j interactions for each plant i

B_P_val <- b_P * (M_P_total / (M_P_total + M_P0))

# Resulting Plant Fitness (W_P) = P_{t+1, i} / P_{t, i}
W_P <- exp(r_P + B_P_val - C_P_val - m_P * P_t)


# --- 4. Compute Animal Fitness ---
# Animal Benefits (Eq 60, 66)
denom_H_anim <- as.vector(t(H_int) %*% A_t) + 1e-9 # Length S_P, competition amongst animals for plant i

# M_A Matrix (S_A x S_P): sweeping along cols (MARGIN=2) for plant-specific denominators
M_A_matrix <- sweep(H_int, 2, denom_H_anim, "/")
M_A_matrix <- sweep(M_A_matrix, 2, P_t, "*")
M_A_total <- rowSums(M_A_matrix)

B_A_val <- b_A * (M_A_total / (M_A_total + M_A0))

# Animal Costs (Eq 148)
C_A_val <- C_A0 * exp((mu_A / mu_A0) + 0.5 * (sigma_A / mu_A0)^2)

# Resulting Animal Fitness (W_A)
W_A <- exp(B_A_val - C_A_val)


# --- 5. Display the Multispecies Evaluation ---
print("----------- PLANT RESULTS -----------")
print(data.frame(
  Species = paste0("P", 1:S_P),
  Corolla_mu = mu_P,
  Benefits = B_P_val,
  Costs = C_P_val,
  Fitness_W = W_P
))

print("---------- ANIMAL RESULTS -----------")
print(data.frame(
  Species = paste0("A", 1:S_A),
  Proboscis_mu = mu_A,
  Benefits = B_A_val,
  Costs = C_A_val,
  Fitness_W = W_A
))
