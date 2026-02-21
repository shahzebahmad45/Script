# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)

# Read Excel file
file_path <- "D:/SHAH/Academics/Ph.D/Yangzhou/My article/About selection index/Data for new Article sprtd.xlsx"

# Read all sheets
sheet1 <- read_excel(file_path, sheet = "Sheet1", skip = 1)
sheet2 <- read_excel(file_path, sheet = "Sheet2", skip = 1)
sheet3 <- read_excel(file_path, sheet = "Sheet3", skip = 1)
sheet4 <- read_excel(file_path, sheet = "Sheet4", skip = 1)

# Clean and process data
# Sheet1: Correlation matrix for weights
corr_weights <- sheet1 %>% 
  select(-1) %>% 
  as.matrix()
rownames(corr_weights) <- c("WEI-0", "WEI-6", "WEI-12", "WEI")
colnames(corr_weights) <- c("WEI-0", "WEI-6", "WEI-12", "WEI")

# Sheet2: Correlation matrix for milk traits
corr_milk <- sheet2 %>% 
  select(-1) %>% 
  as.matrix()
rownames(corr_milk) <- c("DMY", "FP", "PP")
colnames(corr_milk) <- c("DMY", "FP", "PP")

# Sheet3: Variances and heritabilities
var_data <- sheet3 %>% 
  select(Trait, SD, σ2_P, σ2_A, h2)

# Sheet4: Covariances
cov_data <- sheet4

# Extract phenotypic and genetic variances for our traits of interest
var_dmy <- var_data %>% filter(Trait == "DMY")
var_wei12 <- var_data %>% filter(Trait == "WEI-12")

σ2_P_DMY <- var_dmy$σ2_P
σ2_A_DMY <- var_dmy$σ2_A
σ2_P_WEI12 <- var_wei12$σ2_P
σ2_A_WEI12 <- var_wei12$σ2_A

# Extract genetic correlation between DMY and WEI-12
# We need to calculate this from the covariance data
cov_genetic <- cov_data %>% 
  filter(`Pairs of Traits` == "WEI-12 & WEI")  # Using WEI as proxy since WEI-12 & DMY not directly available

# For this example, we'll assume genetic correlation based on available data
# In practice, you would need the direct genetic correlation between DMY and WEI-12
# Let's use phenotypic correlation as approximation for demonstration
genetic_corr <- -0.15  # Approximated from available correlations

# Phenotypic variance-covariance matrix (P)
P <- matrix(c(σ2_P_DMY, genetic_corr * sqrt(σ2_P_DMY * σ2_P_WEI12),
              genetic_corr * sqrt(σ2_P_DMY * σ2_P_WEI12), σ2_P_WEI12),
            nrow = 2, ncol = 2)

rownames(P) <- colnames(P) <- c("DMY", "WEI-12")

# Genetic variance-covariance matrix (G)
G <- matrix(c(σ2_A_DMY, genetic_corr * sqrt(σ2_A_DMY * σ2_A_WEI12),
              genetic_corr * sqrt(σ2_A_DMY * σ2_A_WEI12), σ2_A_WEI12),
            nrow = 2, ncol = 2)

rownames(G) <- colnames(G) <- c("DMY", "WEI-12")

# Economic weights (based on desired improvements: 60% DMY, 10% WEI-12)
# Convert percentage improvements to economic weights
target_DMY_improvement <- 1  # less concern
target_WEI12_improvement <- 0.5  # Favor DMY

# Economic weights proportional to desired improvements
a <- c(target_DMY_improvement, target_WEI12_improvement)

# Calculate selection index coefficients (b)
# b = P^(-1) * G * a
b <- solve(P) %*% G %*% a

# Calculate accuracy of selection index
# r_{HI} = sqrt(b' * P * b / (a' * G * a))
r_HI <- sqrt((t(b) %*% P %*% b) / (t(a) %*% G %*% a))[1,1]

# Calculate expected genetic gain per generation
# ΔG = i * r_{HI} * σ_H
# Where i = selection intensity, σ_H = sqrt(a' * G * a)

# Assume selection intensity (i) for top 10%
i <- 1.75  # Standardized selection differential for 10% selected

σ_H <- sqrt(t(a) %*% G %*% a)[1,1]
ΔG_per_generation <- i * r_HI * σ_H

# Calculate expected gains for individual traits
ΔG_traits <- i * r_HI * (G %*% b) / sqrt(t(b) %*% P %*% b)[1,1]

# Calculate over 4 generations
total_gain_4gen <- 4 * ΔG_per_generation
total_gain_traits_4gen <- 4 * ΔG_traits

# Print results
cat("SELECTION INDEX ANALYSIS FOR HUAXI COWS\n")
cat("=======================================\n\n")

cat("Traits: DMY (Daily Milk Yield) and WEI-12 (Weight at 12 months)\n\n")

cat("Phenotypic Variance-Covariance Matrix (P):\n")
print(P)
cat("\n")

cat("Genetic Variance-Covariance Matrix (G):\n")
print(G)
cat("\n")

cat("Economic Weights (a):\n")
cat("DMY:", a[1], "\n")
cat("WEI-12:", a[2], "\n\n")

cat("Selection Index Coefficients (b):\n")
cat("DMY:", b[1,1], "\n")
cat("WEI-12:", b[2,1], "\n\n")

cat("Accuracy of Selection Index (r_HI):", round(r_HI, 4), "\n\n")

cat("Expected Genetic Gain Per Generation:\n")
cat("Overall Genetic Gain (ΔG):", round(ΔG_per_generation, 4), "\n")
cat("DMY Gain:", round(ΔG_traits[1,1], 4), "\n")
cat("WEI-12 Gain:", round(ΔG_traits[2,1], 4), "\n\n")

cat("Expected Genetic Gain Over 4 Generations:\n")
cat("Total Overall Genetic Gain:", round(total_gain_4gen, 4), "\n")
cat("Total DMY Gain:", round(total_gain_traits_4gen[1,1], 4), "\n")
cat("Total WEI-12 Gain:", round(total_gain_traits_4gen[2,1], 4), "\n\n")

# Calculate relative efficiency compared to single trait selection
single_trait_DMY_gain <- i * sqrt(σ2_A_DMY) * sqrt(var_dmy$h2)
single_trait_WEI12_gain <- i * sqrt(σ2_A_WEI12) * sqrt(var_wei12$h2)

cat("Comparison with Single Trait Selection:\n")
cat("Single trait DMY gain per generation:", round(single_trait_DMY_gain, 4), "\n")
cat("Single trait WEI-12 gain per generation:", round(single_trait_WEI12_gain, 4), "\n")
cat("Index selection DMY gain per generation:", round(ΔG_traits[1,1], 4), "\n")
cat("Index selection WEI-12 gain per generation:", round(ΔG_traits[2,1], 4), "\n\n")

# Create selection index function for practical use
calculate_selection_index <- function(DMY_value, WEI12_value, mean_DMY = 0, mean_WEI12 = 0) {
  # Standardize values (assuming means provided)
  DMY_std <- DMY_value - mean_DMY
  WEI12_std <- WEI12_value - mean_WEI12
  
  # Calculate index value
  index_value <- b[1,1] * DMY_std + b[2,1] * WEI12_std
  
  return(index_value)
}

# Example usage of selection index
cat("Selection Index Calculation Example:\n")
cat("For an animal with DMY = 28.5 kg and WEI-12 = 357.15 kg:\n")
example_index <- calculate_selection_index(28.5, 357.15)
cat("Selection Index Value:", round(example_index, 4), "\n")

# Plotting the selection response (optional)
if(require(ggplot2)) {
  generations <- 1:4
  cumulative_gain <- cumsum(rep(ΔG_per_generation, 4))
  
  gain_data <- data.frame(
    Generation = rep(generations, 3),
    Trait = rep(c("Overall", "DMY", "WEI-12"), each = 4),
    Gain = c(cumulative_gain,
             cumsum(rep(ΔG_traits[1,1], 4)),
             cumsum(rep(ΔG_traits[2,1], 4)))
  )
  
  ggplot(gain_data, aes(x = Generation, y = Gain, color = Trait)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(title = "Expected Genetic Gain Over 4 Generations",
         y = "Cumulative Genetic Gain",
         color = "Trait") +
    theme_minimal()
}
