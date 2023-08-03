library(vegan)

# Suppose you have two matrices representing the PCA and PCoA results
# Example data for illustration purposes (replace with your actual data)
pca_matrix <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)  # PCA data (replace with your PCA data)
pcoa_matrix <- matrix(c(2, 3, 4, 5, 6, 7), ncol = 2)  # PCoA data (replace with your PCoA data)

# Perform Procrustes analysis
procrustes_result <- procrustes(pca_matrix, pcoa_matrix)

# Get the Procrustes correlation coefficient (r-value)
procrustes_corr <- procrustes_result$common

# Perform a permutation test with 1000 permutations
permutation_result <- protest(pca_matrix, pcoa_matrix, permutations = 1000)

# Get the p-value from the permutation test
p_value <- permutation_result$p.value

# Print the results
print(paste("Procrustes correlation coefficient (r-value):", procrustes_corr))
print(paste("Permutation test p-value:", p_value))
