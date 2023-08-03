
#Filter out berries
ps.rar.berry <- subset_samples(ps.phyF.tree, sample.type == "Berry")


#ordination
ord.bray.berry <- ordinate(ps.rar.berry, "PCoA", "bray", weighted=F)


#Axes
bray.1.2.berry <- plot_ordination(ps.rar.berry, 
                               ord.bray.berry, color="vineyard", 
                               axes = c(1, 2))

bray.1.3.berry <- plot_ordination(ps.rar.berry, 
                               ord.bray.berry, color="vineyard", 
                               axes = c(1, 3))

#Plot them
berry.br.1.2 <- bray.1.2.berry +
  scale_x_reverse() +
  geom_point(size=7.5) +
  geom_path() +
  geom_text(aes(label = stage), color = "black", fontface = "bold") +
  scale_colour_manual(values=colours) +
  theme_pcoa +
  theme(legend.position = "none")

berry.br.1.3 <- bray.1.3.berry +
  scale_x_reverse() +
  geom_point(size=7.5) +
  geom_path() +
  geom_text(aes(label = stage), color = "black", fontface = "bold") +
  scale_colour_manual(values=colours) +
  theme_pcoa +
  theme(legend.position = "none")

#PERMANOVA

library(vegan)

ps.rar.berry <- subset_samples(ps.phyF.prevF.rar, sample.type == "Berry")

metadata.berry <- as(sample_data(ps.rar.berry), "data.frame")

adonis2(distance(ps.rar.berry, method="bray") ~ vineyard, data = metadata.berry)


#merge!
ggarrange(berry.br.1.2, berry.br.1.3,
          nrow = 1, ncol = 2, common.legend = TRUE, labels = c("A)", "B)"), 
          font.label = list(size=18, face="bold", color="black"))

ggsave(filename = "R_output/Figure8.png", 
       width = 15, height = 7.5, dpi = 600)


#Make the matrices
pcoa_matrix <- ord.bray.berry$vectors
pcoa_matrix <- pcoa_matrix[, 1:5]
pcoa_matrix <- pcoa_matrix[2:9,]

pca_matrix <- read.csv("pca_matrix.csv", row.names = 1)


## Mantel Test ##

# Calculate the distance matrices from the PCA and PCoA matrices
pca_dist_matrix <- vegdist(pca_matrix)
pcoa_dist_matrix <- vegdist(pcoa_matrix)

# Add a constant value to make all elements of PCA and PCoA matrices non-negative
constant_value <- abs(min(c(min(pca_matrix), min(pcoa_matrix))))
pca_matrix_nonnegative <- pca_matrix + constant_value
pcoa_matrix_nonnegative <- pcoa_matrix + constant_value

# Calculate the distance matrices from the non-negative PCA and PCoA matrices
pca_dist_matrix <- vegdist(pca_matrix_nonnegative)
pcoa_dist_matrix <- vegdist(pcoa_matrix_nonnegative)

# Perform the Mantel test using the distance matrices
mantel_result <- mantel(pca_dist_matrix, pcoa_dist_matrix)

# Print the Mantel correlation coefficient and p-value
print(paste("Mantel correlation coefficient:", mantel_result$statistic))
print(paste("Mantel p-value:", mantel_result$signif))


## Procrustes Test ##

# Perform Procrustes analysis
procrustes_result <- procrustes(pca_matrix, pcoa_matrix)

# Get the Procrustes correlation coefficient
procrustes_corr <- 1 - procrustes_result$ss / sum(procrustes_result$X^2)

# Number of permutations (increase this for more accurate p-value estimation)
n_permutations <- 999

# Permutation testing to estimate p-value
permuted_correlations <- numeric(n_permutations)
for (i in 1:n_permutations) {
  shuffled_pcoa_matrix <- pcoa_matrix[sample(nrow(pcoa_matrix)), ]
  permuted_result <- procrustes(pca_matrix, shuffled_pcoa_matrix)
  permuted_corr <- 1 - permuted_result$ss / sum(permuted_result$X^2)
  permuted_correlations[i] <- permuted_corr
}

# Calculate the p-value (two-tailed test)
observed_correlation <- procrustes_corr
p_value <- sum(abs(permuted_correlations) >= abs(observed_correlation)) / n_permutations

# Print the Procrustes correlation coefficient and p-value
print(paste("Procrustes correlation coefficient:", observed_correlation))
print(paste("Procrustes p-value:", p_value))

