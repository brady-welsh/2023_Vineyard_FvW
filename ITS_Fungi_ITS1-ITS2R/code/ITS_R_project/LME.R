# Load necessary libraries
library(lme4)
library(multcomp)
library(dplyr)

ps.inoc <- subset(ps.phyF.prevF.rar.div.df, sample.type != "Wild")
ps.wild <- subset(ps.phyF.prevF.rar.div.df, sample.type != "Inoculated")

## LME for ASV richness ##
# Inoculated Wines #
# Fit mixed-effects models for each vineyard separately
model_organic.in <- lm(observed ~ stage, data = subset(ps.inoc, vineyard == "Organic"))
model_biodynamic.in <- lm(observed ~ stage, data = subset(ps.inoc, vineyard == "Biodynamic"))
model_conventional.in <- lm(observed ~ stage, data = subset(ps.inoc, vineyard == "Conventional"))

# View the summaries of each model
summary(model_organic.in)
summary(model_biodynamic.in)
summary(model_conventional.in)

# Wild Wines
# Fit mixed-effects models for each vineyard separately
model_organic.wi <- lm(observed ~ stage, data = subset(ps.wild, vineyard == "Organic"))
model_biodynamic.wi <- lm(observed ~ stage, data = subset(ps.wild, vineyard == "Biodynamic"))
model_conventional.wi <- lm(observed ~ stage, data = subset(ps.wild, vineyard == "Conventional"))

# View the summaries of each model
summary(model_organic.wi)
summary(model_biodynamic.wi)
summary(model_conventional.wi)

model.list <- list(summary(model_organic.wi)$coefficients, summary(model_organic.in)$coefficients, 
                   summary(model_biodynamic.wi)$coefficients, summary(model_biodynamic.in)$coefficients, 
                   summary(model_conventional.wi)$coefficients, summary(model_conventional.in)$coefficients)


# Create a data frame from the extracted data
model.list.df <- data.frame(
  Coefficient = rownames(model.list),
  Estimate = coef_data[, "Estimate"],
  Std_Error = coef_data[, "Std. Error"],
  P_Value = coef_data[, "Pr(>|t|)"]
)

write.csv(model.list, file = "model-list_ASV-richness.csv", row.names = FALSE)


## LME for Shannon diversity ##
# Inoculated Wines #
# Fit mixed-effects models for each vineyard separately
model_organic.in <- lm(diversity_shannon ~ stage, data = subset(ps.inoc, vineyard == "Organic"))
model_biodynamic.in <- lm(diversity_shannon ~ stage, data = subset(ps.inoc, vineyard == "Biodynamic"))
model_conventional.in <- lm(diversity_shannon ~ stage, data = subset(ps.inoc, vineyard == "Conventional"))

# View the summaries of each model
summary(model_organic.in)
summary(model_biodynamic.in)
summary(model_conventional.in)

# Wild Wines #
# Fit mixed-effects models for each vineyard separately
model_organic.wi <- lm(diversity_shannon ~ stage, data = subset(ps.wild, vineyard == "Organic"))
model_biodynamic.wi <- lm(diversity_shannon ~ stage, data = subset(ps.wild, vineyard == "Biodynamic"))
model_conventional.wi <- lm(diversity_shannon ~ stage, data = subset(ps.wild, vineyard == "Conventional"))

# View the summaries of each model
summary(model_organic.wi)
summary(model_biodynamic.wi)
summary(model_conventional.wi)

model.list <- list(summary(model_organic.wi)$coefficients, summary(model_organic.in)$coefficients, 
                   summary(model_biodynamic.wi)$coefficients, summary(model_biodynamic.in)$coefficients, 
                   summary(model_conventional.wi)$coefficients, summary(model_conventional.in)$coefficients)


# Create a data frame from the extracted data
model.list.df <- data.frame(
  Coefficient = rownames(model.list),
  Estimate = coef_data[, "Estimate"],
  Std_Error = coef_data[, "Std. Error"],
  P_Value = coef_data[, "Pr(>|t|)"]
)

write.csv(model.list, file = "model-list_Shannon-diversity.csv", row.names = FALSE)