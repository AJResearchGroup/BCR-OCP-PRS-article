#Author: Christina Chatsatourian 

# PRS_Stratification_Statistics_Trend_of_BC_Risk_Across_PRS

# Load necessary libraries
library(dplyr)

# Ensure no missing values in PRS_BC
filtered_data <- filtered_data %>%
  filter(!is.na(prs))

# Exclude NA from Age at Menarche
filtered_data <- filtered_data[!is.na(filtered_data$AgeMenarch), ]

# Exclude NA from BMI
filtered_data <- filtered_data[!is.na(filtered_data$bmi), ]

# Exclude NA from TDI
filtered_data <- filtered_data[!is.na(filtered_data$tdi), ]

# Only run this line when you do the ever-never analysis_Exclude NA from number_children (Don't do this steps for during-after)
filtered_data <- filtered_data[!is.na(filtered_data$number_children), ]

# Load necessary libraries
library(ggplot2)
library(pROC)

# Descriptive Statistics
# Summary statistics for PRS_BC
summary(filtered_data$prs)
mean_PRS <- mean(filtered_data$prs, na.rm = TRUE)
sd_PRS <- sd(filtered_data$prs, na.rm = TRUE)

# Print mean and standard deviation
cat("Mean PRS:", mean_PRS, "\n")
cat("Standard Deviation of PRS:", sd_PRS, "\n")

# Plot PRS distribution
ggplot(filtered_data, aes(x = prs)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(title = "Distribution of PRS for Breast Cancer", x = "PRS_BC", y = "Count")

# Stratification into Deciles
# Calculate decile thresholds
decile_thresholds <- quantile(filtered_data$prs, probs = seq(0, 1, 0.1), na.rm = TRUE)

# Create PRS decile groups
filtered_data <- filtered_data %>%
  mutate(PRS_decile = cut(prs, breaks = decile_thresholds, labels = 1:10, include.lowest = TRUE))

# Check the distribution of participants in each decile
table(filtered_data$PRS_decile)

# Logistic Regression

# Convert PRS_decile to a factor for logistic regression
filtered_data$PRS_decile <- as.factor(filtered_data$PRS_decile)

# Run logistic regression: BC_status is the binary outcome variable
model <- glm(BC ~ PRS_decile, family = binomial(link = "logit"), data = filtered_data)

# Summary of the logistic regression model
summary(model)

# Extract odds ratios and confidence intervals
library(broom)
ORs <- exp(cbind(OR = coef(model), confint(model)))

# Print odds ratios and confidence intervals
print(ORs)

# Assess Predictive Ability (AUC)
# Predicted probabilities from the logistic model
pred <- predict(model, type = "response")

# Compute ROC and AUC
roc_obj <- roc(filtered_data$BC, pred)

# Print AUC
cat("AUC:", auc(roc_obj), "\n")

# Plot the ROC curve
plot(roc_obj, col = "blue", main = "ROC Curve for PRS Prediction")

# 6. Trend Analysis Across Deciles
# Calculate mean breast cancer risk per decile
risk_by_decile <- filtered_data %>%
  group_by(PRS_decile) %>%
  summarize(mean_risk = mean(BC, na.rm = TRUE))

# Plot the trend
ggplot(risk_by_decile, aes(x = as.numeric(PRS_decile), y = mean_risk)) +
  geom_line(color = "blue") +
  geom_point(size = 3, color = "red") +
  labs(x = "PRS Decile", y = "Mean Breast Cancer Risk", 
       title = "Trend of Breast Cancer Risk Across PRS Deciles") +
  theme_minimal()
