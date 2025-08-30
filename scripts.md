# Load required libraries
library(randomForest)
library(pROC)
library(caret)
library(ggplot2)
library(lattice)

set.seed(42)  # For reproducibility

# -------------------------------
# 1. Load and preprocess expression data
# -------------------------------

# Read expression matrix
expr <- read.table("data/final_matrix-old.txt", header = TRUE, sep = "\t", row.names = 1)

# Clean column names (remove file suffix)
colnames(expr) <- gsub("_mapped.bam", "", colnames(expr))

# Remove genes with zero total expression
expr_filtered <- expr[rowSums(expr) > 0, ]

# Transpose matrix: rows = samples, columns = genes
expr_transposed <- t(expr_filtered)

# Load sample labels
labels <- read.csv("data/labels.csv", row.names = 1)

# Merge expression data with condition labels
expr_final <- data.frame(expr_transposed)
expr_final$Condition <- labels[rownames(expr_final), "Condition"]

# Save merged dataset
write.table(expr_final, file = "data/expr_final.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------------
# 2. Feature selection: top 500 variable genes
# -------------------------------

# Due to large data size and memory limits, we filter to most informative genes
variances <- apply(expr_final[, -ncol(expr_final)], 2, var)
top_genes <- names(sort(variances, decreasing = TRUE))[1:500]

# Create reduced dataset
expr_reduced <- expr_final[, c(top_genes, "Condition")]
expr_reduced$Condition <- as.factor(expr_reduced$Condition)

# -------------------------------
# 3. Train Random Forest model
# -------------------------------

rf_model <- randomForest(Condition ~ ., data = expr_reduced, importance = TRUE)
print(rf_model)

# Save model and feature importance
saveRDS(rf_model, file = "results/rf_model.rds")
write.table(importance(rf_model), file = "results/feature_importance.tsv", sep = "\t", quote = FALSE)

# Save confusion matrix
conf_matrix <- table(Predicted = predict(rf_model, expr_reduced),
                     Actual    = expr_reduced$Condition)
write.table(conf_matrix, file = "results/confusion_matrix.tsv", sep = "\t", quote = FALSE)

# -------------------------------
# 4. ROC Curve and AUC
# -------------------------------

# Predict class probabilities
prob_stress <- predict(rf_model, expr_reduced, type = "prob")[, "stress"]

# Compute ROC and AUC
roc_obj <- roc(expr_reduced$Condition, prob_stress, levels = c("control", "stress"))
plot(roc_obj, col = "#1f77b4", lwd = 2, main = "ROC Curve: stress vs control")
auc_val <- auc(roc_obj)
cat("AUC =", round(auc_val, 3), "\n")

# Save ROC curve (base R)
png("figures/roc_curve_base.png", width = 800, height = 600)
plot(roc_obj, col = "#1f77b4", lwd = 2, main = "ROC Curve: stress vs control")
dev.off()

# Save ROC curve (ggplot2)
roc_df <- data.frame(
  tpr = roc_obj$sensitivities,
  fpr = 1 - roc_obj$specificities
)
ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_line(color = "#1f77b4", size = 1.2) +
  geom_abline(linetype = "dashed", color = "grey") +
  labs(title = "ROC Curve", x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal()
ggsave("figures/roc_curve_ggplot.png", width = 6, height = 5)

# Extract optimal threshold
best_coords <- coords(roc_obj, x = "best", best.method = "youden",
                      ret = c("threshold", "sensitivity", "specificity"))
print(best_coords)

# -------------------------------
# 5. Cross-Validation
# -------------------------------

# Ensure factor levels are correct
expr_reduced$Condition <- factor(expr_reduced$Condition, levels = c("control", "stress"))

# Define CV settings
ctrl <- trainControl(
  method          = "cv",
  number          = 3,
  classProbs      = TRUE,
  summaryFunction = defaultSummary
)

# Train model with CV
cv_model <- train(
  Condition ~ ., 
  data      = expr_reduced,
  method    = "rf",
  metric    = "Accuracy",
  trControl = ctrl
)
print(cv_model)

# -------------------------------
# 6. Feature Importance from CV model
# -------------------------------

imp <- varImp(cv_model)
print(imp)

# Save importance table
imp_df <- data.frame(Gene = rownames(imp$importance), imp$importance)
write.table(imp_df, "results/feature_importance_cv.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Plot top 20 important genes
varImpPlot(cv_model$finalModel, n.var = 20, main = "Top 20 Important Genes (CV)")

# -------------------------------
# 7. Mean expression per class
# -------------------------------

# Compute mean expression of each gene in control vs stress
agg <- aggregate(expr_reduced[, -ncol(expr_reduced)],
                 by = list(Class = expr_reduced$Condition),
                 FUN = mean)
write.table(agg, "results/mean_expression_by_class.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
