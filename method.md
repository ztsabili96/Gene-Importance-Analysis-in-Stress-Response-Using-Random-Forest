# Methods: Gene Importance Analysis in Stress Response

## 1. Data Preprocessing

Raw RNA-seq expression data were processed to construct a sample-wise expression matrix. Genes with zero total expression across all samples were removed. To reduce dimensionality and focus on informative features, the top 500 genes with the highest variance were selected. Sample labels were assigned as either "control" or "stress" based on experimental conditions. Expression values were log-transformed and normalized prior to model training.

## 2. Feature Selection

Unsupervised filtering based on variance was applied to retain genes with the most dynamic expression profiles. No supervised selection (e.g., based on class labels) was performed prior to model training to avoid bias and ensure fair evaluation of feature importance.

## 3. Model Training

A Random Forest classifier was implemented using the `randomForest` package in R. The model was trained on the reduced expression matrix containing 500 genes and binary class labels ("control" vs. "stress"). Feature importance scores were extracted using the MeanDecreaseGini metric, which reflects each gene’s contribution to classification accuracy.

## 4. Model Evaluation

Model performance was assessed using multiple metrics:

- A confusion matrix and classification statistics were computed via `caret::confusionMatrix`.  
- ROC curve and Area Under the Curve (AUC) were calculated using the `pROC` package.  
- Stratified 3-fold cross-validation was performed to evaluate model stability and generalizability.

## 5. Gene Importance Analysis

The top 20 genes were ranked based on their MeanDecreaseGini scores. Among them, **Os10g0417600** emerged as the most influential gene in distinguishing stress from control samples. Mean expression analysis revealed higher expression of Os10g0417600 in control samples (~9376.67) compared to stress (~5763.33), suggesting potential downregulation under stress conditions.

## 6. Reproducibility and Documentation

All scripts used in this analysis are fully documented and reproducible. Output files, including feature importance tables, confusion matrices, and visualizations (e.g., ROC curves, gene ranking plots), are saved in the `results/` and `figures/` directories. The project structure and workflow are described in the accompanying `README.md` file.



