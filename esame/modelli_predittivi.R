library(DataExplorer)     # IDA
library(dplyr)
library(ggstatsplot)
library(dlookr)
library(flextable)
library(ggpubr)
library(SmartEDA)
library(performance)
library(tidyr)
if (!require(pwr)) install.packages("pwr")
library(pwr)
install.packages("effsize") # Uncomment and run once if you haven't installed it
library(effsize)
install.packages("caret")
library(caret)
library(pROC)
library(rpart)



injury_df <- read.csv('college_injuries.csv', header = T, sep = ',')
injury_df$Sport_Type <- as.factor(injury_df$Sport_Type)
injury_df$Activity_Type <- as.factor(injury_df$Activity_Type)
injury_df$Injury_Occurred <- as.factor(injury_df$Injury_Occurred)


# imposto la k-fold cross validation
# --- Impostazioni per la K-fold Cross Validation ---
set.seed(123) # Per riproducibilità
k <- 5 # Numero di fold
folds <- createFolds(injury_df$Injury_Occurred, k = k, list = TRUE, returnTrain = FALSE)

# --- Funzione per valutare il modello e raccogliere le metriche ---
evaluate_model <- function(model, test_data, actual_labels) {
  prob <- predict(model, newdata = test_data, type = "response")
  pred_class <- ifelse(prob > 0.5, 1, 0)
  
  conf_matrix <- table(Predicted = pred_class, Actual = actual_labels)
  
  TP <- ifelse(length(conf_matrix[2, 2]) > 0, conf_matrix[2, 2], 0)
  TN <- ifelse(length(conf_matrix[1, 1]) > 0, conf_matrix[1, 1], 0)
  FP <- ifelse(length(conf_matrix[2, 1]) > 0, conf_matrix[2, 1], 0)
  FN <- ifelse(length(conf_matrix[1, 2]) > 0, conf_matrix[1, 2], 0)
  
  accuracy <- (TP + TN) / sum(conf_matrix)
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  # Calcola AUC solo se ci sono entrambe le classi nel test set
  auc_val <- NA
  if(length(unique(actual_labels)) == 2 && length(unique(prob)) > 1) {
    roc_obj <- roc(actual_labels, prob, quiet = TRUE)
    auc_val <- auc(roc_obj)
  }
  
  return(list(accuracy = accuracy, sensitivity = sensitivity, specificity = specificity, auc = auc_val))
}


# --- Funzione per valutare il modello Decision Tree e raccogliere le metriche ---
evaluate_dt_model <- function(model, test_data, actual_labels) {
  # Per rpart, predict(type = "prob") restituisce una matrice con le probabilità per ogni classe.
  # Vogliamo la probabilità della classe "1" (seconda colonna).
  prob <- predict(model, newdata = test_data, type = "prob")[, "1"]
  pred_class <- factor(ifelse(prob > 0.5, "1", "0"), levels = c("0", "1"))
  
  # Usa confusionMatrix dal pacchetto caret per calcolare metriche più complete
  conf_matrix_obj <- confusionMatrix(pred_class, actual_labels, positive = "1")
  
  accuracy <- conf_matrix_obj$overall["Accuracy"]
  # Precision e Recall sono accessibili tramite $byClass
  precision <- conf_matrix_obj$byClass["Pos Pred Value"] # Precision
  recall <- conf_matrix_obj$byClass["Sensitivity"]      # Recall (Sensibilità)
  specificity <- conf_matrix_obj$byClass["Specificity"]  # Specificity
  
  # Calcola AUC solo se ci sono entrambe le classi nel test set e probabilità valide
  auc_val <- NA
  if(length(unique(actual_labels)) == 2 && length(unique(prob)) > 1) {
    roc_obj <- roc(actual_labels, prob, quiet = TRUE, levels = c("0", "1"))
    auc_val <- auc(roc_obj)
  }
  
  return(list(accuracy = accuracy, precision = precision, recall = recall,
              specificity = specificity, auc = auc_val,
              pred_prob = prob, actual_labels = actual_labels))
}

#############################################################################
# LOGISTIC REGRESSION VERSIONE BASE (solo 2 feature)
#############################################################################
all_metrics_base <- data.frame()

for (i in 1:k) {
  cat(paste0("\n--- Fold ", i, " (Modello Base) ---\n"))
  test_indices <- folds[[i]]
  train_data <- injury_df[-test_indices, ]
  test_data <- injury_df[test_indices, ]
  
  # Training del modello
  pred1 <- glm(data = train_data, Injury_Occurred ~ Impact_Force_Newtons + Cumulative_Fatigue_Index, family = binomial)
  
  # Valutazione
  metrics <- evaluate_model(pred1, test_data, test_data$Injury_Occurred)
  all_metrics_base <- rbind(all_metrics_base, metrics)
  
  cat("Accuracy: ", round(metrics$accuracy, 3), "\n")
  cat("Sensitivity (Recall): ", round(metrics$sensitivity, 3), "\n")
  cat("Specificity: ", round(metrics$specificity, 3), "\n")
  cat("AUC: ", round(metrics$auc, 3), "\n")
}

# --- Media delle metriche per il modello base ---
mean_accuracy_base <- mean(all_metrics_base$accuracy, na.rm = TRUE)
mean_sensitivity_base <- mean(all_metrics_base$sensitivity, na.rm = TRUE)
mean_specificity_base <- mean(all_metrics_base$specificity, na.rm = TRUE)
mean_auc_base <- mean(all_metrics_base$auc, na.rm = TRUE)

cat("\n--- Medie delle Metriche (Modello Base) ---\n")
cat("Accuracy Media: ", round(mean_accuracy_base, 3), "\n")
cat("Sensitivity Media: ", round(mean_sensitivity_base, 3), "\n")
cat("Specificity Media: ", round(mean_specificity_base, 3), "\n")
cat("AUC Media: ", round(mean_auc_base, 3), "\n")

# --- Calcolo e stampa della Confusion Matrix e ROC per il Modello Base Finale ---
cat("\n--- Calcolo finale per il Modello Base ---\n")
# Addestramento del modello base sull'intero dataset per le metriche finali
final_model_base <- glm(data = injury_df, Injury_Occurred ~ Impact_Force_Newtons + Cumulative_Fatigue_Index, family = binomial)

# 1. Probabilità predette
prob_base_final <- predict(final_model_base, type = "response")

# 2. Classe predetta (con soglia 0.5)
pred_class_base_final <- ifelse(prob_base_final > 0.5, 1, 0)

# 3. Confusion matrix
cat("\nConfusion Matrix (Modello Base Finale):\n")
conf_matrix_base_final <- table(Predicted = pred_class_base_final, Actual = injury_df$Injury_Occurred)
print(conf_matrix_base_final)

# 4. ROC curve e AUC
roc_obj_base_final <- roc(injury_df$Injury_Occurred, prob_base_final)
auc_val_base_final <- auc(roc_obj_base_final)

cat("AUC (Modello Base Finale): ", round(auc_val_base_final, 3), "\n")

# Plot della ROC curve del modello base finale
plot(roc_obj_base_final, col = "blue", main = "ROC Curve (Modello Base Finale)")
abline(a = 0, b = 1, lty = 2, col = "gray") # linea random

#############################################################################
# LOGISTIC REGRESSION CON INTERAZIONE TRA LE 2 FEATURE
#############################################################################
all_metrics_interaction <- data.frame()

for (i in 1:k) {
  cat(paste0("\n--- Fold ", i, " (Modello con Interazione) ---\n"))
  test_indices <- folds[[i]]
  train_data <- injury_df[-test_indices, ]
  test_data <- injury_df[test_indices, ]
  
  # Training del modello
  pred2 <- glm(data = train_data, Injury_Occurred ~ Impact_Force_Newtons + Cumulative_Fatigue_Index +
                 Impact_Force_Newtons * Cumulative_Fatigue_Index, family = binomial)
  
  # Valutazione
  metrics <- evaluate_model(pred2, test_data, test_data$Injury_Occurred)
  all_metrics_interaction <- rbind(all_metrics_interaction, metrics)
  
  cat("Accuracy: ", round(metrics$accuracy, 3), "\n")
  cat("Sensitivity (Recall): ", round(metrics$sensitivity, 3), "\n")
  cat("Specificity: ", round(metrics$specificity, 3), "\n")
  cat("AUC: ", round(metrics$auc, 3), "\n")
}

# --- Media delle metriche per il modello con interazione ---
mean_accuracy_interaction <- mean(all_metrics_interaction$accuracy, na.rm = TRUE)
mean_sensitivity_interaction <- mean(all_metrics_interaction$sensitivity, na.rm = TRUE)
mean_specificity_interaction <- mean(all_metrics_interaction$specificity, na.rm = TRUE)
mean_auc_interaction <- mean(all_metrics_interaction$auc, na.rm = TRUE)

cat("\n--- Medie delle Metriche (Modello con Interazione) ---\n")
cat("Accuracy Media: ", round(mean_accuracy_interaction, 3), "\n")
cat("Sensitivity Media: ", round(mean_sensitivity_interaction, 3), "\n")
cat("Specificity Media: ", round(mean_specificity_interaction, 3), "\n")
cat("AUC Media: ", round(mean_auc_interaction, 3), "\n")

# --- Calcolo e stampa della Confusion Matrix e ROC per il Modello Base Finale ---
cat("\n--- Calcolo finale per il Modello Base ---\n")
# Addestramento del modello base sull'intero dataset per le metriche finali
final_model_interaction <- glm(data = injury_df, Injury_Occurred ~ Impact_Force_Newtons + Cumulative_Fatigue_Index, family = binomial)

# 1. Probabilità predette
prob_base_final <- predict(final_model_interaction, type = "response")

# 2. Classe predetta (con soglia 0.5)
pred_class_base_final <- ifelse(prob_base_final > 0.5, 1, 0)

# 3. Confusion matrix
cat("\nConfusion Matrix (Modello Base Finale):\n")
conf_matrix_base_final <- table(Predicted = pred_class_base_final, Actual = injury_df$Injury_Occurred)
print(conf_matrix_base_final)

# 4. ROC curve e AUC
roc_obj_base_final <- roc(injury_df$Injury_Occurred, prob_base_final)
auc_val_base_final <- auc(roc_obj_base_final)

cat("AUC (Modello Base Finale): ", round(auc_val_base_final, 3), "\n")

# Plot della ROC curve del modello base finale
plot(roc_obj_base_final, col = "blue", main = "ROC Curve (Modello con Interazione)")
abline(a = 0, b = 1, lty = 2, col = "gray") # linea random


#############################################################################
# DECISION TREE CON 2 FEATURE
#############################################################################
all_metrics_dt_base <- data.frame()

for (i in 1:k) {
  cat(paste0("\n--- Fold ", i, " (Decision Tree Base) ---\n"))
  test_indices <- folds[[i]]
  train_data <- injury_df[-test_indices, ]
  test_data <- injury_df[test_indices, ]
  
  # Training del modello Decision Tree
  # Usiamo le stesse impostazioni di controllo rpart.control
  tree_model_fold <- rpart(Injury_Occurred ~ Impact_Force_Newtons + Cumulative_Fatigue_Index,
                           data = train_data, method = "class", control = rpart.control(minsplit = 30, cp = 0.0001))
  
  # Valutazione
  metrics <- evaluate_dt_model(tree_model_fold, test_data, test_data$Injury_Occurred)
  
  # Aggiungiamo solo le metriche numeriche al dataframe
  all_metrics_dt_base <- rbind(all_metrics_dt_base,
                               data.frame(accuracy = metrics$accuracy,
                                          precision = metrics$precision,
                                          recall = metrics$recall,
                                          specificity = metrics$specificity,
                                          auc = metrics$auc))
  
  cat("Accuracy: ", round(metrics$accuracy, 3), "\n")
  cat("Precision: ", round(metrics$precision, 3), "\n")
  cat("Recall (Sensitivity): ", round(metrics$recall, 3), "\n")
  cat("Specificity: ", round(metrics$specificity, 3), "\n")
  cat("AUC: ", round(metrics$auc, 3), "\n")
}

# --- Media delle metriche per il modello Decision Tree Base ---
mean_accuracy_dt_base <- mean(all_metrics_dt_base$accuracy, na.rm = TRUE)
mean_precision_dt_base <- mean(all_metrics_dt_base$precision, na.rm = TRUE)
mean_recall_dt_base <- mean(all_metrics_dt_base$recall, na.rm = TRUE)
mean_specificity_dt_base <- mean(all_metrics_dt_base$specificity, na.rm = TRUE)
mean_auc_dt_base <- mean(all_metrics_dt_base$auc, na.rm = TRUE)

cat("\n--- Medie delle Metriche (Decision Tree Base) ---\n")
cat("Accuracy Media: ", round(mean_accuracy_dt_base, 3), "\n")
cat("Precision Media: ", round(mean_precision_dt_base, 3), "\n")
cat("Recall Media: ", round(mean_recall_dt_base, 3), "\n")
cat("Specificity Media: ", round(mean_specificity_dt_base, 3), "\n")
cat("AUC Media: ", round(mean_auc_dt_base, 3), "\n")

