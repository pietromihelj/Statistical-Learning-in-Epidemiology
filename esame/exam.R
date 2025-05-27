library(DataExplorer)     # IDA
library(dplyr)
library(ggstatsplot)
library(dlookr)
library(flextable)
library(ggpubr)
library(SmartEDA)
library(performance)
library(tidyr)

injury_df <- read.csv('sports_injury_detection_dataset.csv', header = T, sep = ',')
summary(injury_df)
injury_df$Sport_Type <- as.factor(injury_df$Sport_Type)
injury_df$Activity_Type <- as.factor(injury_df$Activity_Type)
injury_df$Injury_Occurred <- as.factor(injury_df$Injury_Occurred)
levels(injury_df$Injury_Occurred) <- c('No','Yes')


## How roes and cols are ##
introduce(injury_df) %>% t()

## CATEGORICAL VARS FIRST LOOK ##
plot_bar(injury_df)
plot_bar(injury_df, by = "Injury_Occurred") #grouped by target
ggbarstats(data  = injury_df, x = Injury_Occurred, y = Sport_Type, label = "both") #look at category relation with target
ggbarstats(data  = injury_df, x = Injury_Occurred, y = Activity_Type, label = "both")

## NUMERICAL VAR SIRST LOOK ##
ft1 <- dlookr::describe(injury_df) %>% flextable()
ft1 <- colformat_double(x = ft1,digits = 2)
ft1

ft2 <- injury_df %>% 
  group_by(Injury_Occurred) %>% 
  univar_numeric()
knitr::kable(ft2, digits=2) #summary statistics

plot_density(injury_df) #density visualization for normality check

plot_qq(injury_df, by = "Injury_Occurred") #normality check

normality(injury_df) %>% mutate_if(is.numeric, ~round(., 3)) %>% flextable() #another normality check    

# normalization with log to check normality
# Elenco delle variabili da plottare (sostituisci con le tue)
vars_to_plot <- c("Blood_Oxygen_Level_Percent", "Cumulative_Fatigue_Index", "Duration_Minutes", "Heart_Rate_BPM", 
                  "Impact_Force_Newtons", "Injury_Risk_Score", "Respiratory_Rate_BPM", "Skin_Temperature_C")


injury_df_log <- injury_df %>%
  filter(across(all_of(vars_to_plot), ~ . > 0)) %>%
  mutate(across(all_of(vars_to_plot), log, .names = "log_{.col}")) 

# estrai le nuove variabili log e aggiungi Injury_Occurred
log_vars <- paste0("log_", vars_to_plot)
injury_df_log_subset <- injury_df_log %>%
  select(all_of(log_vars), Injury_Occurred)

# plot QQ delle variabili log
plot_qq(injury_df_log_subset, by = "Injury_Occurred") # non c'è nessun miglioramento neanche 
# con trasformazione logaritmica



  
# Applichiamo log, poi creiamo i qqplot per ogni variabile e gruppo
injury_df %>%
  # Filtro per valori positivi (log richiede valori > 0)
  filter(across(all_of(vars_to_plot), ~ . > 0)) %>%
  # Raggruppamento per Injury_Occurred
  group_by(Injury_Occurred) %>%
  # Applichiamo il log e passiamo a formato lungo per ggplot
  mutate(across(all_of(vars_to_plot), log, .names = "log_{.col}")) %>%
  pivot_longer(cols = starts_with("log_"), names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(sample = Value)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  facet_grid(Variable ~ Injury_Occurred, scales = "free") +
  labs(title = "QQ-plot delle variabili log-trasformate",
       subtitle = "Raggruppate per Injury_Occurred",
       x = "Quantili teorici", y = "Quantili campione") +
  theme_minimal()




ggqqplot(injury_df, 'Heart_Rate_BPM', facet.by = 'Injury_Occurred')
ggqqplot(injury_df, 'Respiratory_Rate_BPM', facet.by = 'Injury_Occurred')
ggqqplot(injury_df, 'Skin_Temperature_C', facet.by = 'Injury_Occurred')
ggqqplot(injury_df, 'Blood_Oxygen_Level_Percent', facet.by = 'Injury_Occurred')
ggqqplot(injury_df, 'Impact_Force_Newtons', facet.by = 'Injury_Occurred')
ggqqplot(injury_df, 'Cumulative_Fatigue_Index', facet.by = 'Injury_Occurred')
ggqqplot(injury_df, 'Duration_Minutes', facet.by = 'Injury_Occurred')
ggqqplot(injury_df, 'Injury_Risk_Score', facet.by = 'Injury_Occurred')

#no normality except for injury_risk-> i'm gonna use non parametrical approaches 


## FIRST BIVARIATE ANALISY ##
plot_boxplot(injury_df, by = "Injury_Occurred") #boxplot controll on target
ExpNumViz(injury_df, target = "Injury_Occurred", Page = c(2,4))


#looking in to difference of groups defined by target in continuous var
ggbetweenstats(data = injury_df, x = Injury_Occurred, y = Impact_Force_Newtons, type = "np")
ggbetweenstats(data = injury_df, x = Injury_Occurred, y = Cumulative_Fatigue_Index, type = "np")
ggbetweenstats(data = injury_df, x = Injury_Occurred, y = Injury_Risk_Score, type = "np")


#correlation between predictive vars
out <- check_outliers(injury_df$Cumulative_Fatigue_Index, method = c("zscore", "iqr", "cook", "mahalanobis"))
cc <- correlate(injury_df, method="spearman")
plot(cc) #corr matrix
ft8 <- correlation::correlation(data = injury_df, method= "pearson")
ft8 #table

#correlation between numerical 
var = as.numeric(injury_df$Injury_Occurred)-1
cor(injury_df$Skin_Temperature_C, var, method = 'pearson')
cor(injury_df$Impact_Force_Newtons, var, method = 'pearson')
cor(injury_df$Cumulative_Fatigue_Index, var, method = 'pearson')

## BUILDING A PREDICTOR ## 
all <- glm(data = injury_df, Injury_Occurred ~., family=binomial)
modello_best <- step(all, direction = "both", trace = 1)


# LOGISTIC REGRESSION BASE #
levels(injury_df$Injury_Occurred) <- c(0,1)
pred1 = glm(data = injury_df, Injury_Occurred ~ Impact_Force_Newtons + Cumulative_Fatigue_Index, family = binomial)
summary(pred1)

pred1 = glm(data = injury_df, Injury_Occurred ~ ., family = binomial)
summary(pred2)

pred1 <- glm(Injury_Occurred ~ ., data = injury_df[,-c(1,3,12)], family = binomial(), control = list(maxit = 50))
summary(pred1)

# 1. Probabilità predette
prob <- predict(pred1, type = "response")

# 2. Classe predetta (con soglia 0.5)
pred_class <- ifelse(prob > 0.5, 1, 0)

# 3. Confusion matrix
conf_matrix <- table(Predicted = pred_class, Actual = injury_df$Injury_Occurred)
print(conf_matrix)

# 4. Metriche di classificazione
TP <- conf_matrix[2, 2]
TN <- conf_matrix[1, 1]
FP <- conf_matrix[2, 1]
FN <- conf_matrix[1, 2]

accuracy <- (TP + TN) / sum(conf_matrix)
sensitivity <- TP / (TP + FN)   # recall / true positive rate
specificity <- TN / (TN + FP)   # true negative rate

cat("Accuracy: ", round(accuracy, 3), "\n")
cat("Sensitivity (Recall): ", round(sensitivity, 3), "\n")
cat("Specificity: ", round(specificity, 3), "\n")

# 5. ROC curve e AUC
library(pROC)
roc_obj <- roc(injury_df$Injury_Occurred, prob)
auc_val <- auc(roc_obj)
plot(roc_obj, col = "blue", main = "ROC Curve")
abline(a = 0, b = 1, lty = 2, col = "gray")  # linea random

cat("AUC: ", round(auc_val, 3), "\n")

# il modello ha delle prestazioni appena sufficienti, proviamo a fare qualcosa di più complesso


# LOGISTIC REGRESSION CON ALTRE RELAZIONI TRA I PREDITTORI
pred2 = glm(data = injury_df, Injury_Occurred ~ Impact_Force_Newtons + Cumulative_Fatigue_Index + 
              Impact_Force_Newtons*Cumulative_Fatigue_Index, family = binomial)
summary(pred2)
# 1. Probabilità predette
prob <- predict(pred2, type = "response")

# 2. Classe predetta (con soglia 0.5)
pred_class <- ifelse(prob > 0.5, 1, 0)

# 3. Confusion matrix
conf_matrix <- table(Predicted = pred_class, Actual = injury_df$Injury_Occurred)
print(conf_matrix)

# 4. Metriche di classificazione
TP <- conf_matrix[2, 2]
TN <- conf_matrix[1, 1]
FP <- conf_matrix[2, 1]
FN <- conf_matrix[1, 2]

accuracy <- (TP + TN) / sum(conf_matrix)
sensitivity <- TP / (TP + FN)   # recall / true positive rate
specificity <- TN / (TN + FP)   # true negative rate

cat("Accuracy: ", round(accuracy, 3), "\n")
cat("Sensitivity (Recall): ", round(sensitivity, 3), "\n")
cat("Specificity: ", round(specificity, 3), "\n")

# 5. ROC curve e AUC
library(pROC)
roc_obj <- roc(injury_df$Injury_Occurred, prob)
auc_val <- auc(roc_obj)
plot(roc_obj, col = "blue", main = "ROC Curve")
abline(a = 0, b = 1, lty = 2, col = "gray")  # linea random

cat("AUC: ", round(auc_val, 3), "\n")

# nessun tipo di miglioramento anche aggiungendo l'inteazione tra i 2 predittori

# DECISION TREE #
# 1. Carica pacchetti (installa se serve)
if (!require(rpart)) install.packages("rpart")
if (!require(rpart.plot)) install.packages("rpart.plot")
if (!require(pROC)) install.packages("pROC")

library(rpart)
library(rpart.plot)
library(pROC)

# 2. Modello ad albero decisionale
injury_shuffle =  injury_df[sample(nrow(injury_df)), ]
injury_tr = injury_shuffle[1:800,]
injury_te = injury_shuffle[801:1000,]
tree_model <- rpart(Injury_Occurred ~ Impact_Force_Newtons + Cumulative_Fatigue_Index,
                    data = injury_tr, method = "class", control = rpart.control(minsplit = 20, cp = 0.001))

# cp va lasciato così altrimenti fa cagare, provare a giocare su minsplit nel range 15-20 (più è basso e più l'albero è complesso)


# 3. Visualizzazione dell’albero
rpart.plot(tree_model, type = 3, extra = 102, fallen.leaves = TRUE)

# 4. Predizioni (classe e probabilità)
pred_class <- predict(tree_model, type = "class")  # Classe 0/1
pred_prob <- predict(tree_model)[, 2]              # Probabilità classe "1"

# 5. Confusion matrix
conf_matrix <- table(Predicted = pred_class, Actual = injury_df$Injury_Occurred)
print(conf_matrix)

# 6. Metriche di valutazione
TP <- conf_matrix["1", "1"]
TN <- conf_matrix["0", "0"]
FP <- conf_matrix["1", "0"]
FN <- conf_matrix["0", "1"]

accuracy <- (TP + TN) / sum(conf_matrix)
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)

cat("Accuracy: ", round(accuracy, 3), "\n")
cat("Sensitivity: ", round(sensitivity, 3), "\n")
cat("Specificity: ", round(specificity, 3), "\n")

# 7. ROC Curve e AUC
roc_obj <- roc(injury_df$Injury_Occurred, pred_prob)
auc_val <- auc(roc_obj)

plot(roc_obj, col = "darkred", main = "ROC Curve - Decision Tree")
abline(a = 0, b = 1, lty = 2, col = "gray")
cat("AUC: ", round(auc_val, 3), "\n")

# abbiamo cercato un DT con un trade off tra complessità e prestazioni, va decisamente meglio della logistic regression e siamo riusciti a non andare in overfitting

# PROVO A CREARE LA HEATMAP
library(ggplot2)
library(rpart)

# 1. Costruisci una griglia più fitta
x_seq <- seq(
  from = min(injury_df$Impact_Force_Newtons, na.rm = TRUE),
  to   = max(injury_df$Impact_Force_Newtons, na.rm = TRUE),
  length.out = 300    # risoluzione aumentata
)
y_seq <- seq(
  from = min(injury_df$Cumulative_Fatigue_Index, na.rm = TRUE),
  to   = max(injury_df$Cumulative_Fatigue_Index, na.rm = TRUE),
  length.out = 300
)
grid <- expand.grid(
  Impact_Force_Newtons     = x_seq,
  Cumulative_Fatigue_Index = y_seq
)

# 2. Predict con probabilità
grid$pred_prob <- predict(
  tree_model,
  newdata = grid,
  type   = "prob"
)[,'Yes'] # estrai la colonna di probabilità per la classe "1"

# 3. Plot
ggplot(grid, aes(x = Impact_Force_Newtons, y = Cumulative_Fatigue_Index)) +
  # a) background colorato in base alla probabilità
  geom_tile(aes(fill = pred_prob)) +
  # b) contorni sulle iso-probabilità
  geom_contour(
    aes(z = pred_prob),
    breaks = seq(0, 1, by = 0.1),  # ogni 0.1 di probabilità
    color  = "white",
    size   = 0.3
  ) +
  # 4. scala di colori continua
  scale_fill_gradient(
    name    = "P(Injury)",
    low     = "lightgreen",
    high    = "salmon"
  ) +
  # 5. tema ed etichette
  labs(
    x = "Impact Force (Newtons)",
    y = "Cumulative Fatigue Index",
    title = "Mappa di probabilità da albero decisionale"
  ) +
  theme_minimal(base_size = 14)


# VARIANTE HEATMAP SOLO CON ZONE ROSSE CHIARO E SCURO
library(ggplot2)

# 1. Prepara la griglia (risoluzione alta)
x_seq <- seq(
  from = min(injury_df$Impact_Force_Newtons, na.rm = TRUE),
  to   = max(injury_df$Impact_Force_Newtons, na.rm = TRUE),
  length.out = 300
)
y_seq <- seq(
  from = min(injury_df$Cumulative_Fatigue_Index, na.rm = TRUE),
  to   = max(injury_df$Cumulative_Fatigue_Index, na.rm = TRUE),
  length.out = 300
)
grid <- expand.grid(
  Impact_Force_Newtons     = x_seq,
  Cumulative_Fatigue_Index = y_seq
)

# 2. Ottieni le probabilità
grid$pred_prob <- predict(
  tree_model,
  newdata = grid,
  type   = "prob"
)[, "Yes"]

# 3. Definisci il “fill” in base alle soglie
grid$fill <- with(grid, 
                  ifelse(pred_prob > 0.80, "darkred",
                         ifelse(pred_prob > 0.60, "lightcoral", NA))
)

# 4. Plottaggio
ggplot(grid, aes(x = Impact_Force_Newtons, y = Cumulative_Fatigue_Index)) +
  # a) tiles con solo le aree colorate (NA → bianco)
  geom_tile(aes(fill = fill), na.rm = TRUE) +
  # b) usa i colori così come sono definiti
  scale_fill_identity(
    name = "Probabilità\nInfortunio",
    guide = "legend",
    labels = c("> 0.80" = "darkred", 
               "0.60–0.80" = "lightcoral"),
    breaks = c("lightcoral", "darkred")
  ) +
  labs(
    x = "Impact Force (Newtons)",
    y = "Cumulative Fatigue Index",
    title = "Aree con P(infortunio) > 0.60"
  ) +
  theme_minimal(base_size = 14)


# VARIANTE HEATMAP SOLO CON ZONE ROSSO SCURO
library(ggplot2)

# Ricorda di aver già creato:
#   grid$pred_prob = predict(tree_model, newdata = grid, type = "prob")[, "1"]

ggplot(grid, aes(
  x = Impact_Force_Newtons,
  y = Cumulative_Fatigue_Index,
  z = pred_prob
)) +
  # 1) Riempi soltanto l'area sopra 0.75
  stat_contour_filled(
    breaks     = c(0.80, 1),   # un unico intervallo (0.75,1]
    show.legend = FALSE
  ) +
  # 2) Applica solo il colore darkred al livello estratto
  scale_fill_manual(values = c("darkred")) +
  # 3) Titoli e tema
  labs(
    x     = "Impact Force (Newtons)",
    y     = "Cumulative Fatigue Index",
    title = "Zona a P(infortunio) > 0.75"
  ) +
  theme_minimal(base_size = 14)










library(dplyr)

# Assumiamo che `fill` sia attualmente colori.
# Creiamo una nuova colonna con i livelli testuali (range) corrispondenti:
grid <- grid %>%
  mutate(prob_range = case_when(
    fill == "green"       ~ "0.0-0.25",
    fill == "lightgreen"  ~ "0.25–0.50",
    fill == "red"         ~ "0.50-0.75",
    fill == "darkred"     ~ "> 0.75",
    TRUE                  ~ NA_character_
  ))

# Ora usiamo questa nuova colonna per fill:
ggplot(grid, aes(x = Impact_Force_Newtons, y = Cumulative_Fatigue_Index)) +
  geom_tile(aes(fill = prob_range), na.rm = TRUE) +
  scale_fill_manual(
    name = "Probabilità\nInfortunio",
    values = c(
      "0.0-0.25" = "green",
      "0.25–0.50" = "lightgreen",
      "0.50-0.75" = "red",
      "> 0.75" = "darkred"
    )
  ) +
  labs(
    x = "Impact Force (Newtons)",
    y = "Cumulative Fatigue Index",
    title = "Aree con P(infortunio) > 0.60"
  ) +
  theme_minimal(base_size = 14)














