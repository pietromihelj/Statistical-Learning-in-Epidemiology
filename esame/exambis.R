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
library(effsize)

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

# CALCOLO DEL SAMPLE SIZE PER TEST X^2 (ggbarstats)
alpha_level = 0.05 # Livello di significatività
power_level = 0.80 # Potenza desiderata (80%)

# --- Sport_Type vs. Injury_Occurred ---
cat("\n--- Calcolo Sample Size per Test Chi-quadro (Sport_Type vs. Injury_Occurred) ---\n")

# 1. Crea la tabella di contingenza
table_sport_injury <- table(injury_df$Sport_Type, injury_df$Injury_Occurred)
print("Tabella di Contingenza (Sport_Type vs. Injury_Occurred):")
print(table_sport_injury)

# 2. Calcola la dimensione dell'effetto 'w' da questa tabella
# ES.w2() calcola l'effetto size f^2, che per il chi-quadro si usa per derivare w
# w = sqrt(f^2)
# Puoi anche calcolare w direttamente usando la formula o con altre funzioni se disponibili.
# Per semplicità, pwr.chisq può prendere w direttamente.
# Una stima di w dai dati attuali è la radice quadrata del chi-quadro normalizzato.
chi_sq_test_sport <- chisq.test(table_sport_injury)
cohens_w_sport <- sqrt(chi_sq_test_sport$statistic / sum(table_sport_injury))
cat("Cohen's w stimato per Sport_Type vs. Injury_Occurred:", round(cohens_w_sport, 3), "\n")


# 3. Calcola i gradi di libertà (df)
df_sport_type <- (nrow(table_sport_injury) - 1) * (ncol(table_sport_injury) - 1)
cat("Gradi di Libertà (df):", df_sport_type, "\n")

# 4. Calcola il sample size
sample_size_chisq_sport <- pwr.chisq.test(w = cohens_w_sport, df = df_sport_type,
                                     sig.level = alpha_level, power = power_level)
print(sample_size_chisq_sport)
cat("Dimensione totale del campione richiesta:", ceiling(sample_size_chisq_sport$N), "\n")


# --- Activity_Type vs. Injury_Occurred ---
cat("\n--- Calcolo Sample Size per Test Chi-quadro (Activity_Type vs. Injury_Occurred) ---\n")

# 1. Crea la tabella di contingenza
table_activity_injury <- table(injury_df$Activity_Type, injury_df$Injury_Occurred)
print("Tabella di Contingenza (Activity_Type vs. Injury_Occurred):")
print(table_activity_injury)

# 2. Calcola la dimensione dell'effetto 'w' da questa tabella
chi_sq_test_activity <- chisq.test(table_activity_injury)
cohens_w_activity <- sqrt(chi_sq_test_activity$statistic / sum(table_activity_injury))
cat("Cohen's w stimato per Activity_Type vs. Injury_Occurred:", round(cohens_w_activity, 3), "\n")

# 3. Calcola i gradi di libertà (df)
df_activity_type <- (nrow(table_activity_injury) - 1) * (ncol(table_activity_injury) - 1)
cat("Gradi di Libertà (df):", df_activity_type, "\n")

# 4. Calcola il sample size
sample_size_chisq_activity <- pwr.chisq.test(w = cohens_w_activity, df = df_activity_type,
                                        sig.level = alpha_level, power = power_level)
print(sample_size_chisq_activity)
cat("Dimensione totale del campione richiesta:", ceiling(sample_size_chisq_activity$N), "\n")


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

# SAMPLE SIZE PER TEST non parametrico di MANN-WHITNEY U (ggbetween stats type = np)
# Per Impact_Force_Newtons vs. Injury_Occurred
cat("\n--- Calcolo Sample Size per Test Mann-Whitney U (Impact_Force_Newtons vs. Injury_Occurred) ---\n")
# Calcola Cohen's d dai tuoi dati attuali
# Nota: cohen.d() calcola d per due gruppi indipendenti
# Assicurati che Injury_Occurred sia numerica (0/1) per la funzione se lo richiede,
# oppure che la funzione accetti un fattore a 2 livelli. injury_df$Injury_Occurred è un fattore.
cohens_d_impact <- cohen.d(injury_df$Impact_Force_Newtons, injury_df$Injury_Occurred)$estimate
cat("Cohen's d stimato per Impact_Force_Newtons:", round(cohens_d_impact, 3), "\n")

# Usa il 'd' stimato per il calcolo del sample size
sample_size_mannwhitney_impact <- pwr.t.test(d = abs(cohens_d_impact), sig.level = alpha_level, power = power_level,
                                             type = "two.sample", alternative = "two.sided")
print(sample_size_mannwhitney_impact)
cat("Dimensione del campione richiesta per gruppo:", ceiling(sample_size_mannwhitney_impact$n), "\n")
cat("Dimensione totale del campione richiesta:", ceiling(sample_size_mannwhitney_impact$n) * 2, "\n")


# Per Cumulative_Fatigue_Index vs. Injury_Occurred
cat("\n--- Calcolo Sample Size per Test Mann-Whitney U (Cumulative_Fatigue_Index vs. Injury_Occurred) ---\n")
cohens_d_fatigue <- cohen.d(injury_df$Cumulative_Fatigue_Index, injury_df$Injury_Occurred)$estimate
cat("Cohen's d stimato per Cumulative_Fatigue_Index:", round(cohens_d_fatigue, 3), "\n")

sample_size_mannwhitney_fatigue <- pwr.t.test(d = abs(cohens_d_fatigue), sig.level = alpha_level, power = power_level,
                                              type = "two.sample", alternative = "two.sided")
print(sample_size_mannwhitney_fatigue)
cat("Dimensione del campione richiesta per gruppo:", ceiling(sample_size_mannwhitney_fatigue$n), "\n")
cat("Dimensione totale del campione richiesta:", ceiling(sample_size_mannwhitney_fatigue$n) * 2, "\n")


# Per Injury_Risk_Score vs. Injury_Occurred
cat("\n--- Calcolo Sample Size per Test Mann-Whitney U (Injury_Risk_Score vs. Injury_Occurred) ---\n")
cohens_d_risk <- cohen.d(injury_df$Injury_Risk_Score, injury_df$Injury_Occurred)$estimate
cat("Cohen's d stimato per Injury_Risk_Score:", round(cohens_d_risk, 3), "\n")

sample_size_mannwhitney_risk <- pwr.t.test(d = abs(cohens_d_risk), sig.level = alpha_level, power = power_level,
                                           type = "two.sample", alternative = "two.sided")
print(sample_size_mannwhitney_risk)
cat("Dimensione del campione richiesta per gruppo:", ceiling(sample_size_mannwhitney_risk$n), "\n")
cat("Dimensione totale del campione richiesta:", ceiling(sample_size_mannwhitney_risk$n) * 2, "\n")




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
tree_model <- rpart(Injury_Occurred ~ Impact_Force_Newtons + Cumulative_Fatigue_Index,
                    data = injury_df, method = "class", control = rpart.control(minsplit = 20, cp = 0.001))

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
)[, "1"]  # estrai la colonna di probabilità per la classe "1"

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
)[, "1"]

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



