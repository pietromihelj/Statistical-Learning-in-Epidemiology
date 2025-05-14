library(DataExplorer)     # IDA
library(dplyr)
library(ggstatsplot)
library(dlookr)
library(flextable)
library(ggpubr)
library(SmartEDA)
library(performance)
library(tidyr)
library(neuralnet)

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
out <- check_outliers(injury_df$Cumulative_Fatigue_Index, method = c("zscore"))
plot(out)
cc <- correlate(injury_df[,1:11], method="spearman")
plot(cc) #corr matrix
ft8 <- correlation::correlation(data = injury_df[,1:11], method= "pearson")
ft8 #table

#correlation between numerical 
var = as.numeric(injury_df$Injury_Occurred)-1
cor(injury_df$Skin_Temperature_C, var, method = 'pearson')
cor(injury_df$Impact_Force_Newtons, var, method = 'pearson')
cor(injury_df$Cumulative_Fatigue_Index, var, method = 'pearson')

## BUILDING A PREDICTOR ## 
all <- glm(data = injury_df, Injury_Occurred ~., family=binomial)
modello_best <- stepAIC(all, direction = "both", trace = 1)





