library(epiR)
library(Epi)
library(popEpi)
library(ggplot2)
library(scales)
library(lubridate)
library(tidyverse)
library(OptimalCutpoints)
library(pROC)
library(MASS)
library(here)
library(epitools)

### FIRST BLOCK ###

##PREVALENCE##
ncas <- 4; npop <- 200
tmp <- as.matrix(cbind(ncas, npop))
epi.conf(tmp, ctype = "prevalence", method = "exact", N = 1000, design = 1, conf.level = 0.95) * 100

#now we calculate empirical and teoretichal prevalence
16281/5534738
X <- 16281
N <- 5534738
mp <- glm(cbind(X, N-X)~1, family=binomial) #linear model fitted on a binomial for the teoretical prevalence
round(ci.pred(mp)*100,4)

#another example
data(pr)
prt <- stat.table(index=list(A,sex), contents=ratio(X,N,100), data=pr)[1,,]
round(prt[70:75,],1)
matplot(0:99+0.5, prt, type="l", lwd=2, col=c(4,2),lty=1 )
abline(v=c(76,81)+0.5)



##INCIDENCE##
#incidence rate
ncas <- 136; ntar <- 22050
tmp <- as.matrix(cbind(ncas, ntar))
epi.conf(tmp, ctype = "inc.rate", method = "exact", N = 1000, design = 1, conf.level = 0.95) * 1000

#mortality rate
2648/150000 #empirical rate
2648/300    #mortality rate with pearson year consideration 
#teoretical mortality rate
D <- 2648   
Y <- 300000
mm <- glm(D~1,offset=log(Y), family=poisson)
round(ci.pred(mm,newdata=data.frame(Y=1000)),3)

##exercise##
#1
26/481000 #incidence rate for 1999
113/(mean(c(496000,497000,496000,495000,491000))*5)#IR for 1993-1997
#2
D <- 3   
Y <- 481000
mm <- glm(D~1,offset=log(Y), family=poisson)
round(ci.pred(mm,newdata=data.frame(Y=481)),3)#DR for 1999

D<-22
Y<-mean(c(496000,497000,496000,495000,491000))*5
mm <- glm(D~1,offset=log(Y), family=poisson)
round(ci.pred(mm,newdata=data.frame(Y=Y)),3)#DR for 1993-1997

### SECOND BLOCK ###  
##standardized rates
M <- matrix( c(10,1157,46.0,0.9,22,1109,41.9,2.0,0.44
               ,76,809,32.0,9.4,68,786,29.7,8.6,1.09
               ,305,455,18.0,67,288,524,19.8,55,1.22
               ,201,102,4.0,196,354,229,8.6,155,1.27), nrow=4, byrow=T )
M <- data.frame(M)
names(M) <- c("mca","mpy","mp","mr",
              "fca","fpy","fp","fr","rr")
M

#crude ratio
rates <- with( M, c( sum(mca)/sum(mpy)*100,
                     sum(fca)/sum(fpy)*100 ) )
rates[3] <- rates[1]/rates[2]
names(rates) <- c("M rate","F rate","M/F RR")
round( rates, 2 )

#standard ratios
#direct
wm <- with( M, mpy/sum(mpy) ) #compute the category weigths
rates <-with( M, c( sum(mca/mpy*wm)*100, 
                    sum(fca/fpy*wm)*100 ) ) #compute the wheigted rates
rates[3] <- rates[1]/rates[2]
names(rates) <- c("M rate","F rate","M/F RR")
round( rates, 2 )

#undirect
WSP <- c(120,100,90,90,80,80,60,60,60,60,50,40,40,30,20,10,5,3,2) #world distribution
wt <- sum(WSP[1:7])
wt[2] <- sum(WSP[8:11])
wt[3] <- sum(WSP[12:15])
wt[4] <- sum(WSP[16:19])
wt <- wt/sum(wt)         #world weights for category
rates <- with(M, c(sum(mca/mpy*wt)*100,sum(fca/fpy*wt)*100))
rates[3] <- rates[1]/rates[2]
names(rates) <- c("M rate","F rate","M/F RR")
round( rates, 2 )

##exercise

std.rates <- read.delim(here("datasets","std-rates.txt"))
head(std.rates)
std <- std.rates
std
colon_rates <- with(std,c(sum(D[sex == 'M' & typ == 'Colon'])/2521177 * 100, sum(D[sex == 'F' & typ == 'Colon'])/2596061 * 100))
colon_rates[3] <- colon_rates[1]/colon_rates[2]

rectum_rates <- with(std,c(sum(D[sex == 'M' & typ == 'Rectum'])/2521177 * 100, sum(D[sex == 'F' & typ == 'Rectum'])/2596061 * 100))
rectum_rates[3] <- rectum_rates[1]/rectum_rates[2]

lung_rates <- with(std,c(sum(D[sex == 'M' & typ == 'Lung'])/2521177 * 100, sum(D[sex == 'F' & typ == 'Lung'])/2596061 * 100))
lung_rates[3] <- lung_rates[1]/lung_rates[2]

rates <- rbind(colon_rates, rectum_rates, lung_rates)
colnames(rates) <- c("M rate","F rate","M/F RR") 
round( rates, 2 )

WSP <- c(120,100,90,90,80,80,60,60,60,60,50,40,40,30,20,10,5,3,2) #world distribution
cwt <- sum(WSP[1:7])/sum(WSP)
rcwt <- sum(WSP[8:14])/sum(WSP)
lwt <- sum(WSP[14:length((WSP))])/sum(WSP)

sc_rates <- with(std,c(sum(D[sex == 'M' & typ == 'Colon']/2521177)*cwt*100, sum(D[sex == 'F' & typ == 'Colon']/2521177)*cwt*100))          
sc_rates[3] <- sc_rates[1]/sc_rates[2]

rc_rates <- with(std,c(sum(D[sex == 'M' & typ == 'Rectum']/2521177)*cwt*100, sum(D[sex == 'F' & typ == 'Rectum']/2521177)*cwt*100))          
rc_rates[3] <- rc_rates[1]/rc_rates[2]

lc_rates <- with(std,c(sum(D[sex == 'M' & typ == 'Lung']/2521177)*cwt*100, sum(D[sex == 'F' & typ == 'Lung']/2521177)*cwt*100))          
lc_rates[3] <- lc_rates[1]/lc_rates[2]

srates <- rbind(sc_rates, rc_rates, lc_rates)
colnames(srates) <- c("M rate","F rate","M/F RR") 
round( srates, 4 )

### THIRD BLOCK ###
##diagnostic test##

dat <- as.table(matrix(c(815,115,208,327), nrow = 2, byrow = TRUE))
colnames(dat) <- c("Dis+","Dis-")
rownames(dat) <- c("Test+","Test-")
dat
rval <- epi.tests(dat, conf.level = 0.95) #this calculate all with the confidence intervall
print(rval)     

##example##
load(here("datasets","CK.Rdata"))
CK$ck.100 <- ifelse(CK$ck>80, 'CK > 80', 'CK <= 80') #creating the categories
pp <- table(CK$ck.100,CK$disease) #creating the table
dat <- as.table(matrix(c(26,34,1,58), nrow = 2, byrow = TRUE)) #putting the numbers in the right order
colnames(dat) <- c("AMI+","AMI-")
rownames(dat) <- c("Test+","Test-")
rval <- epi.tests(dat, conf.level = 0.95)
print(rval)

##exercise##
dat <- as.table(matrix(c(800,3200,100,95900), nrow = 2, byrow = TRUE)) #putting the numbers in the right order
colnames(dat) <- c("AMI+","AMI-")
rownames(dat) <- c("Test+","Test-")
rval <- epi.tests(dat, conf.level = 0.95)
print(rval)


### FOURTH BLOCK ###
##rock curve
#example
data(aSAH)
with(aSAH, table(gender, outcome))
roc(aSAH$outcome, aSAH$s100b,percent=TRUE, plot=TRUE, ci=TRUE)
#adding the vonfidence intervall
rocobj <- plot.roc(aSAH$outcome, aSAH$s100b)
ci.sp.obj <- ci.sp(rocobj, sensitivities=seq(0, 1, .01), boot.n=100)
plot(rocobj)
plot(ci.sp.obj, type="shape", col="blue")
#comparing two markers
roc.test(response = aSAH$outcome, predictor1= aSAH$wfns, predictor2 = aSAH$s100)
roc.s100b <- roc(aSAH$outcome, aSAH$s100b)
roc.wfns <- roc(aSAH$outcome, aSAH$wfns)
plot(roc.s100b)
plot(roc.wfns, add=TRUE, col="red")
#finding the optimal cutpoints
perf.s100b <-optimal.cutpoints(X ="s100b", status ="outcome", tag.healthy ="Good",
                               methods = "Youden", data =aSAH, control = control.cutpoints(), 
                               ci.fit = TRUE, 
                               conf.level = 0.95, trace = FALSE)
perf.s100b

##excercise
diabetes <- read.csv(here("datasets","diabetes.csv"))
glimpse(diabetes)
diabetes$Outcome = factor(diabetes$Outcome, labels=c("No", "Yes"))
table(diabetes$Outcome)

roc(diabetes$Outcome, diabetes$Glucose, plot=TRUE, ci=TRUE)
perf.s100b <-optimal.cutpoints(X ="Glucose", status ="Outcome", tag.healthy ="No",
                               methods = "MinValueNPV", data =diabetes, control = control.cutpoints(), 
                               ci.fit = TRUE, 
                               conf.level = 0.95, trace = FALSE)
perf.s100b

diabetes$Glucose <- ifelse(diabetes$Glucose>113, 'G > 124', 'G <= 124')
tab <- table(diabetes$Glucose, diabetes$Outcome)
dat <- as.table(matrix(c(183,128,85,372), nrow = 2, byrow = TRUE)) #putting the numbers in the right order
colnames(dat) <- c("AMI+","AMI-")
rownames(dat) <- c("Test+","Test-")
rval <- epi.tests(dat, conf.level = 0.95) #this calculate all with the confidence intervall
print(rval) 

### FIFTH BLOCK ###
##Relative risk
data.type <- matrix(c(178,79,1141,1486), 2, 2)
dimnames(data.type) <- list("Behavior type" = c("Type A", "Type B"),
                            "Occurrence of CHD Event" = c("Yes", "No"))
epi.2by2(data.type, method = "cohort.count")
##Odds ratio
data.coffee <- matrix(c(347,20,555,88), 2, 2)
dimnames(data.coffee) <- list("Coffee Drinking (cups per day)" = c(">=1", "0"),
                              "Pancreatic cancer" = c("cases", "controls"))
oddsratio(data.coffee, method="wald")

##Attributable risk
bwt <- birthwt
head(bwt)
bwt$low    <- factor(bwt$low, levels = c(1,0))
bwt$smoke  <- factor(bwt$smoke, levels = c(1,0))
bwt$race   <- factor(bwt$race, levels = c(1,2,3))
low.tab <- table(bwt$smoke, bwt$low, dnn = c("Smoke", "Low BW"))
low.tab
epi.2by2(dat = low.tab, method = "cohort.count", conf.level = 0.95, 
         units = 100, outcome = "as.columns")

##Exercises
RNE = 1/100000
risks = c(1/50000,1/10000,1/1000)
RRs = risks/RNE
AREs = (RRs-1)/RRs
tab <- data.frame(risks, AREs)
tab

Pe = 5/100
ARs = (Pe*(RRs-1))/(Pe*(RRs-1)+1)
tab <- data.frame(risks,AREs,ARs)
tab

data <- matrix(c(85,77,1736,2400),2,2)
dimnames(data)<-list("Type of living" = c('Rented', 'Property'), 'Coronary ARtery Disesase' = c('Yes','No'))
epi.2by2(data, method = "cohort.count")
4.67/3.11
(85/1821)/(1-(85/1821))
(77/2477)/2400


diabetes <- read.csv(here("datasets","ugdp.csv"))
glimpse(diabetes)
tab <- table(diabetes$Treatment, diabetes$Status)
mat <- as.matrix(tab)    # Converti in matrice
mat[c(2,1), ] <- mat[c(1,2), ]  # Scambia righe 1 e 2
tab <- as.table(mat)     # Riconverti in tabella
dimnames(tab) <- list(rows = c('Tolbutamide', 'Placebo'), Cols = c('Death', 'Survivor'))
epi.2by2(tab, method = "cohort.count")



























