library(here)
library(yarrr)
### RCT ###
dsdat <- read.csv(here("Exercises block 2-20250415/datasets", "dse-rli-trial-data-archive.csv"))
dsdat <- dsdat[dsdat$included == 1]
dsdat$group <- as.factor(dsdat$group)
levels(dsdat$group) <- c("Intervention","Control")
pirateplot(letter_sound_t2~group,theme=1,cex.lab=1.5,data=dsdat, point.o=1, xlab="", ylab="Post-intervention score")
tab <- psych::describeBy(dsdat$letter_sound_t2, group=dsdat$group, mat=TRUE, digits=2)
tab <- tab[,c(2,4:6)]
colnames(tab) <- c('Group','N','Mean','SD')
tab[1:2,1] <- c("Intervention","Control")
tab <- flextable::flextable(tab)
tab

## T-TEST ##
myt1 <- t.test(dsdat$letter_sound_t2~dsdat$group,var.equal=T,alternative='greater')
tab_T <- c(round(myt1$statistic,3), round())