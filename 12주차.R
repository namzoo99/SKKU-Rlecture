install.packages("survival")
install.packages("survminer")
library(survival) 
library(survminer)

# Drawing survival plot with Kaplan-Meier estimation method
fit <- survfit(Surv(time, status) ~ sex, data = lung) 

ggsurvplot(fit)


## Surv() creates a survival object
Surv(lung$time, lung$status)

## survfit() computes a survival curve for censored data
survfit(Surv(time, status) ~ sex, data = lung)


# Exercise
## 1. Draw survival plot of non-small cell lung cancer by smoking status
survival <- read.csv("/yourpath/survival.csv")


## 2. Draw survival plot of non-small cell lung cancer by subtype
## lecture_survival.csv, data_clinical_sample.txt
## Change color of lines to "#EFBC9B", "#EE92C2", "#A8B4A5", "#725D68"
subtype <- read.csv("/yourpath/data_clinical_sample.txt", sep = '\t')


# https://rpkgs.datanovia.com/survminer/reference/ggsurvplot.html
# https://www.rdocumentation.org/packages/survminer/versions/0.1.1/topics/ggsurvplot