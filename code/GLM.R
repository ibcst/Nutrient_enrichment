#Isabela B. Castro (ibc10@case.edu)
#GLM for soil analysis 

library(dplyr)
library(plyr)
library(easystats)
#load data
setwd("C:/Users/55619/Desktop/D_temp/Cap_Div/solo")
datf <- read.csv("soila.csv", sep= ";", header=T) 
datf
str(datf)



datf$PARC <- factor(datf$parcela)
datf$TRAT <- factor(datf$tratamento)



modelo1 <- glm(NH4 ~ tratamento , family=gaussian, data=datf)
modelo1
summary(modelo1)

#easy stats 
library(see)
report(modelo1)
check_model(modelo1)
model_performance(modelo1)
r2(modelo1)
model_parameters(modelo1)


contrasts <- estimate_contrasts(modelo1)
contrasts









