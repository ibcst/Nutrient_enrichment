

library(funspace)
library(dplyr)
#set the directory
setwd("C:/Users/55619/Desktop/D_temp/Cap_Div/TPD")
traitsdf <- read.table("traits_bruto3.csv",sep=";",header=T)
traitsdf = na.omit(traitsdf)
####################################################
####################################################
####################################################
####################################################
###########POR TRATAMENTO!!!!!!!!!!#################
#CONTROLE
control_matrix = traitsdf[traitsdf$TRAT=='C',] #select control
control_matrix
library(sjmisc)
#put the data in log and scale it
controlmatrix2 = log (control_matrix[,4:ncol(control_matrix)])
controlmatrix2 = scale (controlmatrix2)
controlmatrix2

traitsc = replace_columns(control_matrix, controlmatrix2, add.unique = TRUE)
traitsc

treat = as.factor (control_matrix$PARC)
summary(treat)
#transform it into a matrix
traits.matrix = as.matrix(traitsc[,4:ncol(traitsc)])
traits.matrix = na.omit(traits.matrix)
#discover how many dimensions you should use
funspaceDim(traits.matrix) #two dimensions retained
#run PCA
pca.trait <- princomp(traits.matrix, cor = TRUE)
#build the functional space with the 2 dimensions
trait_space_C <- funspace(x = pca.trait ,n_divisions = 700, group.vec = treat)
summary(trait_space_C)

#1C
plot(trait_space_C, type="groups", which.group= "1C",threshold = 0.95,
     quant.plot=T, quant.col = "black",arrows=T,arrows.length = 0.7, colors=c("#c2e699", "#006837"),
     pnt=T,pnt.pch = 21, pnt.col = "black", xlim= c(-7, 7), ylim=c(-8,8))
#5C
plot(trait_space_C, type="groups", which.group= "5C",
     quant.plot=T, quant.col = "black",arrows=T,arrows.length = 0.7, colors=c("#c2e699", "#006837"),
     pnt=T, pnt.pch = 21,pnt.col = "black", xlim= c(-7, 7), ylim=c(-8,8))
#11C
plot(trait_space_C, type="groups", which.group= "11C",
     quant.plot=T, quant.col = "black",arrows=T,arrows.length = 0.7, colors=c("#c2e699", "#006837"),
     pnt=T, pnt.pch = 21,pnt.col = "black", xlim= c(-7, 7), ylim=c(-8,8))
#21C
plot(trait_space_C, type="groups", which.group= "21C",
     quant.plot=T, quant.col = "black",arrows=T,arrows.length = 0.7, colors=c("#c2e699", "#006837"),
     pnt=T, pnt.pch = 21,pnt.col = "black", xlim= c(-7, 7), ylim=c(-8,8))
#####################################################################################################
#repeat process for each treatment and threshold 



