###Isabela C. 2023
#ibc10@case.edu
#Hypervolume (Blonder et al 2018) construction using Rius et al (2023) as a base and code

############################################################

#Data treatment and PCA 

############################################################
#packages
library(dplyr)
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(hypervolume)
library(ggpubr)
library(dplyr)
#set the directory 
setwd("C:/Users/55619/Desktop/D_temp/Cap_Div/hypervolume")
traits <- read.table("traits_bruto.csv",sep=";",header=T) 

#clean the data 

traits <- na.omit(traits)

traits <- traits %>% select(-SP, -H) #get rid of the species column 
str(traits)

####################PCA
#perform PCA 

set.seed(111)

pca_all <- prcomp(traits[,-1:-2], #get rid of treatment column 
               center = TRUE,
               scale. = TRUE)
attributes(pca_all)
print(pca_all)

fviz_eig(pca_all, addlabels = TRUE) #check the contribution of each axes to variation


#plot the pca 

paleta3 = c("black","#BD3A17","#61C78B", "#19488F" ,"#C763B1")

fviz_pca_biplot(pca_all, geom="point", pointshape=19, col.var = "black", habillage = traits$TRAT,
                addlabel=T,
                palette = paleta3,
                alpha.var=2,
                pointsize = 2,
                ggtheme = theme_minimal())

#get score factors 

mydf=summary(pca_all)
pca_all$rotation
pca_all$center
pca_all$scale

new_variables<-pca_all$x
df_new_variables<-as.data.frame(new_variables)
colnames(df_new_variables)
final_table<-bind_cols(traits,df_new_variables)
write.csv(final_table, file = "C:/Users/55619/Desktop/D_temp/Cap_Div/hypervolume/pcanormalfactor.csv")

############################################################

#Hypervolume construction 

############################################################


#get table 
traits_we <- read.table("pcanormalfactor.csv",sep=",",header=T) 
traits_we


#perform hypervolume for each treatment 
#select data for control

HVD_C = traits_we[traits_we$TRAT=='C',] #select control
HVD_C <- HVD_C %>% select(PC1, PC2, PC3) #select pca axes (in this case 82% explained var)
str(HVD_C)
#hypervolume  

hv_c=hypervolume(HVD_C, method = "gaussian", name="control")
summary(hv_c)

#plot simple
plot(hv_c, color="black")


#select data for N treatment

HVD_N = traits_we[traits_we$TRAT=='N',] #select control
HVD_N <- HVD_N %>% select(PC1, PC2, PC3) #select pca axes (in this case 82% explained var)
str(HVD_N)
#hypervolume 

hv_n=hypervolume(HVD_N, method = "gaussian",name="nitrogen")
summary(hv_n)

#plot simple
plot(hv_n, color="#61C78B")


#select data for P treatment 

HVD_P = traits_we[traits_we$TRAT=='P',] #select control
HVD_P <- HVD_P %>% select(PC1, PC2, PC3) #select pca axes (in this case 82% explained var)
str(HVD_P)
#hypervolume  

hv_p=hypervolume(HVD_P, method = "gaussian", name="phosphorus")
summary(hv_p)

#plot simple
plot(hv_p, color="#C763B1")

#select data for NP treatment 

HVD_NP = traits_we[traits_we$TRAT=='NP',] #select control
HVD_NP <- HVD_NP %>% select(PC1, PC2, PC3) #select pca axes (in this case 82% explained var)
str(HVD_NP)
#hypervolume 

hv_np=hypervolume(HVD_NP, method = "gaussian", name="np")
summary(hv_np)

#plot simple
plot(hv_np, color="#19488F")

#select data for Ca treatment 

HVD_CA = traits_we[traits_we$TRAT=='Ca',] #select control
HVD_CA <- HVD_CA %>% select(PC1, PC2, PC3) #select pca axes (in this case 82% explained var)
str(HVD_CA)
#hypervolume 

hv_ca=hypervolume(HVD_CA, method = "gaussian", name="liming")
summary(hv_ca)

#plot simple
plot(hv_ca, color="#BD3A17")


#########################################
#plot them together 
set.seed(111)
#CONTROL VS N 
cl=c("black","#61C78B")
list=hypervolume_join(hv_c,hv_n)
plotc_n=plot.HypervolumeList(list,show.density=T,show.legend = F, show.random = T,
                             cex.axis = 1, show.data=T,show.centroid = T,show.contour = T,
                             contour.lwd=.5,contour.type ="kde", colors = cl,cex.data = 0,
                             cex.legend = 0,cex.centroid=2,contour.raster.resolution = 1,reshuffle = F,
                             contour.kde.level = 1e-02, verbose=T)
??plot.HypervolumeList 
#CONTROL VS P

cl=c("black","#C763B1")
list=hypervolume_join(hv_c,hv_p)
plotc_p=plot.HypervolumeList(list,show.density=T,show.legend = F, show.random = T,
                             cex.axis = 1, show.data=T,show.centroid = T,show.contour = T,
                             contour.lwd=.5,contour.type ="kde", colors = cl,cex.data = 0.3,
                             cex.legend = 0,cex.centroid=2,contour.raster.resolution = 20,reshuffle = F,
                             contour.kde.level = 1e-02, verbose=T)


#CONTROL VS NP

cl=c("black","#19488F")
list=hypervolume_join(hv_c,hv_np)
plotc_np=plot.HypervolumeList(list,show.density=T,show.legend = F, show.random = T,
                             cex.axis = 1, show.data=T,show.centroid = T,show.contour = T,
                             contour.lwd=.5,contour.type ="kde", colors = cl,cex.data = 0.3,
                             cex.legend = 0,cex.centroid=2,contour.raster.resolution = 20,reshuffle = F,
                             contour.kde.level = 1e-02, verbose=T)

#CONTROL VS CA

cl=c("black","#BD3A17")
list=hypervolume_join(hv_c,hv_ca)
plotc_ca=plot.HypervolumeList(list,show.density=T,show.legend = F, show.random = T,
                             cex.axis = 1, show.data=T,show.centroid = T,show.contour = T,
                             contour.lwd=.5,contour.type ="kde", colors = cl,cex.data = 0.3,
                             cex.legend = 0,cex.centroid=2,contour.raster.resolution = 20,reshuffle = F,
                             contour.kde.level = 1e-02, verbose=T)





ggarrange( plotc_n, plotc_p, plotc_np, plotc_ca, 
           ncol = 3, nrow = 1, widths = c(0.7, 0.7), heights = (4), align = "v",
           font.label = list(size = 5))

#statistics 

#control and nitrogen 
hv_c_n_set<- hypervolume_set(hv_c,hv_n, check.memory=FALSE)

c_n_stats=hypervolume_overlap_statistics(hv_c_n_set)
c_n_vol=get_volume(hv_c_n_set)

c_n_stats <- as.data.frame(c_n_stats)
c_n_vol <- as.data.frame(c_n_vol)

#control and p

hv_c_p_set<- hypervolume_set(hv_c,hv_p, check.memory=FALSE)

c_p_stats=hypervolume_overlap_statistics(hv_c_p_set)
c_p_vol=get_volume(hv_c_p_set)

c_p_stats <- as.data.frame(c_p_stats)
c_p_vol <- as.data.frame(c_p_vol)

#control and np

hv_c_np_set<- hypervolume_set(hv_c,hv_np, check.memory=FALSE)

c_np_stats=hypervolume_overlap_statistics(hv_c_np_set)
c_np_vol=get_volume(hv_c_np_set)

c_np_stats <- as.data.frame(c_np_stats)
c_np_vol <- as.data.frame(c_np_vol)

#control and ca

hv_c_ca_set<- hypervolume_set(hv_c,hv_ca, check.memory=FALSE)

c_ca_stats=hypervolume_overlap_statistics(hv_c_ca_set)
c_ca_vol=get_volume(hv_c_ca_set)

c_ca_stats <- as.data.frame(c_ca_stats)
c_ca_vol <- as.data.frame(c_ca_vol)


#get table 
stats_table<-bind_cols(c_n_stats,c_p_stats,c_np_stats,c_ca_stats)
write.csv(stats_table, file = "C:/Users/55619/Desktop/D_temp/Cap_Div/hypervolume/stats_table_w.csv")

volume_table<-bind_cols(c_n_vol,c_p_vol,c_np_vol,c_ca_vol)
write.csv(volume_table, file = "C:/Users/55619/Desktop/D_temp/Cap_Div/hypervolume/volume_table_w.csv")


####################


############# stats for overlap metrics
 
#between c and ca 

method2_path = hypervolume_permute("c_ca2_permutation", hv_c, hv_ca, n = 999, cores = 5)
result4 = hypervolume_overlap_test(hv_c, hv_ca, method2_path, cores = 4)
resultc_ca = as.data.frame(result4$p_values)
write.csv(resultc_ca, file = "C:/Users/55619/Desktop/D_temp/Cap_Div/hypervolume/STATS_CCA.csv")

# Graphical Results of null sorensen statistic
plot1 <- result4$plots$sorensen +
  xlab("Sorensen distance") +
  ylab("Density") +
  ggtitle("Control x Liming") +
  xlim(0.5, 1)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = c(0,10,16))+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
plot1

#between c and n 

method2_path = hypervolume_permute("c_n2_permutation", hv_c, hv_n, n = 999, cores = 4)
resultcn = hypervolume_overlap_test(hv_c, hv_n, method2_path, cores = 4)
resultc_n = as.data.frame(resultcn$p_values)
write.csv(resultc_n, file = "C:/Users/55619/Desktop/D_temp/Cap_Div/hypervolume/STATS_CN.csv")

# Graphical Results of null sorensen statistic
plotcn = resultcn$plots$sorensen + 
  xlab("Sorensen distance") +
  ylab("Density") +
  ggtitle("Control x N ") +
  xlim(0.5, 1) +
  ylim(0,20)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = c(0,10,19))+
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
plotcn

#between c and p 

method2_path = hypervolume_permute("c_p_permutation", hv_c, hv_p, n = 999, cores = 4)
resultcp = hypervolume_overlap_test(hv_c, hv_p, method2_path, cores = 4)
resultc_p = as.data.frame(resultcp$p_values)
write.csv(resultc_p, file = "C:/Users/55619/Desktop/D_temp/Cap_Div/hypervolume/STATS_CP.csv")

# Graphical Results of null sorensen statistic
plotcp = resultcp$plots$sorensen + 
  xlab("Sorensen distance") +
  ylab("Density") +
  ggtitle("Control x P") +
  xlim(0.5, 1) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = c(0,10,16))+
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

plotcp

#between c and np 

method2_path = hypervolume_permute("c_np_permutation", hv_c, hv_np, n = 999, cores = 4)
resultcnp = hypervolume_overlap_test(hv_c, hv_np, method2_path, cores = 4)
resultc_np = as.data.frame(resultcnp$p_values)
write.csv(resultc_np, file = "C:/Users/55619/Desktop/D_temp/Cap_Div/hypervolume/STATS_CNP.csv")

# Graphical Results of null sorensen statistic
plotcnp = resultcnp$plots$sorensen + 
  xlab("Sorensen distance") +
  ylab("Density") +
  ggtitle("Control x NP") +
  xlim(0.5, 1) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = c(0,10,19))+

 theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

plotcnp



ggarrange(plotcn,plotcp, plotcnp, plot1, ncol=4)


###########################


#plot 

plot(hv_c, contour.type="kde", color="#424C55", cex.axis = 1, show.data=T,show.centroid = T,show.contour = T,
     contour.lwd=1,cex.data = .9,cex.random=.6, cex.centroid=2.7,reshuffle = F,
      verbose=T) 


 plot(hv_c, contour.type="kde", color="#424C55", cex.axis = 1, show.data=T,show.centroid = T,show.contour = T,
     contour.lwd=2,cex.data = 4,cex.random=2, cex.centroid=4,reshuffle = F,
     verbose=T, show.3d=T) 



plot(hv_n, contour.type="kde", color="#9E2B25", cex.axis = 1, show.data=T,show.centroid = T,show.contour = T,
     contour.lwd=1,cex.data = .9,cex.random=.6, cex.centroid=2.7,reshuffle = F,
     verbose=T) 

plot(hv_p, contour.type="kde", color="#104547", cex.axis = 1, show.data=T,show.centroid = T,show.contour = T,
     contour.lwd=1,cex.data = .9,cex.random=.6, cex.centroid=2.7,reshuffle = F,
     verbose=T) 

plot(hv_np, contour.type="kde", color="#B66D0D", cex.axis = 1, show.data=T,show.centroid = T,show.contour = T,
     contour.lwd=1,cex.data = .9,cex.random=.6, cex.centroid=2.7,reshuffle = F,
     verbose=T) 

plot(hv_ca, contour.type="kde", color="#240B36", cex.axis = 1, show.data=T,show.centroid = T,show.contour = T,
     contour.lwd=1,cex.data = .9,cex.random=.6, cex.centroid=2.7,reshuffle = F,
     verbose=T) 

"#61C78B"
"#C763B1"
"#19488F"
"#BD3A17"

"#240B36"
"#424C55"
{
cl=c("#240B36","#BD3A17")
list=hypervolume_join(hv_c,hv_ca)
plotc_ca=plot.HypervolumeList(list,show.density=T,show.legend = F, show.random = T,
                              cex.axis = 0.8, show.data=T,show.centroid = T,show.contour = T,
                              contour.lwd=2.5,contour.type ="kde", colors = cl,cex.data = 0,
                              cex.legend = 0.2, cex.random=0,cex.centroid=3,reshuffle = F,
                               verbose=T,  contour.kde.level = 1e-02)



cl=c("#240B36","#61C78B")
list=hypervolume_join(hv_c,hv_n)
plotc_c_n=plot.HypervolumeList(list,show.density=T,show.legend = F, show.random = T,
                              cex.axis = 0.8, show.data=T,show.centroid = T,show.contour = T,
                              contour.lwd=1.7,contour.type ="kde", colors = cl,cex.data = 0.8,
                              cex.legend = 0.2, cex.random=0.8,cex.centroid=2.7,reshuffle = F,
                              verbose=T)


cl=c("#240B36","#B66D0D")
list=hypervolume_join(hv_c,hv_np)
plotc_np=plot.HypervolumeList(list,show.density=T,show.legend = F, show.random = T,
                               cex.axis = 0.8, show.data=T,show.centroid = T,show.contour = T,
                               contour.lwd=1.7,contour.type ="kde", colors = cl,cex.data = 0.8,
                               cex.legend = 0.2, cex.random=0.8,cex.centroid=2.7,reshuffle = F,
                               verbose=T)

cl=c("#240B36","#C763B1")
list=hypervolume_join(hv_c,hv_p)
plotc_p=plot.HypervolumeList(list,show.density=T,show.legend = F, show.random = T,
                             cex.axis = 0.8, show.data=T,show.centroid = T,show.contour = T,
                             contour.lwd=1.7,contour.type ="kde", colors = cl,cex.data = 0.8,
                             cex.legend = 0.2, cex.random=0.8,cex.centroid=2.7,reshuffle = F,
                             verbose=T)
}
