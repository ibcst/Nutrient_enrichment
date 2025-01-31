###Isabela C. 2023
#ibc10@case.edu

############################################################

#NMDS 

############################################################
library(vegan)
library (ggplot2)
library(permute)
library(lattice)
library(PERMANOVA)


setwd("C:/Users/55619/Desktop/D_temp/Cap_Div/taxonomic_div/Diversity")
d <- read.table("sp_parc_abund.csv",sep=";",header=T, fileEncoding="latin1") 
d

#make community matrix - extract columns with abundance information
input = d[,4:ncol(d)]
str(input)
#turn abundance data frame into a matrix
m_input = as.matrix(input)


any(is.na(m_input))
summary(m_input)

set.seed(123)
nmds = metaMDS(m_input, distance = "bray", trymax=99, k=2 )
nmds
plot(nmds)
scores(nmds)

####obtain scores
scoresp <- read.table("score_sp.csv",sep=";",header=T, fileEncoding="latin1") 

scorespred <- read.table("sp_redox.csv",sep=";",header=T) 

scoresite <- read.table("score_site.csv",sep=";",header=T, fileEncoding="latin1") 


#ggplot plot option

ggplot() + 
  geom_polygon(data=scoresite,aes(x=NMDS1,y=NMDS2,fill=sites,group=sites),alpha=0.75) + # add the convex hulls
  geom_text(data=scorespred,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.7) +
  geom_point(data=scoresite,aes(x=NMDS1,y=NMDS2,shape=sites,colour=sites),size=3, shape=21) + # add the point markers
  scale_colour_manual(values=c("C" = "grey", "N" = "#61C78B", "NP"="#19488F", "P"="#C763B1", "Ca"="#BD3A17")) +
  scale_fill_manual(values=c("C" = "grey", "N" = "#61C78B", "NP"="#19488F", "P"="#C763B1", "Ca"="#BD3A17")) +
  #annotate("text", x = c(0.2), y = c(0.2), label = c("Stress= 0.24") , color="black", size=7 , angle=0, fontface="bold")
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=11), # remove x-axis labels
        axis.title.y = element_text(size=11), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())


#vegan plot option

try = c("C" = "#240B36", "N" = "#BD3A17", "NP"="#61C78B", "P"="#Bb671c", "Ca"="#913D88")

TRAT <- d[,1]

pchv <- 1:7
colv <- 1:7
{
timeLine <- c(-1, +1)

  plot.new()
  plot.window(xlim = timeLine, ylim = timeLine, asp = 1)
  #plot(nmds, display = "sites")
  ordihull(nmds, TRAT, col= try, lwd=2)
  #with(scorespred, text(species, labels = rownames(species),
            #      col = "grey", cex = 0.7))
 # text(nmds, display = "species", cex = 0.8, col = "grey")
  #orditorp(nmds, "species", pcol="grey",cex=0.7)
  ordiellipse(nmds, TRAT, conf=0.95,col=try, kind = "ehull", lwd=2)
  ordiellipse(nmds, TRAT, col=try, draw="polygon")
  ordispider(nmds, TRAT, col=try, label = TRUE)
  ordispider(nmds, TRAT, col=try, w=weights(nmds,display="sites"),label = TRUE)
  points(nmds, disp="sites", pch=21, col="black", bg="white", cex=1.3)
  axis(side = 1)
  axis(side = 2)
  title(xlab = "NMDS 1", ylab = "NMDS 2")
  box()

}



especies <- d[3:ncol(d)]
TRAT <- as.factor(TRAT)

dist.bray <- vegdist(m_input, method = "bray", binary = TRUE)

btadisper = betadisper(dist.bray, TRAT, type = c("median"), bias.adjust = FALSE,
                       sqrt.dist = FALSE, add = FALSE)
tukey = TukeyHSD(btadisper, which = "group", ordered = FALSE,
                 conf.level = 0.95)
tukey


P <- PERMANOVA(dist.bray, TRAT, nperm = 999)


dist.bray


C = matrix(c(0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1), nrow=4, byrow=TRUE)
rownames(C)=c("C - Ca", "C - N", "C - NP", "C - P")


all.sor <- DistContinuous(m_input, coef="Bray_Curtis")

res <- PERMANOVA(all.sor,C=C, TRAT, CoordPrinc=TRUE)
summary(res)
plot(res)
