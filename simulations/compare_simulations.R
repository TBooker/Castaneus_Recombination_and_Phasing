rm(list=ls())
## script to compare the output of the recombination mapping with and without phasing errors

true_map <-read.csv("project/6.Phasing_Castaneus+Recombination_Rates/simulations/chr1_segment/chr1_segment_true.csv")
swit <- read.table("project/6.Phasing_Castaneus+Recombination_Rates/simulations/chr1_segment/with_switches/chr1_segment_N10_with_switches_long.b200.csv",sep=",",head=T)
b100 <- read.csv("project/6.Phasing_Castaneus+Recombination_Rates/simulations/chr1_segment/with_switches/chr1_segment_N10_with_switches.b100.csv")

true_map$sour <-rep("Recombination Map",nrow(true_map))
str(swit)
str(b100)
library(plyr)
b100<-rename(b100,c("p0.050"="p0.025"))
str(swit)
swit$sour <- rep("b200",nrow(swit))
b100$sour <- rep("b100",nrow(b100))
all <- rbind(swit,b100)
all$pos <-(all$right_snp+all$left_snp)/2



library(ggplot2)
str(all)
str(true_map)
ggplot(data=all,aes(x=pos/1000,y=mean,col=sour))+
  geom_line(lwd= 1.5)+
  geom_line(data=true_map,aes(x=true_map$pos/1000,y=true_map$mean),lwd = 1.5)+
  scale_color_discrete("")+
  scale_y_continuous(limits=c(0,0.2))+
    xlab("Position (bp)")+
  ylab("Rho/bp")+
  theme(axis.title.x = element_text(face="bold", colour="black", size=20),
        axis.text.x  = element_text(vjust=0.5, size=16),
        axis.title.y = element_text(face="bold", colour="black",  size=22,angle=90),
        axis.text.y  = element_text(size=16),
        strip.text.x = element_text(face="bold",size=20))

