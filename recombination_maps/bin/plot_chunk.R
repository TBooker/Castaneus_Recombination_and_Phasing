rm(list=ls())
b100 <- read.csv('project/6.Phasing_Castaneus+Recombination_Rates/recombination_maps/chr17/processed/chr17_15457850-15733984.b100.txt', head=F, skip=3, sep =' ') 
b100$pos <- 15457850+b100$V1+(b100$V2-b100$V1)/2
b100$bp <- '100'


b10 <- read.csv('project/6.Phasing_Castaneus+Recombination_Rates/recombination_maps/chr17/processed/chr17_15457850-15733984.b10.txt', head=F, skip=3, sep =' ') 
b10$pos <- 15457850+b10$V1+(b10$V2-b10$V1)/2
b10$bp <- '10'

all<- rbind(b100, b10)

with(
  b100, plot(pos/1e6, V3, type ='l', log ='y')
)
with(
  b10, lines(V3 ~ pos/1e6, type ='l', lty ='dashed')
)
library(ggplot2)

expression('rho')

ggplot(data = all, aes(x = pos/1e6, y= V3, col=bp))+
  geom_line()+
  scale_y_log10(expression(rho^-bp))+
  scale_x_continuous('Position (Mbp)')+
  scale_color_discrete(name="Block\nPenalty")+
  theme_bw()
#  theme(legend.title=element_blank())

