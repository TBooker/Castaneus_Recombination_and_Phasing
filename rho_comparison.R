rm(list=ls())
library(ggplot2)

x<-read.csv("/media//booker//Elements//Mouse_data//6.Phasing_Castaneus//chr5_19_comparison_RHO.csv")

chr5<- x[which(x$chr=="chr5"),]
chr19<- x[which(x$chr=="chr19"),]
plot(chr5$rho~chr5$source,log="y",ylab ="log10(4Ner)", xlab="",main ="Chr5")
plot(chr19$rho~chr19$source,log="y",ylab ="log10(4Ner)", xlab="",main = "Chr19")


summary(mod)
step(mod)
ggplot(data = x, aes(x=rho))+
  geom_density(binwidth = 0.020)+
  scale_x_log10()+
  facet_grid(chr ~ source)
summary(x)
