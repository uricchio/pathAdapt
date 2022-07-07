library(ggplot2)
library(cowplot)
library(wesanderson)

read.table("~/projects/pathAdapt/simData/simDeltaPfix/sim.deltaFix.txt")->pF

plA<-ggplot(data=pF,aes(V1,V2,col=as.factor(V4)))+geom_point()+geom_line()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="none")+scale_color_manual(values=wes_palette("FantasticFox1"))+xlab(expression(italic(x[i]) * " (initial frequency)"))+ylab(expression(italic(x[f]) * " (final frequency)"))

plB<-ggplot(data=pF,aes(V1,V3,col=as.factor(V4)))+geom_point()+geom_line() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_color_manual(values=wes_palette("FantasticFox1"),name="selection coefficient")+xlab(expression(italic(x[i]) * " (initial frequency)"))+ylab(expression(italic(Delta) * italic(P[f]) * " (change in prob."))
 

plot_grid(plA,plB,labels=c("A","B"),rel_widths=c(1,1.5))

ggsave("~/projects/pathAdapt/figures/Fig3.pdf",height=3,width=9)