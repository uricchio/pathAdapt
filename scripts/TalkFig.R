library(ggplot2)
library(cowplot)
library(wesanderson)

read.table("~/projects/pathAdapt/simData/testRatePop.txt")->ratePop

ratePop<-cbind(ratePop,data.frame(model=rep("Infectious disease",length(ratePop$V1))))

adap<-data.frame(V1=c(50,100,200,500,1000,2000,5000,10000),V4=c(0.1,0.15,0.2,0.3,0.5,0.7,0.8,0.85),model=rep("Conventional pop gen"))

allData<-rbind(newF,adap)

ggplot(adap,aes(V1,V4,group=model,col=model))+geom_point(cex=2)+geom_line(lwd=0.8)+xlab("population size")+ylab("adaptation rate")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_x_log10()+scale_color_manual(values=wes_palette("Darjeeling1"))
