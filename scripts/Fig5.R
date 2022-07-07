library(ggplot2)
library(cowplot)
library(wesanderson)
library(viridis)

read.table("~/projects/pathAdapt/simData/Fig5Data.dimReturns.txt")->dimAl
read.table("~/projects/pathAdapt/simData/Fig5Data.nearNeutral.txt")->neutAl
read.table("~/projects/pathAdapt/simData/Fig5Data.severe.txt")->sevAl

dfNeut<-cbind(neutAl,data.frame(label=rep("nearly netural",length(neutAl$V1))))
dfDim<-cbind(dimAl,data.frame(label=rep("diminishing returns",length(dimAl$V1))))
dfSev<-cbind(sevAl,data.frame(label=rep("severe",length(sevAl$V1))))

dfAll<-rbind(dfNeut,dfDim,dfSev)

plA<-ggplot(dfAll,aes(V3,V1,col=label))+geom_line()+geom_point()+scale_x_log10(name=expression( italic(N) ))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")+scale_color_manual(values=wes_palette("Darjeeling1"),name="Model")+scale_y_continuous(name=expression(alpha),limits=c(0,0.6))
        
        
read.table("~/projects/pathAdapt/simData/Fig5Data.mid.dimReturns.txt")->dimAl
read.table("~/projects/pathAdapt/simData/Fig5Data.mid.nearNeutral.txt")->neutAl
read.table("~/projects/pathAdapt/simData/Fig5Data.mid.severe.txt")->sevAl

dfNeut<-cbind(neutAl,data.frame(label=rep("nearly netural",length(neutAl$V1))))
dfDim<-cbind(dimAl,data.frame(label=rep("diminishing returns",length(dimAl$V1))))
dfSev<-cbind(sevAl,data.frame(label=rep("severe",length(sevAl$V1))))

dfAll<-rbind(dfNeut,dfDim,dfSev)

plB<-ggplot(dfAll,aes(V3,V1,col=label))+geom_line()+geom_point()+scale_x_log10(name=expression( italic(N) ))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")+scale_color_manual(values=wes_palette("Darjeeling1"),name="Model")+scale_y_continuous(name=expression(alpha),limits=c(0,0.6))
        
        
read.table("~/projects/pathAdapt/simData/Fig5Data.low.dimReturns.txt")->dimAl
read.table("~/projects/pathAdapt/simData/Fig5Data.low.nearNeutral.txt")->neutAl
read.table("~/projects/pathAdapt/simData/Fig5Data.low.severe.txt")->sevAl

dfNeut<-cbind(neutAl,data.frame(label=rep("nearly netural",length(neutAl$V1))))
dfDim<-cbind(dimAl,data.frame(label=rep("diminishing returns",length(dimAl$V1))))
dfSev<-cbind(sevAl,data.frame(label=rep("severe",length(sevAl$V1))))

dfAll<-rbind(dfNeut,dfDim,dfSev)

plC<-ggplot(dfAll,aes(V3,V1,col=label))+geom_line()+geom_point()+scale_x_log10(name=expression( italic(N) ))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_color_manual(values=wes_palette("Darjeeling1"),name="Model")+scale_y_continuous(name=expression(alpha),limits=c(0,0.6))
        
plot_grid(plA,plB,plC,rel_widths=c(1,1,1.5),ncol=3,labels=c("A","B","C"))
        
ggsave("~/projects/pathAdapt/figures/Fig5.pdf",height=3,width=9)
        
        