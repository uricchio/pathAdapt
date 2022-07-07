library(ggplot2)
library(cowplot)
library(wesanderson)
library(viridis)

# plA, alpha total
read.table("~/projects/pathAdapt/simData/alpha.nearlyNeutral.txt")->neutAl
read.table("~/projects/pathAdapt/simData/alpha.diminishingReturns.txt")->dimAl
read.table("~/projects/pathAdapt/simData/alpha.severe.txt")->sevAl

dfNeut<-cbind(neutAl,data.frame(label=rep("nearly netural",length(neutAl$V1))))
dfDim<-cbind(dimAl,data.frame(label=rep("diminishing returns",length(dimAl$V1))))
dfSev<-cbind(sevAl,data.frame(label=rep("severe",length(sevAl$V1))))

dfAll<-rbind(dfNeut,dfDim,dfSev)

plA<-ggplot(dfAll,aes(V4,V1,col=label))+geom_point()+geom_line()+scale_y_log10(name=expression(italic(alpha)),breaks=c(1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1,3e-1),labels = scales::scientific)+xlab(expression(italic(mu) * " (death rate)"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="none")+scale_color_manual(values=wes_palette("Darjeeling1"),name="Model")

# plB, per s site
read.table("~/projects/pathAdapt/simData/alpha.perS.nearlyNeutral.txt")->neutAl
read.table("~/projects/pathAdapt/simData/alpha.perS.diminishingReturns.txt")->dimAl
read.table("~/projects/pathAdapt/simData/alpha.perS.severe.txt")->sevAl


dfNeut<-cbind(neutAl,data.frame(label=rep("nearly netural",length(neutAl$V1))))
dfDim<-cbind(dimAl,data.frame(label=rep("diminishing returns",length(dimAl$V1))))
dfSev<-cbind(sevAl,data.frame(label=rep("severe",length(sevAl$V1))))

dfAll<-rbind(dfNeut,dfDim,dfSev)

plB<-ggplot(dfAll,aes(abs(V1),V2,col=label))+geom_point()+geom_line()+scale_x_log10(name=expression("|" * italic(s) * "|"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_color_manual(values=wes_palette("Darjeeling1"),name="Model")+ylab("excess substitutions per site")

plot_grid(plA,plB,rel_widths=c(1,1.5),labels=c("A","B"))

ggsave("~/projects/pathAdapt/figures/Fig4.pdf",height=3,width=9)