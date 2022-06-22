# make conceptual figure showing the traj of deleterious allele 
library(ggplot2)
library(cowplot)
library(wesanderson)
library(png)
library("patchwork") 

# get the concept figure
img <- readPNG( "~/projects/pathAdapt/figures/Concept1A.png", native=T)


plA<-ggplot()+ theme_void()
plA<-plA+inset_element(p=img,
                left = 0.,
                bottom = 0.0,
                right =1,
                top = 1)


#first make simulated trajectory
read.table("~/projects/pathAdapt/simData/trajSim.txt")->traj

traj<-data.frame(t=seq(1,length(traj$V1)),f=traj$V1)

plB<-ggplot(traj,aes(t,f))+geom_line(col='#B399D4')+xlab("Time (generations)")+scale_y_continuous(name="Allele frequency")

# next, read in Epi dynamics
read.table("~/projects/pathAdapt/simData/simEpiT.txt")->simEpiT

pal <- wes_palette("Zissou1", 9, type = "continuous")


plC<-ggplot(simEpiT,aes(V2,V1,group=V3,fill=V3))+geom_area()+scale_fill_manual(values = pal,name="Compartments")+ylab("Number of individuals")+scale_x_continuous(name="Time (dimensionless)",limits=c(0,4))


plBC<-plot_grid(plB,plC,labels=c("B","C"),ncol=1)

plABC<-plot_grid(plA,plBC,labels=c("A",""),ncol=2)

plABC

ggsave("~/projects/pathAdapt/figures/conceptFig.pdf",plABC,height=5,width=11)



