library(ggplot2)
library(cowplot)
library(wesanderson)
library(viridis)


# first plot the tradeoff models we will use

tradeOff<-function(s,sHalf,scale,BetaMax) {
	return (BetaMax/(1+exp(scale*(log(s)-log(sHalf)))))
}

xvals= c()
i = 0.0000001
while(i < 0.1) {
	i <- i*1.1
	xvals<-c(xvals,i)
}
df<-data.frame()

Lambda*(mu+r) = (Beta/Lambda)-r


#Severe tradeoff
BetaMax = 0.001
scale = 3
sHalf = 0.0005
df<-data.frame(cbind(x=as.numeric(xvals),y=tradeOff(xvals,sHalf,scale,BetaMax),model=rep("Severe",length(xvals))))

#Diminishing returns
BetaMax = 0.001
scale = 1
sHalf = 0.00005
df<-rbind(df, cbind(x=xvals,y=tradeOff(xvals,sHalf,scale,BetaMax),model=rep("Diminishing Returns",length(xvals))))

#Nearly Neutral
BetaMax = 0.001
scale = 2
sHalf = 0.000005
df<-rbind(df, cbind(x=xvals,y=tradeOff(xvals,sHalf,scale,BetaMax),model=rep("Nearly Neutral",length(xvals))))

pal<-wes_palette("Darjeeling1")

plA<-ggplot(df,aes(as.numeric(x),as.numeric(y),color=model))+geom_line(size=1)+scale_color_manual(values=pal)+scale_x_continuous(trans="log10")+xlab(expression("selection coefficient (|" * italic(s) * "|)"))+ylab(expression(mu))+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# plB, plC will be constraints on trait values

dfC<-data.frame()
Lambda=2;

for (r in seq(0,1,0.01)) {
	for (mu in seq(0,1,0.01)) {
		dfC<-rbind(dfC,data.frame(r=r,mu=mu,beta=Lambda*(r+mu)))
	}
}

pal <- wes_palette("Zissou1", 30, type = "continuous")

plC<-ggplot(dfC,aes(r,mu,fill=beta))+geom_tile()

plC<-ggplot(dfC,aes(r,mu,fill=beta))+geom_tile()+xlab(expression(italic(r)))+ylab(expression(mu))+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+scale_y_continuous(limits=c(0,1),expand=c(0,0)) + labs(fill=expression(beta))+scale_x_continuous(expand = c(0, 0),limits=c(0,1)) +scale_fill_viridis(option="magma")

plC
plot_grid(plA,plC,ncol=2,labels=c("A","B"),rel_widths=c(1.7,1))
ggsave("~/projects/pathAdapt/figures/Fig2.pdf",width=8,height=2.2)


