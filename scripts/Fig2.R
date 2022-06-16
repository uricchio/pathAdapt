library(ggplot2)
library(cowplot)
library(wesanderson)

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

plA<-ggplot(df,aes(as.numeric(x),as.numeric(y),color=model))+geom_line()+scale_color_manual(values=pal)+scale_x_continuous(trans="log10")+xlab(expression("selection coefficient (|" * italic(s) * "|)"))+ylab(expression(mu))+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# plB, plC will be constraints on trait values
pal<-wes_palette("Darjeeling2")

Beta<-function(muPr,Lambda) {
	return(Lambda*(muPr))
}

Mu<-function(r,LambdaTb) {
	return (LambdaTb - r)
}

lVals = c(0.05,0.1,0.2,0.5,1,2)
muPrv = seq(0.001,0.2,0.001)
dfB<-data.frame()
i <-1;
while (i < length(lVals)) {
  dfB<-rbind(dfB,data.frame(muPr = muPrv, b = Beta(muPrv,lVals[i]),Lambda=rep(lVals[i],length(muPrv))))
  i<-i+1
}

plB<-ggplot(dfB,aes(as.numeric(muPr),as.numeric(b),group=Lambda,color=factor(Lambda)))+geom_line()+scale_color_manual(values=pal)+xlab(expression(mu * " + " * italic(r)))+ylab(expression(beta))+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ labs(color=expression(Lambda )) 
plB


lTbVals = c(0.05,0.1,0.2,0.5,1,2)*0.1
rv = seq(-0.1,0.2,0.00001)
dfC<-data.frame()
i <-1;
while (i < length(lVals)) {
  dfC<-rbind(dfC,data.frame(r = rv, mu = Mu(rv,lTbVals[i]),LambdaTb=rep(lTbVals[i],length(rv))))
  i<-i+1
}

plC<-ggplot(dfC,aes(r,mu,group=LambdaTb,color=factor(LambdaTb)))+geom_line(lwd=0.8)+scale_color_manual(values=pal)+xlab(expression(italic(r)))+ylab(expression(mu))+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+scale_y_continuous(limits=c(0,1e-1),expand=c(0,0)) + labs(color=expression(Lambda * beta ))+scale_x_continuous(expand = c(0, 0),limits=c(0,0.1)) 



plC
#plBC<-plot_grid(plB, plC, labels=c("B","C"))
plot_grid(plA,plC,ncol=2,labels=c("A","B"),rel_widths=c(1.8,1))

ggsave("~/projects/pathAdapt/figures/Fig2.pdf",width=11,height=2.5)


