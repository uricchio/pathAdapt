library(ggplot2)
library(cowplot)

read.table("~/projects/boots/pathAdapt/simData/total_by_freq_2.txt")->t

df<-data.frame(x=rep(-t$V1,3),y=c(t$V2,t$V5,t$V6),model=(c(rep("difference",55), rep("no spillover",55),rep("rare spillover",55))))
df


plA<-ggplot(data=df,aes(x,y,color=model))+geom_blank()+geom_point()+scale_x_log10()+geom_line(lty=3)+xlab(expression("|" * italic(S) * "|"))+ylab(expression(italic(D)))+theme_bw() 
plA
plB <-ggplot(data=t,aes(-V1,V4))+geom_point()+scale_x_log10()+ylab(expression(italic(beta)))+xlab(expression("|" * italic(S) * "|"))+geom_line()+geom_hline(yintercept=0.003,lty=3)+theme_bw() 
plB      
       

plot_grid(plA, plB,labels=c("A","B"),rel_widths=c(1.4,1))
