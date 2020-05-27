# This code creates Figure 3 and 4.

## Libraries
library(boot)
library(ggplot2)
library(ggExtra)
library(ggpubr)

library(rpart)
library(rpart.plot)
library(rpart.utils)
library(plyr)


## Load matching algorithm and code
library(optmatch)
library(approxmatch)
library(sensitivitymult)
source("fn_match_summv1.R")

#############################################################################
## Simulate the study
#Sample size
n<-600
# two dependent variable
set.seed(81)#51
d1_1<-data.frame(x1=rnorm(n/4,2,1),x2=runif(n/4,0,4),tr=0)
d1_2<-data.frame(x1=rnorm(n/4,6,1),x2=runif(n/4,4,6),tr=0)
d3<-data.frame(x1=rnorm(n/2,4,1),x2=runif(n/2,1,6),tr=1)
mydata<-rbind(d1_1,d1_2,d3)
mydata$y<-rnorm(n,2*mydata$tr+mydata$x1+mydata$x2^2,1)
mydata$ty<-ifelse(mydata$tr==1,"t","c")
mydata$color<-c(rep("gray70",n/2),rep("gray40",n/2))
xl<-min(mydata$x1)-1
xr<-max(mydata$x1)+1
yl<-min(mydata$x2)-1
yr<-max(mydata$x2)+1

var_list <- c("y","x1","x2")
table(mydata$ty)
#############################################################################


#################################################
############################################
  # propensity score stratification
  
  glmfit<-glm(tr~x1+x2,family=binomial,
              x=TRUE,y=TRUE,data=mydata)
  #summary(glmfit)
  mydata$treatment<-glmfit$y
  propscore<-predict(glmfit,type="response")
  
  # Create subclasses
  no.subclasses<-10
  cutoffs<-quantile(propscore,seq(1/no.subclasses,(no.subclasses-1)/no.subclasses,1/no.subclasses))
  subclass<-rep(1,length(propscore))
  for(i in 1:(no.subclasses-1)){  subclass<-subclass+1*(propscore>=cutoffs[i])}
  mydata$subclass<-subclass
  
  
  mydata$ps <- propscore

  match_summv1(data1 = mydata,covrt = "ps",cls="subclass",trt="treatment")$b


#################################################
############################################
  ##############################################################
  # decision tree matching
  ms<-0
  fit <- rpart(tr~x1+x2,
               method="class", data=mydata,parms=list(split = "gini"),
               control = rpart.control(minsplit =ms,cp = 0.0,maxcompete = 20))

  # choose cp in a semi-autonomous manner
  mincpIdx <- which.min(fit$cptable[fit$cptable[,2]<11,1])	# At least 10 nodes
  
  for(cpidx in (mincpIdx):nrow(fit$cptable)){
   cp_choose = fit$cptable[cpidx,1]
   
   pfit <- prune(fit,cp=cp_choose)
   mydata$subclass_tree<-pfit$where
   
   sd_x1 = match_summv1(data1 = mydata,covrt = "x1",cls="subclass_tree",trt="treatment")$b
   sd_x2 = match_summv1(data1 = mydata,covrt = "x2",cls="subclass_tree",trt="treatment")$b
   print(c(cp_choose, sd_x1, sd_x2))
   if( abs(sd_x1) < .1 &   abs(sd_x2) < .1 ){	
	break;
   }
  }

  table(mydata$subclass_tree, mydata$treatment)

 
 
######### PLOTS : Figure 3 #####################################

#windows()
#dev.set(3)
mydata$x22<-mydata$x2-0.2
mydata$tr1<-as.factor(mydata$tr)
ggscatterhist(
  mydata, x = "x1", y = "x2",
  xlab=expression(paste("\n",x[1])), 
  ylab=expression(paste("\n",x[2])),
  	font.x  = c(16, "plain", "grey5"), font.y = c(16, "plain", "grey5"),
	font.legend  = c(14, "plain", "grey5"),  font.tickslab = c(14, "plain", "grey5"),
  color = "tr1", point=T,size = rep(1.5,n), alpha = 0.6,shape="tr1",
  #label="ty",label.size=0.1,stat="identity",vjust=-0.1,
  palette = c("grey55", "grey5"),
  xlim=c(-.5, 8.15), ylim=c(-.55, 6.15), 
  margin.plot = "hist", 
  margin.params = list(fill = "ty", color = "black", size = 0.0),
  ggtheme = theme_set(theme_bw() + theme(
			axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
			axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
			panel.background = element_rect(fill = "white", colour = "grey5"),
			legend.key=element_blank(),
			legend.margin = margin(-10,0,0,0))),
  legend.title="Treatment status",legend="bottom"#,
  #xticks.by = 2,#c(0,2, 4, 6,8),
  #yticks.by = 2,#c(0,2, 4, 6)
)
#############################################################################




#############################################################################
# FIGURE 4(a)
#####------------------------------
# plot propensity strata 
glmline<-data.frame(b0<-rep(0,no.subclasses-1),b1<-rep(0,no.subclasses-1),
                    xtxt<-rep(0,no.subclasses-1),ytxt<-rep(0,no.subclasses-1))
windows()
xtxt <- c(2.5, 4.3, 5.7, 6.45, 7.2, 7.5, 7.8, 8.0, 8.15)
delta <- c(-1, -.75, -.27, -.29, -.2, -.37, -.4, -.23, -.25)
deltax <- c(-.35, -.26, rep(-.2, 7))
colnames(glmline)<-c("b0","b1","xtxt","ytxt")
for(i in 1:length(cutoffs)){
  glmline$b0[i]<-(logit(cutoffs[i])-glmfit$coefficients[1])/glmfit$coefficients[3]
  glmline$b1[i]<- (glmfit$coefficients[2]/glmfit$coefficients[3])*-1
  xtxt1<-glmline$b0[i]+(xl+1)*glmline$b1[i]
  xtxt2<-ifelse(i==1,glmline$b0[i],glmline$b0[i]+(xl+1)*glmline$b1[i])

  xtxt3<-ifelse(i==1,glmline$b0[i],glmline$b0[i]+(xtxt[i])*glmline$b1[i])
  xtxt3<- glmline$b0[i]+(xtxt[i])*glmline$b1[i]

  xtxt4 <- (-.45-glmline$b0[i])/glmline$b1[i]
  
  print(c(xtxt1, xtxt2))
  glmline$xtxt[i]<- xtxt[i] #+ deltax[i] #xtxt[i] #(xl+1)
  glmline$ytxt[i]<- (xtxt3)+delta[i]
  
}

 rbind(glmline$xtxt, glmline$ytxt)


mydata$tr1<-as.factor(mydata$tr)

dev.set(4)
#windows()
ggplot(mydata, aes(x=x1, y=x2, color=color,shape=tr1)) +geom_point(size=1.5)+
  scale_color_grey(start = 0.05, end = 0.55)+theme_bw()+
  theme(	axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
		axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
		legend.position = "none", axis.title=element_text(colour="grey5", size=rel(1.5)),
				axis.text=element_text(colour="grey5", size=rel(1.25)),
	panel.background = element_rect(colour = "black"))+
  geom_abline(intercept = glmline$b0,slope=glmline$b1,colour="grey25")+
  xlab(expression(x[1])) + ylab(expression(x[2])) + 
  annotate("text",x=glmline$xtxt,y=glmline$ytxt,size=5,label=c(1:length(cutoffs)),color="grey20") + 
  annotate("text",x=8.25,y=4.5,label=10,color="grey25", size=5) + 
  #xlim(-.5, 8.86) + ylim(-.58, 6.15) +
  scale_x_continuous(limits=c(-.65, 9), breaks=c(0,2, 4, 6,8)) + 
  scale_y_continuous(limits=c(-.55, 6.15), breaks=c(0,2, 4, 6))



#dev.off()

#############################################################################
## FIGURE 4(b)

#--------------------------------------------------------
# plot tree leaf as strata 
ule_df <- rpart.rules.table(pfit) %>%
  dplyr::filter(Leaf==TRUE) %>%
  dplyr::group_by(Rule) %>%
  dplyr::summarise(Subrules = paste(Subrule, collapse=","))


sunruledf<-rpart.subrules.table(pfit)

labr<-sort(unique(mydata$subclass_tree))

xleft<-xright<-ybot<-ytop<-rep(0,nrow(ule_df))
for(i in 1:nrow(ule_df)){
  s<-sunruledf[which(sunruledf$Subrule %in% strsplit(ule_df$Subrules[i],split=",")[[1]]),]
  
  xleft1<-max(as.numeric(as.character(s$Greater[s$Variable=="x1"])),na.rm = T)
  ybottom1<-max(as.numeric(as.character(s$Greater[s$Variable=="x2"])),na.rm = T)
  xright1<-min(as.numeric(as.character(s$Less[s$Variable=="x1"])),na.rm = T)
  ytop1<-min(as.numeric(as.character(s$Less[s$Variable=="x2"])),na.rm = T)
  
  xleft1<-ifelse(is.na(xleft1)|is.infinite(xleft1),xl+.84,xleft1)
  xright1<-ifelse(is.na(xright1)|is.infinite(xright1),xr-.45,xright1)
  ybottom1<-ifelse(is.na(ybottom1)|is.infinite(ybottom1),yl+.45,ybottom1)
  ytop1<-ifelse(is.na(ytop1)|is.infinite(ytop1),yr-.85,ytop1)
  
  xleft[i]<-xleft1
  xright[i]<-xright1
  ybot[i]<-ybottom1
  ytop[i]<-ytop1
  
}


thisps <- labr %in% c(16, 17, 25, 26)
dev.set(6)
windows()
ggplot(mydata, aes(x=x1, y=x2, color=color,shape=ty)) +geom_point(size=1.2)+
  scale_x_continuous(limits=c(-.65, 9), breaks=c(0,2, 4, 6,8)) + 
  scale_y_continuous(limits=c(-.55, 6.15), breaks=c(0,2, 4, 6)) +
   theme_classic()+theme(legend.position = "none")+
  annotate("rect",xmin=xleft,ymin=ybot,xmax=xright,ymax=ytop,color="gray40",alpha=0)+
  annotate("text",x=((xleft+xright)/2)[thisps],y=((ybot+ytop)/2)[thisps],label=labr[thisps],
			size=3.2,color="black")+
 scale_color_grey(start = 0.05, end = 0.55)+theme_bw()+
  theme(
		axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
		axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
		legend.position = "none", axis.title=element_text(colour="grey5", size=rel(1.5)),
				axis.text=element_text(colour="grey5", size=rel(1.25)),
	panel.background = element_rect(colour = "black"))+
  xlab(expression(x[1])) + ylab(expression(x[2])) + 
 annotate("text",x=.5,y=2,label=3,size=4.2,color="black")+
 annotate("text",x=1,y=5.1,label=4,size=4.2,color="black")+
 annotate("text",x=2.9,y=.75,label=12,size=4.2,color="black")+
 annotate("text",x=4.2,y=.05,label=20,size=4.2,color="black")+
 annotate("text",x=7,y=1.1,label=8,size=4.2,color="black")+
 annotate("text",x=8,y=4.5,label=7,size=4.2,color="black")+
 annotate("text",x=2.8,y=3.5,label=14,size=4.2,color="black")+
 annotate("text",x=2.9,y=5.75,label=18,size=4.2,color="black")+
 annotate("text",x=4.8,y=3.3,label=28,size=4.2,color="black")+
 annotate("text",x=3.85,y=4.4,label=29,size=4.2,color="black")+
 annotate("text",x=4.75,y=5.75,label=27,size=4.2,color="black")
  

#dev.off()


