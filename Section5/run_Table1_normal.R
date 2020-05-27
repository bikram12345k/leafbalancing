#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


iter=as.numeric(args)[1]
set.seed(iter)
print(iter)

### This code has been set up to run on a cluster.
### Run the code in 900 parallel cores to reproduce table 1.
## Locally this code can be run.  Currently it is set up to run 
##	12 iterations.  This can be increased.

library(rmutil)
library(rpart)
library(plyr)

## Load matching algorithm and code
library(optmatch)
library(approxmatch)
library(sensitivitymult)
library(sensitivitymv)
library(senstrat)
source("fn_match_summv1.R")

library(doSNOW, foreach)	

deprintize<-function(f){
 return(function(...) {capture.output(w<-f(...));return(w);});
}	

cl<-makeCluster(1) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

# simulate from Tukey's lambda distribution	# Not used.
rtld <- function(n, l){	
	sam = sapply(runif(n,0,1), function(p) (p^l - (1-p)^l)/l)
	v = ifelse(l > -1/2, (2/l^2)*( 1/(1+2*l) - gamma(l+1)^2/gamma(2*l+2) ), 1)
	sam/sqrt(v)
}

## Study size
n<-600
nn<- 1

for(iter2 in 1:12){	# 12 iterations
  set.seed(iter*900+iter2)	# set seed for replication.
  print(iter*900+iter2)

  ## GENERATE DATA
  #  Covariates and the treatment assignment
  d1_1<-data.frame(x1=rnorm(n/4,2,1),x2=runif(n/4,0,4),tr=0)
  d1_2<-data.frame(x1=rnorm(n/4,6,1),x2=runif(n/4,4,6),tr=0)
  d3<-data.frame(x1=rnorm(n/2,4,1),x2=runif(n/2,1,6),tr=1)

  mydata<-rbind(d1_1,d1_2,d3)
  mydata$ty<-ifelse(mydata$tr==1,"t","c")  

    
  ##############################################################
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
  

  ##############################################################
  # decision tree matching
  ms<-3
  fit <- rpart(tr~x1+x2,#I(2*x1+x2)+I(x1+2*x2),
               method="class", data=mydata,parms=list(split = "information"),
               control = rpart.control(minsplit =ms,cp = 0.0,maxcompete = 20))
  
  #printcp(fit) # display the results
  

  # Automatic strategy to choose a level of prunning. 
  # start with a reasonable number of leaf nodes.
  mincpIdx <- which.min(fit$cptable[fit$cptable[,2]<1.5*sqrt(n/ms),1])

  cp_choose1 = 10
  best_smd1 <- c(100,100, 100, 100)
	mydata$x2.2 = mydata$x2^2
	mydata$x1.2 = mydata$x1^2

  # Search for a good cp.
  for(cpidx in (mincpIdx-1):nrow(fit$cptable)){
   cp_choose = fit$cptable[cpidx,1]
      
   pfit <- prune(fit,cp=cp_choose)
   mydata$subclass_tree<-pfit$where
   
   # Calculate SMD for the covariates and their squares
   best_smd = c(match_summv1(data1 = mydata,covrt = "x1",cls="subclass_tree",trt="treatment")$b, 
			match_summv1(data1 = mydata,covrt = "x2",cls="subclass_tree",trt="treatment")$b,
			match_summv1(data1 = mydata,covrt = "x1.2",cls="subclass_tree",trt="treatment")$b,
			match_summv1(data1 = mydata,covrt = "x2.2",cls="subclass_tree",trt="treatment")$b)

   
   if(all(abs(best_smd) < abs(best_smd1)+.005)){
		best_smd1 = best_smd
		cp_choose1 = cp_choose
	}

   print(c(cp_choose, best_smd))
   # Stop when reasonable balance is achieved
   if((mean(abs(best_smd)) < .08 & max(abs(best_smd)) < 0.12) | all(abs(best_smd) < .08)){
	cp_choose1 = cp_choose
	break;
   }
  }

  print(cp_choose1)
  pfit <- prune(fit,cp=cp_choose1)
  mydata$subclass_tree<-pfit$where

  ##############################################################
  ## pair matching in the leaf nodes

  mydata$SEQN <- rownames(mydata)
  mydata$newtr<-mydata$tr+1

  tabulate <- table(mydata$newtr, mydata$subclass_tree)
  nonpurenodes = as.numeric(colnames(tabulate)[apply(tabulate, 2, function(x) prod(x)>0)])

  distmat <- multigrp_dist_struc(mydata, 'newtr', 
                                  list(mahal=c("x1","x2")), wgts=1)

  allmatch <- NULL
  for(s in 1:length(nonpurenodes)){
	units <- rownames(mydata)[mydata$subclass_tree == nonpurenodes[s]] 
  	unitstr <- units[mydata[units,'newtr']==1]
  	unitscl <- units[mydata[units,'newtr']==2]
      distmat1 <- distmat
      distmat1[[1]] <- distmat[[1]][unitstr,unitscl,drop=FALSE]
	match <- deprintize(kwaymatching)(distmat1, 'newtr', .data=mydata[units,])
	allmatch <- rbind(allmatch, match$matches)     
  }

  
  matches <- data.frame(SEQN = as.numeric(allmatch), gp = rep(1:nrow(allmatch), 2))
  ## Matched dataset
  mydata_treematch <- join(mydata, matches )
  

  ##############################################################3
  # mahalnobis + caliper

  d5<-mydata
  d5$SEQN<-rownames(d5)
  d5$newtr<-d5$tr+1
  distmat1 <- multigrp_dist_struc(d5, 'newtr', 
                                  list(prop=c("x1","x2")), wgts=1)
  distmat2 <- multigrp_dist_struc(d5, 'newtr', 
                                  list(mahal=c("x1","x2")), wgts=1)
  
  distmat<-distmat2  
  distmat[[1]][distmat1[[1]]>.2] <- 500
  
  match <- deprintize(kwaymatching)(distmat, 'newtr', .data=d5)
    
  d5_match <- d5[rownames(d5) %in% as.vector(match$matches),]
  
  matches <- match$matches
  matches <- as.vector(matches)
  matches <- data.frame(SEQN = as.numeric(matches), gp = rep(1:nrow(match$matches), 2))
  ## Matched dataset
  d5_match <- join(d5_match, matches)
  
  
  #############################################################################
  ## INFERNECE FOR THE AVERAGE TREATMENT EFFECT 
  #############################################################################
 

  # To invert the test and find estimates and the confidence intervals
  #	stratified design.
  ff.2 <- function(x,y,z, cl, x1, x2, dev){	
    newy = y-x*z
     # from senstrat package with average difference in the outcomes as the 
	 #	test statistic.
    senstrat(newy, z, cl, gamma = 1)$Result[2]-dev
  }
  
   
  # To invert the test and find estimates and the confidence intervals
  #	matched pairs design.
  ff1.2<-function(x,y,z, cl, x1, x2, dev){
	 y <- y[!is.na(cl)]
	 x1 <- x1[!is.na(cl)]
	 x2 <- x2[!is.na(cl)]
	 z <- z[!is.na(cl)]
	 gp <- cl[!is.na(cl)]
    
	newy = y-x*z;
	newyd <- sapply(unique(gp), function(g){
			id <- which(gp==g)
			(newy[id[1]] - newy[id[2]])*(z[id[1]]-z[id[2]])
			})

	# From sensitivitymv package with permutation t test.
	senmv(newyd, method='t')$deviate - dev	
  }
  

   # Different treatment effects.
   taus = c(0,  .5, .75, 1, 1.25, 1.5)
   
  res <- foreach(j=1:length(taus), .combine =c, .packages = c('sensitivitymult','sensitivitymv','senstrat', 'MASS', 'rmutil')) %do% {

      taus = c(0,  .5, .75, 1, 1.25, 1.5)
	   
	  taui = taus[j]

	  ## Simulate outcome: 	
	  mydata$y<-rnorm(n,taui*mydata$tr+mydata$x1+mydata$x2^2,1)


	  # Set range of tau where we search for the estimate 
	  #	and the confidence interval limits
	  lw <- max(mydata$y[mydata$treatment==1]) - min(mydata$y[mydata$treatment==0])
	  up <- -lw
      up = max(up, lw)
      lw = -up
	  alpha = 0.05	#level
	  
	  ## ANALYSIS for propensity score stratification
	  res_prop <- c()
	  
		# Estimate
	  est <- try(uniroot(ff.2, lower=lw, upper=up, z=mydata$treatment, 
				y = mydata$y, cl=mydata$subclass, x1 = mydata$x1, 
				x2 = mydata$x2, dev=0)$root, silent=TRUE)
	  
	  if(is.null(attr(est, "class"))){
			res_prop <- c(res_prop, est)
	  } else res_prop <- c(res_prop, NA)
	 
		# Lower confidence limit
	  estL <- try(uniroot(ff.2,lower = lw, upper = up, z=mydata$treatment, 
				y = mydata$y, cl=mydata$subclass, x1 = mydata$x1, x2 = mydata$x2,
				dev=(-stats::qnorm(alpha/2)))$root, silent=TRUE)
		# Upper confidence limit 
	  estH <- try(uniroot(ff.2,lower = lw, upper = up, z=mydata$treatment, 
				y = mydata$y, cl=mydata$subclass, x1 = mydata$x1, x2 = mydata$x2,
				dev=(stats::qnorm(alpha/2)))$root, silent=TRUE)
				
	  if(is.null(attr(estL, "class")) & is.null(attr(estH, "class"))) {
			res_prop <- c(res_prop, estL, estH)
	  } else res_prop <- c(res_prop, NA, NA)				
       
	  res_prop

	  ## ANALYSIS for leaf node balancing
	  res_leaf <- c()
      	  
		  # Estimate
	  est <- try(uniroot(ff.2, lower=lw, upper=up, z=mydata$treatment, 
				y = mydata$y, cl=mydata$subclass_tree, x1 = mydata$x1, 
				x2 = mydata$x2, dev=0)$root, silent=TRUE)
	  
	  if(is.null(attr(est, "class"))){
		res_leaf <- c(res_leaf, est)
	  } else res_leaf <- c(res_leaf, NA)
		
		# Lower confidence limit
      estL = try(uniroot(ff.2,lower = lw, upper = up, z=mydata$treatment, 
				y = mydata$y, cl=mydata$subclass_tree, x1 = mydata$x1, x2 = mydata$x2,
				dev=(-stats::qnorm(alpha/2)))$root, silent=TRUE)
		# Upper confidence limit
	  estH <- try(uniroot(ff.2,lower = lw, upper = up, z=mydata$treatment, 
				y = mydata$y, cl=mydata$subclass_tree, x1 = mydata$x1, x2 = mydata$x2,
				dev=(stats::qnorm(alpha/2)))$root, silent=TRUE)
				
	  if(is.null(attr(estL, "class")) & is.null(attr(estH, "class"))) {
		  res_leaf <- c(res_leaf, estL, estH)
	  } else 	res_leaf <- c(res_leaf, NA, NA)

	  res_leaf

	  ## ANALYSIS for mahalanobis matching
	  d5_match$y = mydata$y

	  res_pair <- c()
	  	  # Estimate
	  est <- try(uniroot(ff1.2, lower=lw, upper=up, z=d5_match$treatment, 
				y = d5_match$y, cl=d5_match$gp, x1 = d5_match$x1, x2 = d5_match$x2, 
				dev=0)$root, silent=TRUE)
				
	  if(is.null(attr(est, "class"))){
		res_pair <- c(res_pair, est)
	  } else res_pair <- c(res_pair, NA)
			  
		# Lower confidence limit  
      estL <- try(uniroot(ff1.2,lower = lw, upper = up, z=d5_match$treatment, 
				y = d5_match$y, cl=d5_match$gp, x1 = d5_match$x1, x2 = d5_match$x2,
				dev=(-stats::qnorm(alpha/2)))$root, silent=TRUE)
		# Upper confidence limit
	  estH <- try(uniroot(ff1.2,lower = lw, upper = up, z=d5_match$treatment, 
				y = d5_match$y, cl=d5_match$gp, x1 = d5_match$x1, x2 = d5_match$x2,
				dev=(stats::qnorm(alpha/2)))$root, silent=TRUE)
	  
	  if(is.null(attr(estL, "class")) & is.null(attr(estH, "class"))) {
			res_pair <- c(res_pair, estL, estH)
	  } else res_pair <- c(res_pair, NA, NA)
	  
	  res_pair
 
	  ## ANALYSIS for tree leaf matching
	  mydata_treematch$y = mydata$y
	  mydata_treematch$treatment <- mydata$treatment

	  res_treepair <- c()

	  if(sum(mydata$subclass_tree %in% nonpurenodes) < n*.05){
		# Don't attempt inference if too few units
		res_treepair <- c(res_treepair, NA)
	  } else {
			# Estimate
	    est <- try(uniroot(ff1.2, lower=lw, upper=up, z=mydata_treematch$treatment, 
				y = mydata_treematch$y, cl=mydata_treematch$gp, x1 = mydata_treematch$x1, x2 = mydata_treematch$x2, dev=0)$root, silent=TRUE) 
				
		if(is.null(attr(est, "class"))){		
			res_treepair <- c(res_treepair, est)
		} else   res_treepair <- c(res_treepair, NA)
	  }
	  
      if(sum(mydata$subclass_tree %in% nonpurenodes) < n*.05){
		res_treepair <- c(res_treepair, NA, NA)
	  } else {
			# Lower confidence limit
	     estL <- try(uniroot(ff1.2,lower = lw, upper = up, z=mydata_treematch$treatment, 
				y = mydata_treematch$y, cl=mydata_treematch$gp, x1 = mydata_treematch$x1, x2 = mydata_treematch$x2, dev=(-stats::qnorm(alpha/2)))$root, silent=TRUE)
			# Upper confidence limit
		 estH <- try(uniroot(ff1.2,lower = lw, upper = up, z=mydata_treematch$treatment, 
				y = mydata_treematch$y, cl=mydata_treematch$gp, x1 = mydata_treematch$x1, x2 = mydata_treematch$x2,
				dev=(stats::qnorm(alpha/2)))$root, silent=TRUE)
		 
		 if(is.null(attr(estL, "class")) & is.null(attr(estH, "class"))) {
			res_treepair <- c(res_treepair, estL, estH)
		  } else res_treepair <- c(res_treepair, NA, NA)
	  }
	  res_treepair 

	  c(res_prop, res_leaf, res_pair, res_treepair)

  }

  
  res_tab<-c(res, length(nonpurenodes), sum(mydata$subclass_tree %in% nonpurenodes))#)
  
  write.table(t(c(iter*900+iter2, res_tab)), 'sim_10k_normal_6.csv', row.names = FALSE, col.names = FALSE, append = TRUE, sep=',') 
}  


stopCluster(cl)

