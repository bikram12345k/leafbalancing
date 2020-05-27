### This function is called by other codes to compute the
#	standardized mean difference.
#	data1 <- the dataset
#	covrt <- the name of the covariate in the data for which the 
#			SMD will be calculated
#	cls <- name of the column in the data that contain the strata identifiers
#	trt <- name of the treatment variable.

match_summv1<-function(data1,covrt,cls,trt){
  
  data1<-data1
  covariate<-covrt
  
  colnum<-which(colnames(data1)==covariate)
  clsn<-which(colnames(data1)==cls)
  trtn<-which(colnames(data1)==trt)
  
  
  no.subclasses<-length(table(data1[,clsn]))
  
  sct<-unique(data1[,clsn])
  
  if(!(covrt %in% colnames(data1) && cls %in% colnames(data1) && trt %in% colnames(data1)))
  {
    stop("Coveriate, treatment and subclass variables should be part of dataset")
  }
  
  no.subclasses<-length(unique(data1[,clsn]))
  
  
  
  row_names1<-c(1:no.subclasses)
  col_names1<-c("classno","nct","mean","sd")
  matrix_names1<-c("control","treatment")
  
  tempdf<-array(rep(NA,no.subclasses*2*4),dim=c(no.subclasses,4,2),
                dimnames = list(row_names1,col_names1,matrix_names1))
  
  for(i in 1:no.subclasses){  
    tempdf[i,1,1]<-i  #class no
    tempdf[i,1,2]<-i
    tempdf[i,2,1]<-nrow(data1[data1[,trtn]==0 & data1[,clsn]==sct[i],]) # total observations
    tempdf[i,2,2]<-nrow(data1[data1[,trtn]==1 & data1[,clsn]==sct[i],])
    tempdf[i,3,1]<-mean(data1[data1[,trtn]==0 & data1[,clsn]==sct[i],colnum],na.rm=TRUE) # mean value
    tempdf[i,3,2]<-mean(data1[data1[,trtn]==1 & data1[,clsn]==sct[i],colnum],na.rm=TRUE)
    tempdf[i,4,1]<-sd(data1[data1[,trtn]==0 & data1[,clsn]==sct[i],colnum],na.rm=TRUE) # sd value
    tempdf[i,4,2]<-sd(data1[data1[,trtn]==1 & data1[,clsn]==sct[i],colnum],na.rm=TRUE)
    
  }
  
  overall_diff_full<-mean(data1[data1[,trtn]==1,colnum],na.rm=TRUE)-
    mean(data1[data1[,trtn]==0,colnum],na.rm=TRUE)
  diff_sub<-tempdf[,3,2]-tempdf[,3,1]
  size<-tempdf[,2,1]+tempdf[,2,2]
  nodrop<-c(tempdf[,2,1]>0)&c(tempdf[,2,2]>0) #rep(NA,no.subclasses)
  
  
  vart<-var(data1[data1[,trtn]==1,colnum],na.rm=TRUE)
  varc<-var(data1[data1[,trtn]==0,colnum],na.rm=TRUE)
  overall_diff_sub<-sum(diff_sub[nodrop]*size[nodrop],na.rm=TRUE)/sum(tempdf[which(nodrop==TRUE),2,1:2])
  overall_diff_substd<-overall_diff_sub/sqrt((vart+varc)/2)
  
  return(list(a=overall_diff_sub,b=overall_diff_substd,c=sum(tempdf[which(nodrop==TRUE),2,1:2])/nrow(data1)))
  
}
