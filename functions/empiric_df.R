empiric_df<-function(data,x)
{  
  #-----------------------------------------------------------------------------#
  #   Transform data to standard uniform d.f. through its empirical d.f.        #
  #
  # < CDD for FMRI >
  #   
  #   Copyright (C) 2018 Namgil Lee & Jong-Min Kim
  #-----------------------------------------------------------------------------#
  
  data<-sort(data)
  
  if(min(data)>0) a<-0 else a<-floor(min(data)/100)*100
  if(max(data)<0) b<-0 else b<-ceiling(max(data)/100)*100
  
  for(j in 1:length(x))
  {
    if(x[j]<a) x[j]<-a
    if(x[j]>b) x[j]<-b
  }
  
  data<-c(a,data,b)
  n<-length(data)
  p<-c(rep(0,(n-1)))
  q<-c(rep(0,(n-1)))
  
  for(i in 2:(n-2))
  {
    p[i]<-(data[i]+data[i+1])/2
    q[i]<-(i-1)/(n-2)
  }
  p[1]<-a
  p[n-1]<-b
  q[1]<-0
  q[n-1]<-1
  approx(p,q,xout=c(x))$y
}