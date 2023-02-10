Fit.it=function(tem.X,tem.y,tem.h,tem.S,tem.kernel)
{
  res=Llinear(tem.h, tem.X, tem.y, tem.S,tem.kernel)
  return(res)
}

CV.it=function(tem.X,tem.Y,tem.ind,tem.kernel,tem.interval)
{
  CV.object=function(tem.h)
  {
    tem.fit=Llinear_cv(tem.h, tem.X, tem.Y, tem.ind, tem.kernel)
    tem.resi=tem.Y[tem.ind]-tem.fit
    return(mean(tem.resi^2))
  }
  res=optimize(CV.object,interval=tem.interval)
  return(res$minimum)
}

# Fit.mi=function(tem.X,tem.n,tem.gY,tem.h,tem.S,tem.kernel)
# {
#   tem.N=length(tem.X)
#   Fake.Y=c(rep(tem.gY[1:(length(tem.gY)-ifelse(tem.N%%tem.n,1,0))],each=tem.n), 
#            rep(tem.gY[length(tem.gY)],tem.N%%tem.n))
#   Mu=mean(Fake.Y)
#   tem.gn=matrix(c(rep(tem.n,ceiling(tem.N/tem.n)-1),ifelse(tem.N%%tem.n, tem.N%%tem.n, tem.n)),tem.N,1)
#   run.Y=tem.gn*Fake.Y-(tem.gn-1)*Mu
#   res=Llinear(tem.h, tem.X, run.Y, tem.S,tem.kernel)
#   return(res)
# }
# 
# CV.mi=function(tem.X,tem.n,tem.gY, tem.ind, tem.kernel,tem.interval)
# {
#   tem.N=length(tem.X)
#   Fake.Y=c(rep(tem.gY[1:(length(tem.gY)-ifelse(tem.N%%tem.n,1,0))],each=tem.n), 
#            rep(tem.gY[length(tem.gY)],tem.N%%tem.n))
#   Mu=mean(Fake.Y)
#   tem.gn=matrix(c(rep(tem.n,ceiling(tem.N/tem.n)-1),ifelse(tem.N%%tem.n, tem.N%%tem.n, tem.n)),tem.N,1)
#   run.Y=tem.gn*Fake.Y-(tem.gn-1)*Mu
#   
#   CV.object=function(tem.h)
#   {
#     tem.fit=Llinear_cv(tem.h, tem.X, run.Y, tem.ind, tem.kernel)
#     tem.resi=run.Y[tem.ind]-tem.fit
#     return(mean(tem.resi^2))
#   }
#   res=optimize(CV.object,interval=tem.interval)
#   return(res$minimum)
# }


Fit.pool=function(tem.X,tem.n,tem.gy,tem.h,tem.S,tem.kernel,tem.pool)
{
  res=Llinear_pool(tem.h, matrix(tem.X,length(tem.X)/tem.n,tem.n,byrow=TRUE), tem.gy, tem.S,tem.kernel,tem.pool)
  return(res)
}

CV.pool=function(tem.X,tem.n,tem.gy,tem.indJ,tem.kernel,tem.pool,tem.interval)
{
  tem.gsize=length(tem.X)/tem.n
  tem.gX=matrix(tem.X,tem.gsize,tem.n,byrow=TRUE)
  CV.object=function(tem.h)
  {
    tem.fit=Llinear_pool_cv(tem.h,tem.gX,tem.gy,1:tem.gsize,tem.kernel,tem.pool)
    temp=rowMeans(tem.fit)
    tem.resi=tem.gy[tem.indJ]-temp[tem.indJ]
    return(mean(tem.resi^2))
  }
  res=optimize(CV.object,interval=tem.interval)
  return(res$minimum)
}

#cv.mi
CV.mi.S1=function(X,c,groupy,tem.kernel=0,tem.interval=c(0.01,2)) #W is the sample weight
{
  
  n=length(X)
  start_ind=seq(1,n,by=c)
  end_ind=seq(c,n,by=c)
  
  # mu_star_hat=mean(groupy) #use to generate similar results as the old code - old method to calculate mu
  
  R1=rep(0,length(X))
  for (j in 1:(n/c)){
    index=start_ind[j]:end_ind[j]
    mu_star_hat=c*sum(groupy[-j])/(n-c)
    R1[index]=c*groupy[j]-(c-1)*mu_star_hat
  }
  
  
  W=rep(1,n)
  W_cv1=W^{-1}
  W_sample1=W^{-1}
  U=W*X
  X1=cbind(rep(1,n))
  
  #20200524 for revision comments
  Z_concise=groupy
  W_dotj_time_Z_concise=c*groupy
  pool_size=end_ind-start_ind+1
  CV.object.S1<-function(h) #objective function of cross-validation selection of h
  {
    #res=VCM_cv_cpp(X1,U,R1,h,W_sample1,W_cv1,start_ind,end_ind,kernel_type)
    #kernel_type equals tem.kernel+1#0:guasian, 1:EP -> 1:guasian, 2:EP
    res=VCM_cv_cpp_RP(X1,U,R1,h,W_sample1,W_cv1,start_ind,end_ind,tem.kernel+1,W_dotj_time_Z_concise,pool_size)
    return(res)
  }
  h1=optimize(CV.object.S1,interval=tem.interval)$minimum 
  return (h1)
}




Fit.mi.S1<-function(X,c,groupy,h.m.S1,nx,tem.kernel=0){
  
  n=length(X)
  start_ind=seq(1,n,by=c)
  end_ind=seq(c,n,by=c)
  
  # mu_star_hat=mean(groupy) #use to generate similar results as the old code - old method to calculate mu
  
  R1=rep(0,length(X))
  for (j in 1:(n/c)){
    index=start_ind[j]:end_ind[j]
    mu_star_hat=c*sum(groupy[-j])/(n-c)
    R1[index]=c*groupy[j]-(c-1)*mu_star_hat
  }
  
  
  W=rep(1,n)
  W_cv1=W^{-1}
  W_sample1=W^{-1}
  U=W*X
  X1=cbind(rep(1,n))
  
  
  fit1=VCM_cpp(X1,U,R1,h.m.S1,nx,W_sample1,tem.kernel+1) #u_grid=U
  return (fit1[1,])
}



#cv.mi
CV.mi.S2=function(X,c,groupy,tem.kernel=0,tem.interval=c(0.01,2),fit1_U) #W is the sample weight
{
  
  n=length(X)
  start_ind=seq(1,n,by=c)
  end_ind=seq(c,n,by=c)
  
  R2=rep(0,n)
  for (j in 1:length(start_ind)){
    index=start_ind[j]:end_ind[j]
    for (i in 1:length(index)){
      R2[index[i]]=c*groupy[j]-sum(fit1_U[index[-i]])
    }
  }
  
  W=rep(1,n)
  W_cv2=W^{-1}
  W_sample2=W^{-1}
  U=W*X
  X2=cbind(rep(1,n))
  
  CV.object.S2<-function(h) #objective function of cross-validation selection of h
  {
    res=VCM_cv_cpp(X2,U,R2,h,W_sample2,W_cv2,start_ind,end_ind,tem.kernel+1)
    return(res)
  }
  h2=optimize(CV.object.S2,interval=tem.interval)$minimum 
  return (h2)
}



Fit.mi.S2<-function(X,c,groupy,h.m.S2,nx,tem.kernel=0,fit1_U){
  n=length(X)
  start_ind=seq(1,n,by=c)
  end_ind=seq(c,n,by=c)
  

  R2=rep(0,n)
  for (j in 1:length(start_ind)){
    index=start_ind[j]:end_ind[j]
    for (i in 1:length(index)){
      R2[index[i]]=c*groupy[j]-sum(fit1_U[index[-i]])
    }
  }
  
  W=rep(1,n)
  W_cv2=W^{-1}
  W_sample2=W^{-1}
  U=W*X
  X2=cbind(rep(1,n))
  
  fit2=VCM_cpp(X2,U,R2,h.m.S2,nx,W_sample2,tem.kernel+1) 
  
  return (fit2[1,])
}



