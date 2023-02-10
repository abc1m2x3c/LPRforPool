packages <- c("Rcpp","RcppArmadillo")
install.packages(setdiff(packages, rownames(installed.packages())))  

rm(list=ls(all=TRUE))
library(maxLik)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("pool.cpp")
source("Pool.R")

########################################
# Data generation
########################################
set.seed(12345)
n = 600 # sample size
eta=function(x){
  return(x^3*exp(x^4/1000))
}
x=rbinom(n,1,0.8)
x1=(runif(n,0,1)-0.5)/0.0625
x1=sign(x1)*(sign(x1)*x1)^(1/3)
x2=runif(n,-1,1)
x=x*x1+(1-x)*x2
nx=seq(-1.9,1.9,by=0.1); nx = round(nx,2)
ystd=0.6
y=eta(x)+rnorm(n,0,ystd)
truth=eta(nx)



#input: 
# x: independent variable
# y: dependent variable
# tem.ind: indexes of individuals (for individual data) or pools (for pooled data) that are used in the bandwidth selection
# tem.kernel: kernel, 0: gaussian; 1: Ep kernel
# tem.interval: range of bandwidth for optimization
# tem.pool: 0 = average weighted estimator; 1 = product weight estimator
# nx: a vector at which the fitted value will be calculated
#output:
# function h.* returns the optimal bandwidth based on the proposed cross validatiion method
# function Fit.* returns the fitted value at the vector nx


########################################
# Fitting individual data
#########################################
h.it=CV.it(x,y,tem.ind=1:length(x),tem.kernel=0,tem.interval=c(0.01,2))
IT.res=Fit.it(x,y,h.it,nx,tem.kernel=0)
plot(x[x>=nx[1]&x<=nx[length(nx)]],y[x>=nx[1]&x<=nx[length(nx)]],ylab='y',xlab='x');lines(nx,truth,col="black",lwd=2.5);lines(nx,IT.res,col="yellow",lwd=2.5)
legend('topleft',legend=c('True','Individual est','average est','product est','marginal intergration est'),col=c('black','yellow','red','green','blue'),lty=1,lwd=2.5)


########################################
# Data generation for random pooling
########################################
c=2 # pool size
gsize = n/c # number of groups
groupy=colMeans(matrix(y,c,gsize)) # grouped dependent variable by pooling every c individual y's

########################################
# random pooling estimation
########################################
h.s=CV.pool(x,c,groupy,tem.ind=1:gsize,tem.kernel=0,tem.pool=0,tem.interval=c(0.01,2))
lps_hat=Fit.pool(x,c,groupy,h.s,nx,tem.kernel=0,tem.pool=0)
lines(nx,lps_hat,col="red",lwd=2.5)

# product weighted estimator 
h.p=CV.pool(x,c,groupy,1:gsize,tem.kernel=0,tem.pool=1,tem.interval=c(0.01,2))
lpp_hat=Fit.pool(x,c,groupy,h.p,nx,tem.kernel=0,tem.pool=1)
lines(nx,lpp_hat,col="green",lwd=2.5)

# marginal integration estimator; 
h.m.S1=CV.mi.S1(x,c,groupy,tem.kernel=0,tem.interval=c(0.01,2)) #W is the sample weight
lpm_hat=Fit.mi.S1(x,c,groupy,h.m.S1,nx,tem.kernel=0)
lines(nx,lpm_hat,col="blue",lwd=2.5)

########################################
# Data generation for homogeneous pooling
########################################
yh=y[order(x)]
xh=x[order(x)]
indJ=1:gsize
grouphy=rowMeans(matrix(yh,gsize,c,byrow=TRUE))

plot(x[x>=nx[1]&x<=nx[length(nx)]],y[x>=nx[1]&x<=nx[length(nx)]],ylab='y',xlab='x');lines(nx,truth,col="black",lwd=2.5);lines(nx,IT.res,col="yellow",lwd=2.5)
legend('topleft',legend=c('True','Individual est','average est','product est','marginal intergration est'),col=c('black','yellow','red','green','blue'),lty=1,lwd=2.5)

########################################
# homogeneous pooling estimation
########################################
h.s=CV.pool(xh,c,grouphy,indJ,tem.kernel=1,tem.pool=0,tem.interval=c(0.01,2))
lps_hat=Fit.pool(xh,c,grouphy,h.s,nx,tem.kernel=0,tem.pool=0)
lines(nx,lps_hat,col="red",lwd=2.5)

h.p=CV.pool(xh,c,grouphy,indJ,tem.kernel=1,tem.pool=1,tem.interval=c(0.01,2))
lpp_hat=Fit.pool(xh,c,grouphy,h.p,nx,tem.kernel=0,tem.pool=1)
lines(nx,lpp_hat,col="green",lwd=2.5)

h.m.S1=CV.mi.S1(xh,c,grouphy,tem.kernel=0,tem.interval=c(0.01,2)) #W is the sample weight
lpm_hat=Fit.mi.S1(xh,c,grouphy,h.m.S1,nx,tem.kernel=0)
lines(nx,lpm_hat,col="blue",lwd=2.5)

