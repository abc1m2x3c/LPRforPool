#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Gaussian kernel 
// [[Rcpp::export]]
vec GKernel(vec x){
  vec res(x.size());
  res=exp(-pow(x,2)/2);
  return(res);
}

// Ep kernel 
// [[Rcpp::export]]
vec EKernel(vec x){
  vec res(x.size());
  res=(1-pow(x,2))*3/4;
  res.elem(find(pow(x,2)>1)).zeros();
  return(res);
}

// local linear estimator
// Arguments:
// h = bandwidth
// X
// Y
// nx
// kernel: 0 Gaussian; 1 Ep kernel
// Output:
// Local linear fix at nx

// [[Rcpp::export]]
vec Llinear(double h, vec X, vec Y, vec nx , int kernel){
  vec res(nx.size());
  vec Xt(X.size());
  vec Wt(X.size());
  double a,b,d,e,f;
  for(int k=0; k<nx.size(); k++){
    double t=nx[k];
    Xt=X-t;
    if(kernel==0) //Gaussian kernel
    {
      Wt=GKernel(Xt/h);
    }
    if(kernel==1) //Ep kernel
    {
      Wt=EKernel(Xt/h);
    }
    a=sum(Wt);
    b=sum(Xt%Wt);
    d=sum(Xt%Xt%Wt);
    e=sum(Wt%Y);
    f=sum(Wt%Xt%Y);
    if(a*d-b*b==0){
      res[k]=10^6;
    }else{
      res[k]=(d*e-b*f)/(a*d-b*b);
    }
    //res[k]=(a*f-b*e)/g;
  }
  return(res);
}

// [[Rcpp::export]]
vec dele(vec x, int k){
  x.shed_row(k-1);
  return(x);
}

// cv local linear estimator
// Arguments:
// h = bandwidth
// X
// Y
// nx
// kernel: 0 Gaussian; 1 Ep kernel
// Output:
// Local linear fix at nx

// [[Rcpp::export]]
vec Llinear_cv(double h, vec X, vec Y, IntegerVector ind, int kernel){
  vec res(ind.size());
  vec Xt(X.size()-1);
  vec Wt(X.size()-1);
  vec Yt(Y.size()-1);
  double a,b,d,e,f;
  for(int k=0; k<ind.size(); k++){
    double t=X[ind[k]-1];
    Xt=dele(X-t,ind[k]);
    Yt=dele(Y,ind[k]);
    if(kernel==0) //Gaussian kernel
    {
      Wt=GKernel(Xt/h);
    }
    if(kernel==1) //Ep kernel
    {
      Wt=EKernel(Xt/h);
    }
    a=sum(Wt);
    b=sum(Xt%Wt);
    d=sum(Xt%Xt%Wt);
    e=sum(Wt%Yt);
    f=sum(Wt%Xt%Yt);
    if(a*d-b*b==0){
      res[k]=10^6;
    }else{
      res[k]=(d*e-b*f)/(a*d-b*b);
    }
    //res[k]=(a*f-b*e)/g;
  }
  return(res);
}



// Gaussian kernel 

// [[Rcpp::export]]
mat GGKernel(mat x){
  mat res(1,1,fill::zeros);
  res=exp(-pow(x,2)/2);
  return(res);
}

// Ep kernel 

// [[Rcpp::export]]
mat EEKernel(mat x){
  mat res(1,1,fill::zeros);
  //res=exp(-pow(x,2)/2);
  res=(1-pow(x,2))*3/4;
  res.elem(find(pow(x,2)>1)).zeros();
  return(res);
}

// // [[Rcpp::export]]
// vec test(mat x){
//   vec res(x.n_rows);
//   res=mean(x,1);
//   return(res);
// }

// local linear estimator based sum of kernels
// Arguments:
// h = bandwidth
// gX = matrix, the jth row is the X of the jth group: J times c
// gY = vec, group respons
// nx
// kernel: 0 Gaussian; 1 Ep kernel
// pool: 0 sum of kernels, 1 product of kernels
// Output:
// Local linear fix using sum of kernels at nx

// [[Rcpp::export]]
vec Llinear_pool(double h, mat gX, vec gY, vec nx , int kernel, int pool){
  vec res(nx.size());
  vec gXt(gY.size());
  vec gWt(gY.size());
  double a,b,d,e,f;
  for(int k=0; k<nx.size(); k++){
    double t=nx[k];
    gXt=mean(gX-t,1);
    
    if(kernel==0) //Gaussian kernel
    {
      if(pool==0){
        gWt=mean(GGKernel((gX-t)/h),1);
      }
      if(pool==1){
        gWt=prod(GGKernel((gX-t)/h),1);
      }
    }
    if(kernel==1) //Ep kernel
    {
      if(pool==0){
        gWt=mean(EEKernel((gX-t)/h),1);
      }
      if(pool==1){
        gWt=prod(EEKernel((gX-t)/h),1);
      }
    }
    a=sum(gWt);
    b=sum(gXt%gWt);
    d=sum(gXt%gXt%gWt);
    e=sum(gWt%gY);
    f=sum(gWt%gXt%gY);
    if(a*d-b*b==0){
      res[k]=10^6;
    }else{
    res[k]=(d*e-b*f)/(a*d-b*b);
    }
    //res[k]=(a*f-b*e)/g;
  }
  return(res);
}


// [[Rcpp::export]]
vec Lconstant_pool(double h, mat gX, vec gY, vec nx , int kernel, int pool){
    vec res(nx.size());
    vec gXt(gY.size());
    vec gWt(gY.size());
    double a,e;
    for(int k=0; k<nx.size(); k++){
        double t=nx[k];
        gXt=mean(gX-t,1);
        
        if(kernel==0) //Gaussian kernel
        {
            if(pool==0){
                gWt=mean(GGKernel((gX-t)/h),1);
            }
            if(pool==1){
                gWt=prod(GGKernel((gX-t)/h),1);
            }
        }
        if(kernel==1) //Ep kernel
        {
            if(pool==0){
                gWt=mean(EEKernel((gX-t)/h),1);
            }
            if(pool==1){
                gWt=prod(EEKernel((gX-t)/h),1);
            }
        }
        a=sum(gWt);
        e=sum(gWt%gY);
        if(a==0){
            res[k]=10^6;
        }else{
            res[k]=e/a;
        }
        //res[k]=(a*f-b*e)/g;
    }
    return(res);
}

// // [[Rcpp::export]]
// mat deleRow(mat x, int k1, int k2){
//   x.shed_rows(k1-1,k2-1);
//   return(x);
// }

// [[Rcpp::export]]
mat Llinear_pool_cv(double h, mat gX, vec gY, IntegerVector indJ , int kernel, int pool){
  int c=gX.n_cols;
  mat res(indJ.size(), c);
  vec gXt(gY.size());
  vec gWt(gY.size());
  vec gYt(gY.size()-1);
  double a,b,d,e,f;
  for(int j=0; j<indJ.size();j++){
    for(int i=0;i<c;i++){
      double t=gX(j,i);
      
      gXt=mean(gX-t,1);
      gXt=dele(gXt,j+1);
      if(kernel==0) //Gaussian kernel
      {
        if(pool==0){
          gWt=mean(GGKernel((gX-t)/h),1);
        }
        if(pool==1){
          gWt=prod(GGKernel((gX-t)/h),1);
        }
      }
      if(kernel==1) //Ep kernel
      {
        if(pool==0){
          gWt=mean(EEKernel((gX-t)/h),1);
        }
        if(pool==1){
          gWt=prod(EEKernel((gX-t)/h),1);
        }
      }
      gWt=dele(gWt,j+1);
      gYt=dele(gY,j+1);
      a=sum(gWt);
      b=sum(gXt%gWt);
      d=sum(gXt%gXt%gWt);
      e=sum(gWt%gYt);
      f=sum(gWt%gXt%gYt);
      if(a*d-b*b==0){
        res(j,i)=10^6;
      }else{
      res(j,i)=(d*e-b*f)/(a*d-b*b);
      }
    }    
  }
  return(res);
}



// [[Rcpp::export]]
mat VCM_cpp(mat X, vec U, vec Y, double h, vec nx, vec W_sample ,int kernel){
  int p = X.n_cols;
  mat res(2*p, nx.size());
  vec Uu(U.size());
  vec W_u(U.size());
  mat Uu_X(size(X)); 
  mat Gamma_u(size(X)[0],2*size(X)[1]);
  Gamma_u.cols(0,X.n_cols-1)=X;
  mat W_u_Gamma_u(size(Gamma_u));
  for(int k=0; k<nx.size(); k++){
    double u=nx[k];
    Uu=U-u;
    if(kernel==1) //Gaussian kernel
    {
      W_u=1/h*GKernel(Uu/h);
    }
    if(kernel==2) //Ep kernel
    {
      W_u=1/h*EKernel(Uu/h);
    }
    W_u = W_u%W_sample;
    Uu_X=X;
    Uu_X.each_col() %= Uu;  
    Gamma_u.cols(X.n_cols,2*X.n_cols-1)=Uu_X;
    // Gamma_u=join_cols(X,Uu_X);
    W_u_Gamma_u=Gamma_u;
    W_u_Gamma_u.each_col() %= W_u;
    // res.col(k) = (Gamma_u.t()*W_u_Gamma_u).i()*(W_u_Gamma_u.t()*Y);
    res.col(k) = solve(Gamma_u.t()*W_u_Gamma_u,W_u_Gamma_u.t()*Y);
  }
  return (res);
}



// [[Rcpp::export]]
mat deleRows(mat x, int k1, int k2){
  x.shed_rows(k1-1,k2-1);
  return(x);
}



// [[Rcpp::export]]
double VCM_cv_cpp(mat X, vec U, vec Y, double h, vec W_sample, vec W_cv,IntegerVector start_ind, IntegerVector end_ind, int kernel){
  int p = X.n_cols;
  vec res(start_ind.size());
  for (int index=0; index<start_ind.size(); index++){
    vec u_grid(end_ind(index)-start_ind(index)+1);
    u_grid=U.rows(start_ind[index]-1,end_ind[index]-1);
    mat X_rest=deleRows(X,start_ind[index],end_ind[index]);
    vec U_rest=deleRows(U,start_ind[index],end_ind[index]);
    vec Y_rest=deleRows(Y,start_ind[index],end_ind[index]);
    vec W_sample_rest=deleRows(W_sample,start_ind[index],end_ind[index]);
    mat beta = VCM_cpp(X_rest, U_rest, Y_rest, h, u_grid, W_sample_rest, kernel);
    beta.shed_rows(p,2*p-1);
    mat X_keep=X.rows(start_ind[index]-1,end_ind[index]-1);
    vec Y_keep=Y.rows(start_ind[index]-1,end_ind[index]-1);
    vec W_cv_keep=W_cv.rows(start_ind[index]-1,end_ind[index]-1);
    res[index]=sum(W_cv_keep%pow(Y_keep-sum(X_keep%beta.t(),1),2));
  }
  return (mean(res));
}


//this is a special version to acommadate reviewer's suggestion of removing correlation when calculating the overall mean (Wang, Mou, Li, Huang's paper)
//add two extra input in the end: vec Z of total number of individuals,  vec Z_concise, W_dotj_time_Z_concise of length number of pools
// [[Rcpp::export]]
double VCM_cv_cpp_RP(mat X, vec U, vec Y, double h,vec W_sample, vec W_cv,IntegerVector start_ind, IntegerVector end_ind, int kernel, vec W_dotj_time_Z_concise, vec pool_size){
  int p = X.n_cols;
  vec res(start_ind.size());
  for (int index=0; index<start_ind.size(); index++){
    vec u_grid(end_ind(index)-start_ind(index)+1);
    u_grid=U.rows(start_ind[index]-1,end_ind[index]-1);
    mat X_rest=deleRows(X,start_ind[index],end_ind[index]);
    vec U_rest=deleRows(U,start_ind[index],end_ind[index]);
    //new from 2020/05/24
    vec W_dotj_time_Z_concise_rest=deleRows(W_dotj_time_Z_concise,index+1,index+1);
    vec pool_size_rest=deleRows(pool_size,index+1,index+1);
    vec Z_concise_rest(pool_size_rest.size());
    for (int i=0; i<Z_concise_rest.size(); i++){
      vec W_dotj_time_Z_concise_rest_rest=deleRows(W_dotj_time_Z_concise_rest,i+1,i+1);
      Z_concise_rest[i]=W_dotj_time_Z_concise_rest[i]-(pool_size_rest[i]-1)*sum(W_dotj_time_Z_concise_rest_rest)/(U_rest.size()-pool_size_rest[i]);
    }
    vec Z_rest;
    for (int i=0; i<pool_size_rest.size(); i++){
      vec temp(pool_size_rest[i]);
      temp.fill(Z_concise_rest[i]);
      Z_rest=join_cols(Z_rest,temp);
    }
    vec Y_rest=deleRows(Y,start_ind[index],end_ind[index]);
    vec W_sample_rest=deleRows(W_sample,start_ind[index],end_ind[index]);
    mat beta = VCM_cpp(X_rest, U_rest, Z_rest, h, u_grid, W_sample_rest, kernel);
    beta.shed_rows(p,2*p-1);
    mat X_keep=X.rows(start_ind[index]-1,end_ind[index]-1);
    vec Y_keep=Y.rows(start_ind[index]-1,end_ind[index]-1);
    vec W_cv_keep=W_cv.rows(start_ind[index]-1,end_ind[index]-1);
    res[index]=sum(W_cv_keep%pow(Y_keep-sum(X_keep%beta.t(),1),2));
  }
  return (mean(res));
}



