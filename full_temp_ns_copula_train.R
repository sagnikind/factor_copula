data = read.csv(file ="/Users/mondals/Documents/Research/Non_Gaussian_Random_Field/model_5/data_try/Raphael_data/refl_Uscale_ave_train.csv",sep = ",",header = TRUE)

library("geosphere")
library("tlrmvnmvt")
library("pracma")
library("parallel")
library("mnormt")

coord = matrix(data = NA, nrow = 12,ncol = 2)
coord[,1] = as.vector(data$Lon) # longitude 
coord[,2] = as.vector(data$Lat) # latitude 

################################################################################################################
# dist_mat = distm(coord, fun = distHaversine)*10^(-3) # in KM
dist_mat = distm(coord, fun = distHaversine)*10^(-5) # in KM * 10^-2
loc = coord

U = as.matrix(t(as.matrix(data[,3:154])))
colnames(U)=NULL
rownames(U)=NULL

############################################
############################################
############################################


N = nrow(U) #number of replicates
n = ncol(U) #number of locations










invfunc = function(q,func,th,l.b=1e-12,u.b=1-1e-12,tol=1e-8){ 
  lq = length(q);
  lth= length(th); 
  v0 = rep(l.b,lq);
  v1 = rep(u.b,lq);
  out = rep(0, lq);
  iconv = rep(F, lq); 
  func0 = function(v,q,th){return(func(v,th)-q)};
  f0 = func0(v0,q,th); f1 = func0(v1,q,th);
  l0 = (f0 > 0 & f1 > 0); 
  l1 = (f0 < 0 & f1 < 0); 
  iconv[l0 | l1] = T; 
  tol0 = rep(1,lq); tol0[iconv] = 0; 
  while(prod(iconv) < 0.5) {  
    v0a = v0[!iconv]; v1a = v1[!iconv]; 
    v2a = (v0a+v1a)/2; qa = q[!iconv]; 
    f0a = f0[!iconv]; f1a = f1[!iconv];
    tha=th; if(lth > 1) tha=th[!iconv];
    f2a = func0(v2a,qa,tha);                      
    i1 = (f0a*f2a < 0); i2 = (f1a*f2a < 0);
    v1a[i1] = v2a[i1]; v1[!iconv] = v1a; 
    v0a[i2] = v2a[i2]; v0[!iconv] = v0a;
    
    f0[!iconv][i1] = f0a[i1];
    f0[!iconv][i2] = f2a[i2];
    f1[!iconv][i1] = f2a[i1];
    f1[!iconv][i2] = f1a[i2];
    
    tol0[!iconv] = abs(f2a); 
    iconv = (abs(tol0) < tol)     
  }
  
  out=(v0+v1)/2
  out[l0] = l.b; out[l1] = u.b
  out
}

nq = 151
#cc = gaussLegendre(nq, 0,300)
cc = gaussLaguerre(nq, a = 0)
ccw = cc$w
ccx = cc$x

f_0 = function(x){
  if(x<10^-8){
    return(0)
  }else{
    return(1/x)
  }
}
f_0 = Vectorize(f_0)

Inv = function(A){
  U = eigen(A)$vectors
  l = eigen(A)$values
  if(min(l)<10^-8){
    l_inv = f_0(l)
    return((U)%*%diag(l_inv)%*%t(U))
  }else{solve(A)}
}

llh = function(a){
  phi = exp(a[1]); nu = exp(a[2]) # Sigma_Z
  b = a[3]; a_1 = a[4]; a_2 = a[5] # alpha_0
  alpha = exp(a[6]) #alpha
  
  
  alpha_0 = vector(length = n)
  for(i in 1:n){
    alpha_0[i] = exp(b+a_1*0.01*loc[i,1]+a_2*0.01*loc[i,2])
  }
  
  if(phi > 5000) return(-1e100)
  if(phi < 0.01) return(-1e100)
  
  if(nu > 5) return(-1e100)
  if(nu < 0.01) return(-1e100)
  
  if(max(alpha,alpha_0) > 10) return(-1e100)
  if(min(alpha,alpha_0) < 0.01) return(-1e100)
  
  rho_Z = function(h){
    if(h==0){
      return(1)
    }else{
      temp_1 = 1/(gamma(nu)*2^(nu-1))
      temp_2 = h/phi
      temp_3 = besselK(temp_2, nu, expon.scaled = FALSE)
      return( temp_1 * temp_2^nu * temp_3)
    }
  } # correlation function
  rho_Z = Vectorize(rho_Z)
  Sigma_Z = matrix(nrow = n,ncol = n)
  for(j in 1:n){
    Sigma_Z[,j] = rho_Z(dist_mat[,j])
  }
  
  W = matrix(nrow = N,ncol = n)
  for(j in 1:n){
    Dist_F = function(z,th=0){   ####add th parameter to make it compartible with invfunc()
      temp_1_1 = log(alpha) + ((0.5/alpha^2) - (z/alpha)) + pnorm(z - (1/alpha),log.p = TRUE)
      term_1 = exp(temp_1_1)/(alpha - alpha_0[j])
      
      temp_1_2 = log(alpha_0[j]) + ((0.5/alpha_0[j]^2) - (z/alpha_0[j])) + pnorm(z - (1/alpha_0[j]),log.p = TRUE)
      term_2 = exp(temp_1_2)/(alpha - alpha_0[j])
      
      temp = pnorm(z) - term_1 + term_2
      return(temp)
    }
    
    W[,j] = invfunc(U[,j],Dist_F,0,l.b=-1000,u.b=1000,tol=1e-10)
  }
  
  den_f = function(z,j){
    temp_1 = exp(((0.5/alpha^2) - (z/alpha)) + pnorm((z- (1/alpha)),log.p = TRUE))/(alpha - alpha_0[j])
    temp_2 = exp(log(alpha) + ((0.5/alpha^2) - (z/alpha)) + dnorm((z- (1/alpha)),log = TRUE))/(alpha - alpha_0[j])
    temp_3 = exp(((0.5/alpha_0[j]^2) - (z/alpha_0[j])) + pnorm((z- (1/alpha_0[j])),log.p = TRUE))/(alpha - alpha_0[j])
    temp_4 = exp( log(alpha_0[j]) + ((0.5/alpha_0[j]^2) - (z/alpha_0[j])) + dnorm((z- (1/alpha_0[j])),log = TRUE))/(alpha - alpha_0[j])
    
    temp = dnorm(z) + temp_1 - temp_2 - temp_3 + temp_4
    return(temp)
  }
  
  H = diag(rep(1/alpha,n)) %*% Sigma_Z %*% diag(rep(1/alpha,n))
  Inv_Sigma_Z = Inv(Sigma_Z)
  Inv_H = diag(rep(alpha,n)) %*% Inv_Sigma_Z %*% diag(rep(alpha,n))
  
  llh_rep = function(i){
    W_v = matrix(data = NA,nrow = nq,ncol = n)
    for(k in 1:nq){
      W_v[k,] = W[i,] - alpha_0*ccx[k]
    }
    S_v = W_v%*%Inv_Sigma_Z
    T_v = alpha*S_v - matrix(1,nrow = nq,ncol = n)
    Y_v = T_v%*%H
    Sum = vector(length = nq)
    for(k in 1:nq){
      term_1 = 0.5*(as.numeric(t(Y_v[k,])%*%Inv_H%*%(Y_v[k,]))-as.numeric(t(W_v[k,])%*%Inv_Sigma_Z%*%(W_v[k,])))
      term_2 = -n*log(alpha)
      term_3 = log(abs(pmnorm(rep(0,n),mean = -Y_v[k,], varcov = H)[1]))
      Sum[k] = ccw[k]*exp(term_1+term_2+term_3)
    }
    
    univ_log_den = vector(length = n)
    for(j_0 in 1:n){
      univ_log_den[j_0] = log(den_f(W[i,j_0],j_0)) 
    }
    return(log(sum(Sum))- sum(univ_log_den))
    #return(sum(Sum) - sum(univ_log_den))
  }
  llh_rep = Vectorize(llh_rep)
  return(sum(as.numeric(unlist(mclapply(seq(1:N),llh_rep,mc.cores = 15)))))
}

negative_llh = function(a){
  return(-llh(a))
}

ini_theta = c(5.5410899, -1.9298124, -0.7240428, 40.3519924, -3.8351316, -0.6038322)

a = Sys.time()
llh(ini_theta)
b = Sys.time()
b - a


H = optim(ini_theta,negative_llh,method = "BFGS")

phi_est = exp(H$par[1])
nu_est = exp(H$par[2])
b_est = H$par[3]
a_1_est = H$par[4]
a_2_est = H$par[5]
alpha_est = exp(H$par[6])