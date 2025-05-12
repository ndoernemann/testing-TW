library(foreach)
library(doParallel)
library(MASS)
library(RMTstat)
library(nlshrink)

kerne = 20
#alpha = 0.05 # significance level 
# see quantile_lim.R.out
quantile = 4.102893
quantile_onatski_1 = 6.981246 # one gap
quantile_onatski_10 = 19.9251  # ten gap 

psi_deriv <-function (alpha, y) { 
  # if LSD of H = delta_1  
  return (1-y/(alpha-1)^2)
}


sim <- function (n,p,d,l,normal,epsilon) {
  cl <- makeCluster(kerne)
  registerDoParallel(cl) #, outline="")
  y =p/n
  QUEST = 1 #use population eigenvalues (then 0) or quest estimator(then 1)? 
  
  lambda_bulk = rep(1, p-1)
  #polynomial decreasing eigenvalues 
  #j <- 1:(p-1)  # Sequence of integers from 1 to 10
  #c <-  1    # Exponent value
  #lambda_bulk <- (j)^(-c)
  #lambda_bulk[1:150] = rep(1,150)
  
  lambda = c(l, lambda_bulk)
  
  #### find true xi
  f <- function (x) {
    return ( 1/p* sum( (lambda*x)^2 / (1 - lambda*x)^2 ) - n/(p)  )
     }
  xi = uniroot (f, lower=0, upper = 1/max(lambda), tol = 0.0001)$root
  xil = l*xi 
  
  Sigma =diag(lambda)
  
  ### function for calculation of onatski's test statistic
  max_gap <- function (lambda, i) {
    gap = rep(0,i)
    for (k in 1:i){
      gap[k] =  ( lambda[k] - lambda[k+1] ) / ( lambda[k+1] - lambda[k+2] )
    }
    return (max(gap)) 
  }
  
  result <- foreach(k = 1:d, .packages = c("MASS", "RMTstat", "nlshrink"), .combine = "cbind") %dopar% {
    # Generate data
    print(k)
    T = rep(0,d) 
    if(normal) {
    X = mvrnorm(n, rep(0,p), Sigma)
    }
    else {
      deg = 10 # freiheitsgrade der t verteilung
      var_T = deg / (deg - 2) # varianz der t verteilung
      X <- matrix(rt(n*p, df=deg)/sqrt(var_T), nrow = n)
      X = t(X)
      X <- sqrt(Sigma) %*% X
      X=t(X)
    }
    #Estimation of population eigenvalues
    if (QUEST == 1)  
    {
      lambda_quest = tau_estimate(X, k = 0, method = "nlminb", control = list())
      lambda = lambda_quest 
    }
    lambda_hat = eigen(1/n*X%*%t(X), only.values = TRUE)$values 
  
    #Estimation of xi
    xi_hat = -1/(n)*sum( 1/( lambda_hat[2:n] - lambda_hat[1]  ) )
  
    ###Estimation of the variance 
    lambda = pmin(lambda, 1/( xi_hat * (1 + epsilon) ))
    sigma_3 = 1/xi_hat^3 * (1 + p/n*1/p*sum( ( ( lambda * xi_hat ) / ( 1 - lambda *xi_hat) )^3 ) )
    sigma_hat = abs(sigma_3)^(1/3) 
    
    #Computation of test statistics 
    T = (n)^(2/3)/sigma_hat * ( lambda_hat[1] - lambda_hat[2])
    T_onatski_1 = max_gap(lambda_hat,1)
    T_onatski_10 = max_gap(lambda_hat,10)
    c(T, sigma_hat , xi_hat, T_onatski_1, T_onatski_10) 
     }
  T = result[1,]
  sigma_hat= result[2,]
  xi_hat = result[3,]
  T_onatski_1 = result[4,]
  T_onatski_10 = result[5,]
  y = p/n
  
  indices <-  T>quantile
  count <- sum(indices) /d
  count2 = 0
  count3=0
  count_onatski_1 = sum (T_onatski_1 > quantile_onatski_1) /d
  count_onatski_10 = sum (T_onatski_10 > quantile_onatski_10) /d
  
  stopCluster(cl)
  return (c(count, count2, mean(T), var(T), count3, xil, mean(xi_hat), var(xi_hat), mean(sigma_hat), var(sigma_hat), count_onatski_1, count_onatski_10 ))
} 

n=600
p=400
d=500
l=c(1,1.25,1.3,1.5,1.7,2, 2.25, 2.5,2.75,3,3.5) 
normal = TRUE 
epsilon = 0.2
count = rep(0,length(l))
count2= rep(0,length(l))
mean = rep(0,length(l))
var = rep(0,length(l))
count3 = rep(0,length(l))
count_onatski_1 = rep(0,length(l))
count_onatski_10 = rep(0,length(l))
xil = rep(0,length(l))
mean_xi_hat =  rep(0,length(l))
var_xi_hat = rep(0, length(l))
mean_sigma_hat = rep(0, length(l))
var_sigma_hat = rep(0, length(l))
for (i in 1:length(l)) {
  hilf = sim(n,p,d,l[i],normal, epsilon)
  count[i] =hilf[1]
  count2[i] =hilf[2]
  mean[i] = hilf[3]
  var[i] = hilf[4]
  count3[i] = hilf[5]
  xil[i] = hilf[6]
  mean_xi_hat[i] = hilf[7]
  var_xi_hat[i] = hilf[8]
  mean_sigma_hat[i] = hilf[9]
  var_sigma_hat[i] = hilf[10]
  count_onatski_1[i] = hilf[11]
  count_onatski_10[i] = hilf[12]
}
xil
count
count_onatski_1
count_onatski_10
mean
var
mean_xi_hat
var_xi_hat
mean_sigma_hat
var_sigma_hat
#simulated mean for GOE: 2.068558 
# simulated var for GOE: var 1.260833 