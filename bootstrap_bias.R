kerne = 10
library(foreach)
library(doParallel)
library(MASS)
library(RMTstat)
library(nlshrink)

sim <- function (n,p,d,B,l) { # d: realizations of X_n; B: size of each bootstrap sample; d_b: number of bootstrap samples
  cl <- makeCluster(kerne)
  registerDoParallel(cl)
  
  normal=TRUE
  lambda = c( l, rep(1,p-length(l)) )
  #polynomial decreasing eigenvalues 
  #j <- 1:(p-length(l))  
  #c <-  1   
  #lambda_bulk <- (j)^(-c)
  #lambda_bulk[1:150] = rep(1,150)
  #lambda = c(l, lambda_bulk)
  
  Sigma=diag(lambda,p,p)
  epsilon = 0.2
  
  #### find true xi
  f <- function (x) {
    return ( 1/p* sum( (lambda*x)^2 / (1 - lambda*x)^2 ) - n/(p)  )
  }
  xi = uniroot (f, lower=0, upper = 1/max(lambda), tol = 0.0001)$root
  
  ###Calculation of r 
  r = 1/xi * (1 + p/n*1/p*sum( ( ( lambda * xi ) / ( 1 - lambda *xi) ) ) )
  
  #Computation of test statistics 
  T = lambda[1] - r
  
  ###intialize for bootstrap statistics
  g <- function (x) {
    return ( 1/p* sum( (lambda_quest*x)^2 / (1 - lambda_quest*x)^2 ) - n/(p)  )
  }

  
  result <- foreach(k = 1:d, .packages = c("MASS", "RMTstat", "nlshrink" ), .combine = "rbind") %dopar% {
    # Generate data
    lambda_star = rep(0,B)
    if(normal) {
      X = mvrnorm(n, rep(0,p), Sigma)
    }
    else {
      generateRandomPDMatrix <- function(p,eigenvalues) { #returns matrix with sqrt(eigenvalues), like sqrtm(sigma) rather than Sigma 
        random_matrix <- matrix(rnorm(p^2), ncol = p)
        random_orthogonal_matrix <- qr.Q(qr(random_matrix))
         D <- diag(eigenvalues)
        pd_matrix <- random_orthogonal_matrix %*% sqrt(D) %*% t(random_orthogonal_matrix)
        return(pd_matrix)
      }
      Sigma_sqrt <- generateRandomPDMatrix(p,lambda) 
      
      ### t-distribution
      deg = 10 
      var_T = deg / (deg - 2) 
      X <- matrix(rt(n*p, df=deg)/sqrt(var_T), nrow = n)
      
      X = t(X)
      X <- Sigma_sqrt %*% X 
      X=t(X)
    }
    lambda_quest = tau_estimate(X, k = 0, method = "nlminb", control = list())
    lambda_quest = sort(lambda_quest, decreasing=TRUE)
    
    lambda_hat = eigen(1/n*X%*%t(X), only.values = TRUE)$values 
    xi_hat = -1/(n)*sum( 1/( lambda_hat[2:n] - lambda_hat[1]  ) )
    lambda_quest = pmin(lambda_quest, 1/(xi_hat*(1+epsilon) ) )
    
    for(j in 1:B) {
      X_star = mvrnorm(n, rep(0,p), diag(lambda_quest,p,p))
      lambda_hat_star = eigen(1/n*X_star%*%t(X_star), only.values = TRUE)$values 
      lambda_star[j] =  lambda_hat_star[1]
        }
    T_star_1 = lambda_quest[1] 
    T_star_2 = mean(lambda_star) 
    c( r,  T_star_1, T_star_2, lambda_hat[1])
  }
 
  stopCluster(cl)
  return(result)
} 
d = 600
B = 500
n=600
p=400
#l=list(1, 1.4,1.7,c(1.5,1.1),c(1.7,1.1)) #für n=400,p=600
#l=list(1, 1.2,1.4,c(1.3,1.1),c(1.4,1.1)) # für n=500, p=300
#l=list(1, 1.1, 1.2, 1.3, 1.4, c(1.4,1.2), c(1.4,1.1)) #für n=600, p=400
l=list(1, 1.1, 1.2,  c(1.2,1.1))
mean_bootstrap = rep(0,length(l))
mean = rep(0,length(l)) 
sd_bootstrap= rep(0,length(l))
for (i in 1:length(l)) {
hilf=sim(n,p,d,B,unlist(l[i]))
r= hilf[1,1] 
T_star_1 = hilf[,2]
T_star_2 = hilf[,3]
lambda_hat=hilf[,4]

mean_bootstrap[i] = mean(T_star_1 - T_star_2)
mean[i] = unlist(l[i])[1] - mean(lambda_hat)
sd_bootstrap[i] = sd(T_star_1 - T_star_2)
} 

rbind(mean_bootstrap, mean, sd_bootstrap)
 
