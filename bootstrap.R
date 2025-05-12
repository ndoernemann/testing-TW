kerne = 10
library(foreach)
library(doParallel)
library(MASS)
library(RMTstat)
library(nlshrink)

sim <- function (n,p,d,B) { # d: realizations of X_n; B: size of each bootstrap sample; d_b: number of bootstrap samples
  cl <- makeCluster(kerne)
  registerDoParallel(cl)
  
  normal=TRUE
  l = rep(1.3,5) 
  lambda = c( l, rep(1,p-length(l)) )
  Sigma=diag(lambda,p,p)
  epsilon = 0.2
  
  #### find true xi
  f <- function (x) {
    return ( 1/p* sum( (lambda*x)^2 / (1 - lambda*x)^2 ) - n/(p)  )
    }
  xi = uniroot (f, lower=0, upper = 1/max(lambda), tol = 0.0001)$root
  sigma_3 = 1/xi^3 * (1 + p/n*1/p*sum( ( ( lambda * xi ) / ( 1 - lambda *xi) )^3 ) )
  sigma = abs(sigma_3)^(1/3) 
  
  ###Calculation of r 
  r = 1/xi * (1 + p/n*1/p*sum( ( ( lambda * xi ) / ( 1 - lambda *xi) ) ) )
  
  result <- foreach(k = 1:d, .packages = c("MASS", "RMTstat", "nlshrink", "expm"), .combine = "cbind") %dopar% {
    # Generate data
    print(k)
    T = rep(0,B)
    T_gap=rep(0,B)
    T_star = rep(0,B)
    T_star_gap = rep(0,B)
    if(normal) {
      X = mvrnorm(n, rep(0,p), Sigma)
    }
    else {
      generateRandomPDMatrix <- function(p,eigenvalues) {
        random_matrix <- matrix(rnorm(p^2), ncol = p)
        random_orthogonal_matrix <- qr.Q(qr(random_matrix))
       D <- diag(eigenvalues)
        pd_matrix <- random_orthogonal_matrix %*% D %*% t(random_orthogonal_matrix)
        return(pd_matrix)
      }
      Sigma <- generateRandomPDMatrix(p,lambda) 
      
      deg = 10
      var_T = deg / (deg - 2) 
      X <- matrix(rt(n*p, df=deg)/sqrt(var_T), nrow = n)
      X = t(X)
      X <- sqrtm(Sigma) %*% X 
      X=t(X)
    }
    lambda_quest = tau_estimate(X, k = 0, method = "nlminb", control = list())
      
    lambda_hat = eigen(1/n*X%*%t(X), only.values = TRUE)$values 
    xi_hat = -1/(n)*sum( 1/( lambda_hat[2:n] - lambda_hat[1]  ) )
    lambda_quest = pmin(lambda_quest, 1/(xi_hat*(1+epsilon) ) )
   
    #Computation of test statistics 
    T[1] = (n)^(2/3)/sigma * ( lambda_hat[1] - r)
    T_gap[1] =  (n)^(2/3)/sigma * ( lambda_hat[1] - lambda_hat[2])
    
    ###jetzt bootstrap statistics
      g <- function (x) {
        return ( 1/p* sum( (lambda_quest*x)^2 / (1 - lambda_quest*x)^2 ) - n/(p)  )
      }
      xi_hat_Q = uniroot (g, lower=0, upper = 1/max(lambda_quest), tol = 0.0001)$root
      sigma_3_Q = 1/xi_hat_Q^3 * (1 + p/n*1/p*sum( ( ( lambda_quest * xi_hat_Q ) / ( 1 - lambda_quest *xi_hat_Q) )^3 ) )
      sigma_Q = abs(sigma_3_Q)^(1/3)
      r_hat_Q = 1/xi_hat_Q * (1 + p/n*1/p*sum( ( ( lambda_quest * xi_hat_Q ) / ( 1 - lambda_quest *xi_hat_Q) ) ) )
      for(j in 1:B) {
      X_star = mvrnorm(n, rep(0,p), diag(lambda_quest,p,p))
      lambda_hat_star = eigen(1/n*X_star%*%t(X_star), only.values = TRUE)$values 
      T_star[j] = (n)^(2/3)/sigma_Q * ( lambda_hat_star[1] - r_hat_Q)
      T_star_gap[j] = (n)^(2/3)/sigma_Q * ( lambda_hat_star[1] - lambda_hat_star[2])
      if (j >= 2) { ### more simulations for non-bootstrap world
        X = mvrnorm(n, rep(0,p), Sigma)
        #Computation of test statistics 
        T[j] = (n)^(2/3)/sigma * ( lambda_hat[1] - r)
        T_gap[j] =  (n)^(2/3)/sigma * ( lambda_hat[1] - lambda_hat[2])
      }
      }
      c(T, T_gap, T_star, T_star_gap)
  }
  return(result)
} 
d = 600
B = 500
hilf=sim(500,300,d,B)
 T= hilf[1:B,] 
 T_gap = hilf[(B+1):(2*B),]
 T_star = hilf[(2*B+1):(3*B),] 
 T_star_gap = hilf[(3*B+1):(4*B),]

 mean(T)
 mean_star <- apply(T_star, MARGIN = 2, mean)
 mean(mean_star)
 sd(mean_star)
 
 sd(T)
 sd_star <- apply(T_star, MARGIN = 2, sd)
 mean(sd_star)
 sd(sd_star)
 
 quantile(T,0.95)
 quantile_star <- apply(T_star, MARGIN = 2, quantile, 0.95 )
 quantile_hat= mean(quantile_star)
 quantile_hat
 sd(quantile_star)
 sum(quantile_hat > T)/length(T)
 
 
 mean(T_gap)
 mean_star_gap <- apply(T_star_gap, MARGIN = 2, mean)
 mean(mean_star_gap)
 sd(mean_star_gap) 
 
 sd(T_gap)
 sd_star_gap <- apply(T_star_gap, MARGIN = 2, sd)
 mean(sd_star_gap)
 sd(sd_star_gap)
 
 quantile(T_gap,0.95)
 quantile_star_gap <- apply(T_star_gap, MARGIN = 2, quantile, 0.95 )
 quantile_gap_hat = mean(quantile_star_gap)
 quantile_gap_hat
 sd(quantile_star_gap)
 sum(quantile_gap_hat > T_gap)/length(T_gap)
 
 # output 
 result1 <- rbind(
   round(mean(T), 2),
   round(mean(mean_star), 2),
   round(sd(mean_star), 2),
   round(sd(T), 2),
   round(mean(sd_star), 2),
   round(sd(sd_star), 2),
   round(sum(quantile_hat > T) / length(T), 2)
 )
 
 result2 <- rbind(
   round(mean(T_gap), 2),
   round(mean(mean_star_gap), 2),
   round(sd(mean_star_gap), 2),
   round(sd(T_gap), 2),
   round(mean(sd_star_gap), 2),
   round(sd(sd_star_gap), 2),
   round(sum(quantile_gap_hat > T_gap) / length(T_gap), 2)
 )
 
 print(result1)
 print(result2)