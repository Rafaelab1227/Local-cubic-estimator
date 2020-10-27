#The objective is to create a function to estimate the regression function m(x) non-parametrically
#considering a local cubic estimator.
lc <- function(x, data, h, K){
  
  #Arguments:
  #x ---> vector of evaluation points
  #y ---> data frame or list containing the sample data for X as the 
  #       predictor variable and Y the response variable
  #h ---> the bandwith
  #K ---> kernel function
  
  # We need to compute for each evaluation point of x the estimated m, then we can induce a local
  # parametrization of m using a p order Taylors expansion and finally, turn it into a 
  # a linear regression problem where the unknown parameters will be the betas, and which will need
  # to minimized the difference between Y for the x and the m estimated. Additionally, these needs to
  # be weighted based on the proximity of each sample of X of the evaluated point what will be given by
  # calculating kernels. The estimate for m(x) can be computed as a weighted least squares problem
  # solved in matrix form:
  reg_mat <- sapply(x, function(x) {
    X <- cbind(1, data$X - x, (data$X - x)^2, (data$X - x)^3) #Matrix of n rows or the number of observations
    #of x in the sample data and p+1 columns which correspond to the difference between each sample X and the
    #evaluation point elevated to a power p.
    W <- diag(K((data$X - x) / h) / h) #Matrix of diagonal entries of the kernel density estimator.
    MASS::ginv(t(X)%*%W%*%X)%*%(t(X)%*%W%*%data$Y) #we could also obtained the coefficient Bo by introducing a linear 
    #weighted regression, but for ilustrative construction of the function we are computing manually
    #using matrix form.
    
  })
  
  t(reg_mat)[,1] #we extract the Bo, that will give the estimated for m(x)
}

# b) Test of the implementation
set.seed(12345) #set seed to make it replicable
m <- function(x)(x-1)^2 #set real function
n <- 500 #number of samples
X <- rnorm(n, mean=1, sd=1)
e <- rnorm(n, mean=0, sd = 0.5)

Y <- m(X) + e

data <- data.frame(X,Y)

x <- seq(-2, 4, by = 0.01)
h <- 0.5
K <- dnorm

mTrue <- m(x)# True regression
mEst <- lc(x= x, data= data, h=0.5, K=dnorm) #our function

#We are comparing the real regression with the local cubic estimator and the 
#version Locpoly function of the package KernSmooth
plot(data$X, data$Y, xlab = "x", ylab = "y", xlim = c(-2, 4)) #plot sample
lines(x, mTrue, lwd = 3, col = 3)#real m of the evaluation points
lines(x, mEst, lwd = 3,col = 4) #estimated m of the evaluation points
# Prediction by KernSmooth::locpoly
lp <- KernSmooth::locpoly(x = data$X, y= data$Y, bandwidth = h, degree = 3, range.x = c(-2,4),gridsize = 601L)
lines(x, lp$y, lwd = 1,col = 5) #estimated m of the evaluation points
legend("top", legend=c("True regression", "Local cubic estimator","Locpoly estimator"), lwd=2, col=c(3,4,5))
