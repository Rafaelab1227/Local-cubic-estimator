# Local-cubic-estimator
Function to generate a local cubic estimator

The objective was to create a function to estimate the regression function m(x) non-parametricall, considering a local cubic estimator.

lc(x, data, h, K)
  
Where:
  x ---> vector of evaluation points
  y ---> data frame or list containing the sample data for X as the 
         predictor variable and Y the response variable
  h ---> the bandwith
  K ---> kernel function
