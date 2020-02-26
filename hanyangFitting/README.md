# FIT RESULTS

## Method
The independent variables are x1 and int1. The dependant variable is z.  
x1 --> tns1 column in the input data  
int1 --> Intens1 column in the input data  
z --> Intens column in the input data  
The tns column in the input data was unused.  


The convolution result was chopped off once it reached the input vector length (taking the first 780 entries).  


## Results
```
General model:
     f(x1,int1) = y0+subsref(conv(a*exp(-1/tau*x1),int1),struct('type','()','subs',{{1:length(x1),1}}))
Coefficients (with 95% confidence bounds):
       a =       2.043  (2.039, 2.048)
       tau =      0.1688  (0.1685, 0.169)
       y0 =       11.43  (11.11, 11.76)

Goodness of fit:
  SSE: 1.62e+04
  R-square: 0.9999
  Adjusted R-square: 0.9999
  RMSE: 4.567
```
