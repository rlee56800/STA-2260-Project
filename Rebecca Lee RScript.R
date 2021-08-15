####################### OPEN FILE #######################
### removed file read for privacy

#qsar <- read.csv("")
head(qsar)

n <-dim(qsar)[1] # number of observations
d <-dim(qsar)[2] # number of variables




####################### SCATTERPLOTS #######################

install.packages("leaps")
library(leaps)

panel.hist <- function(x, ...)
{
  usr <-par("usr"); on.exit(par(usr))
  par(usr =c(usr[1:2], 0, 1.5) )
  h <-hist(x, plot =FALSE)
  breaks <-h$breaks; nB <-length(breaks)
  y <-h$counts; y <-y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col ="cyan", ...)
}
panel.cor <-function(x, y, digits =2, prefix ="", cex.cor, ...)
{
  usr <-par("usr"); on.exit(par(usr))
  par(usr =c(0, 1, 0, 1))
  r <-abs(cor(x, y))
  txt <-format(c(r, 0.123456789), digits = digits)[1]
  txt <-paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <-0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex =cex.cor )
}
pairs(qsar_aquatic_toxicity, pch=19, lower.panel=panel.cor, diag.panel=panel.hist,
      upper.panel=panel.smooth)

#plot(qsar) # [don't post]

y <-qsar$V9 #pulls the V9 variable
par(mfrow=c(1,1))
  ## received
  ##  > Error in plot.new() : figure margins too large
  ## set margins to 1,1
for(i in 1:d)
{
  x <-qsar[,i] #predictor i
  plot(x,y,
       xlab = names(qsar)[i], # name of predictor i
       ylab = "V9")
}
#par(mfrow=c(1,1)) # reset the plot window size back to normal.
  ## never changed the plot window size




####################### M1 #######################

# V4
fit.V4 <-lm(V9 ~V4, data =qsar) #linear
fit.V4_2 <-lm(V9 ~V4 +I(V4^2), data =qsar) #quadratic
AIC(fit.V4, fit.V4_2)
##         df      AIC
##fit.V4    3 1926.827
##fit.V4_2  4 1920.820

# V5
fit.V5 <-lm(V9 ~V5, data =qsar) #linear
fit.V5_2 <-lm(V9 ~V5 +I(V5^2), data =qsar) #quadratic
AIC(fit.V5, fit.V5_2)
##         df      AIC
##fit.V5    3 2020.952
##fit.V5_2  4 1976.263

# V6
fit.V6 <-lm(V9 ~V6, data =qsar) #linear
fit.V6_2 <-lm(V9 ~V6 +I(V6^2), data =qsar) #quadratic
AIC(fit.V6, fit.V6_2)
##         df      AIC
##fit.V6    3 2055.500
##fit.V6_2  4 2049.264


summary(fit.V4)
## Call:
##   lm(formula = V9 ~ V4, data = qsar)
## 
## Residuals:
##   Min      1Q  Median      3Q     Max 
## -2.9044 -0.9382 -0.2406  0.6309  5.0898 
## 
## Coefficients:
##   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  3.47399    0.10021   34.67   <2e-16 ***
##   V4           0.51197    0.03462   14.79   <2e-16 ***
##   ---
##   Signif. codes:  
##   0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.408 on 544 degrees of freedom
## Multiple R-squared:  0.2868,	Adjusted R-squared:  0.2855 
## F-statistic: 218.7 on 1 and 544 DF,  p-value: < 2.2e-16


summary(fit.V4_2)
## Call:
##   lm(formula = V9 ~ V4 + I(V4^2), data = qsar)
## 
## Residuals:
##   Min      1Q  Median      3Q     Max 
## -3.6316 -0.8938 -0.2169  0.6507  5.1558 
## 
## Coefficients:
##   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  3.51673    0.10071  34.920  < 2e-16 ***
##   V4           0.38167    0.05744   6.644 7.42e-11 ***
##   I(V4^2)      0.03087    0.01090   2.832   0.0048 ** 
##   ---
##   Signif. codes:  
##   0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.399 on 543 degrees of freedom
## Multiple R-squared:  0.2972,	Adjusted R-squared:  0.2946 
## F-statistic: 114.8 on 2 and 543 DF,  p-value: < 2.2e-16

# to plot the fitted line for a linear one-predictor model
plot(qsar$V4, qsar$V9)
abline(fit.V4, col ="green")

# to plot the fitted line for a quadratic, we need to predict it:
x <-qsar$V4
xmesh <-seq(0.5*min(x), 2*max(x), by=0.1)
yhat <-predict(fit.V4_2,newdata =data.frame(V4=xmesh))

lines(xmesh, yhat, col ="purple")
legend("topleft", #location of legend
       c("Linear", "Quadratic"), # names
       lty =c(1,1),  # just leave this as 1 for each line you're drawing
       col =c("green", "purple"), #colors
       cex = 0.6 ## legend box was enormous
)




####################### M2 #######################

# 5 highest correlation values
m2.1<-lm(V9~V2+I(V2^2)+V5+I(V5^2)+V6+I(V6^2)+V3+V7, data=qsar)
summary(m2.1)
## Call:
##   lm(formula = V9 ~ V2 + I(V2^2) + V5 + I(V5^2) + V6 + I(V6^2) + 
##        V3 + V7, data = qsar)
## 
## Residuals:
##   Min      1Q  Median      3Q     Max 
## -3.1202 -0.8869 -0.1678  0.6859  5.3275 
## 
## Coefficients:
##   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  1.720e+00  5.323e-01   3.232  0.00131 ** 
##   V2          -9.174e-03  2.033e-03  -4.512 7.89e-06 ***
##   I(V2^2)      1.350e-05  5.246e-06   2.574  0.01033 *  
##   V5           2.712e+00  3.273e-01   8.286 9.38e-16 ***
##   I(V5^2)     -2.543e-01  5.560e-02  -4.574 5.94e-06 ***
##   V6          -1.746e+00  6.775e-01  -2.577  0.01024 *  
##   I(V6^2)      3.207e-01  2.848e-01   1.126  0.26066    
## V3          -1.499e-01  6.461e-02  -2.320  0.02070 *  
##   V7          -7.965e-02  4.821e-02  -1.652  0.09904 .  
## ---
##   Signif. codes:  
##   0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.296 on 537 degrees of freedom
## Multiple R-squared:  0.4033,	Adjusted R-squared:  0.3944 
## F-statistic: 45.36 on 8 and 537 DF,  p-value: < 2.2e-16
AIC(m2.1)
## 1843.467


# 3 highest correlation values
m2.2<-lm(V9~V5+I(V5^2)+V6+I(V6^2)+V3, data=qsar)
summary(m2.2)
## Call:
##   lm(formula = V9 ~ V5 + I(V5^2) + V6 + I(V6^2) + V3, data = qsar)
## 
## Residuals:
##   Min      1Q  Median      3Q     Max 
## -3.1171 -0.9005 -0.1600  0.7261  5.4029 
## 
## Coefficients:
##   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  2.36923    0.52892   4.479 9.14e-06 ***
##   V5           2.35034    0.32036   7.337 8.08e-13 ***
##   I(V5^2)     -0.22683    0.05365  -4.228 2.77e-05 ***
##   V6          -2.17544    0.68697  -3.167  0.00163 ** 
##   I(V6^2)      0.42384    0.29000   1.462  0.14445    
## V3          -0.28525    0.03755  -7.598 1.34e-13 ***
##   ---
##   Signif. codes:  
##   0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.325 on 540 degrees of freedom
## Multiple R-squared:  0.3726,	Adjusted R-squared:  0.3668 
## F-statistic: 64.14 on 5 and 540 DF,  p-value: < 2.2e-16
AIC(m2.2)
## 1864.803


# similar shapes
m2.3<-lm(V9~V1+I(V1^2)+V2+I(V2^2)+V5+I(V5^2)+V3, data=qsar)
summary(m2.3)
## Call:
##   lm(formula = V9 ~ V1 + I(V1^2) + V2 + I(V2^2) + V5 + I(V5^2) + 
##        V3, data = qsar)
## 
## Residuals:
##   Min      1Q  Median      3Q     Max 
## -3.9671 -0.9112 -0.0758  0.8172  4.4273 
## 
## Coefficients:
##   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  8.716e-02  4.367e-01   0.200  0.84190    
## V1           1.504e-02  5.324e-03   2.824  0.00491 ** 
##   I(V1^2)     -1.533e-05  3.365e-05  -0.456  0.64883    
## V2          -2.238e-02  3.437e-03  -6.512 1.70e-10 ***
##   I(V2^2)      2.446e-05  1.274e-05   1.921  0.05531 .  
## V5           2.819e+00  3.210e-01   8.781  < 2e-16 ***
##   I(V5^2)     -2.741e-01  5.551e-02  -4.938 1.06e-06 ***
##   V3          -1.183e-01  6.285e-02  -1.883  0.06029 .  
## ---
##   Signif. codes:  
##   0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.32 on 538 degrees of freedom
## Multiple R-squared:  0.3799,	Adjusted R-squared:  0.3718 
## F-statistic: 47.09 on 7 and 538 DF,  p-value: < 2.2e-16
AIC(m2.3)
## 1862.419
