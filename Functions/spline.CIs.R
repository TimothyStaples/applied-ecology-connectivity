spline.CIs<-function(model, num.data, spline.point.n){

  # this function accepts a GAM model with a single spline, and predicts significance
  # and confidence intervals across the spline at splint.point.n points using
  # finite differencing. This won't work for splines fit to different group levels
  # use the specific group function for those models.
  
  x.mesh<-data.frame(PY=seq(min(num.data),max(num.data),length=spline.point.n))
  
  # predict points for slope
gam.predict<-data.frame(fit=predict(model, newdata=x.mesh, 
                                    se.fit=FALSE, type="link"),
                        se.fit=predict(model, newdata=x.mesh, 
                                    se.fit=TRUE, type="link")$se.fit,
                        x.mesh)
                        
# get spline derivative significance #

A<-model


X0<-predict(A, newdata=x.mesh, type="lpmatrix")

eps<- (max(num.data)-min(num.data))*1e-5 # finite difference interval

x.mesh$PY<-as.numeric(x.mesh$PY) + eps

X1<-predict(A, 
            newdata=x.mesh, 
            type="lpmatrix")

Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives


# get slope estimates and CIs across spline, and determine if 95% CIs cross 0.

  Xi <- Xp*0 
  Xi[,-1] <- Xp[,-1] ## Xi%*%coef(b) = smooth deriv i
  df <- Xi%*%coef(A)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%A$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
  CIs<-data.frame(upper=df + 2*df.sd,
                  lower=df - 2*df.sd)
  
  CIs$diff0<-ifelse(CIs[,1]<0 & CIs[,2]<0 |
                      CIs[,1]>0 & CIs[,2]>0, TRUE, FALSE)
  
  return(cbind(gam.predict,
               diff0=CIs[,"diff0"]))
}