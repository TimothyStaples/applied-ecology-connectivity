grouped.spline.CIs<-function(model, num.data, group.data, spline.point.n){

  # this function accepts a GAM model with grouped splines, and predicts significance
  # and confidence intervals across the spline at splint.point.n points using
  # finite differencing
group.no<-length(unique(group.data))

# generate data-frame of x-values to predict off
x.mesh<-data.frame(PY=rep(seq(min(num.data),
                              max(num.data),
                              length=spline.point.n), 
                          group.no),
                   
                   group=rep(levels(group.data), 
                             each=spline.point.n))

# predict points for slope
gam.predict<-data.frame(fit=predict(model, newdata=x.mesh, 
                                    se.fit=TRUE, type="link")$fit,
                        se.fit=predict(model, newdata=x.mesh, 
                                       se.fit=TRUE, type="link")$se.fit,
                        x.mesh)

# get spline derivative significance #

A<-model

X0<-predict(A, newdata=x.mesh,
            type="lpmatrix")

eps<- (max(num.data)-min(num.data))*1e-5 # finite difference interval

x.mesh$PY<-as.numeric(x.mesh$PY) + eps

X1<-predict(A, 
            newdata=x.mesh, 
            type="lpmatrix")

Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives


# get slope estimates and CIs across spline, and determine if 95% CIs cross 0.
CIs<-lapply(levels(group.data), function(x){

  group.cols<-which(grepl(paste0(x), colnames(Xp)))

  Xi <- Xp*0 
  Xi[,group.cols] <- Xp[,group.cols] ## Xi%*%coef(b) = smooth deriv i
  df <- Xi%*%coef(A)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%A$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
  cbind(gam.predict$fit[1:200], df[1:200])
  
  CIs<-data.frame(upper=df + 2*df.sd,
                  lower=df - 2*df.sd)

  CIs$diff0<-ifelse(CIs[,1]<0 & CIs[,2]<0 |
                    CIs[,1]>0 & CIs[,2]>0, TRUE, FALSE)
  
  # remove zeros for where other groups weren't tested
  rows<-which(CIs[,1]!=0)

  cbind(gam.predict[rows,],
        diff0=CIs[rows,"diff0"])
  
  })

return(CIs)
}