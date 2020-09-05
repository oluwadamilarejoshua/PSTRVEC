#########################################################################
##
##  IVTAR.R
##  
##  This is a R procedure.  It computes estimates and confidence
##  intervals for threshold models with endogenous regressors.
##  
##  It computes the estimator described in
##  "Instrumental Variable Estimation of a Threshold Model"
##  written by Mehmet Caner and Bruce E. Hansen
##  
##  Be sure to have the Gauss graphics library active.
##  
##  If you run this program, it is currently set up to generate some simulated data 
##  and estimate the model on this data.  
##  For your own application, you can substitute your own data.
##  
#########################################################################

# Functions #

##************************************************************
##  REGRESS
##  
##  Computes a linear regression. Uses generalized inverse if X'X is singular
##  
##  Format
##  beta <- regress(y,x)
##  
##  Inputs
##  y	      nxm	dependent variable(s)
##  x	      nxk	independent variables (should include constant)
##  
##  Output
##  beta    kxm	Regression slope estimates
##  
##************************************************************

regress <- function(y,x){
  if (qr(x)$rank==ncol(x)) beta <- qr.solve(x,y)
  if (qr(x)$rank<ncol(x)) beta <- (qr(t(x)%*%x)$qr)%*%(t(x)%*%y)
  beta
}

##************************************************************
##  JOINT_THRESH
##  
##  Estimates the threshold in a multivariate threshold model
##  
##  Format
##  output <- joint_thresh(y,x,q)
##  yhat <- output$yhat
##  qhat <- output$qhat
##  
##  Inputs
##  y	      nxm	dependent variable(s)
##  x	      nxk	independent variables (should include constant)
##  q	      nx1	threshold variable
##  
##  Outputs
##  yhat	nxm	Predicted values for y
##  qhat	1x1	Threshold Estimate 
##  
##************************************************************

joint_thresh <- function(y,x,q){
  n <- nrow(y)
  k <- ncol(x)
  e <- y-x%*%regress(y,x)
  s0 <- det(t(e)%*%e)    
  n1 <- round(.05*n)+k
  n2 <- round(.95*n)-k
  qs <- sort(q)
  qs <- qs[n1:n2]
  qs <- as.matrix(unique(qs))
  qn <- nrow(qs)
  sn <- matrix(0,qn,1)
  for (r in 1:qn){
    d <- (q<=qs[r])
    xx <- x*(d%*%matrix(1,1,k))
    xx <- xx-x%*%regress(xx,x)
    ex <- e-xx%*%regress(e,xx)
    sn[r] <- det(t(ex)%*%ex)   
  }
  r <- which.min(sn)
  smin <- sn[r]
  qhat <- qs[r]
  d <- (q<=qhat)
  x1 <- x*(d%*%matrix(1,1,k))
  x2 <- x*((1-d)%*%matrix(1,1,k))
  beta1 <- regress(y,x1)
  beta2 <- regress(y,x2)
  yhat <- x1%*%beta1+x2%*%beta2
  list(yhat=yhat,qhat=qhat)
}

##************************************************************
##  JOINT_THRESH_CI
##  
##  Estimates the threshold in a multivariate threshold model 
##  and computes asymptotic confidence intervals
##  
##  Format
##  output <- joint_thresh(y,x,q,_conf,_graph)
##  yhat <- output$yhat
##  
##  Inputs
##  y	      nxm	dependent variable(s)
##  x	      nxk	independent variables (should include constant)
##  q	      nx1	threshold variable
##  conf	1x1	level for confidence interval, e.g. conf=.9
##  graph	1x1	Set graph=1 to view graph of likelihood ratio
##  		      Set graph=0 to not display graph
##  
##  Outputs
##  yhat	nxm	Predicted values for y
##  qhat	1x1	Threshold Estimate 
##  qcf_0	1x2	Confidence interval for threshold, no heteroskedasticity correction
##  qcf_h1	1x2	Confidence interval for threshold, heteroskedasticity correction 
##  		      using quadratic variance estimate 
##  		      (only if m=1)
##  qcf_h2	1x2	Confidence interval for threshold, heteroskedasticity correction 
##  		      using variance estimated by Epanechnikov kernel with automatic bandwidth
##  		      (only if m=1)
##  
##************************************************************

joint_thresh_ci <- function(y,x,q,conf,graph){
  n <- nrow(y)
  k <- ncol(x)
  e <- y-x%*%regress(y,x)
  s0 <- det(t(e)%*%e)    
  n1 <- round(.05*n)+k
  n2 <- round(.95*n)-k
  qs <- sort(q)
  qs <- qs[n1:n2]
  qs <- as.matrix(unique(qs))
  qn <- nrow(qs)
  sn <- matrix(0,qn,1)
  for (r in 1:qn){
    d <- (q<=qs[r])
    xx <- x*(d%*%matrix(1,1,k))
    xx <- xx-x%*%regress(xx,x)
    ex <- e-xx%*%regress(e,xx)
    sn[r] <- det(t(ex)%*%ex)   
  }
  r <- which.min(sn)
  smin <- sn[r]
  qhat <- qs[r]
  d <- (q<=qhat)
  x1 <- x*(d%*%matrix(1,1,k))
  x2 <- x*((1-d)%*%matrix(1,1,k))
  beta1 <- regress(y,x1)
  beta2 <- regress(y,x2)
  yhat <- x1%*%beta1+x2%*%beta2
  e <- y-yhat
  lr <- n*(sn/smin-1)
  sig2 <- smin/n
  
  if (ncol(y)> 1){
    eta1 <- 1 
    eta2 <- 1
  }else{
    r1 <- (x%*%(beta1-beta2))^2
    r2 <- r1*(e^2)
    qx <- cbind(q^0,q^1,q^2)
    qh <- cbind(qhat^0,qhat^1,qhat^2)  
    m1 <- qr.solve(qx,r1)  
    m2 <- qr.solve(qx,r2)  
    g1 <- qh%*%m1
    g2 <- qh%*%m2
    eta1 <- as.vector((g2/g1)/sig2)
    sigq <- sqrt(mean((q-mean(q))^2))
    hband <- 2.344*sigq/(n^(.2))
    u <- (qhat-q)/hband
    u2 <- u^2
    f <- mean((1-u2)*(u2<=1))*(.75/hband)
    df <- -mean(-u*(u2<=1))*(1.5/(hband^2))
    eps <- r1 - qx%*%m1
    sige <- (t(eps)%*%eps)/(n-3)
    hband <- as.vector(sige/(4*f*((m1[3]+(m1[2]+2*m1[3]*qhat)*df/f)^2)))
    u2 <- ((qhat-q)/hband)^2
    kh <- ((1-u2)*.75/hband)*(u2<=1)
    g1 <- mean(kh*r1)
    g2 <- mean(kh*r2)
    eta2 <- as.vector((g2/g1)/sig2)
  }
  c1 <- -2*log(1-sqrt(conf))
  lr0 <- (lr >= c1)
  lr1 <- (lr >= (c1*eta1))
  lr2 <- (lr >= (c1*eta2))
  if (max(lr0)==1){
    qcf_0 <- cbind(qs[which.min(lr0)],qs[qn+1-which.min(rev(lr0))])
  }else{
    qcf_0 <- cbind(qs[1],qs[qn])
  }
  if (max(lr1)==1){
    qcf_h1 <- cbind(qs[which.min(lr1)],qs[qn+1-which.min(rev(lr1))])
  }else{
    qcf_h1 <- cbind(qs[1],qs[qn])
  }
  if (max(lr2)==1){
    qcf_h2 <- cbind(qs[which.min(lr2)],qs[qn+1-which.min(rev(lr2))])
  }else{
    qcf_h2 <- cbind(qs[1],qs[qn])
  }
  
  if (graph==1){
    x11()
    mtit <- "Confidence Interval Construction for Threshold"
    ytit <- "Likelihood Ratio Sequence in gamma" 
    xtit <- "Threshold Variable"
    clr <- matrix(1,qn,1)*c1
    if (ncol(y) == 1) clr <- cbind(clr,(clr*eta1),(clr*eta2))
    plot(qs,lr,lty=1,col=1,type="l",ann=0)
    lines(qs,clr[,1],lty=2,col=2)
    if (ncol(y) == 1){
      lines(qs,clr[,2],lty=3,col=3)
      lines(qs,clr[,3],lty=4,col=4)
    }
    title(main=mtit,ylab=ytit,xlab=xtit) 
    tit1 <- "LRn(gamma)"
    tit2 <- "90% Critical"
    tit3 <- "Hetero Corrected - 1"
    tit4 <- "Hetero Corrected - 2" 
    if (ncol(y) != 1) legend("topleft",c(tit1,tit2),lty=c(1,2),col=c(1,2))
    if (ncol(y) == 1) legend("topleft",c(tit1,tit2,tit3,tit4),lty=c(1,2,3,4),col=c(1,2,3,4))
  }
  list(yhat=yhat,qhat=qhat,qcf_0=qcf_0,qcf_h1=qcf_h1,qcf_h2=qcf_h2)
}

##************************************************************
##  GMM_LINEAR
##  
##  Computes the GMM estimator of a linear model
##  
##  Format
##  output <- gmm_linear(y,z,x)
##  beta <- output$beta
##  
##  Inputs
##  y	      nx1	dependent variable
##  z	      nxk	rhs variables
##  x	      nxl	instruments variables (should include constant and exogenous parts of z), l>=k
##  
##  Outputs
##  beta	kx1	Regression slope estimates
##  se	kx1	standard errors
##  jstat	1x1	J Statistic
##  
##************************************************************

gmm_linear <- function(y,z,x){
  pihat <- regress(z,x)
  xz <- t(x)%*%z 
  xy <- t(x)%*%y
  beta <- qr.solve((t(pihat)%*%xz),(t(pihat)%*%xy))
  e <- y-z%*%beta
  xe <- x*(e%*%matrix(1,1,ncol(x)))
  g <- solve(t(xe)%*%xe)
  v <- solve(t(xz)%*%g%*%xz)
  beta <- v%*%(t(xz)%*%g%*%xy)
  se <- as.matrix(sqrt(diag(v)))
  m <- t(x)%*%(y-z%*%beta)
  jstat <- t(m)%*%g%*%m
  list(beta=beta,se=se,jstat=jstat)
}

##************************************************************
##  GMM_THRESH
##  
##  Computes threshold model with endogenous variables
##  
##  Format
##  
##  qhat <- gmm_thresh(y,z1,z2,x,q,conf,conf1,conf2,reduced)
##  
##  Inputs
##  y	      nx1	dependent variable
##  z1	nxk1	endogenous rhs variables
##  z2	nxk2	exogenous rhs variables
##  x	      nxm	instrumental variables not included in z2, m>=k1
##  q	      nx1	threshold variable
##  conf	1x1	level for confidence interval for threshold, e.g. conf=.9
##  conf1	1x1	confidence level for threshold interval used as input to slope intervals, e.g. conf2=.8 
##  conf2	1x1	interval dummy - to determine method for threshold interval used as input to slope intervals
##  		      set to 0 for uncorrected (homoskedastic) interval
##  		      set to 1 for heteroskedastic correction by quadratic
##  		      set to 2 for heteroskedastic correction by nonparametric kernel
##  reduced	1x1	reduced form dummy
##  		      set to 0 if reduced form is linear (no threshold)
##  		      set to 1 if reduced form is a linear threshold model
##  
##  Output
##  qhat	1x1	threshold estimate
##  
##  Most of the output is written to the screen.
##  
##************************************************************

gmm_thresh <- function(y,z1,z2,x,q,conf,conf1,conf2,reduced){
  xx <- cbind(x,z2) 
  if (reduced==0) z1hat <- t(xx)%*%regress(z1,xx)
  if (reduced==1){ 
    out <-  joint_thresh(z1,xx,q)
    z1hat <- out$yhat
    rhohat <- out$qhat
  }
  zhat <- cbind(z1hat,z2)
  
  out <- joint_thresh_ci(y,zhat,q,conf,1)
  yhat <- out$yhat
  qhat <- out$qhat
  qcf0 <- out$qcf_0
  qcf1 <- out$qcf_h1 
  qcf2 <- out$qcf_h1 
  
  z <- cbind(z1,z2)
  da <- (q<=qhat)
  db <- 1-da
  ya <- y[da%*%matrix(1,1,ncol(y))>0]
  ya <- matrix(ya,length(ya)/ncol(y),ncol(y))
  xa <- xx[da%*%matrix(1,1,ncol(xx))>0]
  xa <- matrix(xa,length(xa)/ncol(xx),ncol(xx))
  za <- z[da%*%matrix(1,1,ncol(z))>0]
  za <- matrix(za,length(za)/ncol(z),ncol(z))
  yb <- y[db%*%matrix(1,1,ncol(y))>0]
  yb <- matrix(yb,length(yb)/ncol(y),ncol(y))
  xb <- xx[db%*%matrix(1,1,ncol(xx))>0]
  xb <- matrix(xb,length(xb)/ncol(xx),ncol(xx))
  zb <- z[db%*%matrix(1,1,ncol(z))>0]
  zb <- matrix(zb,length(zb)/ncol(z),ncol(z))
  out1 <- gmm_linear(ya,za,xa)
  out2 <- gmm_linear(yb,zb,xb)
  beta1 <- out1$beta
  se1 <- out1$se 
  jstat1 <- out1$jstat
  beta2 <- out2$beta
  se2 <- out2$se 
  jstat2 <- out2$jstat
  
  beta1l <- beta1-se1*1.96
  beta1u <- beta1+se1*1.96
  beta2l <- beta2-se2*1.96
  beta2u <- beta2+se2*1.96
  
  out <- joint_thresh_ci(y,zhat,q,conf1,0)
  yhat <- out$yhat
  qhat <- out$qhat
  qcf0i <- out$qcf_0
  qcf1i <- out$qcf_h1 
  qcf2i <- out$qcf_h1 
  
  if (conf2==0) qcf <- qcf0i 
  if (conf2==1) qcf <- qcf1i 
  if (conf2==2) qcf <- qcf2i 
  qq <- unique(q)
  qq <- as.matrix(sort(qq))
  qq <- qq[qq<=qcf[2]]
  qq <- as.matrix(qq[qq>=qcf[1]])
  
  for (i in 1:nrow(qq)){
    qi <- qq[i]
    dai <- (q<=qi)
    dbi <- 1-dai
    ya <- y[dai%*%matrix(1,1,ncol(y))>0]
    ya <- matrix(ya,length(ya)/ncol(y),ncol(y))
    xa <- xx[dai%*%matrix(1,1,ncol(xx))>0] 
    xa <- matrix(xa,length(xa)/ncol(xx),ncol(xx))
    za <- z[dai%*%matrix(1,1,ncol(z))>0]
    za <- matrix(za,length(za)/ncol(z),ncol(z))
    yb <- y[dbi%*%matrix(1,1,ncol(y))>0]
    yb <- matrix(yb,length(yb)/ncol(y),ncol(y))
    xb <- xx[dbi%*%matrix(1,1,ncol(xx))>0]
    xb <- matrix(xb,length(xb)/ncol(xx),ncol(xx))
    zb <- z[dbi%*%matrix(1,1,ncol(z))>0]
    zb <- matrix(zb,length(zb)/ncol(z),ncol(z))
    out1 <- gmm_linear(ya,za,xa)
    out2 <- gmm_linear(yb,zb,xb)
    beta1i <- out1$beta
    se1i <- out1$se 
    jstat1i <- out1$jstat 
    beta2i <- out2$beta
    se2i <- out2$se 
    jstat2i <- out2$jstat
    beta1l <- apply(t(cbind((beta1i-se1i*1.96),beta1l)),2,min)
    beta1u <- apply(t(cbind((beta1i+se1i*1.96),beta1u)),2,max)
    beta2l <- apply(t(cbind((beta2i-se2i*1.96),beta2l)),2,min)
    beta2u <- apply(t(cbind((beta2i+se2i*1.96),beta2u)),2,max)
  }
  
  cat ("\n")
  cat ("\n")
  cat ("Threshold Estimate                        ", qhat, "\n")
  cat ("Confidence Interval - Uncorrected         ", qcf0, "\n")
  cat ("Confidence Interval - Het Corrected Quad  ", qcf1, "\n")
  cat ("Confidence Interval - Het Corrected NP    ", qcf2, "\n")
  cat ("\n")
  cat ("\n")
  cat ("Regime 1 : Threshold variable less than   ", qhat, "\n")
  cat ("Number of observations                    ", sum(da),"\n")
  cat ("\n")
  cat ("Estimates     S.E.          Lower         Upper", "\n")
  for (i in 1:nrow(beta1)) cat (beta1[i],"  ",se1[i],"  ",beta1l[i],"  ",beta1u[i],"\n")
  cat ("\n")
  cat ("Regime 2 : Threshold variable greater than", qhat, "\n")
  cat ("Number of observations                    ", sum(db),"\n")
  cat ("\n")
  cat ("Estimates     S.E.          Lower         Upper", "\n")
  for (i in 1:nrow(beta2)) cat (beta2[i],"  ",se2[i],"  ",beta2l[i],"  ",beta2u[i],"\n")
  cat ("\n")
  qhat
}

#************************************************************#

# Example using Simulated Data #

n <- 100
kx <- 4 
sig <- matrix(c(1,0.5,0.5,1),2,2)  
e <- matrix(rnorm(n*2),n,2)%*%chol(sig)
q <- matrix(rnorm(n),n,1)  
x <- matrix(rnorm(n*kx),n,kx)  
z2 <- cbind(matrix(1,n,1),q)  
z1 <- x%*%matrix(1,kx,1)*3+e[,2]
zz <- (z1+z2%*%matrix(1,2,1))*.1
y <- zz*(q<0)-zz*(q>=0)+e[,1]

qhat <- gmm_thresh(y,z1,z2,x,q,.9,.8,1,1)


