###################
############
#Fisher information Matrix 4x4
fisherinform<-function(theta){
  finform= array(0,c(4,4))
  
  finform[1,1]<-(1-theta[2]+0.01)/(theta[1]*(1-theta[1]-theta[2])+0.01)
  finform[1,2]<-1/(1-theta[1]-theta[2]) 
  finform[1,3]<-0
  finform[1,4]<-0
  finform[2,1]<-finform[1,2] 
  finform[2,2]<-(1-theta[1]+0.01)/(theta[2]*(1-theta[1]-theta[2])+0.01)
  finform[2,3]<-0
  finform[2,4]<-0  
  finform[3,1]<-0
  finform[3,2]<-0  
  finform[3,3]<-(1-theta[1]-theta[2])*((trigamma(theta[3])-trigamma(theta[3]+theta[4]))+0.01)
  finform[3,4]<-(1-theta[1]-theta[2])*(-trigamma(theta[3]+theta[4])+0.01) 
  finform[4,1]<-0
  finform[4,2]<-0 
  finform[4,3]<-finform[3,4] 
  finform[4,4]<-(1-theta[1]-theta[2])*((trigamma(theta[4])-trigamma(theta[3]+theta[4]))+0.01) 
  return(finform)
}
############
#Fisher information Matrix 4x4
fisherinform1<-function(theta){
  finform= array(0,c(4,4))
  aa<-10^(-5)
  if(theta[1]<10^(-6)){
    finform[1,1]<-(1-theta[2]+aa)/(theta[1]*(1-theta[1]-theta[2])+aa)
  }else{
    finform[1,1]<-(1-theta[2])/(theta[1]*(1-theta[1]-theta[2]))
  }
  
  finform[1,2]<-1/(1-theta[1]-theta[2]) 
  finform[1,3]<-0
  finform[1,4]<-0
  finform[2,1]<-finform[1,2]
  if(theta[2]<10^(-6)){
    finform[2,2]<-(1-theta[1]+aa)/(theta[2]*(1-theta[1]-theta[2])+aa)
  }else{
    finform[2,2]<-(1-theta[1])/(theta[2]*(1-theta[1]-theta[2]))
  }
  
  finform[2,3]<-0
  finform[2,4]<-0  
  finform[3,1]<-0
  finform[3,2]<-0 
  if(beta(theta[3],theta[4])<10^(-6)){
    finform[3,3]<-(1-theta[1]-theta[2])*((trigamma(theta[3])-trigamma(theta[3]+theta[4]))+aa)
    finform[3,4]<-(1-theta[1]-theta[2])*(-trigamma(theta[3]+theta[4])+aa)  
    finform[4,4]<-(1-theta[1]-theta[2])*((trigamma(theta[4])-trigamma(theta[3]+theta[4]))+aa)
    
  }else{
    finform[3,3]<-(1-theta[1]-theta[2])*(trigamma(theta[3])-trigamma(theta[3]+theta[4]))
    finform[3,4]<-(1-theta[1]-theta[2])*(-trigamma(theta[3]+theta[4]))  
    finform[4,4]<-(1-theta[1]-theta[2])*(trigamma(theta[4])-trigamma(theta[3]+theta[4]))
    
  }
  finform[4,1]<-0
  finform[4,2]<-0 
  finform[4,3]<-finform[3,4] 
  return(finform)
}
###equation set
fun<-function(alph, beta, c1, c2){
  f<-array(0,c(2,1))
  f[1]<-digamma(exp(alph))-digamma(exp(alph)+exp(beta))-c1
  f[2]<-digamma(exp(beta))-digamma(exp(alph)+exp(beta))-c2
  return(f)
}
###################
############
#derivative values,Jacobian Matrix 2x2
Jac<-function(alph, beta){
  jacb = array(0,c(2,2))
  jacb[1,1]<-exp(alph)*(trigamma(exp(alph))-trigamma(exp(alph)+exp(beta)))
  jacb[1,2]<--exp(beta)*(trigamma(exp(alph)+exp(beta)))
  jacb[2,1]<--exp(alph)*(trigamma(exp(alph)+exp(beta)))
  jacb[2,2]<-exp(beta)*(trigamma(exp(beta))-trigamma(exp(alph)+exp(beta)))
  return(jacb)
}
#################
NR<-function(alph,beta,c1, c2, tol, max.it){
  xold =c(alph,beta)
  xnew = xold-solve(Jac(alph, beta))%*%fun(alph, beta, c1, c2)   
  for (k in 1:max.it){
    niter = k
    xold = xnew
    #print(fun(xold[1], xold[2], c1, c2))
    #print(xold)
    xnew = xold -solve(Jac(xold[1], xold[2]))%*%fun(xold[1], xold[2], c1, c2)  
    d<-sum(abs(xold-xnew))
    #print(d)
    if((d < tol) | (sum(abs(fun(xnew[1], xnew[2], c1, c2)))<0.001) | max(abs(xold-xnew))/max(abs(xnew))<tol | (niter==max.it) ) break     
    xold<-xnew 
  }                                                                                                                        
  #list(r=xnew,niter=niter) 
  #print("NR OK")
  return(xnew)
}   

#solve(Jac(xold[1], xold[2]))%*%Jac(xold[1], xold[2])

###############
#function of parameters estimate with EM algorithm
#The algorithm for estimating the four-parameter model which is a random effects model with
#random coefficeint

logMLEestimate<-function(binn,binx,intPi0,intPi1,intalph,intbeta,toll=0.00001,iter=1000){
  nbin<-NULL
  xbin<-NULL
  #if(is.matrix(binn)){
  binn<-as.matrix(binn)
  binx<-as.matrix(binx)
  for(j in 1:dim(binn)[2]){
    nbin<-c(nbin, binn[!is.na(binn[,j]),j])
    xbin<-c(xbin, binx[!is.na(binx[,j]),j])
  }  
  #}else{
  #  nbin<-binn
  #  xbin<-binx
  # }
  ntotal<-length(nbin)
  
  buzeronum48c<-sum(xbin==0)/ntotal
  buonenum48c<-sum(xbin==nbin)/ntotal
  intPi0<-buzeronum48c
  intPi1<-buonenum48c
  rate<-xbin/nbin
  rate1<-rate[which(rate!=0 & rate!=1 )]
  if(sum(rate!=0 &rate!=1)>2){
    bumuem48c<-mean(rate1)
    buvarem48c<-var(rate1)
    bualfa48c<-bumuem48c*((bumuem48c-bumuem48c^2)/buvarem48c-1)
    bubeta48c<-(1-bumuem48c)*((bumuem48c-bumuem48c^2)/buvarem48c-1)    
    intalph<-bualfa48c
    intbeta<-bubeta48c
  }
  
  itPi01<-0.1#intPi0
  itPi11<-0.1#intPi1
  italph1<-0.1#intalph
  italbeta1<-0.1#intbeta 
  itPi02<--1
  itPi12<--1
  italph2<--1
  italbeta2<--1
  for(i in 1:iter) {
    
    ##E step: calculate di1,di2,di3 and di4
    if(beta(italph1,italbeta1)>10^(-5)){
      di1<-(itPi01/(itPi01+(1-itPi01-itPi11)*(beta(italph1,nbin+italbeta1)/beta(italph1,italbeta1))))
      sumd1<-sum(di1[which(xbin==0)]) 
      di2<-(itPi11/(itPi11+(1-itPi01-itPi11)*(beta(nbin+italph1,italbeta1)/beta(italph1,italbeta1))))
      sumd2<-sum(di2[which(xbin==nbin)]) 
      di31<-((1-itPi01-itPi11)*(beta(italph1,nbin+italbeta1)/beta(italph1,italbeta1))*(digamma(italph1)-digamma(italph1+nbin+italbeta1)))/(itPi01+(1-itPi01-itPi11)*(beta(italph1,nbin+italbeta1)/beta(italph1,italbeta1)))
      di32<-((1-itPi01-itPi11)*(beta(nbin+italph1,italbeta1)/beta(italph1,italbeta1))*(digamma(nbin+italph1)-digamma(italph1+nbin+italbeta1)))/(itPi11+(1-itPi01-itPi11)*(beta(nbin+italph1,italbeta1)/beta(italph1,italbeta1)))
      di33<-digamma(xbin+italph1)-digamma(italph1+nbin+italbeta1)
      sumd3<-sum(di31[which(xbin==0)])+sum(di32[which(xbin==nbin)])+sum(di33[which(xbin<nbin & xbin>0)])
      di41<-((1-itPi01-itPi11)*(beta(italph1,nbin+italbeta1)/beta(italph1,italbeta1))*(digamma(nbin+italbeta1)-digamma(italph1+nbin+italbeta1)))/(itPi01+(1-itPi01-itPi11)*(beta(italph1,nbin+italbeta1)/beta(italph1,italbeta1)))
      di42<-((1-itPi01-itPi11)*(beta(nbin+italph1,italbeta1)/beta(italph1,italbeta1))*(digamma(italbeta1)-digamma(italph1+nbin+italbeta1)))/(itPi11+(1-itPi01-itPi11)*(beta(nbin+italph1,italbeta1)/beta(italph1,italbeta1)))
      di43<-digamma(nbin-xbin+italbeta1)-digamma(italph1+nbin+italbeta1)
      sumd4<-sum(di41[which(xbin==0)])+sum(di42[which(xbin==nbin)])+sum(di43[which(xbin<nbin &xbin>0)])
    }else{
      di1<-(itPi01/(itPi01+(1-itPi01-itPi11)*((beta(italph1,nbin+italbeta1)+10^(-5))/(beta(italph1,italbeta1)+10^(-5)))))
      sumd1<-sum(di1[which(xbin==0)]) 
      di2<-(itPi11/(itPi11+(1-itPi01-itPi11)*((beta(nbin+italph1,italbeta1)+10^(-5))/(beta(italph1,italbeta1)+10^(-5)))))
      sumd2<-sum(di2[which(xbin==nbin)]) 
      di31<-((1-itPi01-itPi11)*((beta(italph1,nbin+italbeta1)+10^(-5))/(beta(italph1,italbeta1)+10^(-5)))*(digamma(italph1)-digamma(italph1+nbin+italbeta1)))/(itPi01+(1-itPi01-itPi11)*((beta(italph1,nbin+italbeta1)+10^(-5))/(beta(italph1,italbeta1)+10^(-5))))
      di32<-((1-itPi01-itPi11)*((beta(nbin+italph1,italbeta1)+10^(-5))/(beta(italph1,italbeta1)+10^(-5)))*(digamma(nbin+italph1)-digamma(italph1+nbin+italbeta1)))/(itPi11+(1-itPi01-itPi11)*((beta(nbin+italph1,italbeta1)+10^(-5))/(beta(italph1,italbeta1)+10^(-5))))
      di33<-digamma(xbin+italph1)-digamma(italph1+nbin+italbeta1)
      sumd3<-sum(di31[which(xbin==0)])+sum(di32[which(xbin==nbin)])+sum(di33[which(xbin<nbin &xbin>0)])
      di41<-((1-itPi01-itPi11)*((beta(italph1,nbin+italbeta1)+10^(-5))/(beta(italph1,italbeta1)+10^(-5)))*(digamma(nbin+italbeta1)-digamma(italph1+nbin+italbeta1)))/(itPi01+(1-itPi01-itPi11)*((beta(italph1,nbin+italbeta1)+10^(-5))/(beta(italph1,italbeta1)+10^(-5))))
      di42<-((1-itPi01-itPi11)*((beta(nbin+italph1,italbeta1)+10^(-5))/(beta(italph1,italbeta1)+10^(-5)))*(digamma(italbeta1)-digamma(italph1+nbin+italbeta1)))/(itPi11+(1-itPi01-itPi11)*((beta(nbin+italph1,italbeta1)+10^(-5))/(beta(italph1,italbeta1)+10^(-5))))
      di43<-digamma(nbin-xbin+italbeta1)-digamma(italph1+nbin+italbeta1)
      sumd4<-sum(di41[which(xbin==0)])+sum(di42[which(xbin==nbin)])+sum(di43[which(xbin<nbin &xbin>0)])
      
    }
    
    ddd<-c(sumd1,sumd2,sumd3,sumd4)
    # print(ddd)
    ###M step: maximum likelihood estimate Pi01, Pi11, alph and beta
    itPi02<-sumd1/ntotal
    itPi12<-sumd2/ntotal
    fc1<-sumd3/(ntotal-sumd1-sumd2)
    fc2<-sumd4/(ntotal-sumd1-sumd2)
    #niudun-raphsor iterate
    #fc1<--1.2878
    #fc2<--0.576
    #print(i)
    #print(itPi02)
    #print(itPi12)
    alphbetaest<-NR(alph=log(italph1),beta=log(italbeta1),c1=fc1, c2=fc2,tol=0.00001,max.it=500)
    
    #alphbetaest[alphbetaest<0]<-0.1
    italph2<-exp(alphbetaest[1])
    italbeta2<-exp(alphbetaest[2]) 
    
    #print(italph2)
    #print(italbeta2)
    likehold<-sumd1*log(itPi02)+sumd2*log(itPi12)+(ntotal-sumd1-sumd2)*log(1-itPi02-itPi12)-(ntotal-sumd1-sumd2)*log(beta(italph2,italbeta2))+sumd3*italph2+sumd4*italbeta2
    #print(likehold)
    dd<-abs(itPi02-itPi01)+abs(itPi12-itPi11)+abs(italph2-italph1)+abs(italbeta2-italbeta1)
    nitera = i
    #dd<-na.omit(dd)
    #print("dd去除NA")
    #print(dd)
    #print(c(itPi02,itPi12,italph2,italbeta2))
    if((dd<toll) | (nitera==iter) ) break   
    #print("OK")
    itPi01<-itPi02
    itPi11<-itPi12
    italph1<-italph2
    italbeta1<-italbeta2
  }
  canshu<-c(itPi02,itPi12,italph2,italbeta2)
  list(r=canshu,niter=nitera) 
}


###############
#function of parameters estimate with EM algorithm
#The algorithm for estimating the four-parameter model which is a random effects model with
#random coefficeint
e_DMV <- function(e_beta, p0, p1, a, b, n) {
  
  o_iter = 1
  while (o_iter > 1e-5)
  {
    p0_ = p0
    p1_ = p1
    a_ = a
    b_ = b
    
    # E Step
    d1 <- c()
    d2 <- c()
    d3 <- c()
    d4 <- c()
    for (i in 1:n) {
      d1 <- c(d1, (e_beta[i, 2] == 0)*p0/(p0 + (1 - p0 - p1)*beta(a, b + e_beta[i, 1])/beta(a, b)))
      d2 <- c(d2, (e_beta[i, 2] == e_beta[i, 1])*p1/(p1 + (1 - p0 - p1)*beta(a + e_beta[i, 1], b)/beta(a, b)))
      d3 <- c(d3, (e_beta[i, 2] == 0)*(1 - p1 - p0)*beta(a, b + e_beta[i, 1])/beta(a, b)*(digamma(a) - digamma(a + b + e_beta[i, 1]))/
                (p0 + (1 - p0 - p1)*beta(a, b + e_beta[i, 1])/beta(a, b)) +
                (e_beta[i, 2] > 0)*(e_beta[i, 2] < e_beta[i, 1])*(digamma(e_beta[i, 2] + a) - digamma(a + b + e_beta[i, 1])) +
                (e_beta[i, 2] == e_beta[i, 1])*(1 - p1 - p0)*beta(a + e_beta[i, 1], b)/beta(a, b)*(digamma(a + e_beta[i, 1]) - digamma(a + b + e_beta[i, 1]))/
                (p1 + (1 - p0 - p1)*beta(a + e_beta[i, 1], b)/beta(a, b)))
      d4 <- c(d4, (e_beta[i, 2] == 0)*(1 - p1 - p0)*beta(a, b + e_beta[i, 1])/beta(a, b)*(digamma(b + e_beta[i, 1]) - digamma(a + b + e_beta[i, 1]))/
                (p0 + (1 - p0 - p1)*beta(a, b + e_beta[i, 1])/beta(a, b)) +
                (e_beta[i, 2] > 0)*(e_beta[i, 2] < e_beta[i, 1])*(digamma(b + e_beta[i, 1] - e_beta[i, 2]) - digamma(a + b + e_beta[i, 1])) +
                (e_beta[i, 2] == e_beta[i, 1])*(1 - p1 - p0)*beta(a + e_beta[i, 1], b)/beta(a, b)*(digamma(b) - digamma(a + b + e_beta[i, 1]))/
                (p1 + (1 - p0 - p1)*beta(a + e_beta[i, 1], b)/beta(a, b)))
    }
    
    
    # M Step
    p0 = sum(d1)/n
    p1 = sum(d2)/n
    
    a1 = log(a)
    b1 = log(b)
    
    ab = matrix(c(a1, b1))
    fa = matrix(0, nrow = 2, ncol = 1)
    A = matrix(0, nrow = 2, ncol = 2)
    
    k1 = sum(d3)/sum(1 - d1 - d2)
    k2 = sum(d4)/sum(1 - d1 - d2)
    
    iter = 1
    while (iter > 1e-5) {
      a = exp(ab[1, 1])
      b = exp(ab[2, 1])
      A[1, 1] = a*(trigamma(a) - trigamma(a + b))
      A[1, 2] = -b*(trigamma(a + b))
      A[2, 1] = -a*(trigamma(a + b))
      A[2, 2] = b*(trigamma(b) - trigamma(a + b))
      fa[1, 1] = digamma(a) - digamma(a + b) - k1
      fa[2, 1] = digamma(b) - digamma(a + b) - k2
      ab = ab - solve(A) %*% fa
      iter = sqrt((ab[1, 1] - log(a))^2 + (ab[2, 1] - log(b))^2)
    }
    o_iter = sqrt((p0 - p0_)^2 + (p1 - p1_)^2 + (a - a_)^2 + (b - b_)^2)
    #cat(o_iter, '\n')
  }
  return(matrix(c(p0, p1, a, b), nrow = 1))
}


logMLEest<-function(binn,binx,intPi0,intPi1,intalph,intbeta){
  nbin<-NULL
  xbin<-NULL
  #if(is.matrix(binn)){
  binn<-as.matrix(binn)
  for(j in 1:dim(binn)[2]){
    nbin<-c(nbin, binn[!is.na(binn[,j]),j])
    xbin<-c(xbin, binx[!is.na(binx[,j]),j])
  }  
  #}else{
  #  nbin<-binn
  #  xbin<-binx
  # }
  ntotal<-length(nbin)
  
  buzeronum48c<-sum(xbin==0)/ntotal
  buonenum48c<-sum(xbin==nbin)/ntotal
  intPi0<-buzeronum48c
  intPi1<-buonenum48c
  rate<-xbin/nbin
  rate1<-rate[which(rate!=0 & rate!=1 )]
  if(sum(rate!=0 &rate!=1)>2){
    bumuem48c<-mean(rate1)
    buvarem48c<-var(rate1)
    bualfa48c<-bumuem48c*((bumuem48c-bumuem48c^2)/buvarem48c-1)
    bubeta48c<-(1-bumuem48c)*((bumuem48c-bumuem48c^2)/buvarem48c-1)    
    intalph<-bualfa48c
    intbeta<-bubeta48c
  }
  
  out<-cbind(nbin,xbin)
  dim(out)
  para1 = e_DMV(e_beta=out, p0=intPi0, p1=intPi1,a=intalph, b=intbeta, n=ntotal)
  
  return(para1)
}


