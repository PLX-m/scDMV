#simulation
#generate simululation data
setwd("D:/Files/Github/")
source("scDMV/data/functionscDMVm.r")
## 有差异实验
treadn<- read.delim("raw_data/cpg_chr21_all_n.txt")
dim(treadn) # 645756 76

scSMVpvalue<-NULL
absd = NULL
testRegion=vector(mode="list")

for (sim in 1:1000) {
  #################3
  #generate simulation data
  qujian<-10
  n1<-48

  eachbinn48c<-treadn[sample(dim(treadn)[1],qujian,replace =FALSE),4:76]

  for (j in 1:73) {
    nnumb<-c(treadn[!is.na(treadn[,(j+3)]),(j+3)])
    eachbinn48c[,j]<-sample(nnumb,qujian,replace =TRUE)
  }

  eachbinx48c<-eachbinn48c
  for (i in 1:qujian) {
    ###8 cell
    initialPi08c<-0.3
    initialPi18c<-0
    initialalfa8c<-1.8
    initialbeta8c<-1.5
    binn8c<-eachbinn48c[i,1:n1]
    binn8c1<-binn8c[!is.na(binn8c)]
    pio1<-rbinom(length(binn8c1),1,initialPi08c+initialPi18c)
    pio12<-pio1
    pio12[which(pio1==1)]<-rbinom(sum(pio1),1,initialPi18c/(initialPi08c+initialPi18c))
    pio12[which(pio1==0)]<-rbeta(sum(pio1==0),initialalfa8c, initialbeta8c)
    binx8c<-binn8c1
    binx8c[pio12==0]<-0
    nn<-binx8c[pio12>0 & pio12<1]
    binx8c[pio12>0 & pio12<1]<-mapply(rbinom, n = rep(1, length(nn)), size =nn, prob =pio12[pio12>0 & pio12<1])#rbinom(length(nn),nn,pio12[pio12>0 & pio12<1])
    binn8c[!is.na(binn8c)]<-binx8c
    eachbinx48c[i,1:n1]<-binn8c

    ###4 cell
    initialPi04c<-0.3
    initialPi14c<-0
    initialalfa4c<-1.8
    initialbeta4c<-1.5
    binn4c<-eachbinn48c[i,(n1+1):73]
    binn4c1<-binn4c[!is.na(binn4c)]
    pio1<-rbinom(length(binn4c1),1,initialPi04c+initialPi14c)
    pio12<-pio1
    pio12[which(pio1==1)]<-rbinom(sum(pio1),1,initialPi14c/(initialPi04c+initialPi14c))
    pio12[which(pio1==0)]<-rbeta(sum(pio1==0),initialalfa4c, initialbeta4c)
    binx4c<-binn4c1
    binx4c[pio12==0]<-0
    nn<-binx4c[pio12>0 & pio12<1]
    binx4c[pio12>0 & pio12<1]<-mapply(rbinom, n = rep(1, length(nn)), size =nn, prob =pio12[pio12>0 & pio12<1])#rbinom(length(nn),nn,pio12[pio12>0 & pio12<1])
    binn4c[!is.na(binn4c)]<-binx4c
    eachbinx48c[i,(n1+1):73]<-binn4c

  }

  testRegion[[sim]]=vector("list",2)
  names(testRegion[[sim]])=c("x","n")
  testRegion[[sim]][[1]]=eachbinx48c
  testRegion[[sim]][[2]]=eachbinn48c
  
  meansamplep8c<-NULL
  meansamplep4c<-NULL
  for (j in 1:n1) {
    n8c<-sum(eachbinn48c[!is.na(eachbinn48c[,j]),j])
    x8c<-sum(eachbinx48c[!is.na(eachbinx48c[,j]),j])
    meansamplep8c<- c(meansamplep8c, (x8c+0.00001)/n8c)
  }
  for (j in (n1+1):73) {
    n4c<-sum(eachbinn48c[!is.na(eachbinn48c[,j]),j])
    x4c<-sum(eachbinx48c[!is.na(eachbinx48c[,j]),j])
    meansamplep4c<- c(meansamplep4c, (x4c+0.00001)/n4c)
  }
  
  #generate n and x
  ntotal8c<-dim(eachbinn48c)[1]*n1
  ntotal4c<-dim(eachbinn48c)[1]*(73-n1)
  
  #Δ
  absd<-c(absd, abs(x8c/n8c-x4c/n4c))
  
  ########## scSMV test
  pmest8c<-logMLEestimate(binn=eachbinn48c[,1:n1],binx=eachbinx48c[,1:n1],intPi0=initialPi0,intPi1=initialPi1,intalph=initialalfa,intbeta=initialbeta,toll=0.00001,iter=1000)
  pramt8c<-pmest8c$r
  pmest8c
  
  pmest4c<-logMLEestimate(binn=eachbinn48c[,(n1+1):73],binx=eachbinx48c[,(n1+1):73],intPi0=initialPi0,intPi1=initialPi1,intalph=initialalfa,intbeta=initialbeta,toll=0.00001,iter=1000)
  pramt4c<-pmest4c$r
  pmest4c
  
  
  g8c<-pramt8c[2]+(1-pramt8c[1]-pramt8c[2])*(pramt8c[3]/(pramt8c[3]+pramt8c[4]))
  g4c<-pramt4c[2]+(1-pramt4c[1]-pramt4c[2])*(pramt4c[3]/(pramt4c[3]+pramt4c[4]))
  
  #derivative of the function
  derg8c<-c(-(pramt8c[3]/(pramt8c[3]+pramt8c[4])),(pramt8c[4]/(pramt8c[3]+pramt8c[4])),(1-pramt8c[1]-pramt8c[2])*(pramt8c[4]/(pramt8c[3]+pramt8c[4])^2),-(1-pramt8c[1]-pramt8c[2])*(pramt8c[3]/(pramt8c[3]+pramt8c[4])^2))
  derg4c<-c(-(pramt4c[3]/(pramt4c[3]+pramt4c[4])),(pramt4c[4]/(pramt4c[3]+pramt4c[4])),(1-pramt4c[1]-pramt4c[2])*(pramt4c[4]/(pramt4c[3]+pramt4c[4])^2),-(1-pramt4c[1]-pramt4c[2])*(pramt4c[3]/(pramt4c[3]+pramt4c[4])^2))
  
  #Fisher information
  fisher8c<-fisherinform1(theta=pramt8c)
  fisher4c<-fisherinform1(theta=pramt4c)
  var8c<-t(derg8c)%*%solve(fisher8c)%*%derg8c/ntotal8c
  var4c<-t(derg4c)%*%solve(fisher4c)%*%derg4c/ntotal4c
  #ward test statistics
  Tward<-(g8c-g4c)/sqrt((var8c+var4c)/4*qujian) #+0.001
  pvaluesmv<-2*pnorm(-abs(Tward))
  scSMVpvalue<-c(scSMVpvalue,pvaluesmv)
  if(sim%%100==0)  print(sim)

}

simulation<-list(testRegion=testRegion, scSMVpvalue=scSMVpvalue, dlt=absd)
save(simulation, file="data/simulationdata8c8c_1.Rdata")

cutoff<-10^(-2)
sum(scSMVpvalue<cutoff)

sum((scSMVpvalue<cutoff)&(absd>=0.15))








