# main code
#generate simululation data

## 无差异实验
treadn<- read.delim("data/cpg_chr21_all_n.txt")
dim(treadn)
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

}

save(testRegion, file="data/simulationdata8c8c_1.Rdata.Rdata")

