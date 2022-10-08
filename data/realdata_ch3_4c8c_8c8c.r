##本代码作用：数据处理

##输入数据：mat_n,mat_x(分别表示总读数数据矩阵和甲基化读数矩阵)
##输出数据：染色体号_4cVS8c_data.Rdata，染色体号_8cVS8c_data.Rdata （例如chr3_4c8c_data.Rdata）
##输出数据形式均为列表形式，列表中存放了过滤位点后并划分好区域后生成的列表testRegion和将mat_n和mat_x按照位点位置排好序且过滤位点后的矩阵treadn和treadx

##get 50% data:对于一个位点，若所有样本在该位点的读数数据缺失值个数超过样本个数的50%则该位点不参与运算，过滤掉该位点

rm(list=ls())
library("limma")

#过滤函数
#filter function：最终得到一个列表，其中含有生成的testRegion列表、treadn和treadx
filter<-function(mat_n,mat_x,sample8c,sample4c,per){
  #mat_n:total reads matrix
  #mat_x:methylation reads matrix
  #sample8c:Column of 8 cell samples 8细胞样本所在列
  #sample4c:Column of 4 cell samples 4细胞样本所在列
  #per:过滤位点时所设置的缺失值占比值最小值，占比超过per则过滤掉此位点

  ##对原数据mat_n和mat_x按照位点位置先进行排序
  mat_x=mat_x[order(as.numeric(rownames(mat_x))),]
  mat_n=mat_n[order(as.numeric(rownames(mat_n))),]

  ##将4cell和8cell样本的x和n数据分别表示出来
  cell4_x=mat_x[,sample4c]
  cell4_n=mat_n[,sample4c]
  cell8_x=mat_x[,sample8c]
  cell8_n=mat_n[,sample8c]

  ##生成testRegion
  index=c() ###记录过滤的位点有哪些
  precite=as.numeric(rownames(mat_n[1,])) ###记录某区域的第一个位点的位置，这里初始化precite为输入数据中第一个位点的位置
  j=1 ###j记录
  testRegion=vector(mode="list")
  groupindex=c()
  testRegion1=vector(mode="list")
  m=1
  for (i in 1:dim(mat_n)[1]) {
    if((sum(is.na(cell8_n[i,]))/dim(cell8_n)[2])>per | (sum(is.na(cell4_n[i,]))/dim(cell4_n)[2])>per){
	  ###NA值大于样本数50%时记录该位点是第几个
      index=c(index,i)
    }else{
	  ###如果位点被保留，则进行区域划分，以一个区域最大长度为300bp且一个区域至少包括3个位点为标准划分区域
      if((as.numeric(rownames(mat_n[i,]))-precite)<=300){
        groupindex=c(groupindex,i)
        testRegion[[j]]=vector("list",2)
        names(testRegion[[j]])=c("x","n")
        testRegion[[j]][[1]]=mat_x[groupindex,]
        testRegion[[j]][[2]]=mat_n[groupindex,]
      }else{
        if(j==length(testRegion)){
          if(dim(testRegion[[j]]$n)[1]>=3){
            testRegion1[[m]]=vector("list",2)
            names(testRegion1[[m]])=c("x","n")
            testRegion1[[m]]=testRegion[[j]]
            m=m+1
          }
        }
        precite=as.numeric(rownames(mat_n[i,]))
        j=j+1
        groupindex=c(i)
      }
    }
    if(i%%1000==0)  cat(i," ")
  }

  if(length(index)!=0){
    ###在原mat_x和mat_n中删除过滤掉的位点信息，生成新的mat_n和mat_x
    cell8_n=cell8_n[-index,]
    cell8_x=cell8_x[-index,]
    cell4_n=cell4_n[-index,]
    cell4_x=cell4_x[-index,]
    mat_x_sort=mat_x[-index,]
    mat_n_sort=mat_n[-index,]
  }

  return(list(treadn=mat_n_sort,treadx=mat_x_sort,testRegion=testRegion1,index=index))
}


#读取数据
mat_n=read.delim("mat_n")
mat_x=read.delim("mat_x")

#8cell,4cell样本所在列号
sample8c=c(2:12,16,17,19,21:25,27,29,33,35,36,38:40,43:45,47:49,52:55,57:60,62:64,66,69,70,73)
sample4c=c(1,13:15,18,20,26,28,30:32,34,37,41,42,46,50,51,56,61,65,67,68,71,72)
per=0.5 #设置占比值

#得到8cVS4c实验所用的输入数据
read=filter(mat_n,mat_x,sample8c,sample4c,per)
chr3_data=read
save("chr3_data", file="data/chr3_4cVS8c_data.Rdata")


###############################################################################################################################################

#生成8c8c的mat_n和mat_x数据
#从48个8细胞样本中选取40个样本作为无差异实验样本
mat_n_8c8c=data.frame(GSM2986370=mat_n$GSM2986370, GSM2481623=mat_n$GSM2481623, GSM2481642=mat_n$GSM2481642, GSM2986374=mat_n$GSM2986374,
                      GSM2481641=mat_n$GSM2481641, GSM2481622=mat_n$GSM2481622, GSM2481631=mat_n$GSM2481631, GSM2481637=mat_n$GSM2481637,
                      GSM2481619=mat_n$GSM2481619, GSM2986372=mat_n$GSM2986372, GSM2481643=mat_n$GSM2481643, GSM2481618=mat_n$GSM2481618,
                      GSM2986380=mat_n$GSM2986380, GSM2986384=mat_n$GSM2986384, GSM2986368=mat_n$GSM2986368, GSM2481626=mat_n$GSM2481626,
                      GSM2986389=mat_n$GSM2986389, GSM2481630=mat_n$GSM2481630, GSM2986373=mat_n$GSM2986373, GSM2986378=mat_n$GSM2986378,
                      GSM2986369=mat_n$GSM2986369, GSM2986382=mat_n$GSM2986382, GSM2481628=mat_n$GSM2481628, GSM2481621=mat_n$GSM2481621,
                      GSM2481634=mat_n$GSM2481634, GSM2481629=mat_n$GSM2481629, GSM2481624=mat_n$GSM2481624, GSM2481632=mat_n$GSM2481632,
                      GSM2481640=mat_n$GSM2481640, GSM2481627=mat_n$GSM2481627, GSM2986377=mat_n$GSM2986377, GSM2986387=mat_n$GSM2986387,
                      GSM2481633=mat_n$GSM2481633, GSM2986385=mat_n$GSM2986385, GSM2481635=mat_n$GSM2481635, GSM2481636=mat_n$GSM2481636,
                      GSM2986379=mat_n$GSM2986379, GSM2986388=mat_n$GSM2986388, GSM2986383=mat_n$GSM2986383, GSM2986381=mat_n$GSM2986381
)
rownames(mat_n_8c8c)=rownames(mat_n)

mat_x_8c8c=data.frame(GSM2986370=mat_x$GSM2986370, GSM2481623=mat_x$GSM2481623, GSM2481642=mat_x$GSM2481642, GSM2986374=mat_x$GSM2986374,
                      GSM2481641=mat_x$GSM2481641, GSM2481622=mat_x$GSM2481622, GSM2481631=mat_x$GSM2481631, GSM2481637=mat_x$GSM2481637,
                      GSM2481619=mat_x$GSM2481619, GSM2986372=mat_x$GSM2986372, GSM2481643=mat_x$GSM2481643, GSM2481618=mat_x$GSM2481618,
                      GSM2986380=mat_x$GSM2986380, GSM2986384=mat_x$GSM2986384, GSM2986368=mat_x$GSM2986368, GSM2481626=mat_x$GSM2481626,
                      GSM2986389=mat_x$GSM2986389, GSM2481630=mat_x$GSM2481630, GSM2986373=mat_x$GSM2986373, GSM2986378=mat_x$GSM2986378,
                      GSM2986369=mat_x$GSM2986369, GSM2986382=mat_x$GSM2986382, GSM2481628=mat_x$GSM2481628, GSM2481621=mat_x$GSM2481621,
                      GSM2481634=mat_x$GSM2481634, GSM2481629=mat_x$GSM2481629, GSM2481624=mat_x$GSM2481624, GSM2481632=mat_x$GSM2481632,
                      GSM2481640=mat_x$GSM2481640, GSM2481627=mat_x$GSM2481627, GSM2986377=mat_x$GSM2986377, GSM2986387=mat_x$GSM2986387,
                      GSM2481633=mat_x$GSM2481633, GSM2986385=mat_x$GSM2986385, GSM2481635=mat_x$GSM2481635, GSM2481636=mat_x$GSM2481636,
                      GSM2986379=mat_x$GSM2986379, GSM2986388=mat_x$GSM2986388, GSM2986383=mat_x$GSM2986383, GSM2986381=mat_x$GSM2986381
)
rownames(mat_x_8c8c)=rownames(mat_x)

# 8c vs 8c的样本平均分成两份当成两类样本
sample8c1=c(1:20)
sample8c2=c(21:40)

#8c VS 8c得到输入数据
read8c8c=filter(mat_n_8c8c,mat_x_8c8c,sample8c1,sample8c2,per)
chr3_8c8c_data=read8c8c
save("chr3_8c8c_data", file="data/chr3_8cVS8c_data.Rdata")

