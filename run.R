# index <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
# Input_path = Input2
# Output_path = Input3


# index <- 2;core_num <- 36;one_batch <- 3000000
# index <- 2;core_num <- 36;one_batch <- 3000000
index <- as.numeric(index)
core_num <- 5 # as.numeric(commandArgs(trailingOnly=TRUE)[2])
one_batch <- 50000 # as.numeric(commandArgs(trailingOnly=TRUE)[3])
sub_batch_size <- 10000
node_list <- c("Bark.thickness","Conduit.diam.","Crown.diameter",
                "Crown.height","Leaf.area","Leaf.density","Leaf.K","Leaf.N",         
                "Leaf.P","Leaf.thickness","Leaf.Vcmax","Root.depth","Specific.leaf.area",
                "Stem.diameter","Stomatal.conduct.","Tree.height","Wood.density","Seed.dry.mass")

# if(!"Matrix" %in% rownames(installed.packages())){install.packages("Matrix")}
if(!"dplyr" %in% rownames(installed.packages())){install.packages("dplyr")}
if(!"foreach" %in% rownames(installed.packages())){install.packages("foreach")}
if(!"doParallel" %in% rownames(installed.packages())){install.packages("doParallel")}
if(!"tseries" %in% rownames(installed.packages())){install.packages("tseries")}

print("all packages installed.")

# Most centrality measures are for nodes but not all 
library(dplyr)
library(foreach)
library(doParallel)
library(tseries)
library(igraph)
library('Matrix')



ecoregion_trait <- read.csv(paste0(Input_path,'Global_tree_trait/Imputed_trait_data.csv'))
ecoregion_trait <- ecoregion_trait %>% na.omit() %>% filter(ECO_ID>0)
data <- ecoregion_trait[node_list]

batch_list <- seq(1,2^length(node_list),one_batch)


#####     source("../Parallel_packages.R")     #####


#library('MatrixStats')

# Shannon entropy

# Shannon entropy

entropia<-function(a)
{
  a<-a[which(a>0)];
  -sum(a*log(a));
}

#function

alpha<-function(g){
  
  N<-length(V(g))
  
  r<-sort(alpha.centrality(g,exo=degree(g)/(N-1),alpha=1/N))/((N^2))
  
  return(c(r,max(c(0,1-sum(r)))))
  
}

#returns the node distance matrix

node_distance<-function(g){
  n<-length(V(g))
  if(n==1){
    retorno=1
  }
  if(n>1){
    a<-Matrix(0,nrow=n,ncol=n,sparse=TRUE)
    m<-shortest.paths(g,algorithm=c("unweighted"))
    m[which(m=="Inf")]<-n
    quem<-setdiff(intersect(m,m),0)
    for(j in (1:length(quem))){
      l<-which(m==quem[j])/n
      linhas<-floor(l)+1
      posicoesm1<-which(l==floor(l))
      if(length(posicoesm1)>0){
        linhas[posicoesm1]<-linhas[posicoesm1]-1
      }
      a[1:n,quem[j]]<-hist(linhas,plot=FALSE,breaks=(0:n))$counts
    }
    #m<-c()
    retorno=(a/(n-1))
  }
  return(retorno)
}


# nnd

nnd<-function(g){
  
  N<-length(V(g))
  
  nd<-node_distance(g)
  
  pdfm<-colMeans(nd)
  
  norm<-log(max(c(2,length(which(pdfm[1:(N-1)]>0))+1)))
  
  return(c(pdfm,max(c(0,entropia(pdfm)-entropia(nd)/N))/norm))
  
  
}


#function

d<-function(g,h,w1,w2,w3){
  
  first<-0
  
  second<-0
  
  third<-0
  
  # g<-read.graph(g,format=c("edgelist"),directed=FALSE)
  # 
  # h<-read.graph(h,format=c("edgelist"),directed=FALSE)
  
  N<-length(V(g))
  
  M<-length(V(h))
  
  PM<-matrix(0,ncol=max(c(M,N)))
  
  if(w1+w2>0){
    
    pg=nnd(g)
    
    PM[1:(N-1)]=pg[1:(N-1)]
    
    PM[length(PM)]<-pg[N]
    
    ph=nnd(h)
    
    PM[1:(M-1)]=PM[1:(M-1)]+ph[1:(M-1)]
    
    PM[length(PM)]<-PM[length(PM)]+ph[M]
    
    PM<-PM/2
    
    first<-sqrt(max(c((entropia(PM)-(entropia(pg[1:N])+entropia(ph[1:M]))/2)/log(2),0)))
    
    second<-abs(sqrt(pg[N+1])-sqrt(ph[M+1]))
    
    
  }
  
  if(w3>0){
    
    pg<-alpha(g)
    
    ph<-alpha(h)
    
    m<-max(c(length(pg),length(ph)))
    
    Pg<-matrix(0,ncol=m)
    
    Ph<-matrix(0,ncol=m)
    
    Pg[(m-length(pg)+1):m]<-pg
    Ph[(m-length(ph)+1):m]<-ph
    third<-third+sqrt((entropia((Pg+Ph)/2)-(entropia(pg)+entropia(ph))/2)/log(2))
  }
  
  return(w1*first+w2*second+w3*third)
  
  
}


node_distance_w<-function(g){
  n<-length(V(g))
  if(n==1){
    retorno=1
  }
  if(n>1){
    a<-Matrix(0,nrow=n,ncol=n,sparse=TRUE)
    
    m0 <-shortest.paths(g,algorithm=c("unweighted"))
    bar_l <- mean(m0[is.finite(m0) & m0 != 0],na.rm = TRUE)

    m1 <- distances(g,weights = 1/E(g)$weight)
    bar_lw <- mean(m1[is.finite(m1) & m1 != 0],na.rm = TRUE)
    
    if (is.na(bar_l)){
      m = m0
      m[which(m=="Inf")] <- n
    }else{
      m = ceiling(m1*bar_l/bar_lw)
      # m_add = m
      # m_add[m_add==0] = 0
      # m_add[m_add!=0] = 1
      # m = m + m_add
      m[which(m=="Inf")]<- max(m[is.finite(m)]) + 1
    }
    
    quem<-setdiff(intersect(m,m),0)
    for(j in (1:length(quem))){
      l<-which(m==quem[j])/n
      linhas<-floor(l)+1
      posicoesm1<-which(l==floor(l))
      if(length(posicoesm1)>0){
        linhas[posicoesm1]<-linhas[posicoesm1]-1
      }
      a[1:n,quem[j]]<-hist(linhas,plot=FALSE,breaks=(0:n))$counts
    }
    #m<-c()
    retorno=(a/(n-1))
  }
  return(retorno)
}

# nnd
nnd_w<-function(g){
  
  N<-length(V(g))
  
  nd<-node_distance_w(g)
  
  pdfm<-colMeans(nd)
  
  norm<-log(max(c(2,length(which(pdfm[1:(N-1)]>0))+1)))
  
  return(c(pdfm,max(c(0,entropia(pdfm)-entropia(nd)/N))/norm))
}
#function
wd<-function(g,h,w1,w2,w3){
  
  first<-0
  
  second<-0
  
  third<-0
  
  # g<-read.graph(g,format=c("edgelist"),directed=FALSE)
  # 
  # h<-read.graph(h,format=c("edgelist"),directed=FALSE)
  
  N<-length(V(g))
  
  M<-length(V(h))
  
  PM<-matrix(0,ncol=max(c(M,N)))
  
  if(w1+w2>0){
    
    pg=nnd_w(g)
    
    PM[1:(N-1)]=pg[1:(N-1)]
    
    PM[length(PM)]<-pg[N]
    
    ph=nnd_w(h)
    
    PM[1:(M-1)]=PM[1:(M-1)]+ph[1:(M-1)]
    
    PM[length(PM)]<-PM[length(PM)]+ph[M]
    
    PM<-PM/2
    
    first<-sqrt(max(c((entropia(PM)-(entropia(pg[1:N])+entropia(ph[1:M]))/2)/log(2),0)))
    
    second<-abs(sqrt(pg[N+1])-sqrt(ph[M+1]))
    
    
  }
  
  if(w3>0){
    
    pg<-alpha(g)
    
    ph<-alpha(h)
    
    m<-max(c(length(pg),length(ph)))
    
    Pg<-matrix(0,ncol=m)
    
    Ph<-matrix(0,ncol=m)
    
    Pg[(m-length(pg)+1):m]<-pg
    Ph[(m-length(ph)+1):m]<-ph
    third<-third+sqrt((entropia((Pg+Ph)/2)-(entropia(pg)+entropia(ph))/2)/log(2))
  }
  
  return(w1*first+w2*second+w3*third)
  
  
}

log_normalize <- function(x){
  x_star <- BBmisc::normalize(log(x - min(x,na.rm = T) + 1),method = "standardize")
  return(x_star)
}
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  return(data.frame( row = rownames(cormat)[row(cormat)[ut]],
                     column = rownames(cormat)[col(cormat)[ut]],
                     cor  =(cormat)[ut],
                     p = pmat[ut]))
}

network_metrics_with_combination <- function(number,items,all_edges,Net_full){
  # 根据数字确定需要计算的网络
  PowerSetsBinary <- function(number,items=items){
    one_network_nodes = c()
    for (j in c(0:(length(items)-1))) {
      # print(bitwShiftR(number,j))
      if (bitwShiftR(number,j) %% 2 == 1){
        # print(items[j+1])
        one_network_nodes <- c(one_network_nodes,items[j+1])
        # print(one_network_nodes)
      }
    }
    if (length(one_network_nodes)>1){
      return(one_network_nodes)
    } else{
      return(c())
    }
  }
  node_list <- PowerSetsBinary(number,items)
  if (length(node_list)>1){
    # 根据给定的性状组合trait_combination，从数据data中筛选对应的性状列进行网络构建

    ### node mean sd calculate
    node_value_sd <- as.data.frame(node_list)
    colnames(node_value_sd) <- 'Nodes'
    
    ### calculate the edges
    edges <- subset(all_edges,from %in% node_list & to %in% node_list) 
    ### build the network
    Net <- graph_from_data_frame(d=edges,vertices=node_value_sd, directed = FALSE)
    d_nc <- NA
    dw_sr <- NA
    if (length(node_list)>4){
       d_nc  <-  d(Net,Net_full,w1=0.45,w2=0.45,w3=0.1)*10000
       dw_sr <- wd(Net,Net_full,w1=0.45,w2=0.45,w3=0.1)*10000 
    }
    return(c(number,length(node_list),round(d_nc,0),round(dw_sr,0)))
  }
}

# ecoregion_trait <- subset(ecoregion_trait,GrowthForm=='herb')


temp_corr_data <- Hmisc::rcorr(as.matrix(data,type="pearson"))
temp_corr_data <-flattenCorrMatrix(temp_corr_data$r, temp_corr_data$P)
temp_corr_data <-temp_corr_data[(temp_corr_data$p<=0.05),]
colnames(temp_corr_data) <- c("to","from","estimate","p.value")
all_edges <- temp_corr_data[,c("to","from","estimate")]
all_edges$estimate <- abs(all_edges$estimate)
all_edges <-all_edges[(all_edges$estimate>=0.2),]
colnames(all_edges)<- c("from","to","weight")

### build the network based on all edges
Net_full <- graph_from_data_frame(d=all_edges,vertices=as.data.frame(node_list), directed = FALSE)

if (index<=length(batch_list)){
  start_num = batch_list[index]
  # 此处batch_list[index + 1]-1，是因为最后一个数与下一个循环的开头重复。
  end_num = ifelse(index+1<=length(batch_list),batch_list[index + 1]-1,2^length(node_list))
  sub_batch_size <- min(sub_batch_size,one_batch)
  sub_batch_list =seq(start_num,end_num,sub_batch_size)
  log_file = paste0('./network_modularity_log_',index,'.csv')
  if (!file.exists(log_file)){
    log_data <- data.frame(finished = 0)
    write.csv(log_data,log_file,row.names = FALSE)
  }
  # 使用length(sub_batch_list)总体长度是为了防止当不能整除时造成的尾部数据缺失
  for (i in c(1:(length(sub_batch_list)))) {
    sub_start = sub_batch_list[i]
    # 此处sub_batch_list[i+1]-1，是因为最后一个数与下一个循环的开头重复。
    # 最后一个阶段的最终数值为end_num，是一个整数或最终的2^length(node_list)
    sub_end = ifelse(i+1<=length(sub_batch_list),sub_batch_list[i+1]-1,end_num)
    
    log_csv <- read.csv(log_file,header=TRUE)
    if (!(sub_start %in% log_csv$finished)){
          system.time({
            # 网络计算
            cores <- min(detectCores(), core_num)
            print(paste0(Sys.time(),' Start computing ',as.character(sub_start),' --------- ',as.character(sub_end),'with cores:',as.character(cores)))
            cl<- makeCluster(cores)      
            registerDoParallel(cl)       #进程注册
            mydata1 <- foreach(
              # trait_combination = trait_combination,          #输入等待请求的参数
              number_list = c(sub_start:sub_end),          #输入等待请求的参数
              .combine=rbind,  #返回结果的整合
              .packages = c("dplyr","tidyverse","igraph","magrittr","Matrix") #多个进程共享的系统环境
            ) %dopar% network_metrics_with_combination(number_list,items=node_list,all_edges = all_edges,Net_full=Net_full)
            stopCluster(cl)
            colnames(mydata1) <- c('name_number','size','d','dw')
                                   
            write.csv(mydata1,paste0(Output_path,'D_nc_',as.character(sub_start),'_',as.character(sub_end),'.csv'),row.names = FALSE)
            print(paste0(Sys.time(),' Computing finished ',as.character(sub_start),' --------- ',as.character(sub_end)))
            log_csv[nrow(log_csv)+1,]=sub_start
            write.csv(log_csv,log_file,row.names = FALSE)
          })
    }  
  }
} else {
  print('out of list!!!')
}