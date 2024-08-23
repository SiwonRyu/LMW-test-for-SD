#install.packages("foreach")
#install.packages("doParallel")
library("foreach")
library("doParallel")




rm(list=ls(all=TRUE)) 

# Set parameters
alpha <- 0.05
ngrid <- 40
R <- 1
n1 <- 100
n2 <- 100

# Design
design <- c(0,1,0,1)
sample1 <- matrix(rnorm(n1,design[1],design[2]), ncol=1)
sample2 <- matrix(rnorm(n2,design[3],design[4]), ncol=1)
grid <- matrix(seq(min(rbind(sample1,sample2)),max(rbind(sample1,sample2)),length=ngrid))
s <- 1
b1 <-50
b2 <-50

# LMW-test for "stochastic dominance"
# H0 : X1 s-th order SD F2
# Input parameter : sample1, sample2, grid, s(SD order), subsamplesize1, subsamplesize2
# Output : 1x5 vector (LMW-teststatistic, pvalue from subsampling)

lmwtest <-function(sample1,sample2,grid,s,b1,b2,type){
  n1 = dim(sample1)[1]
  n2 = dim(sample2)[1]
  size = n1*n2/(n1+n2)
  operator <- function(x,grid,s){
    return(t(apply(x,1,'<',grid))*apply(grid,1,'-',x)/factorial(s))
  }
  ecdf <- function(x,grid,s){
    return(matrix(colMeans(operator(x,grid,s)),nrow=1))
  }
  D <- function(sample1,sample2,grid,s){
    return(ecdf(sample1,grid,s) - ecdf(sample2,grid,s))
  }
  # For SD
  if (type == "SD"){
    stat <-function(size,sample1,sample2,grid,s){
      stat = sqrt(size) * max(D(sample1,sample2,grid,s))
      return(stat)
    }
  } else if(type == "MSD") {
    stat <-function(size,sample1,sample2,grid,s){
      stat = sqrt(size) * min(max(D(sample1,sample2,grid,s),D(sample2,sample1,grid,s)))
      return(stat)
    }
  } else if(type == "PSD") {
    stat <-function(size,sample1,sample2,grid,s){
      grid_plus = matrix(grid[grid>=0])
      grid_minus = matrix(grid[grid<0])
      stat = sqrt(size) * min(
        max(
          apply(D(sample1,sample2,grid_plus,2),2,'-',D(sample1,sample2,grid_minus,2))
        )
      )
      return(stat)
    }
  }
  
  subsample <-function(sample,subsize,nsub){
    subindex = apply(array(0:nsub,dim=c(nsub,1)),1,'+',array(1:subsize,dim=c(1,subsize)))
    return(subindex)
  }
  substat <- function(sample1,sample2,grid,s,b1,b2){
    n1 = dim(sample1)[1]
    n2 = dim(sample2)[1]
    n = n1*n2/(n1+n2)
    lbd = n2/(n1+n2)
    b = lbd * b1
    nsub = min(n1-b1+1, n2-b2+1)
    
    subindex1 = subsample(sample1,b1,nsub)
    subindex2 = subsample(sample2,b2,nsub)
    substat <- array(0, dim=c(1,ngrid,nsub))
    for(i in 1:nsub){
      subsample1 <- array(sample1[subindex1[,i]],dim=c(b1,1))
      subsample2 <- array(sample2[subindex2[,i]],dim=c(b2,1))
      substat[,,i] <- stat(b,subsample1,subsample2,grid,s)
    }
    return(array(apply(substat,3,max),dim=c(1,1,nsub)))
  }
  lmw <- stat(size,sample1,sample2,grid,1)
  pval <- mean(substat(sample1,sample2,grid,s,b1,b2)>lmw)
  return(c(lmw,pval))
}


result = lmwtest(sample1,sample2,grid,s,b1,b2,"PSD")
result


# Simulation part
# init <- Sys.time()
# R = 100
# Cl <- makeCluster(10,type = "PSOCK")
# registerDoParallel(Cl)
# rej_pre <- foreach(i = 1:R) %dopar%{
#   sample1 <- matrix(rnorm(n1,design[1],design[2]), ncol=1)
#   sample2 <- matrix(rnorm(n2,design[3],design[4]), ncol=1)
#   grid <- matrix(seq(min(rbind(sample1,sample2)),max(rbind(sample1,sample2)),length=ngrid))
#   result <- lmwtest(sample1,sample2,grid,s,b1,b2)
#   result[2]
# }
# stopCluster(Cl)
# rej = matrix(unlist(rej_pre), nrow = R, byrow = TRUE)
# sprintf("Rejection probability = %1.3f", mean(rej[,1]))
# Sys.time() - init


