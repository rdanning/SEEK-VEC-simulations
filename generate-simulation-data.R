# input parameters:
# sim: simulation number
# parameters.fname: filename of a parameters csv file containing the following parameter values:
#   p: number of vocabulary words
#   K: number of topics
#   n: number of documents
#   N.factor: n*N.factor = number of words per document
#   i.pct: proportion of vocabulary that is identifying (stronger loading for 1 topic)
#   m.pct: proportion of vocabulary that is metainformative (stronger loading for 2 topics)
#   scale.factor: factor by which identifying/metainformative loadings are stronger on their corresponding topics (higher scale.factor = stronger signal)
#   wfhh: word frequency heterogeneity/homogeneity factor: smaller values correspond to greater heterogeneity

# filename of parameters for simulation
parameters.fname <- "parameters.csv"

# change 'sim' to correspond to the simulation number
sim <- 1

out.fname <- paste0("data_",sim,".RData")
set.seed(sim)



# function to obtain the topic matrix
get.A <- function(p,K,i.pct,m.pct,scale.factor,wfhh){
  
  i.n <- round(p*i.pct)
  m.n <- round(p*m.pct)
  u.n <- p-i.n-m.n
  
  
  A.i <- get.A.imu(i.n,K,wfhh,scale.factor)
  A.m <- get.A.imu(m.n,K,wfhh,scale.factor,inflate.sample = 2)
  A.u <- get.A.imu(u.n,K,wfhh,do.inflate=FALSE)
  A.all <- rbind(A.i,A.m,A.u)
  A <- prop.table(A.all,2)
  
  return(A)

}

# helper function for obtaining the topic matrix
# performs scaling based on word category
inflate.A.row <- function(i,A,scale.factor,inflate.sample){
  row <- A[i,]
  n.idx <- inflate.sample
  inflate.idx <- sample(1:ncol(A),n.idx)
  inflation <- rep(1,ncol(A))
  inflation[inflate.idx] <- scale.factor
  inflated.row <- row*inflation
  return(inflated.row)
}

# helper function for obtaining the topic matrix
# generates probabilities based on the wfhh value and then scales the probabilities as needed based on the word category (I/M/U)
get.A.imu <- function(p,K,wfhh,scale.factor, do.inflate = TRUE, inflate.sample = 1){
  probs <- (1:p + wfhh)^-1.07
  A <- matrix(rep(probs,each=K),ncol=K,byrow = TRUE)
  if(do.inflate){
    A.scaled <- t(sapply(1:nrow(A),inflate.A.row,A,scale.factor,inflate.sample))
    return(A.scaled)
  } else{return(A)}
}

# function to return a random document-topic matrix
get.W <- function(K,n){
  W <- matrix(runif(n*K), K, n)
  return(prop.table(W,2))
}

# function to obtain the coprus matrix based on the underlying multinomial probabilites
get.D.vec <- function(i,D0,N){
  col <- D0[,i]
  new.col <- rmultinom(1,N,col)
  return(new.col)
}


# function to get the true hallmark matrix
get.true.H <- function(B,i.pct,m.pct){
  Hi <- t(sapply(1:(nrow(B)*i.pct),get.H.row,B,1))
  Hm <- t(sapply((nrow(B)*i.pct)+1:(nrow(B)*m.pct),get.H.row,B,2))
  Hu <- matrix(0,nrow = nrow(B) * (1-i.pct-m.pct),ncol = ncol(B))
  H <- rbind(Hi,Hm,Hu)
  return(H)
}

# helper function for getting the true hallmark matrix
get.H.row <- function(i,B,n){
  row <- B[i,]
  new.row <- rep(0,length(row))
  pick.idx <- order(row,decreasing = TRUE)[1:n]
  new.row[pick.idx] <- 1
  return(new.row)
}


# function to generate simulated data from parameters given in file
generate_data <- function(p,K,n,N.factor,i.pct,m.pct,scale.factor,wfhh,seed){

  
  A <- get.A(p,K,i.pct,m.pct,scale.factor,wfhh)
  B <- prop.table(A,1)
  
  W <- get.W(K,n)
  D0 <- A %*% W
  D <- sapply(1:ncol(D0),get.D.vec,D0,N.factor*p)
  
  H <- get.true.H(B,i.pct,m.pct)
  O.true <- H %*% t(H)
  
  return(list(A = A, # topic matrix
              B = B, # loadings matrix
              H = H, # hallmark matrix
              O.true = O.true, # joint hallmark projection matrix
              D = D # corpus matrix
              ))
  
}



parameters.df <- read.csv(parameters.fname)
parameters.row <- (seed %% nrow(parameters.df)) + 1
parameters <- parameters.df[parameters.row,]


sim.data <- generate_data(parameters$p,
                               parameters$K,
                               parameters$n,
                               parameters$N.factor,
                               parameters$i.pct,
                               parameters$m.pct,
                               parameters$scale.factor,
                               parameters$wfhh)

save(sim.data, file = out.fname)




