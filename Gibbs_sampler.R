source("C:/Users/MERA/Desktop/My_R_Codes/Library.R")

## Inputs to the gibbs samplingr###
mu <- MleCluster$par[1] 
sigma <- MleCluster$par[2] 
df <- MleCluster$par[3]
nDraws <- 5000 #number of iterations
cluster <- cluster6
nc <- ncol(cluster)
parameters <- nc + 3
gibbsdraws3_3 <- matrix(0,nDraws,parameters)
v=seq(1,20,0.1)     # sqence values for degrees of freedom.

## Gibbs sampling starts##
for(i in 1:nDraws)
{
  ## update c given all other parameters
  alpha <- vector("numeric",length = ncol(cluster))
  for(j in 1:length(alpha))
  {
    alpha[j] <- sum(cluster[,j])
  }
  beta <- sum(dt.scaled(year,df = df,mean =mu,sd = sigma))
  cc <- vector("numeric",length = nc)
  cc <- rgamma(nc,alpha,beta)
  gibbsdraws3_3[i,1:nc] <- cc
  
  ## update df given all other parameters using inverse CDF
 
  df <- postSampleDF(data = cluster,year = year,mu=mu,sigma = sigma,c=cc,v=v)
  
  #df<- normalize_df(unnormilized_log_posterior)
  gibbsdraws3_3[i,nc+1] <- df
  
  
  # update sigma slice sampling
  sigma.smp <- MfU.Sample.Run(sigma,lliksigma,nsmp = 1,data=cluster,year=year,mu=mu,c=cc,df=df)
  sigma <- sigma.smp[1]
  gibbsdraws3_3[i,nc+2] <- sigma
  
  #update mu slice sampling
  mu.smp <- MfU.Sample.Run(mu,mu_post, nsmp = 1,sigma=sigma,data = cluster,year=year,c=cc,df=df,mu_prior = 40,sigma_prior = 2.55)
  mu <- mu.smp[1]
  gibbsdraws3_3[i,nc+3] <- mu
  
}


