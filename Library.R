## all the needed packages ##
library(metRology)   # for the density of t distribuion
library(reshape2)    # plotting 
library(ggplot2)     # ploting
library(MfUSampler)  # for the slice sampler
library(MCMCpack)    #  for MCMC of (slice sampling)
library(latex2exp)   # latex functions

### function which returns number of patents of a tech in a given time(lambda_it).###

averagepatents <- function(c,year,mu,sigma,df) 
{
  
  return(c*dt.scaled(year,df=df,mean = mu,sd = sigma))
  
}

datanormalize <- function(x) #normilizes data in the interval0 to 1
{
  y<- c()
  for(i in 1:length(x))
  {
    y[i] <- (x[i]-min(x))/(max(x)-min(x))
  }
  return(y)
}
#### Likelihood functions ####

llikTechs<- function(par,data,year) # logliklihod function for student-t distribution of spec. techs 
{
  lambda <- par[4]*dt.scaled(year,df=par[3],mean = par[1], sd = par[2])
  loglik <- sum(dpois(data,lambda=lambda,log = T))
  return(loglik)
}

llikCluster <- function(par,data,year)         # loglik.  function  of a given clusterfor Student-t 
{
  nc <- ncol(data)                             #number of techs with in a cluster
  c = par[4:(nc+3)]                             # the constant terms of all tecs in the cluster
  lambda<- outer(dt.scaled(year, df=par[3], mean = par[1], sd = par[2]),c)
  llog<- numeric(length = nc)
  for(i in 1:nc)
  {
    llog[i] <- sum(dpois(data[,i],lambda = lambda[,i],log = T))
  }
  llk <- sum(llog)
  return(llk)
}

### optimization functions of the loglikelihood functions####

Optim_tech <- function(time,tech) # there are 4 parameters to be optimized for each tech.
{
  optim(par = rep(1,4),data = tech,year=time,fn = llikTechs,method = "L-BFGS-B",lower = c(1,1,1,1),upper = c(Inf,Inf,20,Inf),hessian = T,control = list(fnscale = -1))
}

Optim_cluster <- function(time,cluster)
{
  num_par <- ncol(cluster) + 3 #number of parmeters in each cluster
  optim(par = rep(1,num_par),data = cluster,year=time,fn = llikCluster,method = "L-BFGS-B",lower = 1,upper = c(Inf,Inf,20,Inf),hessian = T,control = list(fnscale = -1))
}


##### Functions for Gibbs sampling ######
llikmu <- function(sigma,data,year,par,df,c)         # loglik. function for mean of a Student-t 
{
  nc <- ncol(data)                                   #number of techs with in a cluster
  lambda<-outer(dt.scaled(year, df=df, mean = par, sd = sigma),c)
  llog <- numeric(length = length(c))
  for(i in 1:length(c))
  {
    llog[i] <- sum(dpois(data[,i],lambda = lambda[,i],log = T)/log(c[i]))
  }
  
  llk <- sum(llog)
  return(llk)
}


lliksigma<- function(par,data,year,mu,df,c)         # loglik. function of a sigma in Student-t 
{
  nc <- ncol(data)                                  #number of techs with in a cluster
  lambda<-outer(dt.scaled(year, df=df, mean = mu, sd = par),c)
  llog <- numeric(length = length(c))
  for(i in 1:length(c))
  {
    llog[i] <- sum(dpois(data[,i],lambda = lambda[,i],log = T)/log(c[i]))
  }

  llk <- sum(llog)
  return(llk)
}

mu_prior <- function(par,mu_prior,sigma_prior) # prior distribution of mean
{
  lprior <- sum(dnorm(par,mean=mu_prior,sd=sigma_prior,log = T))
  return(lprior)
}


mu_post <- function(data,year,par,sigma,c,df,mu_prior,sigma_prior) # posterior of mean given all other par.
{
  return(llikmu(sigma=sigma,data=data,year=year,par=par,c=c,df=df) + mu_prior(par=par,mu_prior,sigma_prior))
}


postSampleDF<-function(data,year,mu,sigma,c,v) # unormilized posterior of df in a grid of points(1 to 20)
  # Purpose: Simulates on draw from the full conditional posterior of DF for a given cluster with nTechPerClust technologies
  #
  # Input:
  #     data  (nYears-by-nTechPerClust)     Time series of counts over time for the nTechPerClust technologies
  #     year  (nYears-dim vector)           Vector of years for the data
  #     mu    (scalar)                      Mean in the cluster
  #     sigma (scalar)                      Standard deviation in the cluster
  #     c     (nTechPerClust-dim vector)    Constants for the technologies in the cluster
  #     v     (nGridPoints-dim vector)      Vector with gridpoints for the degrees of freedom
  #
  #     Output: (scalar)                    A posterior draw of the degrees of freedom
#    
#     Author: Mehrawi (merastat@gmail.com)
{
  nTechPerClust<-ncol(data)
  llog <- rep(0,length(v))
  for(i in 1:length(v))
  {
    lambda <-  outer(dt.scaled(year, df=v[i], mean = mu, sd = sigma),c)
    for(j in 1:nTechPerClust)
    {
      llog[i] <- llog[i] + sum(dpois(data[,j],lambda = lambda[,j],log = T))/log(c[j])
      
    }
  }
  a<- datanormalize(llog)                            # the variable a will prvent from overfloating as the numebrs are ig.
  unormailez_posterior <- exp(a)        # a is added both in the numerator and denumenator.
  normilized_constant <- sum(exp(a))
  normalized_posterior <- (unormailez_posterior)/normilized_constant
  cdf_df <- cumsum(normalized_posterior)
  inv_cdf <- which((runif(1) < cdf_df)==TRUE)
  df <- v[inv_cdf[1]]
  return(df)
}
###### Functions for the mixture model######

divide_data<-function(data,I,ncomp) # this function updates data based on the allocation indexes.
{
  d<-list()
  for(i in 1:ncomp)
  {
    d[[i]]<- data[which(I==i),]
  }
  return(d)
}


salloc <- function(s,data) #updates the allocation components based on the udated S(probabilties)
{
  n <- nrow(data) #number of techs
  alloc <- rep(0,n)
  for (i in 1:n)
  {
    alloc[i] <- which(s[i,] == 1)
  }
  return(alloc)
}

pi_sampler <- function(s,comp,alpha) # updates the pi values
{
  n<-colSums(s)
  result<-numeric(length=ncomp)
  for(i in 1:ncomp)
  {
  result[i]<-n[i]+alpha[i]
  }
  return(as.vector(rdirichlet(1,result)))
  
}

c_sampler <- function(year,data,ncomp,df,mu,sigma) # updates the C
{
  c<- matrix(0,nrow(data),ncomp)
  
  for(i in 1:nrow(data))
  {
    alpha <- sum(data[i,])
    for(j in 1:ncomp)
    {
      beta  <- sum(dt.scaled(year,df = df[j],mean=mu[j],sd=sigma[j]))
      c[i,j] <- rgamma(1,alpha,beta)
      
    }
    
  }
  return(c)
}

# update the s
s_sampler <- function(data,year,pi,c,df,mu,sigma,ncomp) #Updates the probabilties for component allocation
{
  
  S <- matrix(0,nrow= nrow(data),ncol= ncomp)
  theta <- vector("numeric",length = ncomp)
  a <- vector("numeric",length = ncomp)
  lambda <- matrix(0,nrow=ncol(data),ncol=ncomp)
  for(i in 1:nrow(data)) 
  {
    for(j in 1:length(pi))
    { 
      lambda[,j] <- c[i]*(dt.scaled(year,df=df[j],mean=mu[j],sd=sigma[j]))
      theta[j] <- log(pi[j])+sum(dpois(data[i,],lambda = lambda[,j],log = T))
      
    }
    a<- -max(theta)
    #theta <- theta/a
    S[i,] <- t(rmultinom(1, size = 1 , prob = exp(theta+a)/sum(exp(theta+a))))
  }
  return(S)
}
