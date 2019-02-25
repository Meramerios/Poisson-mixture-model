source("C:/Users/MERA/Desktop/My_R_Codes/Library.R")

#### Input to the mixture Model#####
ncomp <- 3 # number of components
theta <- vector("numeric",length = ncomp) # pobalitiy of observations clustered to each component
mixdata <- t(data)
iter <- 2000
S <- t(rmultinom(ncol(data), size = 1 , prob = rep(1/ncomp,ncomp)))
par<-nrow(mixdata)+(ncomp*4) # Total number of parameters (the MUs,Dfs,SIGMAs and PIs)
parms3_2 <- matrix(0,par,nrow = iter)
mu <- rep(MleCluster$par[1],ncomp) 
sigma <- rep(MleCluster$par[2],ncomp) 
df <- rep(MleCluster$par[3],ncomp)
alphas <- rep(10,ncomp)           # prior on dirichilet distribution.
muprior <- c(10,15,20,20,25,30)   # priors on mean
sprior <- rep(2.55,ncomp)         # prior  on the sigma.
##End of all inputs### 

### Mixture Model ####

   for(i in 1:iter)
  {
    I <- salloc(S,mixdata)
    pi <- pi_sampler(S,ncomp,alphas)
    c <- c_sampler(year,mixdata,ncomp,df,mu,sigma)
    cc <- vector("numeric",length = nrow(mixdata))
     for(j in 1:nrow(mixdata))
    { 
      k <- which(S[j, ]==1)
      cc[j] <- c[j,k]
    }

    newdata <- divide_data(data=mixdata,I=I,ncomp = ncomp)
    
 
    for(d in 1:ncomp) ### update degrees of freedom
    {
      data_df<-unlist(newdata[[d]])
      if(length(cc[which(I==d)])==1)
        data_df<-as.matrix(data_df)
      else
        data_df<-as.matrix(t(data_df))
      df[d]<- postSampleDF(data = data_df,year = year,mu=mu[d],sigma = sigma[d],c=cc[which(I==d)],v=v)
      
    }
    
    ## update means
    
    for(m in 1:ncomp)
    {
      data_df<-unlist(newdata[[m]])
      if(length(cc[which(I==m)])==1)
        data_df<-as.matrix(data_df)
      else
        data_df<-as.matrix(t(data_df))
      mu_s<- (MfU.Sample.Run(mu[m],mu_post,nsmp = 1,sigma=sigma[m],data=t(newdata[[m]]),year=year,df=df[m],c=cc[which(I==m)],mu_prior=muprior[m],sigma_prior=sprior[m]))[1]
      mu[m]<-mu_s
    }
   
    ## Update sigmas
    for(s in 1:ncomp)
    {
      data_df<-unlist(newdata[[s]])
      if(length(cc[which(I==s)])==1)
        data_df<-as.matrix(data_df)
      else
        data_df<-as.matrix(t(data_df))
      sigma_s<-(MfU.Sample.Run(sigma[s],lliksigma,nsmp = 1,data=t(newdata[[s]]),year = year,mu=mu[s],c=cc[which(I==s)],df=df[s]))[1]
      sigma[s]<-sigma_s
    }
    
    #update S 
    S <- s_sampler(data=mixdata,year=year,pi=pi,c=cc,df=df,mu=mu,sigma=sigma,ncomp=ncomp)
    parms3_2[i,]<-c(cc,mu,sigma,df,pi)
  }






