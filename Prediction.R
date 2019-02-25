# This file contains a code for prediction purposes as an
# example i will predict for technology number 60 taking the first 10 years as base

test1<-data[,60]
first_few_data <- test1[1:10] # data of tech 1 at its merging stage
n<- length(first_few_data)
draws <- 2000
estimated_C <- matrix(NA,draws,ncomp)
predictions <- matrix(NA,(47-n),draws)
estimates <- matrix(NA,(47-n),ncomp)
lambda <- matrix(0,nrow=n,ncol=ncomp)
theta <- vector("numeric",length = ncomp)
pis <- vector("numeric",length = ncomp)

  for(j in 1:draws)
  {
    for (i in 1:ncomp) {
     
      estimated_C[j,i] <- rgamma(1,sum(first_few_data),sum(dt.scaled(year[1:n],df=parms6_2[j,135+i],
                                                                   mean = parms6_2[j,123+i],sd=parms6_2[j,129+i])))
    
  }


  pis <- abs(parms6_2[j,142:147])
  


for(t in 1:length(pis))
{ 
  lambda[,t] <- mean(estimated_C[j,t])*(dt.scaled(year[1:n],df=mean(parms6_2[j,135+t]),
                                                 mean=mean(parms6_2[j,123+t]),sd=mean(parms6_2[j,t+129])))
  theta[t] <- log(pis[t])+sum(dpois(first_few_data,lambda = lambda[,t],log = T))
  
}
a<--max(theta)
theta<-exp(theta+a)/sum(exp(theta+a))
theta2<-sample(1:ncomp,size = 1,prob = theta)

lambda2<- estimated_C[j,theta2]*(dt.scaled(year[(n+1):47],df=parms6_2[1,135+theta2],
                                              mean=parms6_2[1,123+theta2],sd=parms6_2[1,theta2+129]))

  predictions[,j]<-rpois((47-n),lambda = lambda2)
}
##### End ########

### Sample plots####
plot(0,0,xlim = c(0,50),ylim = c(300,2500),ylab = "number of patents",xlab = "year")
lines(first_few_data)
lines(((n+1):47),rowMeans(predictions)*0.8,col=3)
points(11:47,test1[11:47],col=2)
#lines(11:47,lambda2,col=3)
