source("C:/Users/MERA/Desktop/My_R_Codes/Library.R")

# This file is written for EDA purposes
data <- read.csv("techs_data.csv",sep = ",",header = T) #data containing year in years (1 to 26) and number of patentsn of tec. i in each yeYear
year <- data[,1] 
data <- data[,-1]
year_for_plotting <- seq(1,70,1)
cluster <- cluster3 
tech <- data[,1]
MleTechs <- Optim_tech(time=year,tech=tech)
MleCluster <- Optim_cluster(time=year,cluster=cluster)
MleCluster$par

##A line plot for data an initial observation of data trend
plot_two_tecs<-as.data.frame(cbind(year,tech))
ggplot()+geom_line(data = plot_two_tecs,aes(x = year,y = data[,7])) + geom_point(data = plot_two_tecs,aes(x = year,y = data[,7]))+
  theme_classic() + labs(title = "Trend of patent documents for a technology(1970-216)") + ylab("number of patents")

# plot true vs fitted of MLE of techs i.

year_for_plotting<-seq(1,70,1)
fitted<-MleTechs$par[4]*dt.scaled(year_for_plotting,MleTechs$par[3],MleTechs$par[1],MleTechs$par[2])
plot(year_for_plotting,fitted,type="l")

#plot for bayesian(gibbs) estimates
fitted_g <- mean(gibbstec60[100:1000,1])*dt.scaled(year_for_plotting,df = mean(gibbstec60[100:1000,2]),mean = mean(gibbstec60[100:1000,4]), sd = mean(gibbstec60[100:1000,3]))
data8<-data[,22]
fitted<-cbind(data8,year)
p <- figure(title = "number of patents of one technoloy from 1970-2016",tools = NULL) %>%
  ly_lines(year,data8, data = fitted,color = "black") %>%
  ly_points(year,data8, data = fitted,color = "black",size = 4)%>%
  x_axis(label = "year")%>%y_axis(label = "number of patents")
p
# plot true vs fitted of MLE of clusters

fitted_cluster <- dt.scaled(year_for_plotting,df = MleCluster$par[3],mean = MleCluster$par[1], sd = MleCluster$par[2])
plot(year_for_plotting,fitted_cluster,type = "l",main = "cluster6",xlab="year",ylab="pdf",col=2,lwd=2)


## Scatter plots of all estimates using both bayesian and MLE
scatter<-matrix(NA,nrow = nrow(data),ncol = 3)
for (i in 1:nrow(data)) {
  scatter[i,1]<-Optim_tech(year,data[,i])$par[1]
  scatter[i,2]<-Optim_tech(year,data[,i])$par[2]
  scatter[i,3]<-Optim_tech(year,data[,i])$par[3]
  

}
scatter<-cbind(scatter,clusters)
scatter<- as.data.frame(scatter)
names(scatter) <- c("Mean","Sigma","df","cluster")
ggplot()+geom_point(data = scatter,aes(x=Sigma,y=Mean))+labs(x=expression(sigma),y=expression(mu))
         

