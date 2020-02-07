rm(list=ls())
fig_add <- "./ncov2019/fig/"
file_add <- "./ncov2019/"
source(paste0(file_add,"nCov2019_source.R"))
date <- Sys.Date()
casetype <- "Others"
begin <-  chron(dates.= "01/22/2020")

beta_ls <- seq(0.05,0.5,length.out = 10)
gamma_ls <- seq(0.05,0.5,length.out = 10)

pair_expand <- data.frame(expand.grid(beta_ls,gamma_ls))
names(pair_expand) <- c("beta","gamma")

nid <- nrow(pair_expand)
id_ls <- 1:nid

#server_num=1
server_num=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))# 1-100 
casename <- paste0(casetype,"_",server_num,"_",date)
beta_mean <- pair_expand[server_num,1]
beta_var <- 0.01
(lognorm_beta_parm <- lognorm.parm(beta_mean,beta_var))

gamma_mean <- pair_expand[server_num,2]
gamma_var <- 0.01
(lognorm_gamma_parm <- lognorm.parm(gamma_mean,gamma_var))

x <- get_nCov2019()
time(x)

N=1.386e9-58.5e6  ## Hubei population size: 58.5e6
#NI <- c(198,270,549,729,1052,1423,2714)
#earlyNI <- rep(c(1,4,6,7,8,9,14,18,21,29,30,33,35,37,40,41),c(9,5,2,1,1,1,1,1,1,1,1,1,1,4,1,12)) #Dec 1-Jan12
#NI <- as.numeric(c(summary(x)$confirm))[1:5] #Jan 13-17
#impute <- seq(NI[5],270,length.out = 4)[2:3] #Jan 18-19
#NI2 <- c(270,375,444,549,729,1052,1423,2714,3554,4903,5806,7153,9074,11177,13522,	16678) #Jan 20- Feb 4
#NI_complete <- c(NI,impute,NI2)
NI_complete <- c(21,65,155,332,601,953,1435,1919,2541,3260,4006,4748,5416,6164,7008,7769)
Y <- NI_complete/N #Jan13->Feb 4

T_prime <-length(Y)

#date_vec <- summary(x)$date

T_fin <- 200 # let the total time be 200

chron_ls <- chron(begin:(begin+T_fin))
################ MCMC ##########
model.string <-　paste0("
                       model{
                       for(t in 2:(T_prime+1)){
                       Km[t-1,1] <- -beta*theta[t-1,1]*theta[t-1,2]
                       Km[t-1,9] <- gamma*theta[t-1,2] 
                       Km[t-1,5] <- -Km[t-1,1]-Km[t-1,9]
                       
                       Km[t-1,2] <- -beta*(theta[t-1,1]+0.5*Km[t-1,1])*(theta[t-1,2]+0.5*Km[t-1,5]) 
                       Km[t-1,10] <- gamma*(theta[t-1,2]+0.5*Km[t-1,5])
                       Km[t-1,6] <- -Km[t-1,2]-Km[t-1,10]
                       
                       Km[t-1,3] <- -beta*(theta[t-1,1]+0.5*Km[t-1,2])*(theta[t-1,2]+0.5*Km[t-1,6]) 
                       Km[t-1,11] <- gamma*(theta[t-1,2]+0.5*Km[t-1,6])
                       Km[t-1,7] <- -Km[t-1,3]-Km[t-1,11]
                       
                       Km[t-1,4] <- -beta*(theta[t-1,1]+Km[t-1,3])*(theta[t-1,2]+Km[t-1,7])  
                       Km[t-1,12] <- gamma*(theta[t-1,2]+Km[t-1,7])
                       Km[t-1,8] <- -Km[t-1,4]-Km[t-1,12]
                       
                       alpha[t-1,1] <- theta[t-1,1]+(Km[t-1,1]+2*Km[t-1,2]+2*Km[t-1,3]+Km[t-1,4])/6
                       alpha[t-1,2] <- theta[t-1,2]+(Km[t-1,5]+2*Km[t-1,6]+2*Km[t-1,7]+Km[t-1,8])/6
                       alpha[t-1,3] <- theta[t-1,3]+(Km[t-1,9]+2*Km[t-1,10]+2*Km[t-1,11]+Km[t-1,12])/6
                       
                       theta[t,1:3] ~ ddirch(k*alpha[t-1,1:3])
                       Y[t-1] ~ dbeta(lambda*theta[t,2],lambda*(1-theta[t,2]))
                       }
                       theta[1,1] ~ dnorm(",0.95,",",1e4,")
                       theta[1,2] ~ dbeta(",1,",",N/NI_complete[1],")
                       theta[1,3] <- 1- theta[1,1]- theta[1,2]
                       gamma ~  dlnorm(",lognorm_gamma_parm[1],",",1/lognorm_gamma_parm[2],")
                       beta ~ dlnorm(",lognorm_beta_parm[1],",",1/lognorm_beta_parm[2],")
                       k ~  dgamma(2,0.0001)
                       lambda ~ dgamma(2,0.0001)
                       
                       }
                       ")

model.spec <- textConnection(model.string)

posterior <- jags.model(model.spec,data=list('Y'=Y,'T_prime'=T_prime),n.chains =nchain, n.adapt = nadapt)

update(posterior,nburnin) #burn-in

jags_sample <-jags.samples(posterior,c('alpha','theta','gamma','beta','Y','lambda','k'),n.iter=M*nchain,thin=thn)


####################

#dim(as.mcmc.list(jags2_sample$theta)[[1]])
png(paste0(fig_add,casename,"theta_p.png"), width = 700, height = 900)
plot(as.mcmc.list(jags_sample$theta)[[1]][,(1:3)*(T_prime+1)]) # posterior true porbabilities
dev.off()

theta_p<-array(as.mcmc.list(jags_sample$theta)[[1]],dim=c(len,T_prime+1,3))
#array(rbind(1:12,13:24),dim=c(2,4,3))
#array(rbind(1:12,13:24),dim=c(2,4,3))[,,1]
#nrow(theta_p)-len # should be 0
#ncol(theta_p)

png(paste0(fig_add,casename,"beta_p.png"), width = 700, height = 350)
plot(beta_p<-as.mcmc.list(jags_sample$beta)[[1]])
dev.off()

png(paste0(fig_add,casename,"gamma_p.png"), width = 700, height = 350)
plot(gamma_p<-as.mcmc.list(jags_sample$gamma)[[1]])
dev.off()

png(paste0(fig_add,casename,"lambda_p.png"), width = 700, height = 350)
plot(lambda_p<-as.mcmc.list(jags_sample$lambda)[[1]])
dev.off()

png(paste0(fig_add,casename,"k_p.png"), width = 700, height = 350)
plot(k_p<-as.mcmc.list(jags_sample$k)[[1]])
dev.off()

theta_p_mean <- apply(theta_p[,T_prime+1,],2,mean)
theta_p_ci <- as.vector(apply(theta_p[,T_prime+1,],2,quantile,c(0.025,0.5,0.975)))

beta_p_mean <- mean(beta_p)
beta_p_ci <- quantile(beta_p,c(0.025,0.5,0.975))

gamma_p_mean <- mean(gamma_p)
gamma_p_ci <- quantile(gamma_p,c(0.025,0.5,0.975))

lambda_p_mean <- mean(lambda_p)
lambda_p_ci <- quantile(lambda_p,c(0.025,0.5,0.975))

k_p_mean <- mean(k_p)
k_p_ci <- quantile(k_p,c(0.025,0.5,0.975))

#### Forecast ####
theta_pp <- array(NA,dim=c(len,T_fin-T_prime,3))
alpha_pp <- array(NA,dim=c(len,T_fin-T_prime,3))
Y_pp <- matrix(NA,nrow=len,ncol=T_fin-T_prime)

for(l in 1:len){
  thetat1 <- theta_p[l,T_prime+1,1]
  thetat2 <- theta_p[l,T_prime+1,2]
  thetat3 <- theta_p[l,T_prime+1,3]
  betat <- c(beta_p)[l]
  gammat <- c(gamma_p)[l]
  kt <- c(k_p)[l]
  lambdat <- c(lambda_p)[l]
  if(betat<0 |gammat<0 |thetat1<0 |thetat2<0 |thetat3<0) next
  for(t in 1:(T_fin-T_prime)){
    Km<-NULL
    Km[1] <- -betat*thetat1*thetat2
    Km[9] <- gammat*thetat2 
    Km[5] <- -Km[1]-Km[9]
    
    Km[2] <- -betat*(thetat1+0.5*Km[1])*(thetat2+0.5*Km[5]) 
    Km[10] <- gammat*(thetat2+0.5*Km[5])
    Km[6] <- -Km[2]-Km[10]
    
    Km[3] <- -betat*(thetat1+0.5*Km[2])*(thetat2+0.5*Km[6]) 
    Km[11] <- gammat*(thetat2+0.5*Km[6])
    Km[7] <- -Km[3]-Km[11]
    
    Km[4] <- -betat*(thetat1+Km[3])*(thetat2+Km[7])  
    Km[12] <- gammat*(thetat2+Km[7])
    Km[8] <- -Km[4]-Km[12]
    #if(is.na(thetat1)|is.na(thetat2)) stop("NA1")
    alpha_pp[l,t,1] <- thetat1+(Km[1]+2*Km[2]+2*Km[3]+Km[4])/6
    alpha_pp[l,t,2] <- thetat2+(Km[5]+2*Km[6]+2*Km[7]+Km[8])/6
    alpha_pp[l,t,3] <- thetat3+(Km[9]+2*Km[10]+2*Km[11]+Km[12])/6
    
    thetat_tmp <- rdirichlet(1,kt*c(alpha_pp[l,t,]))
    thetat1<-theta_pp[l,t,1] <- thetat_tmp[1]
    thetat2<-theta_pp[l,t,2] <- thetat_tmp[2]
    thetat3<-theta_pp[l,t,3] <- thetat_tmp[3]
    #if(is.na(thetat1)|is.na(thetat2)) stop("NA2")
    Y_pp[l,t] <- rbeta(1,lambdat*thetat2,lambdat*(1-thetat2))
    
  }
}

par(mfrow=c(1,1))

Y_band <- data.frame(t(apply(Y_pp,2,quantile,probs=c(0.025,0.5,0.975),na.rm=T)))
thetaI_band <- data.frame(t(apply(theta_p[,-1,2],2,quantile,probs=c(0.025,0.5,0.975),na.rm=T)))
Y_mean <- c(colMeans(Y_pp,na.rm = T))
thetaI_mean <- c(colMeans(theta_p[,-1,2],na.rm = T))

colnames(Y_band)<- c("lower", "median", "upper")
colnames(thetaI_band)<- c("lower", "median", "upper")
data_pre <- data.frame(time=1:T_prime,Y)
data_post <-data.frame(time=1:T_prime,thetaI_band)
data_fore <- data.frame(time=(T_prime+1):T_fin,Y_band,Y_mean)

data_comp<-data.frame(time=1:T_fin,rbind(thetaI_band ,Y_band), phase=c(rep('pre',nrow(thetaI_band)),rep('post',nrow(Y_band))),mean=c(thetaI_mean,Y_mean))

data_poly<-data.frame(y=c(thetaI_band$upper,rev(thetaI_band$lower),Y_band$upper,rev(Y_band$lower)),x=c(1:T_prime,T_prime:1,(T_prime+1):T_fin,T_fin:(T_prime+1)),phase=c(rep('pre',T_prime*2),rep('post',(T_fin-T_prime)*2)),value=c(rep(col2[1],T_prime*2),rep(col2[2],(T_fin-T_prime)*2)))

## First-order derivative check
thetaI_mat <- cbind(theta_p[,-1,2],theta_pp[,,2])
thetaS_mat <- cbind(theta_p[,-1,1],theta_pp[,,1])
first_dthetaI <- colMeans((thetaS_mat*thetaI_mat)*replicate(T_fin,c(beta_p))-thetaI_mat*replicate(T_fin,c(gamma_p)))

png(paste0(fig_add,casename,"deriv.png"), width = 700, height = 350)
plot(y=first_dthetaI,x=chron_ls,type='l',ylab="1st order derivative",main="Infection Proportion")
dev.off()

(first_order_stationary <- (1:T_fin)[which.min(first_dthetaI>0)]) # first order derivative=0
(second_order_stationary <- (1:T_fin)[which.max(first_dthetaI)])# first second order derivative=0
R0 <- beta_mean/gamma_mean
plot <- ggplot(data = data_poly, aes(x = x, y = y)) +geom_polygon(alpha = 0.5,aes(fill=value, group=phase)) +labs(title=substitute(paste("Infection Forecast of Case ",casety,"#",num,": ",beta[0],"=",v1,",",gamma[0], "=",v2," and R0=",v3), list(casety=casetype,num=server_num,v1=beta_mean,v2=gamma_mean,v3=R0)),x = "time", y = "infection %")+geom_line(data=data_comp,aes(x=time,y=median),color=2,linetype=2)+geom_vline(xintercept = T_prime,color="blue")+geom_vline(xintercept = second_order_stationary,color="orange")+geom_vline(xintercept = first_order_stationary,color="purple")+geom_line(data=data_comp,aes(x=time,y=mean),color="darkgray")+geom_point(data=data_pre,aes(x=time,y=Y))+theme_bw()+theme( plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")+scale_x_continuous(labels= as.character(chron_ls)[seq(1,T_fin,50)],breaks=seq(1,T_fin,50))


#plot

ggsave(paste0(fig_add,casename,"_forecast.png"))



out_table<-c(theta_p_mean,theta_p_ci,beta_p_mean,beta_p_ci,gamma_p_mean,gamma_p_ci,k_p_mean,k_p_ci,lambda_p_mean,lambda_p_ci,first_order_stationary,second_order_stationary)

#names(out_table)<-c("thetaS_p_mean","thetaI_p_mean","thetaR_p_mean","thetaS_p_ci_low","thetaS_p_ci_med","thetaS_p_ci_up","thetaI_p_ci_low","thetaI_p_ci_med","thetaI_p_ci_up","thetaR_p_ci_low","thetaR_p_ci_med","thetaR_p_ci_up","beta_p_mean","beta_p_ci_low","beta_p_ci_med","beta_p_ci_up","gamma_p_mean","gamma_p_ci_low","gamma_p_ci_med","gamma_p_ci_up","k_p_mean","k_p_ci_low","k_p_ci_med","k_p_ci_up","lambda_p_mean","lambda_p_ci_low","lambda_p_ci_med","lambda_p_ci_up","first_order_stationary","second_order_stationary")

write.csv(out_table,file=paste0(file_add,casename,"_summary.csv"))