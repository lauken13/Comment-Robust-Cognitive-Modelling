########################################################################
# Prior predictive check for non-hierachical model #####################
########################################################################
#Simulate the number of pumps for different probs of popping
set.seed(9458)
prior_1 <- data.frame()
for (p in seq(.05,.5,.05)){
  print(p)
  for (ps in 1:200){
    print(ps)
    gplus <- runif(1,0,10)
    beta <- runif(1,0,10)
    omega <- -gplus/log(1-p)
    for(j in 1:30){
      d=0
      d_tmp =1
      k=0
      pop=0
      while(k <2000 & (d_tmp==1 & pop==0)){
        theta <- 1/(1 + exp(beta*(k-omega)))
        iter <- j
        d_tmp <- rbinom(1,1,theta)
        d <- d+d_tmp
        k <- k+1
        pop=rbinom(1,1,p)
      }
      prior_1 <- rbind(prior_1,data.frame(npump = d, iter=iter,prob_burst=p,ps=ps, upper=u))
    }
  }
}

summary_prior_1 <- prior_1 %>%
  group_by(prob_burst,ps,upper)%>%
  summarize(mean_npump = mean(npump)) %>%
  ungroup()

saveRDS(prior_1,"results/ppc_nh_adjppc.RDS")

plot_1 <- ggplot(summary_prior_1, aes(x=as.factor(prob_burst),y=mean_npump))+ geom_violin(fill="navy")+
  theme_bw() +geom_hline(yintercept=10,colour="#bd0026")+ylab("Mean number of pumps per ps")+xlab("Probability of Bursting")
plot_1
ggsave("images/Meanpumps_Probburst_ppcadj.png",width=15,height=12,units="cm")

###Simulate using usual ppcs
prior_1b <- data.frame()
for (p in seq(.05,.5,.05)){
  print(p)
  for (ps in 1:200){
    print(ps)
    gplus <- runif(1,0,10)
    beta <- runif(1,0,10)
    omega <- -gplus/log(1-p)
    for(j in 1:30){
      d=0
      d_tmp =1
      k=0
      while(k <2000 & (d_tmp==1)){
        theta <- 1/(1 + exp(beta*(k-omega)))
        iter <- j
        d_tmp <- rbinom(1,1,theta)
        d <- d+d_tmp
        k <- k+1
      }
      prior_1b <- rbind(prior_1b,data.frame(npump = d, iter=iter,prob_burst=p,ps=ps, upper=u))
    }
  }
}

summary_prior_1b <- prior_1b %>%
  group_by(prob_burst,ps,upper)%>%
  summarize(mean_npump = mean(npump)) %>%
  ungroup()

saveRDS(prior_1b,"results/ppc_nh_ppc.RDS")

plot_1b <- ggplot(summary_prior_1b, aes(x=as.factor(prob_burst),y=mean_npump))+ geom_violin(fill="navy")+
  theme_bw() +geom_hline(yintercept=10,colour="#bd0026")+ylab("Mean number of pumps per ps")+xlab("Probability of Bursting")
plot_1b
ggsave("images/Meanpumps_Probburst_ppc.png",width=15,height=12,units="cm")




#######################################################################################
#####                   Hierarchical model                                        #####
#######################################################################################

#Now simulate for a hierachical model. with adapted PPCs
prior_2 <- data.frame()
alc_tmp <- c('Sober','Tipsy',"Drunk")
for (p in seq(0.05,0.50,.05)){  
    u = 10
    for (ps in 1:200){
      print(ps)
      mu_gplus <- abs(rnorm(1,0,u))
      sigma_gplus <- abs(rnorm(1,0,u))
      mu_beta <- abs(rnorm(1,0,u))
      sigma_beta <- abs(rnorm(1,0,u))
      for (i in 1:3){
        gplus <- abs(rnorm(1,mu_gplus,sigma_gplus))
        beta <- abs(rnorm(1,mu_beta,sigma_beta))
        omega <- -gplus/log(1-p)
        for(j in 1:30){
          d=0
          d_tmp =1
          k=0
          pop=0
          while(k <2000 & (d_tmp==1 & pop==0)){
            theta <- 1/(1 + exp(beta*(k-omega)))
            iter <- j
            d_tmp <- rbinom(1,1,theta)
            d <- d+d_tmp
            k <- k+1
            pop=rbinom(1,1,p)
          }
          prior_2 <- rbind(prior_2,data.frame(npump = d, iter=iter,prob_burst=p,Alc=alc_tmp[i],ps=ps, upper=u))
        }
      }  
    }
  }
}
prior_2$prob_burst = prior_2$prob_burst*100



ggplot(prior_2, aes(x=as.factor(prob_burst),y=npump,fill=Alc))+ geom_violin(alpha=.3)+
  theme_bw() +geom_hline(yintercept=10,colour="#bd0026")+ylab("Number of pumps")+xlab("Upper bound of Uniform")

saveRDS(prior_2,"results/ppc_h_ppcadj.RDS")

summary_pumps_within <- prior_2 %>%
  group_by(ps,prob_burst,Alc,upper)%>%
  mutate(sd_cond=sd(npump),mean_cond=mean(npump))%>%
  ungroup()%>%
  group_by(ps,prob_burst,upper)%>%
  mutate(sd_cond_comp=sd_cond[Alc=="Sober"]-sd_cond,
         mean_cond_comp=mean_cond[Alc=="Sober"]-mean_cond)


ggplot(summary_pumps_within[summary_pumps_within$Alc!="Sober",], aes(x=as.factor(prob_burst/100),y=mean_cond_comp,fill=Alc))+ geom_violin(alpha=.9)+
  theme_bw() +ylab("Difference in mean")+xlab("Probability of balloon popping")+
  scale_fill_manual(values=c("#c7e9b4","#41b6c4","#225ea8"))+
  theme(legend.position = "bottom")
ggsave('images/hier_difmean_ppcadj.png',width=15,height=12,units="cm")
ggplot(summary_pumps_within[summary_pumps_within$Alc!="Sober",], aes(x=as.factor(prob_burst/100),y=sd_cond_comp,fill=Alc))+ geom_violin(alpha=.9)+
  theme_bw() +ylab("Difference in sd")+xlab("Probability of balloon popping")+
  scale_fill_manual(values=c("#c7e9b4","#41b6c4","#225ea8"))+
  theme(legend.position = "none")
ggsave('images/hier_difvar_ppcadj.png',width=15,height=12,units="cm")

#Now simulate for a hierachical model with PPCs
prior_2b <- data.frame()
alc_tmp <- c('Sober','Tipsy',"Drunk")
for (p in seq(0.05,0.50,.05)){  
  u = 10
  for (ps in 1:200){
    print(ps)
    mu_gplus <- abs(rnorm(1,0,u))
    sigma_gplus <- abs(rnorm(1,0,u))
    mu_beta <- abs(rnorm(1,0,u))
    sigma_beta <- abs(rnorm(1,0,u))
    for (i in 1:3){
      gplus <- abs(rnorm(1,mu_gplus,sigma_gplus))
      beta <- abs(rnorm(1,mu_beta,sigma_beta))
      omega <- -gplus/log(1-p)
      for(j in 1:30){
        d=0
        d_tmp =1
        k=0
        while(k <2000 & (d_tmp==1)){
          theta <- 1/(1 + exp(beta*(k-omega)))
          iter <- j
          d_tmp <- rbinom(1,1,theta)
          d <- d+d_tmp
          k <- k+1
        }
        prior_2b <- rbind(prior_2b,data.frame(npump = d, iter=iter,prob_burst=p,Alc=alc_tmp[i],ps=ps, upper=u))
      }
    }  
  }
}

prior_2b$prob_burst = prior_2b$prob_burst*100

ggplot(prior_2b, aes(x=as.factor(prob_burst),y=npump,fill=Alc))+ geom_violin(alpha=.3)+
  theme_bw() +geom_hline(yintercept=10,colour="#bd0026")+ylab("Number of pumps")+xlab("Upper bound of Uniform")

saveRDS(prior_2b,"results/ppc_h_ppc.RDS")

summary_pumps_within <- prior_2b %>%
  group_by(ps,prob_burst,Alc,upper)%>%
  mutate(sd_cond=sd(npump),mean_cond=mean(npump))%>%
  ungroup()%>%
  group_by(ps,prob_burst,upper)%>%
  mutate(sd_cond_comp=sd_cond[Alc=="Sober"]-sd_cond,
         mean_cond_comp=mean_cond[Alc=="Sober"]-mean_cond)


ggplot(summary_pumps_within[summary_pumps_within$Alc!="Sober",], aes(x=as.factor(prob_burst),y=mean_cond_comp,fill=Alc))+ geom_violin(alpha=.9)+
  theme_bw() +ylab("Difference in mean")+xlab("Probability of balloon popping")+
  scale_fill_manual(values=c("#d4b9da","#980043"))+
  theme(legend.position = "bottom")
ggsave('images/hier_difmean_ppc.png',width=15,height=12,units="cm")
ggplot(summary_pumps_within[summary_pumps_within$Alc!="Sober",], aes(x=as.factor(prob_burst),y=sd_cond_comp,fill=Alc))+ geom_violin(alpha=.9)+
  theme_bw() +ylab("Difference in sd")+xlab("Probability of balloon popping")+
  scale_fill_manual(values=c("#d4b9da","#980043"))+
  theme(legend.position = "none")
ggsave('images/hier_difvar_ppc.png',width=15,height=12,units="cm")

#######################################################################################
####                      Thinking about model comparison                         #####
#######################################################################################

ntrials <- 90  # Number of trials for the BART

Data <- list(
  matrix(data=as.numeric(as.matrix(read.table("Explosive BART/data/GeorgeSober.txt"))[-1, ]),
         ntrials, 8),
  matrix(data=as.numeric(as.matrix(read.table("Explosive BART/data/GeorgeTipsy.txt"))[-1, ]),
         ntrials, 8),
  matrix(data=as.numeric(as.matrix(read.table("Explosive BART/data/GeorgeDrunk.txt"))[-1, ]),
         ntrials, 8)
)

nconds <- length(Data)                         
cash   <- npumps <- matrix (, nconds, ntrials)  # Cashes and nr. of pumps
d      <- array (-99, c(nconds, ntrials, 30))      # Data in binary format

for (i in 1:nconds) {
  cash[i, ]   <- (Data[[i]][, 7] != 0) * 1  # Cash or burst?
  npumps[i, ] <- Data[[i]][, 6]             # Nr. of pumps
  
  for (j in 1:ntrials) {
    if (npumps[i, j] > 0) {d[i, j, 1:npumps[i, j]] <- rep(0, npumps[i, j])}
    if (cash[i, j] == 1) {d[i, j, (npumps[i, j] + 1)] <- 1}
  }
}
options <- cash + npumps  # Nr. of decision possibilities

p_balloon <- Data[[1]][,4]/100
p <- as.numeric(levels(as.factor(Data[[1]][,4]/100)))
p_balloon <- as.numeric(as.factor(Data[[1]][,4]/100))
# To be passed on to Stan:
library(rstan)
library(bridgesampling)
library(loo)
df_bf <- expand.grid(u = 1:10, log_bf =0,loo10=0,loo_sd=0, k=1:10,time_1=-1,warm_1=-1,time_0=-1,warm_0=-1)
hier_model <- stan_model(file= "Explosive BART/code/BART_hierachical.stan")
non_hier_model <- stan_model(file= "Explosive BART/code/BART_not_hierachical.stan")
for (reps in 1:nrow(df_bf)){
  print(reps)
  dat <- list(nconds=nconds, ntrials=ntrials, p=p, options=options, d=d,
              up_prior=df_bf$u[reps],p_balloon = p_balloon) 
  
  fit_1 <- sampling(data=dat, 
                    object = hier_model,
                    chains=4, cores=4, control = list(adapt_delta = .99,max_treedepth=12),iter=500, init=1)
  
  fit_0 <- sampling(data=dat, 
                    object = non_hier_model,
                    chains=4, cores=4, control = list(adapt_delta = .99,max_treedepth=12),iter=500 )
  #BF 
  bs_1 <- bridge_sampler(fit_1)
  bs_0 <- bridge_sampler(fit_0)
  
  #LOO
  log_lik_1 <- extract_log_lik(fit_1, merge_chains = FALSE)
  log_lik_1narm <- log_lik_1[,,!is.na(log_lik_1[1,1,])]
  r_eff <- relative_eff(exp(log_lik_1narm)) 
  loo_1 <- loo(log_lik_1narm, r_eff = r_eff, cores = 2)
  
  log_lik_0 <- extract_log_lik(fit_0, merge_chains = FALSE)
  log_lik_0narm <- log_lik_0[,,!is.na(log_lik_0[1,1,])]
  r_eff <- relative_eff(exp(log_lik_0narm)) 
  loo_0 <- loo(log_lik_0narm, r_eff = r_eff, cores = 2)
  comp <- compare(loo_1, loo_0)
  
  df_bf$log_bf[reps] <- log(bf(bs_1,bs_0)$bf)
  df_bf$loo10[reps] <- comp[1]
  df_bf$loo_sd[reps] <- comp[2]
  df_bf$time_1[reps] <- max(apply(get_elapsed_time(fit_1),1,sum))
  df_bf$warm_1[reps] <- max(get_elapsed_time(fit_1)[1,])
  df_bf$time_0[reps] <- max(apply(get_elapsed_time(fit_0),1,sum))
  df_bf$warm_0[reps] <- max(get_elapsed_time(fit_0)[1,])
  print(df_bf[reps,])
}

ggplot(df_bf, aes(x=u, y=-loo10, group = k))+
  geom_ribbon(data = df_bf,mapping = aes(ymin = -loo10 - loo_sd, ymax = -loo10 + loo_sd,group = k), fill = "navy",alpha=.2)+
  geom_line()+theme_bw() + ylab("Preference for H1 by LOO")+xlab("Upper bound of uniform")
ggsave("Explosive BART/Images/Gdat_LOO.png",width=15, height=12,units="cm")
ggplot(df_bf, aes(x=u, y=loo_sd, group = k))+xlab("Upper bound of uniform")+
  geom_line()+theme_bw() + ylab("Variance of LOO")+ggtitle("George data")
ggsave("Explosive BART/Images/Gdat_LOO_var.png",width=15, height=12,units="cm")
ggplot(df_bf, aes(x=u, y=log_bf, group = k))+xlab("Upper bound of uniform")+
  geom_line()+theme_bw() + ylab("BF in favour of H1")
ggsave("Explosive BART/Images/Gdat_BF.png",width=15, height=12,units="cm")

ggplot()+
  geom_line(df_bf,mapping=aes(x=u, y=time_1, group = k),colour="red")+
  geom_line(df_bf,mapping=aes(x=u, y=time_0, group = k),colour="navy blue")+
  theme_bw() + ylab("Time to run")+ggtitle("George data")

ggsave("Explosive BART/Images/Gdat_runtime.png")

saveRDS(df_bf,'Explosive BART/Results/BF_Georgedata.RDS')

#################################################################################################################
#  What if we permuted the data set?                                                                        ####
#################################################################################################################

p_balloon <- Data[[1]][,4]/100
p <- as.numeric(levels(as.factor(Data[[1]][,4]/100)))
p_balloon <- as.numeric(as.factor(Data[[1]][,4]/100))
ntrials <- 90  # Number of trials for the BART

df_bf_null <- expand.grid(u = 1:10, log_bf =0,loo10=0,loo_sd=0, k=1:10,time_1=-1,warm_1=-1,time_0=-1,warm_0=-1)
hier_model <- stan_model(file= "Explosive BART/code/BART_hierachical.stan")
non_hier_model <- stan_model(file= "Explosive BART/code/BART_not_hierachical.stan")
for (reps in 1:nrow(df_bf_null)){
  print(reps)
  if(df_bf_null$u[reps] == 1){
    Data <- list(
      matrix(data=as.numeric(as.matrix(read.table("Explosive BART/data/GeorgeSober.txt"))[-1, ]),
             ntrials, 8),
      matrix(data=as.numeric(as.matrix(read.table("Explosive BART/data/GeorgeTipsy.txt"))[-1, ]),
             ntrials, 8),
      matrix(data=as.numeric(as.matrix(read.table("Explosive BART/data/GeorgeDrunk.txt"))[-1, ]),
             ntrials, 8)
    )
    #randomly permute the data:
    Data_alt <- Data
    permute1 <- c(sample(1:30,30),sample(31:60,30),sample(61:90),30)
    permute2 <- c(sample(1:30,30),sample(31:60,30),sample(61:90),30)
    permute3 <- c(sample(1:30,30),sample(31:60,30),sample(61:90),30)
    for(i in 1:90){
      samp_cond <- sample(1:3,3)
      Data_alt[[samp_cond[1]]][i,] <- Data_alt[[1]][permute1[i],]
      Data_alt[[samp_cond[2]]][i,] <- Data_alt[[2]][permute2[i],]
      Data_alt[[samp_cond[3]]][i,] <- Data_alt[[3]][permute3[i],]
    }
    nconds <- length(Data_alt)                         
    cash   <- npumps <- matrix (, nconds, ntrials)  # Cashes and nr. of pumps
    d      <- array (-99, c(nconds, ntrials, 30))      # Data in binary format
    
    for (i in 1:nconds) {
      cash[i, ]   <- (Data_alt[[i]][, 7] != 0) * 1  # Cash or burst?
      npumps[i, ] <- Data_alt[[i]][, 6]             # Nr. of pumps
      
      for (j in 1:ntrials) {
        if (npumps[i, j] > 0) {d[i, j, 1:npumps[i, j]] <- rep(0, npumps[i, j])}
        if (cash[i, j] == 1) {d[i, j, (npumps[i, j] + 1)] <- 1}
      }
    }
    options <- cash + npumps  # Nr. of decision possibilities
    
  }
  dat <- list(nconds=nconds, ntrials=ntrials, p=p, options=options, d=d,
              up_prior=df_bf_null$u[reps],p_balloon = p_balloon) 
  
  fit_1 <- sampling(data=dat, 
                    object = hier_model,
                    chains=4, cores=4, control = list(adapt_delta = .99,max_treedepth=12),iter=500, init=1)
  
  fit_0 <- sampling(data=dat, 
                    object = non_hier_model,
                    chains=4, cores=4, control = list(adapt_delta = .99,max_treedepth=12),iter=500 )
  #BF 
  bs_1 <- bridge_sampler(fit_1)
  bs_0 <- bridge_sampler(fit_0)
  
  #LOO
  log_lik_1 <- extract_log_lik(fit_1, merge_chains = FALSE)
  log_lik_1narm <- log_lik_1[,,!is.na(log_lik_1[1,1,])]
  r_eff <- relative_eff(exp(log_lik_1narm)) 
  loo_1 <- loo(log_lik_1narm, r_eff = r_eff, cores = 2)
  
  log_lik_0 <- extract_log_lik(fit_0, merge_chains = FALSE)
  log_lik_0narm <- log_lik_0[,,!is.na(log_lik_0[1,1,])]
  r_eff <- relative_eff(exp(log_lik_0narm)) 
  loo_0 <- loo(log_lik_0narm, r_eff = r_eff, cores = 2)
  comp <- compare(loo_1, loo_0)
  
  df_bf_null$log_bf[reps] <- log(bf(bs_1,bs_0)$bf)
  df_bf_null$loo10[reps] <- comp[1]
  df_bf_null$loo_sd[reps] <- comp[2]
  df_bf_null$time_1[reps] <- max(apply(get_elapsed_time(fit_1),1,sum))
  df_bf_null$warm_1[reps] <- max(get_elapsed_time(fit_1)[1,])
  df_bf_null$time_0[reps] <- max(apply(get_elapsed_time(fit_0),1,sum))
  df_bf_null$warm_0[reps] <- max(get_elapsed_time(fit_0)[1,])
  print(df_bf_null[reps,])
}

ggplot(df_bf_null, aes(x=u, y=-loo10, group = k))+
  geom_ribbon(data = df_bf_null,mapping = aes(ymin = -loo10 - loo_sd, ymax = -loo10 + loo_sd,group = k), fill = "navy",alpha=.2)+
  geom_line()+theme_bw() + ylab("Preference for H1 by LOO")+xlab("Upper bound of uniform")
ggsave("Explosive BART/Images/Permuted_Gdat_LOO.png",width=15, height=12,units="cm")
ggplot(df_bf_null, aes(x=u, y=loo_sd, group = k))+
  geom_line()+theme_bw() + ylab("Variance of LOO")+xlab("Upper bound of uniform")
ggsave("Explosive BART/Images/Permuted_Gdat_LOO_var.png",width=15, height=12,units="cm")
ggplot(df_bf_null, aes(x=u, y=log_bf, group = k))+
  geom_line()+theme_bw() + ylab("BF in favour of H1")+xlab("Upper bound of uniform")
ggsave("Explosive BART/Images/Perumted_Gdat_BF.png",width=15, height=12,units="cm")

ggplot()+
  geom_line(df_bf_null,mapping=aes(x=u, y=time_1, group = k),colour="red")+
  geom_line(df_bf_null,mapping=aes(x=u, y=time_0, group = k),colour="navy blue")+
  theme_bw() + ylab("Time to run")+ggtitle("Permuted George data")

ggsave("Explosive BART/Images/Permuted_Gdat_runtime.png")

saveRDS(df_bf_null,'Explosive BART/Results/BF_GeorgedataPermute.RDS')


################################################################################
# Prior predictive check for non-hierachical model using different priors ######
################################################################################
#Simulate the number of pumps for different probs of popping
prior_3 <- data.frame()
p=.1
for (u in 1:10){
  print(u)
  for (ps in 1:200){
    print(ps)
    gplus <- runif(1,0,u)
    beta <- runif(1,0,u)
    omega <- -gplus/log(1-p)
    for(j in 1:30){
      d=0
      d_tmp =1
      k=0
      pop=0
      while(k <2000 & (d_tmp==1 & pop==0)){
        theta <- 1/(1 + exp(beta*(k-omega)))
        iter <- j
        d_tmp <- rbinom(1,1,theta)
        d <- d+d_tmp
        k <- k+1
        pop=rbinom(1,1,p)
      }
      prior_3 <- rbind(prior_3,data.frame(npump = d, iter=iter,prob_burst=p,ps=ps, upper=u))
    }
  }
}

summary_prior_3 <- prior_3 %>%
  group_by(prob_burst,ps,upper)%>%
  summarize(mean_npump = mean(npump)) %>%
  ungroup()

saveRDS(prior_3,"results/ppc_nh_adjppc_upper.RDS")

plot_3 <- ggplot(prior_3, aes(x=as.factor(upper),y=npump))+ geom_violin(fill="navy")+
  theme_bw() +geom_hline(yintercept=10,colour="#bd0026")+ylab("Number of pumps per trial")+xlab("Upper bound of uniform prior")
plot_3
ggsave("images/Meanpumps_Probburst_ppcadj_upper.png",width=15,height=12,units="cm")

###Simulate using usual ppcs
prior_3b <- data.frame()
for (u in 1:10){
  print(u)
  for (ps in 1:200){
    print(ps)
    gplus <- runif(1,0,u)
    beta <- runif(1,0,u)
    omega <- -gplus/log(1-p)
    for(j in 1:30){
      d=0
      d_tmp =1
      k=0
      while(k <2000 & (d_tmp==1)){
        theta <- 1/(1 + exp(beta*(k-omega)))
        iter <- j
        d_tmp <- rbinom(1,1,theta)
        d <- d+d_tmp
        k <- k+1
      }
      prior_3b <- rbind(prior_3b,data.frame(npump = d, iter=iter,prob_burst=p,ps=ps, upper=u))
    }
  }
}

saveRDS(prior_3b,"results/ppc_nh_ppc_upper.RDS")

plot_3b <- ggplot(prior_3b, aes(x=as.factor(upper),y=npump))+ geom_violin(fill="navy")+
  theme_bw() +geom_hline(yintercept=10,colour="#bd0026")+ylab("Number of pumps per trial")+xlab("Upper bound of uniform prior")
plot_3b
ggsave("images/Meanpumps_Probburst_ppc_upper.png",width=15,height=12,units="cm")

