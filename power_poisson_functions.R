#### frequentist power ####

# internal function #
power_freq<-function(r.alloc,n,lamE,theta)
{
  nP=r.alloc[1]*n # sample size in the placebo arm
  nR=r.alloc[2]*n # sample size in the reference arm
  nE=r.alloc[3]*n # sample size in the experimental arm
  
  lamE_null<-lamP + theta*(lamR-lamP) # value of lamE under null.
  
  mu_null=lamE_null-lamP
  
  mu=lamE-lamP
  mv=lamR-lamP
  
  sigma2u_null<-lamE_null/nE + lamP/nP
  sigma2u<-lamE/nE + lamP/nP
  
  sigma2v<-lamR/nR + lamP/nP
  
  rho_null<-(lamP/nP)/(sqrt(sigma2u_null*sigma2v))
  rho<-(lamP/nP)/(sqrt(sigma2u*sigma2v))
  
  alpha<-(-mv/sqrt(sigma2v))
  
  phi=dnorm(alpha,0,1)
  
  c<-1-pnorm(alpha,0,1)
  
  E1_null<-mu_null+sqrt(sigma2u_null)*rho_null*phi/c
  
  E1<-mu+sqrt(sigma2u)*rho*phi/c
  
  E2<-mv+sqrt(sigma2v)*phi/c
  
  V1_null<-sigma2u_null*(1+((rho_null^2)*alpha*phi/c)-(rho_null*phi/c)^2)
  
  V1<-sigma2u*(1+((rho^2)*alpha*phi/c)-(rho*phi/c)^2)
  
  V2<-sigma2v*(1-(phi/c)*((phi/c)-alpha))
  
  E12_null<-(sqrt(sigma2u_null*sigma2v)*rho_null*(alpha*phi+c)/c) + (sqrt(sigma2u_null)*(mv*rho_null*phi/c)) +
    (sqrt(sigma2v)*(mu_null*phi/c)) + mu_null*mv
  
  E12<-(sqrt(sigma2u*sigma2v)*rho*(alpha*phi+c)/c) + (sqrt(sigma2u)*(mv*rho*phi/c)) +
    (sqrt(sigma2v)*(mu*phi/c)) + mu*mv
  
  cov_null<-E12_null-E1_null*E2
  
  cov<-E12-E1*E2
  
  mw_th_null<-E1_null-theta*E2
  
  mw_th_alt<-E1-theta*E2
  
  vw_null<-V1_null+(theta^2)*V2-2*theta*cov_null
  
  vw_alt<-V1+(theta^2)*V2-2*theta*cov
  
  z<-qnorm(0.975,0,1)
  
  k<-mw_th_null + z*sqrt(vw_null)
  
  power<-1 - pnorm(((k-mw_th_alt)/sqrt(vw_alt)),0,1)
  
  return(power)
}

# Power curve for a specific theta #
power_theta_fn_freq<-function(theta)
{
  power1f<-NULL
  
  for(i in lamE)
  {
    power1f<-c(power1f,power_freq(r.alloc=c(1,1,1),n=n,lamE=i,theta=theta))
  }
  
  power1f<-round(power1f,4)
  
  return(power1f)
}

# sample size calculation #
samplesize_fn_freq<-function(r.alloc,lamE,theta)
{
  a<-1:1000
  
  power<-NULL
  
  for(i in a)
  {
    power<-c(power,power_freq(r.alloc=r.alloc,n=i,lamE=lamE,theta=theta))
  }
  n<-a[min(which(power>=0.8))]
  
  return(n)
}


#### fully bayesian power ####

# internal function #
power_fbayes<-function(r.alloc,n,lamE,theta)
{
  set.seed(123)
  
  aE<-0.5
  bE<-0.00001
  
  aR<-0.5
  bR<-0.00001
  
  aP<-0.5
  bP<-0.00001
  
  Est_Prob<-rep(0,n_star)
  
  nP=r.alloc[1]*n
  nR=r.alloc[2]*n
  nE=r.alloc[3]*n
  
  tE<-tR<-tP<-1
  
  count<-0
  
  ########## Test Procedure ##############
  for(i in 1:n_star)
  {
    count1=0
    count2=0
    
    xP=rpois(1,nP*lamP*tP)
    xR=rpois(1,nR*lamR*tR)
    xE=rpois(1,nE*lamE*tE)
    
    while (count1<T)
    {
      
      lamEgdata=rgamma(1,aE+xE,bE+nE*tE)
      lamRgdata=rgamma(1,aR+xR,bR+nR*tR)
      lamPgdata=rgamma(1,aP+xP,bP+nP*tP)
      
      dEPn=lamEgdata-lamPgdata
      dRPn=lamRgdata-lamPgdata
      
      if(dRPn>0)
      { 
        count1=count1+1
        
        if (dEPn/dRPn>theta)
          count2=count2+1
      }
    }
    
    Est_Prob[i]=count2/T
    
    if(Est_Prob[i]>p_star)
      count=count+1
  }
  
  power<-count/n_star
  return(power)
}   

# Power curve for a specific theta #
power_theta_fn_fbayes<-function(theta)
{
  power1f<-NULL
  
  for(i in lamE)
  {
    power1f<-c(power1f,power_fbayes(r.alloc=c(1,1,1),n=n,lamE=i,theta=theta))
  }
  
  power1f<-round(power1f,4)
  
  return(power1f)
}

# internal function#
#### fully bayesian power ####

# internal function#

power_fbayes<-function(r.alloc,n,lamE,theta)
{
  set.seed(123)
  
  Est_Prob<-rep(0,n_star)
  
  nP=r.alloc[1]*n
  nR=r.alloc[2]*n
  nE=r.alloc[3]*n
  
  tE<-tR<-tP<-1
  
  count<-0
  
  ########## Test Procedure ##############
  for(i in 1:n_star)
  {
    count1=0
    count2=0
    
    xP=rpois(1,nP*lamP*tP)
    xR=rpois(1,nR*lamR*tR)
    xE=rpois(1,nE*lamE*tE)
    
    while (count1<T)
    {
      aE<-21; bE<-1; aR<-21; bR<-1; aP<-1; bP<-1
      lamEgdata=rgamma(1,aE+xE,bE+nE*tE)
      lamRgdata=rgamma(1,aR+xR,bR+nR*tR)
      lamPgdata=rgamma(1,aP+xP,bP+nP*tP)
      
      dEPn=lamEgdata-lamPgdata
      dRPn=lamRgdata-lamPgdata
      
      if(dRPn>0)
      { 
        count1=count1+1
        
        if (dEPn/dRPn>theta)
          count2=count2+1
      }
    }
    
    Est_Prob[i]=count2/T
    
    if(Est_Prob[i]>p_star)
      count=count+1
  }
  
  power<-count/n_star
  return(power)
} 
# sample size calculation #
samplesize_fn_fbayes<-function(r.alloc,lamE,theta,a_max)
{
  a_min<-max(1,(a_max-40))
  a_max<- a_max+10
  a<-a_min:a_max
  
  power<-NULL
  
  for(i in a)
  {
    power<-c(power,power_fbayes(r.alloc=r.alloc,n=i,lamE=lamE,theta=theta))
  }
  n<-try(a[min(which(power>=0.8))],silent=T)
  if(class(n)=="try-error") n<-NA
  
  return(n)
}


#### fully bayesian informative power ####

# internal function #
powerfbayes.info<-function(r.alloc,n,lamE,theta)
{
  set.seed(123)
  
  aE<-70
  bE<-10
  
  aR<-70
  bR<-10
  
  aP<-10
  bP<-10
  
  Est_Prob<-rep(0,n_star)
  
  nP=r.alloc[1]*n
  nR=r.alloc[2]*n
  nE=r.alloc[3]*n
  
  tE<-tR<-tP<-1
  
  count<-0
  
  ########## Test Procedure ##############
  for(i in 1:n_star)
  {
    count1=0
    count2=0
    
    xP=rpois(1,nP*lamP*tP)
    xR=rpois(1,nR*lamR*tR)
    xE=rpois(1,nE*lamE*tE)
    
    while (count1<T)
    {
      
      lamEgdata=rgamma(1,aE+xE,bE+nE*tE)
      lamRgdata=rgamma(1,aR+xR,bR+nR*tR)
      lamPgdata=rgamma(1,aP+xP,bP+nP*tP)
      
      dEPn=lamEgdata-lamPgdata
      dRPn=lamRgdata-lamPgdata
      
      if(dRPn>0)
      { 
        count1=count1+1
        
        if (dEPn/dRPn>theta)
          count2=count2+1
      }
    }
    
    Est_Prob[i]=count2/T
    
    if(Est_Prob[i]>p_star)
      count=count+1
  }
  
  power<-count/n_star
  return(power)
} 

# Power curve for a specific theta for info bayes #
power.info_theta_fn<-function(theta)
{
  power1f<-NULL
  
  for(i in lamE)
  {
    power1f<-c(power1f,powerfbayes.info(r.alloc=c(1,1,1),n=n,lamE=i,theta=theta))
  }
  
  power1f<-round(power1f,4)
  
  return(power1f)
}

#### fully bayesian non-informative power ####

# internal function #
powerfbayes.noninfo<-function(r.alloc,n,lamE,theta)
{
  set.seed(123)
  
  aE<-7
  bE<-1
  
  aR<-7
  bR<-1
  
  aP<-1
  bP<-1
  
  Est_Prob<-rep(0,n_star)
  
  nP=r.alloc[1]*n
  nR=r.alloc[2]*n
  nE=r.alloc[3]*n
  
  tE<-tR<-tP<-1
  
  count<-0
  
  # Test Procedure #
  for(i in 1:n_star)
  {
    count1=0
    count2=0
    
    xP=rpois(1,nP*lamP*tP)
    xR=rpois(1,nR*lamR*tR)
    xE=rpois(1,nE*lamE*tE)
    
    while (count1<T)
    {
      
      lamEgdata=rgamma(1,aE+xE,bE+nE*tE)
      lamRgdata=rgamma(1,aR+xR,bR+nR*tR)
      lamPgdata=rgamma(1,aP+xP,bP+nP*tP)
      
      dEPn=lamEgdata-lamPgdata
      dRPn=lamRgdata-lamPgdata
      
      if(dRPn>0)
      { 
        count1=count1+1
        
        if (dEPn/dRPn>theta)
          count2=count2+1
      }
    }
    
    Est_Prob[i]=count2/T
    
    if(Est_Prob[i]>p_star)
      count=count+1
  }
  
  power<-count/n_star
  return(power)
} 

# Power curve for a specific theta for non-info bayes #
power.noninfo_theta_fn<-function(theta)
{
  power1f<-NULL
  
  for(i in lamE)
  {
    power1f<-c(power1f,powerfbayes.noninfo(r.alloc=c(1,1,1),n=n,lamE=i,theta=theta))
  }
  
  power1f<-round(power1f,4)
  
  return(power1f)
}

#### Approximate bayesian power ####

# internal function #
power_approxbayes<-function(r.alloc=r.alloc,n=i,lamE=lamE,theta=theta)
{
  nP=r.alloc[1]*n
  nR=r.alloc[2]*n
  nE=r.alloc[3]*n
  aE<-21; bE<-1; aR<-21; bR<-1; aP<-1; bP<-1
  me=aE/bE
  mr=aR/bR
  mp=aP/bP
  
  sigma2e=aE/(bE^2)
  sigma2r=aR/(bR^2)
  sigma2p=aP/(bP^2)
  
  sigma2_T<-lamE/nE + (theta^2)*lamR/nR + ((1-theta)^2)*lamP/nP 
  
  
  m_etaEP=me-mp
  m_etaRP=mr-mp
  
  sigma2_etaEP<-sigma2e + sigma2p
  
  sigma2_etaRP<-sigma2r + sigma2p 
  
  rho<-sigma2p/(sqrt(sigma2_etaRP*sigma2_etaEP)) 
  
  alpha<-(-m_etaRP/sqrt(sigma2_etaRP))
  
  phi=dnorm(alpha,0,1)
  
  c<-1-pnorm(alpha,0,1)
  
  delta_alpha=(phi/c)*((phi/c)-alpha)
  
  E1<-m_etaEP+sqrt(sigma2_etaEP)*rho*phi/c
  E2<-m_etaRP+sqrt(sigma2_etaRP)*phi/c
  
  V1<-sigma2_etaEP*(1+(rho^2*alpha*phi/c)-(rho*phi/c)^2)
  V2<-sigma2_etaRP*(1-delta_alpha)
  
  E12<-(sqrt(sigma2_etaEP*sigma2_etaRP)*rho*(alpha*phi+c)/c) + (sqrt(sigma2_etaEP)*m_etaRP*rho*phi/c) +
    (sqrt(sigma2_etaRP)*m_etaEP*phi/c) + m_etaEP*m_etaRP
  
  cov<-E12-E1*E2
  
  mu_star_th<-E1-theta*E2
  
  sigma2_star_nu<-V1+(theta^2)*V2-2*theta*cov
  
  mu_T_theta<-lamE - theta*lamR -(1-theta)*lamP
  
  z<-qnorm(1-p_star,0,1)
  
  point<-z*sqrt((1/sigma2_T) + (1/sigma2_star_nu))*sqrt(sigma2_T) + (mu_star_th/sigma2_star_nu)*(sqrt(sigma2_T)) + (mu_T_theta/sqrt(sigma2_T))
  
  power<-pnorm(point,0,1)
  
  return(power)
}

# sample size calculation #
samplesize_fn_approxbayes<-function(r.alloc,lamE,theta)
{
  a<-1:1000
  
  power<-NULL
  
  for(i in a)
  {
    power<-c(power,sample_pwr_fn_approxbayes(r.alloc=r.alloc,n=i,lamE=lamE,theta=theta))
  }
  power<-round(power,3)
  n<-a[min(which(power>=0.8))]
  
  return(n)
}


#### Fully bayesian power under non-conjugate prior ####

# internal function #
power_nonconjug<-function(N,lamE,theta,r0,p0,aE,bE)
{
  Est_Prob<-rep(0,n_star)
  
  nE<-N/3
  nR<-N/3
  nP<-N/3
  
  count<-0
  tE<-tR<-tP<-1
  
  for(i in 1:n_star)
  {
    count2=0
    
    xe<-rpois(1,nE*lamE*tE)
    xr<-rpois(1,nR*lamR*tR)
    xp<-rpois(1,nP*lamP*tP)
    
    model.string <-"
    
    model{
    
    xr~dpois(nR*lamR)
    xp~dpois(nP*lamP)
    
    u1~dbeta(alpha,beta)
    u2~dgamma(r,p)
    
    alpha<-10
    beta<-60
    
    lamP<-u1*u2
    lamR<-u2
    }"

    model.spec<-textConnection(model.string)
    
    load.module("bugs")
    load.module("mix")
    dat<-list(xp=xp,xr=xr,nP=nP,nR=nR,r=r0,p=p0)

    jags.inits<-list(u1=1/7,u2=7)
    parameters<-c("lamP","lamR")
    
    beta_binom<-jags.model(model.spec,data=dat,inits=jags.inits,n.chains=1,n.adapt=1000)
    
    update(beta_binom,100000)
    samples<-coda.samples(beta_binom,variable.names=parameters,n.iter=1000,thin=1)
    post_samples<-unlist(samples)
    
    lamp_post<-post_samples[1:1000]
    lamr_post<-post_samples[1001:2000]
    
    ae_post<-aE+xe
    be_post<-bE+nE*tE
    
    lame_post<-rgamma(T,ae_post,be_post)
    
    dEPn<-rep(0,T)
    dRPn<-rep(0,T)
    
    for(t in 1:T)
    {
      dEPn[t]=lame_post[t]-lamp_post[t]
      dRPn[t]=lamr_post[t]-lamp_post[t]
      
      if (dEPn[t]/dRPn[t]>theta)
        count2=count2+1
    }
    
    Est_Prob[i]=count2/T
    
    if(Est_Prob[i]>p_star)
      count=count+1
  }
  power<-count/n_star
}

# Power curve for different N and theta #
power_theta_fn_nonconjug<-function(N,theta)
{
  if(N==300 && theta==0.7)
  {
    aE=r0=210;bE=p0=30;
  }else if(N==300 && theta==0.6)
  {
    aE=r0=140;bE=p0=20;
  }else if(N==150 && theta==0.7)
  {
    aE=r0=140;bE=p0=20;
  }
  else{
    aE=r0=70;bE=p0=10;
  }
  
  power1f<-NULL
  
  for(i in lamE)
  {
    power1f<-c(power1f,power_nonconjug(N=N,lamE=i,theta=theta,r0,p0,aE,bE))
  }
  
  power1f<-round(power1f,4)
  
  return(power1f)
}

###################### Real data: Fully bayesian conjugate prior #########################

###################### 1 year #########################

# non-informative #
data.noninfo1<-function(theta)
{
  count1=0
  
  dEPn.vec<-NULL
  dRPn.vec<-NULL
  
  dataP<-rep(0:4,round(c(27.5/2,13.7/2,40/2,13.7/2,3.9/2)))
  dataR<-rep(0:3,round(c(73.9*46/100,17.4*46/100,6.5*46/100,2.2*46/100)))
  dataE<-rep(0:3,round(c(45.8*48/100,33.3*48/100,16.7*48/100,4.2*48/100)))
  
  xP=sum(dataP)
  xR=sum(dataR)
  xE=sum(dataE)
  tE<-tR<-tP<-1
  
  nP=length(dataP)
  nR=length(dataR)
  nE=length(dataE)
  while (count1<T)
  {
    aE<-0.5; bE<-0.00001; aR<-0.5; bR<-0.00001; aP<-0.5; bP<-0.00001
    
    lamEgdata=rgamma(1,aE+xE,bE+nE)
    lamRgdata=rgamma(1,aR+xR,bR+nR)
    lamPgdata=rgamma(1,aP+xP,bP+nP)
    
    dEPn<-lamPgdata-lamEgdata
    dRPn<-lamPgdata-lamRgdata
    
    if(dRPn>0)
    {
      count1=count1+1
      dEPn.vec=rbind(dEPn.vec,lamPgdata-lamEgdata)
      dRPn.vec=rbind(dRPn.vec,lamPgdata-lamRgdata)
    }
  }
  
  count2=0
  
  for(t in 1:T)
  {
    if (dEPn.vec[t]/dRPn.vec[t]>theta)
      count2=count2+1
  }
  
  Est_Prob=count2/T
  Est_Prob
  
  return(Est_Prob)
} 

# informative #
data.info1<-function(theta)
{
  count1=0
  
  dEPn.vec<-NULL
  dRPn.vec<-NULL
  
  dataP<-rep(0:4,round(c(27.5/2,13.7/2,40/2,13.7/2,3.9/2)))
  dataR<-rep(0:3,round(c(73.9*46/100,17.4*46/100,6.5*46/100,2.2*46/100)))
  dataE<-rep(0:3,round(c(45.8*48/100,33.3*48/100,16.7*48/100,4.2*48/100)))
  
  xP=sum(dataP)
  xR=sum(dataR)
  xE=sum(dataE)
  tE<-tR<-tP<-1
  
  nP=length(dataP)
  nR=length(dataR)
  nE=length(dataE)
  
  while (count1<T)
  {
    aE<-8; bE<-10; aR<-20; bR<-5; aP<-12; bP<-8
    
    lamEgdata=rgamma(1,aE+xE,bE+nE)
    lamRgdata=rgamma(1,aR+xR,bR+nR)
    lamPgdata=rgamma(1,aP+xP,bP+nP)
    
    dEPn<-lamPgdata-lamEgdata
    dRPn<-lamPgdata-lamRgdata
    
    if(dRPn>0)
    {
      count1=count1+1
      dEPn.vec=rbind(dEPn.vec,lamPgdata-lamEgdata)
      dRPn.vec=rbind(dRPn.vec,lamPgdata-lamRgdata)
    }
  }
  
  count2=0
  
  for(t in 1:T)
  {
    if (dEPn.vec[t]/dRPn.vec[t]>theta)
      count2=count2+1
  }
  
  Est_Prob=count2/T
  Est_Prob
  
  return(Est_Prob)
}

###################### 2 year #########################

# non-informative #
data.noninfo2<-function(theta)
{
  count1=0
  
  dEPn.vec<-NULL
  dRPn.vec<-NULL
  
  dataP<-rep(0:8,round(c(17.6/2,11.8/2,15.7/2,9.8/2,18/2,13.7/2,7.8/2,2/2,2/2)))
  dataR<-rep(0:4,round(c(47.8*46/100,39.1*46/100,8.7*46/100,2.2*46/100,2.2*46/100)))
  dataE<-rep(0:5,round(c(33.3*48/100,27.1*48/100,20.8*48/100,16.7*48/100,0*48/100,2.1*48/100)))
  
  xP=sum(dataP)
  xR=sum(dataR)
  xE=sum(dataE)
  tE<-tR<-tP<-1
  
  nP=length(dataP)
  nR=length(dataR)
  nE=length(dataE)
  
  while (count1<T)
  {
    aE<-0.5; bE<-0.00001; aR<-0.5; bR<-0.00001; aP<-0.5; bP<-0.00001
    
    lamEgdata=rgamma(1,aE+xE,bE+nE)
    lamRgdata=rgamma(1,aR+xR,bR+nR)
    lamPgdata=rgamma(1,aP+xP,bP+nP)
    
    dEPn<-lamPgdata-lamEgdata
    dRPn<-lamPgdata-lamRgdata
    
    if(dRPn>0)
    {
      count1=count1+1
      dEPn.vec=rbind(dEPn.vec,lamPgdata-lamEgdata)
      dRPn.vec=rbind(dRPn.vec,lamPgdata-lamRgdata)
    }
  }
  
  count2=0
  
  for(t in 1:T)
  {
    if (dEPn.vec[t]/dRPn.vec[t]>theta)
      count2=count2+1
  }
  
  Est_Prob=count2/T
  Est_Prob
  
  return(Est_Prob)
} 

# informative #
data.info2<-function(theta)
{
  count1=0
  
  dEPn.vec<-NULL
  dRPn.vec<-NULL
  
  dataP<-rep(0:8,round(c(17.6/2,11.8/2,15.7/2,9.8/2,18/2,13.7/2,7.8/2,2/2,2/2)))
  dataR<-rep(0:4,round(c(47.8*46/100,39.1*46/100,8.7*46/100,2.2*46/100,2.2*46/100)))
  dataE<-rep(0:5,round(c(33.3*48/100,27.1*48/100,20.8*48/100,16.7*48/100,0*48/100,2.1*48/100)))
  
  xP=sum(dataP)
  xR=sum(dataR)
  xE=sum(dataE)
  tE<-tR<-tP<-1
  
  nP=length(dataP)
  nR=length(dataR)
  nE=length(dataE)
  
  while (count1<T)
  {
    aE<-64; bE<-49.6; aR<-12; bR<-17; aP<-60; bP<-20.4
    
    lamEgdata=rgamma(1,aE+xE,bE+nE)
    lamRgdata=rgamma(1,aR+xR,bR+nR)
    lamPgdata=rgamma(1,aP+xP,bP+nP)
    
    dEPn<-lamPgdata-lamEgdata
    dRPn<-lamPgdata-lamRgdata
    
    if(dRPn>0)
    {
      count1=count1+1
      dEPn.vec=rbind(dEPn.vec,lamPgdata-lamEgdata)
      dRPn.vec=rbind(dRPn.vec,lamPgdata-lamRgdata)
    }
  }
  
  count2=0
  
  for(t in 1:T)
  {
    if (dEPn.vec[t]/dRPn.vec[t]>theta)
      count2=count2+1
  }
  
  Est_Prob=count2/T
  Est_Prob
  
  return(Est_Prob)
} 

###################### Real data: Fully bayesian non-conjugate prior #########################

###################### 1 year #########################

# non-informative #
data.nc.noninfo1<-function(theta)
{
  dataP<-rep(0:4,round(c(27.5/2,13.7/2,40/2,13.7/2,3.9/2)))
  dataR<-rep(0:3,round(c(73.9*46/100,17.4*46/100,6.5*46/100,2.2*46/100)))
  dataE<-rep(0:3,round(c(45.8*48/100,33.3*48/100,16.7*48/100,4.2*48/100)))
  
  xP=sum(dataP)
  xR=sum(dataR)
  xE=sum(dataE)
  tE<-tR<-tP<-1
  
  nP=length(dataP)
  nR=length(dataR)
  nE=length(dataE)
  
  aE<-0.8
  bE<-1
  
  model.string <-"
  model{
  
  xr~dpois(46*lamR)
  xp~dpois(50*lamP)
  
  u1~dbeta(alpha,beta)
  u2~dgamma(r,p)
  
  alpha<-1
  beta<-3.1
  
  r<-1.5
  p<-1
  
  lamR<-u1*u2
  lamP<-u2
  }"

  model.spec<-textConnection(model.string)
  
  load.module("bugs")
  load.module("mix")
  
  dat<-list(xP,xR)
  jags.inits<-list(u1=.244,u2=1.5)
  parameters<-c("lamP","lamR")
  
  beta_binom<-jags.model(model.spec,data=dat,inits=jags.inits,n.chains=1,n.adapt=1000)
  
  update(beta_binom,100000)
  samples<-coda.samples(beta_binom,variable.names=parameters,n.iter=50000,thin=50)
  post_samples<-unlist(samples)
  
  lamp_post<-post_samples[1:1000]
  lamr_post<-post_samples[1001:2000]
  
  ae_post<-aE+xE
  be_post<-bE+nE*tE
  
  lame_post<-rgamma(T,ae_post,be_post)
  
  dEPn<-rep(0,T)
  dRPn<-rep(0,T)
  
  for(t in 1:T)
  {
    
    dEPn[t]=lame_post[t]-lamp_post[t]
    dRPn[t]=lamr_post[t]-lamp_post[t]
  }
  
  count2=0
  
  for(t in 1:T)
  {
    if (dEPn[t]/dRPn[t]>theta)
      count2=count2+1
  }
  
  Est_Prob=count2/T
  return(Est_Prob)
}

# informative #
data.nc.info1<-function(theta)
{
  dataP<-rep(0:4,round(c(27.5/2,13.7/2,40/2,13.7/2,3.9/2)))
  dataR<-rep(0:3,round(c(73.9*46/100,17.4*46/100,6.5*46/100,2.2*46/100)))
  dataE<-rep(0:3,round(c(45.8*48/100,33.3*48/100,16.7*48/100,4.2*48/100)))
  
  xP=sum(dataP)
  xR=sum(dataR)
  xE=sum(dataE)
  tE<-tR<-tP<-1
  
  nP=length(dataP)
  nR=length(dataR)
  nE=length(dataE)
  
  aE<-160
  bE<-200
  
  model.string <-"
  model{
  
  xr~dpois(46*lamR)
  xp~dpois(50*lamP)
  
  u1~dbeta(alpha,beta)
  u2~dgamma(r,p)
  
  alpha<-200
  beta<-630
  
  r<-600
  p<-400
  
  lamR<-u1*u2
  lamP<-u2
  }"

  model.spec<-textConnection(model.string)
  
  load.module("bugs")
  load.module("mix")
  
  dat<-list(xP,xR)
  jags.inits<-list(u1=.241,u2=1.5)
  parameters<-c("lamP","lamR")
  
  beta_binom<-jags.model(model.spec,data=dat,inits=jags.inits,n.chains=1,n.adapt=1000)
  
  update(beta_binom,100000)
  samples<-coda.samples(beta_binom,variable.names=parameters,n.iter=50000,thin=50)
  post_samples<-unlist(samples)
  
  lamp_post<-post_samples[1:1000]
  lamr_post<-post_samples[1001:2000]
  
  ae_post<-aE+xE
  be_post<-bE+nE*tE
  
  lame_post<-rgamma(T,ae_post,be_post)
  
  dEPn<-rep(0,T)
  dRPn<-rep(0,T)
  
  for(t in 1:T)
  {
    
    dEPn[t]=lame_post[t]-lamp_post[t]
    dRPn[t]=lamr_post[t]-lamp_post[t]
  }
  
  count2=0
  
  for(t in 1:T)
  {
    if (dEPn[t]/dRPn[t]>theta)
      count2=count2+1
  }
  
  Est_Prob=count2/T
  return(Est_Prob)
}
###################### 2 year #########################

# non-informative #
data.nc.noninfo2<-function(theta)
{
  dataP<-rep(0:8,round(c(17.6/2,11.8/2,15.7/2,9.8/2,18/2,13.7/2,7.8/2,2/2,2/2)))
  dataR<-rep(0:4,round(c(47.8*46/100,39.1*46/100,8.7*46/100,2.2*46/100,2.2*46/100)))
  dataE<-rep(0:5,round(c(33.3*48/100,27.1*48/100,20.8*48/100,16.7*48/100,0*48/100,2.1*48/100)))
  
  xP=sum(dataP)
  xR=sum(dataR)
  xE=sum(dataE)
  tE<-tR<-tP<-1
  
  nP=length(dataP)
  nR=length(dataR)
  nE=length(dataE)
  
  aE<-0.8
  bE<-0.62
  
  model.string <-"
  model{
  
  xr~dpois(46*lamR)
  xp~dpois(50*lamP)
  
  u1~dbeta(alpha,beta)
  u2~dgamma(r,p)
  
  alpha<-1.2
  beta<-3.7
  
  r<-.75
  p<-.255
  
  lamR<-u1*u2
  lamP<-u2
  }"

  model.spec<-textConnection(model.string)
  
  load.module("bugs")
  load.module("mix")
  
  dat<-list(xP,xR)
  jags.inits<-list(u1=.245,u2=2.94)
  parameters<-c("lamP","lamR")
  
  beta_binom<-jags.model(model.spec,data=dat,inits=jags.inits,n.chains=1,n.adapt=1000)
  
  update(beta_binom,100000)
  samples<-coda.samples(beta_binom,variable.names=parameters,n.iter=50000,thin=50)
  post_samples<-unlist(samples)
  
  lamp_post<-post_samples[1:1000]
  lamr_post<-post_samples[1001:2000]
  
  ae_post<-aE+xE
  be_post<-bE+nE*tE
  
  lame_post<-rgamma(T,ae_post,be_post)
  
  dEPn<-rep(0,T)
  dRPn<-rep(0,T)
  
  for(t in 1:T)
  {
    
    dEPn[t]=lame_post[t]-lamp_post[t]
    dRPn[t]=lamr_post[t]-lamp_post[t]
  }
  
  count2=0
  
  for(t in 1:T)
  {
    if (dEPn[t]/dRPn[t]>theta)
      count2=count2+1
  }
  
  Est_Prob=count2/T
  return(Est_Prob)
}

# informative #
data.nc.info2<-function(theta)
{
  dataP<-rep(0:8,round(c(17.6/2,11.8/2,15.7/2,9.8/2,18/2,13.7/2,7.8/2,2/2,2/2)))
  dataR<-rep(0:4,round(c(47.8*46/100,39.1*46/100,8.7*46/100,2.2*46/100,2.2*46/100)))
  dataE<-rep(0:5,round(c(33.3*48/100,27.1*48/100,20.8*48/100,16.7*48/100,0*48/100,2.1*48/100)))
  
  xP=sum(dataP)
  xR=sum(dataR)
  xE=sum(dataE)
  tE<-tR<-tP<-1
  
  nP=length(dataP)
  nR=length(dataR)
  nE=length(dataE)
  
  aE<-80
  bE<-62
  
  model.string <-"
  model{
  
  xr~dpois(46*lamR)
  xp~dpois(50*lamP)
  
  u1~dbeta(alpha,beta)
  u2~dgamma(r,p)
  
  alpha<-12
  beta<-37
  
  r<-75
  p<-25.5
  
  lamR<-u1*u2
  lamP<-u2
  }"

  model.spec<-textConnection(model.string)
  
  load.module("bugs")
  load.module("mix")
  
  dat<-list(xP,xR)
  jags.inits<-list(u1=.246,u2=1.5)
  parameters<-c("lamP","lamR")
  
  beta_binom<-jags.model(model.spec,data=dat,inits=jags.inits,n.chains=1,n.adapt=1000)
  
  update(beta_binom,100000)
  samples<-coda.samples(beta_binom,variable.names=parameters,n.iter=50000,thin=50)
  post_samples<-unlist(samples)
  
  lamp_post<-post_samples[1:1000]
  lamr_post<-post_samples[1001:2000]
  
  ae_post<-aE+xE
  be_post<-bE+nE*tE
  
  lame_post<-rgamma(T,ae_post,be_post)
  
  dEPn<-rep(0,T)
  dRPn<-rep(0,T)
  
  for(t in 1:T)
  {
    
    dEPn[t]=lame_post[t]-lamp_post[t]
    dRPn[t]=lamr_post[t]-lamp_post[t]
  }
  
  count2=0
  
  for(t in 1:T)
  {
    if (dEPn[t]/dRPn[t]>theta)
      count2=count2+1
  }
  
  Est_Prob=count2/T
  return(Est_Prob)
}
