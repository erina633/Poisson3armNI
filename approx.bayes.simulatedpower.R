source("./power_poisson_functions.R")

ptm<-proc.time()
set.seed(123)

n_star<- 10000
p_star<-0.975

n1<- c(78,112,175,310,697)
n2<- c(38,48,65,91,138)
n3<- c(40,57,89,157,352)
n4<- c(19,25,33,46,70)
n5<- c(32,46,71,127,284)
n6<- c(15,20,27,37,56)

n<-c(n1,n2,n3,n4,n5,n6)
r.alloc<-matrix(c(1,1,2,2,2,3),2,3)

power<-rep(0,30)
nc<-1

lamRt<-21
lamPt<-7

p_star<-0.975

lamEt<-c(20,19.7,19.4,19.1,18.8)

for(l in c(1:3))
{
  lr=r.alloc[1,l]
  le=r.alloc[2,l]
  
  for(theta in c(0.8,0.75))
  {
    for(j in 1:5)
    {
      count<- 0
      nP=n[nc]
      nR=lr*n[nc]
      nE=le*n[nc]
      
      for(i in 1:n_star)
      {
        xP=rpois(1,nP*lamPt)
        xR=rpois(1,nR*lamRt)
        xE=rpois(1,nE*lamEt[j])
        
        lamP<-xP/nP
        lamR<-xR/nR
        lamE<-xE/nE
        
        lamE_null<-lamPt + theta*(lamRt-lamPt)
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
        
        mu_tilt<- (mu_T_theta/sigma2_T + mu_star_th/sigma2_star_nu)
        
        sigma2_tilt<- 1/((1/sigma2_T) + (1/sigma2_star_nu))
        
        point<- mu_tilt*sqrt(sigma2_tilt)
        
        if(pnorm(point,0,1)>p_star)
        {
          count=count+1
        }
      }
      
      power[nc]<- count/n_star
      nc=nc+1
    }
  }
}

time.taken<-(proc.time()-ptm)[3]
time.taken 

power
