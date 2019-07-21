source('./power_poisson_functions.R')

time<-proc.time()
set.seed(123)

n_star<- 10000
T<- 1000
p_star<-0.975


n1<- c(79,112,175,302,685)
n2<- c(38,48,65,88,133)
n3<- c(37,52,81,153,351)
n4<- c(18,24,33,45,64)
n5<- c(31,44,71,125,277)
n6<- c(15,18,26,37,54)

n<-c(n1,n2,n3,n4,n5,n6)
r.alloc<-matrix(c(1,1,2,2,2,3),2,3)

power<-rep(0,30)
nc<-1

lamR<-21
lamP<-7

p_star<-0.975

lamE<-c(20,19.7,19.4,19.1,18.8)

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
      Est_Prob<-rep(0,n_star)
      
      for(i in 1:n_star)
      {
        count1=0
        count2=0
        
        tP=tR=tE<-1
        
        xP=rpois(1,nP*lamP*tP)
        xR=rpois(1,nR*lamR*tR)
        xE=rpois(1,nE*lamE[j]*tE)
        
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
      power[nc]<- count/n_star
      nc=nc+1
    }
  }
}

(proc.time()-time)[3]

power
