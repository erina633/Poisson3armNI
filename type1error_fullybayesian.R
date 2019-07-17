#### fully bayesian power ####
power_fbayes<-function(r.alloc, n, lamE, theta)
{
  set.seed(123)
  n_star<- 10000
  p_star<-0.975
  T_val = 1000
  lamR<-21
  lamP<-7
  
  Est_Prob = rep(0, n_star)
  
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
    
    while (count1<T_val)
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
        
        if (dEPn/dRPn<=((lamE - lamP)/(lamR - lamP)))
          count2=count2+1
      }
    }
    
    Est_Prob[i]=count2/T_val
    
    if(Est_Prob[i]>p_star)
      count=count+1
  }
  
  power<-count/n_star
  return(power)
} 

n1 = c(78,112,175,310,697)
n2 = c(38,48,65,91,138)
n3 = c(40,57,89,157,352)
n4 = c(19,25,33,46,70)
n5 = c(32,46,71,127,284)
n6 = c(15,20,27,37,56)

r.allocmat = matrix(c(1,1,1,2,2,1,3, 2,1),3,3)

lamE = c(20,19.7,19.4,19.1,18.8)
theta = .75

fb_type1 = rep(NA, length(lamE))
for(i in 1: length(lamE)){
    fb_type1[i] = power_fbayes(r.allocmat[,1], n2, lamE[i], theta)
}
fb_type1
#write.csv(fb_type1,file="C:\\Erina\\Paper\\Analysis\\Shrabanti_Poisson_Paper\\type1.csv")
