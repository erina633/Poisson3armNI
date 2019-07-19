#### fully bayesian type I error ####
type1_fbayes = function(alloc, n, l_E, l_R, l_P, theta, aE, bE, aR, bR, aP, bP){
  set.seed(123)
  n_star<- 10000
  
  Est_Prob = rep(0, n_star)
  
  #Sample size correponding to experimental, reference, and placebo arm
  n_arm = rep(NA, 3)
  for(i in 1:3){
    n_arm[i] = n*alloc[i]
  }
  
  nP=n_arm[1]
  nR=n_arm[2]
  nE=n_arm[3]
  
  tE<-tR<-tP<-1
  
  count<-0
  
  ########## Test Procedure ##############
  for(i in 1:n_star)
  {
    count1=0
    count2=0
    
    xP=rpois(1,nP*l_P*tP)
    xR=rpois(1,nR*l_R*tR)
    xE=rpois(1,nE*l_E*tE)
    
    while (count1<T_val)
    {
      
      l_Egdata=rgamma(1,aE+xE,bE+nE*tE)
      l_Rgdata=rgamma(1,aR+xR,bR+nR*tR)
      l_Pgdata=rgamma(1,aP+xP,bP+nP*tP)
      
      dEPn=l_Egdata-l_Pgdata
      dRPn=l_Rgdata-l_Pgdata
      
      if(dRPn>0)
      { 
        count1=count1+1
        
        if (dEPn/dRPn<=((l_E - l_P)/(l_R - l_P)))
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

alloc = matrix(c(1,1,1,2,2,1,3, 2,1),3,3)

l_E = c(20,19.7,19.4,19.1,18.8)
theta = .75
p_star = 0.975
T_val = 1000
l_R = 21
l_P = 7

aE = 21; bE = 1; aR = 21; bR = 1; aP = 1; bP = 1

fb_type1 = rep(NA, length(l_E))
for(i in 1: length(l_E)){
  fb_type1[i] = type1_fbayes(alloc[,1], n = n1, l_E[i], l_R, l_P, theta, aE, bE, aR, bR, aP, bP)
}
fb_type1


