### Sample size calculation for power=0.8 ###

lamR<-21
lamP<-7

p_star<-0.975

### sample size table ###
n.table.freq<-NULL
r.alloc.list<-matrix(c(1,1,1,1,2,2,1,2,3),3,3)

for (i.r.alloc in c(1:3))
{
  r.alloc<-r.alloc.list[i.r.alloc,]
  for (theta in c(0.8,0.75))
  {
    for (lamE in c(20,19.7,19.4,19.1,18.8))
    {
      n<-samplesize_fn_freq(r.alloc=r.alloc,lamE=lamE,theta=theta)
      n.table.freq<-rbind(n.table.freq,c(r.alloc=r.alloc,lamE=lamE,theta=theta,n.freq=n))
    }
  }
}

save(n.table.freq,file="n.table.freq.RData")

############## full Bayes ###################
n_star<-1000 ### No. of simulations ###
T<-1000 ### samples ###

p_star<-0.975

lamR<-21
lamP<-7

r.alloc.list<-matrix(c(1,1,1,1,2,2,1,2,3),3,3)
load("n.table.freq.RData")
fullbayes.table<-NULL
counter=1
for (i.r.alloc in c(1:3))
{
  r.alloc<-r.alloc.list[i.r.alloc,]
  for (theta in c(0.8,0.75))
  {
    for (lamE in c(20,19.7,19.4,19.1,18.8))
    {
      n<-samplesize_fn_fbayes(r.alloc=r.alloc,lamE=lamE,theta=theta,a_max=n.table.freq[counter,6])
      counter=counter+1 
      fullbayes.table<-rbind(fullbayes.table,c(r.alloc=r.alloc,lamE=lamE,theta=theta,n=n))
    }
  }
}

#################### Approx Bayes ###################

p_star<-0.975
lamR<-21
lamP<-7

n.table.ab<-NULL
r.alloc.list<-matrix(c(1,1,1,1,2,2,1,2,3),3,3)

for (i.r.alloc in c(1:3))
{
  r.alloc<-r.alloc.list[i.r.alloc,]
  for (theta in c(0.8,0.75))
  {
    for (lamE in c(20,19.7,19.4,19.1,18.8))
    {
      n<-samplesize_fn_approxbayes(r.alloc=r.alloc,lamE=lamE,theta=theta)
      n.table.ab<-rbind(n.table.ab,c(r.alloc=r.alloc,lamE=lamE,theta=theta,n.ab=n))
    }
  }
}
