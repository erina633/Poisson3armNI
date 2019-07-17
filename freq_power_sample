#Sample size calculation under frequentist setup

frequentist_power = function(n_total, alloc, theta, l_E, l_R, l_P, alpha_val){
  
  #Sample size correponding to experimental, reference, and placebo arm
  n_arm = rep(NA, 3)
  for(i in 1:3){
    n_arm[i] = n_total*alloc[i]
  }
  
  #Mean and variance of W = (U - theta*V)|V>0 = w1 - theta*w2, where U = l_E - l_P and V = l_R - l_P
  #See appendix A.1 for the calculation of mean and variance of W
  
  #Null value of l_E
  l_E0 = l_P + theta*(l_R - l_P) 
  
  mu_u = l_E - l_P
  mu_u0 = l_E0 - l_P
  
  var_u  = (l_E/n_arm[1]) + (l_P/n_arm[3])
  var_u0 = (l_E0/n_arm[1]) + (l_P/n_arm[3])
  
  mu_v = l_R - l_P
  var_v = (l_R/n_arm[2]) + (l_P/n_arm[3])
  
  d_val = -mu_v/sqrt(var_v)
  phi_val = dnorm(d_val, 0, 1)
  c_val = 1 - pnorm(d_val, 0, 1)
  
  rho_val = (l_P/n_arm[3])/(sqrt(var_u*var_v))
  rho_val0 = (l_P/n_arm[3])/(sqrt(var_u0*var_v))

  mu_w1 = mu_u + sqrt(var_u)*rho_val*phi_val/c_val
  mu_w1_0 = mu_u0 + sqrt(var_u0)*rho_val0*phi_val/c_val
  mu_w2 = mu_v + sqrt(var_v)*phi_val/c_val
  
  var_w1 = var_u*(1 + ((rho_val^2)*d_val*phi_val/c_val) - (rho_val*phi_val/c_val)^2)
  var_w1_0 = var_u0*(1 + ((rho_val0^2)*d_val*phi_val/c_val) - (rho_val0*phi_val/c_val)^2)
  var_w2 = var_v*(1 - (phi_val/c_val)*((phi_val/c_val) - d_val))
  
  cov_w1w2 = (sqrt(var_u*var_v)*rho_val*(c_val + d_val*phi_val)/c_val) + 
             (sqrt(var_u)*mu_v*rho_val*phi_val/c_val) +
             (sqrt(var_v)*mu_u*phi_val/c_val) + mu_u*mu_v - mu_w1*mu_w2
  cov_w1w2_0 = (sqrt(var_u0*var_v)*rho_val0*(c_val + d_val*phi_val)/c_val) + 
             (sqrt(var_u0)*mu_v*rho_val0*phi_val/c_val) +
             (sqrt(var_v)*mu_u0*phi_val/c_val) + mu_u0*mu_v - mu_w1_0*mu_w2
  
  mu_w_0 = mu_w1_0 - theta*mu_w2
  mu_w_1 = mu_w1 - theta*mu_w2
  
  var_w_0 = var_w1_0 + (theta^2)*var_w2 - 2*theta*cov_w1w2_0
  var_w_1 = var_w1 + (theta^2)*var_w2 - 2*theta*cov_w1w2
  
  z_val = qnorm(1-alpha_val, 0, 1)
  
  k_alpha = mu_w_0 + z_val*sqrt(var_w_0)
  
  power_freq = 1 - pnorm(((k_alpha - mu_w_1)/sqrt(var_w_1)), 0, 1)
  
  return(power_freq)
}


# Power curve 
alpha_val = 0.025
beta_val = 0.8
l_R = 21
l_P = 7
l_E = seq(14, 22, by = 0.1)
n_total = 300
theta = .8
alloc = matrix(c(1,1,1,2,2,1,3,2,1), byrow = T, 3, 3)

power_curve = matrix(NA, length(l_E), nrow(alloc))
for(i in 1: nrow(alloc)){
  for(j in 1:length(l_E)){
    power_curve[j, i] =  frequentist_power(n_total, alloc = alloc[i,], theta, l_E[j], l_R, l_P, alpha_val)
  }
}

library(ggplot2)
library(reshape2)
df = melt(power_curve)
colnames(df) = c("l_E", "allocation", "value")
df[,1] = rep(l_E, nrow(alloc))
df[,2]= as.factor(df$allocation)

ggplot(df, aes(x = l_E, y = value, color = allocation)) +
  geom_line(aes(linetype = allocation, color = allocation)) + 
  labs(x = expression(lambda[E]), y = "Power") +
  xlab(expression(paste("N = 300, ", theta," = 0.8"))) +
  scale_linetype_manual(values = 1:3, labels = c("1:1:1", "2:2:1", "3:2:1")) +
  scale_color_manual(values = c("green", "red", "blue"), labels = c("1:1:1", "2:2:1", "3:2:1"))+
  ggtitle("Power curves for different allocations under Frequentist approach 
            with allocation ratio nE:nR:nP") + theme(plot.title = element_text(hjust = 0.5))


# Power curve 
alpha_val = 0.025
beta_val = 0.8
l_R = 7#21
l_P = 1#7
l_E = seq(4, 9, by = 0.1)#seq(14, 22, by = 0.1)
n_total = 300
theta = seq(.65, .8, by = .05)
alloc = c(1, 1, 1)

power_curve = matrix(NA, length(l_E), length(theta))
for(i in 1: length(theta)){
  for(j in 1:length(l_E)){
    power_curve[j, i] =  frequentist_power(n_total, alloc = alloc, theta[i], l_E[j], l_R, l_P, alpha_val)
  }
}

library(ggplot2)
library(reshape2)
df = melt(power_curve)
colnames(df) = c("l_E", "theta", "value")
df[,1] = rep(l_E, length(theta))
df[,2]= as.factor(df$theta)

ggplot(df, aes(x = l_E, y = value, color = theta)) +
  geom_line(aes(linetype = theta, color = theta)) + 
  labs(x = expression(lambda[E]), y = "Power") +
  xlab("N = 300") +
  scale_linetype_manual(values = 1:4, labels = theta) +
  scale_color_manual(values = c("green", "red", "blue", "black"), labels = theta)+
  ggtitle("Power curves for different theta under Frequentist approach") + theme(plot.title = element_text(hjust = 0.5))



# sample size calculation 
alpha_val = 0.025
beta_val = 0.8
l_E = c(20,19.7,19.4,19.1,18.8)
l_R = 21
l_P = 7
theta = c(0.8,0.75)
max_n = 1000
alloc = matrix(c(1,1,1,2,2,1,3,2,1), byrow = T, 3, 3)


freq_table = NULL
for(i in 1:nrow(alloc)){
  for(j in 1:length(theta)){
    for(k in 1:length(l_E)){
      
      power = NULL
      for(l in 1:max_n)
      {
        power = c(power, frequentist_power(n_total = l, alloc[i,], theta[j], l_E[k], l_R, l_P, alpha_val))
      }
      
      n = min(which(power>= beta_val))
      power_n = round(power[n], 3)
      
      freq_table = rbind(freq_table,c(alloc = alloc[i,],l_E = l_E[k], theta = theta[j],
                                      n = n, power = power_n))
    }
  }
}



