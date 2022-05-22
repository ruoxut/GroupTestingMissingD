#-------------------------------------------------------------------------------
#This file performs the simulation for \hat{p}_1.

#Reference: Delaigle, A. and Tan, R. (2022+) Group testing regression analysis 
#with missing data and imperfect test. Statistica Sinica.

#Author: Ruoxu Tan; last updated: 2022/May/22.
#-------------------------------------------------------------------------------

library("parallel")

print(Sys.time())

#-------------------------------------------------------------------------------
#Simulation Models 
#-------------------------------------------------------------------------------

#Number of groups J
J_sample = c(250,500,1000,2000)

#Models for p(x); uncomment one part and comment other ones to select the model as needed.

#Model (i)----------------------------------------------------------------------
p = function(x){
  if (x > 2.828 | x< -2.828) {res = 1}
  else {res = x^2/8}
  return(res)
}
#-------------------------------------------------------------------------------

#Model (ii)---------------------------------------------------------------------
#p = function(x){
#  if (x > 3.08) {res = 0}
#  else if (x < -3) {res = 1}
#  else {res = 1/(1+exp(2*x+4))+(x-0.4)^2*sin(pi*x)/20+0.1}
#  return(res)
#}
#-------------------------------------------------------------------------------

#Model (iii)--------------------------------------------------------------------
#p = function(x){
#  1/(1+exp(2*x+3))
#}
#-------------------------------------------------------------------------------

pv = function(x){sapply(x,p)}

#Models for the missing data mechanism; uncomment one part and comment other ones 
#to select the model as needed.

#Model (1)----------------------------------------------------------------------
pmis = function(x){
  0.7+0.3*sin((x-1)^2)
}
#-------------------------------------------------------------------------------

#Model (2)----------------------------------------------------------------------
#pmis = function(x){
#  exp(sin(x)+0.5)/(1+exp(sin(x)+0.5))
#}
#-------------------------------------------------------------------------------

pmisv = function(x){sapply(x,pmis)}

#Sepcify sp and se.
s_1 = 0.01# 1 - specificity = false positive
s_2 = 0.15# 1 - sensitivity = false negative

#Define the x-range of interest
x_range = seq(-1.5,1.5,length.out = 100)

#-------------------------------------------------------------------------------
#Simulation 
#-------------------------------------------------------------------------------
for (s in 1:4)
{
  
  #Group setting, J is group number, n is the vector of group sizes
  J = J_sample[s]
  
  #Grouping (A)-----------------------------------------------------------------
  #N = J*6
  #-----------------------------------------------------------------------------
  
  #Grouping (B)-----------------------------------------------------------------
  #N = J*5
  #-----------------------------------------------------------------------------
  
  #Grouping (C)-----------------------------------------------------------------
  N = J*12
  #-----------------------------------------------------------------------------
  
  sim = function(i)
  { 
    set.seed(3*i)
    
    #Generate X
    X = rnorm(N,0,0.75)
    
    #Generate D
    D = rbinom(N,1,pv(X))
    
    #Generate R^D
    R = rbinom(N,1,pmisv(X))
    
    #Generate grouped outcomes
    N_obs = sum(R)
    
    #Grouping (A)---------------------------------------------------------------
    #if (N_obs %% 12 == 0){
    #  J = N_obs/6
    #  n = rep(0,J)
    #  for (i in 1:J/2){n[i]=4}
    #  for (i in (J/2+1):J){n[i]=8}
    #} else
    #{
    #  J = N_obs %/% 12 * 2 +1
    #  n = rep(0,J)
    #  for (i in 1:((J-1)/2)){n[i]=4}
    #  for (i in ((J-1)/2+1):(J-1)){n[i]=8}
    #  n[J] = N_obs %% 12
    #}
    #---------------------------------------------------------------------------
    
    #Grouping (B)---------------------------------------------------------------
    #if (N_obs %% 5 == 0){
    #  J = N_obs/5
    #  n = rep(5,J)
    #} else
    #{
    #  J = N_obs %/% 5 + 1
    #  n = rep(5,J)
    #  n[J] = N_obs %% 5
    #}
    #---------------------------------------------------------------------------
    
    #Grouping (C)---------------------------------------------------------------
    if (N_obs %% 12 == 0){
      J = N_obs/12
      n = rep(12,J)
    } else
    {
      J = N_obs %/% 12 + 1
      n = rep(12,J)
      n[J] = N_obs %% 12
    }
    #---------------------------------------------------------------------------
    
    n1 = c(0,n)
    D = D[R==1]
    X = X[R==1]
    
    D_star = rep(0,J)
    for (j in 1:J)
    {  D_star[j] =  max(D[(sum(n1[1:j])+1):sum(n1[1:(j+1)])],na.rm=TRUE)
    }
    
    Y_star = rep(0,J)
    Y_star[D_star==0] = rbinom(sum(D_star==0),1,s_1)
    Y_star[D_star==1] = rbinom(sum(D_star==1),1,1-s_2)
    
    Z_star = 1-Y_star
    ZZ = rep(Z_star,n)
    
    #Compute \hat{p}_1
    ISE_1 = function(x_eva)
    {
      source("loclin of p_tilde.R",local=TRUE)
      
      sqer = function(x){(hatpvec_1(x)-pv(x))^2}
      
      
      res = list(hatp_1_res = hatpvec_1(x_eva), 
                 ISE_1_res = integrate(sqer,lower=min(x_eva),upper=max(x_eva))$value)
      
      return(res)
    }
    
    
    #Store results
    res = ISE_1(x_range)
    
    return(res)
  }
  
  
  #Apply the simulation for iter times with parallel computing.
  iter = 200
  
  cl = makeCluster(16)
  
  clusterExport(cl, c('N','J','p','pv','pmis','pmisv','x_range','s_1','s_2','sim'))
  
  simulation = parLapply(cl, 1:iter, function(i){try(sim(i),T)})
  
  file_name = paste0("201120_i_1_ungr_J", J, ".RData")
  save(simulation,file=file_name)
  
  stopCluster(cl)
}  

print(Sys.time())

#-------------------------------------------------------------------------------
#Print interquartile ISEs.
#-------------------------------------------------------------------------------

file_names = list.files(pattern="201120_i_1_ungr_J")

for(k in 1:length(file_names))
{
  eval(parse(text=paste('load("',file_names[k],'")', sep = "")))
  ISE_1 = rep(0,200)
  
  for(i in 1:200)
  {
    ISE_1[i] = simulation[[i]]$ISE_1_res
  }
  
  m_1 = format(round(median(ISE_1)*10^3, 2), nsmall = 2)
  o_1 = order(ISE_1)
  iqr_1 = format(round((ISE_1[o_1[150]]-ISE_1[o_1[50]])*10^3,2), nsmall = 2)
  
  
  print(paste0(file_names[k], " ISE", ": ",  m_1," ", "(", iqr_1, ")"))
}