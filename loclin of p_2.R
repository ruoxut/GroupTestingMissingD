#-------------------------------------------------------------------------------
#Compute \hat{p}_3
#-------------------------------------------------------------------------------

#Define |I_j|
W = rep(0,J)
for (j in 1:J)
{W[j]=sum(R[(sum(n1[1:j])+1):sum(n1[1:(j+1)])])}

WW = rep(W,n)

#-------------------------------------------------------------------------------
#Point estimations
#-------------------------------------------------------------------------------

#hatER is the sample mean of R^D
hatER = sum(W)/N

#Define -log L (E(R^DD)|Z^*_1,...,Z^*_J) 
fobject = function(x)
{    
  ERD = x
  
  L = 0
  
  for (j in 1:J)
  {    if (Z_star[j]==1) 
  {
    if(s_2 + (1-s_1-s_2)*(1-ERD)^n[j] - (1-s_1)*(1-hatER)^n[j] <= 0)
    { L = L - 1e3 } else
    { L = L + log(s_2 + (1-s_1-s_2)*(1-ERD)^n[j] - (1-s_1)*(1-hatER)^n[j]) }
  } else if (Z_star[j]==0) 
  {
    if(1 - s_2 - (1-s_1-s_2)*(1-ERD)^n[j] - s_1*(1-hatER)^n[j] <= 0)
    { L = L - 1e3 } else
    { L = L + log(1 - s_2 - (1-s_1-s_2)*(1-ERD)^n[j] - s_1*(1-hatER)^n[j]) }
  } else 
  {
    if(((1-hatER)^n[j]) <= 0)
    { L = L - 1e3 } else
    { L = L + log((1-hatER)^n[j])}
  }
  }
  
  res = - L
  return(res)
}


#Find MLE of E(R^DD) and calculate q_{RD}^{n_j-1}
MLE = optim(par=hatER/2, fn=fobject, method="Brent", lower=0.001, upper=hatER)
hatERD = MLE$par

q = (1-hatERD)^(n-1)
qq = rep(q,n)

tau = (1-hatER)^(n-1)
tauu = rep(tau,n)


V = 1*(Y_star==0) + (1-s_1)*(Y_star==-1)  
VV = rep(V,n) #W_j
VV_qs = 1 - (VV-s_2)/((1-s_1-s_2)*qq) #U_{b,j}

WR = W-(n-1)*hatER
WRR = rep(WR,n) #U_{d,j}


#-------------------------------------------------------------------------------
#psi and bandwidth selection
#-------------------------------------------------------------------------------
#Local constant estimator at a single point x
fcons = function(x,X,YY,h)
{    
  a0 = dnorm((X-x)/h)

  S0 = sum(a0)
  T0 = sum(YY*a0)
  
  res = T0/S0
  
  if(res<0) {
    res = 0
  } else if(res>1) {
    res = 1
  }

  return(res)
}


#Local constant estimator at a vector x
fconsvec =function(x,X,YY,h)
{    
  sapply(x,fcons,X,YY,h)
}

#Crossvalidation bandwidth with local constant estimator
CVcons = function(X,YY,h)
{    
  index1 = cumsum(n1)+1
  index2 = cumsum(n)

  index1 = rep(index1[1:J],n)
  index2 = rep(index2,n)

  qantX = quantile(X,c(0.1,0.9))
  We = dunif(X,qantX[1],qantX[2]) * (qantX[2]-qantX[1])

  RS = 0
  for(i in 1:N)
  {   
    out = index1[i]:index2[i]

    if(i %in% which(We>0))
    {   
      RS = RS + (YY[i] - fcons(X[i],X[-out],YY[-out],h))^2
    }
  }

  return(RS)
}

CVobjectb = function(h){CVcons(X=X,YY=VV_qs,h)}

h1 = optim(par=sd(X)*length(X)^(-0.2), fn=CVobjectb, method="Brent", lower=(max(x_eva)-min(x_eva))*length(X)^(-0.5), upper=(max(x_eva)-min(x_eva))/2, control=list(maxit=10))
h1 = h1$par#CV bandwidth of b(x)

CVobjectd = function(h){CVcons(X=X,YY=WRR,h)}

h2 = optim(par=sd(X)*length(X)^(-0.2), fn=CVobjectd, method="Brent", lower=(max(x_eva)-min(x_eva))*length(X)^(-0.5), upper=(max(x_eva)-min(x_eva))/2, control=list(maxit=10))
h2 = h2$par#CV bandwidth of d(x)



#Define the local constant estimator of b_{PILOT}(x) and d_{PILOT}(x)
bPIL0 = function(x){fconsvec(x,X=X,YY=VV_qs,h=h1)}
dPIL0 = function(x){fconsvec(x,X=X,YY=WRR,h=h2)}


#Calculate psi
qantX = quantile(X,c(0.1,0.9))

binsize = (qantX[2]-qantX[1])/100
bPIL0_X = t(as.matrix(bPIL0(seq(qantX[1],qantX[2],binsize))))
dPIL0_X = t(as.matrix(dPIL0(seq(qantX[1],qantX[2],binsize))))
q_col = as.matrix(q)
n_col = as.matrix(n)
tau_col = as.matrix(tau)
ones_col = as.matrix(rep(1,J))
ones_row = t(as.matrix(rep(1,length(bPIL0_X))))

V_1 = ((((1-s_1-s_2)*q_col)^(-1)) %*% (1-bPIL0_X)) * (1-2*s_2-(1-s_1-s_2)*(q_col %*% (1-bPIL0_X))) -
      ((((1-s_1-s_2)*q_col)^(-2)) %*% ones_row) * (s_1*(1-s_1)*tau_col %*% (1-dPIL0_X) - s_2*(1-s_2))
V_2 = ((n_col-1)*hatER) %*% (1-bPIL0_X) + ones_col %*% (bPIL0_X*(1-bPIL0_X)) - 
      (q_col^(-1)*(n_col-1)*(hatER-hatERD)*(1-hatERD)^(n_col-2)) %*% (1-bPIL0_X)
V_3 = ((n_col-1)*hatER*(1-hatER)) %*% ones_row + ones_col %*% (dPIL0_X*(1-dPIL0_X))
V_F = V_1/(ones_col %*% (dPIL0_X)) - V_2*(ones_col %*% (2*bPIL0_X/(dPIL0_X^2))) + V_3*(ones_col %*% (bPIL0_X^2/(dPIL0_X^3)))

psi = (rowSums(abs(V_F)*binsize)+1e-4)^(-1)
psi = rep(psi,n)


#Global cubic estimator of m^{(2)}(x) at a point x
fcub = function(x,X,YY)
{    
  S0 = sum(N)
  S1 = sum(X)
  S2 = sum(X^2)
  S3 = sum(X^3)
  S4 = sum(X^4)
  S5 = sum(X^5)
  S6 = sum(X^6) 

  T0 = sum(YY)
  T1 = sum(YY*X)
  T2 = sum(YY*X^2)
  T3 = sum(YY*X^3)

  a3 = (S3^3*T0 + S2^2*S5*T0 + S2^2*S4*T1 - S0*S4^2*T1 - S0*S2*S5*T2 + S1*(S4^2*T0 - S2*S5*T1 - S2*S4*T2) - S2^3*T3 + S0*S2*S4*T3 - S3^2*(S2*T1 + S1*T2 + S0*T3) + S3*(-2*S2*S4*T0 - S1*S5*T0 + S1*S4*T1 + S0*S5*T1 + S2^2*T2 + S0*S4*T2 + 2*S1*S2*T3) + S1^2*(S5*T2 - S4*T3))/
      (S3^4 + S2^2*S4^2 - S0*S4^3 + S1^2*S5^2 - S2^3*S6 - S1^2*S4*S6 - S3^2*(3*S2*S4 + 2*S1*S5 + S0*S6) + S2*(-2*S1*S4*S5 - S0*S5^2 + S0*S4*S6) + 2*S3*((S2^2 + S0*S4)*S5 + S1*(S4^2 + S2*S6)))
  a2 = (-S1*S4*S5*T0 - S2^2*S6*T0 + S3^3*T1 + S0*S4*S5*T1 - S0*S4^2*T2 - S1^2*S6*T2 + S1^2*S5*T3 - S3^2*(S4*T0 + S2*T2 + S1*T3) + S3*(S2*S5*T0 + S1*S6*T0 - S2*S4*T1 - S1*S5*T1 - S0*S6*T1 + 2*S1*S4*T2 + S2^2*T3 + S0*S4*T3) + S2*(S4^2*T0 + S1*S6*T1 + S0*S6*T2 - S1*S4*T3 - S0*S5*T3))/
    (S3^4 + S2^2*S4^2 - S0*S4^3 + S1^2*S5^2 - S2^3*S6 - S1^2*S4*S6 - S3^2*(3*S2*S4 + 2*S1*S5 + S0*S6) + S2*(-2*S1*S4*S5 - S0*S5^2 + S0*S4*S6) + 2*S3*((S2^2 + S0*S4)*S5 + S1*(S4^2 + S2*S6)))

  return(6*a3*x+2*a2)
}

#Global cubic estimator of m^{(2)}(x) at a vector x
fcubvec =function(x,X,YY)
{    
  sapply(x,fcub,X,YY)
}

#Calculate estimator of Theta^{(2)}
X_we = X[X>=qantX[1] & X<=qantX[2]]
Theta2 = sum((dPIL0(X_we)*fcubvec(X_we,X,VV_qs)-bPIL0(X_we)*fcubvec(X_we,X,WRR))^2 / (dPIL0(X_we)^{3}+1e-4))/(N*hatER)

#Calculate the plug-in estimator of h_{PI}
hPI = (1/(2*sqrt(pi)*hatER*Theta2*sum(psi)))^(1/5)


#-------------------------------------------------------------------------------
#Local linear estimation
#-------------------------------------------------------------------------------
#Ratio of two local linear estimator of at a single point x
hatp = function(x,X,YY,ZZ,psi,h)
{    
  a0 = dnorm((X-x)/h)
  a1 = a0*(X-x)/h
  a2 = a1*(X-x)/h

  S1 = sum(psi*a1) 
  S2 = sum(psi*a2) 
  T0_n = sum(YY*psi*a0) 
  T1_n = sum(YY*psi*a1) 
  T0_d = sum(ZZ*psi*a0)
  T1_d = sum(ZZ*psi*a1)

  res = (T0_n*S2-T1_n*S1)/(T0_d*S2-T1_d*S1)
  if(res<0){
    res = 0
  }else if(res>1){
    res = 1
  }
  
  return(res)
}

#Calculate the estimator at a vector x
hatpvec_2 = function(x)
{    
  sapply(x,hatp,X,VV_qs,WRR,psi,hPI)
}
 