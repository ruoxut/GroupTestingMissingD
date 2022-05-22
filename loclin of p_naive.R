#-------------------------------------------------------------------------------
#Compute \hat{p}_{nai}
#-------------------------------------------------------------------------------

#Delete empty groups
if (2 %in% Z_star)
{
  n = n[-which(Z_star==2)]
  n1 = c(0,n)
  X = X[-which(ZZ==2)]
  
  Z_star = Z_star[-which(Z_star==2)]
  ZZ = ZZ[-which(ZZ==2)]
  
  J = length(Z_star)
  N = length(ZZ)
}

#-------------------------------------------------------------------------------
#Point estimations
#-------------------------------------------------------------------------------

#Define -log L (ED|Z^*_1,...,Z^*_J) and find MLE of ED
fobject = function(ED)
{    
  f_1 = 0
  f_0 = 0

  for (j in 1:J)
  { if (Z_star[j]==1) {
    if ((1-s_1-s_2)*(1-ED)^n[j]+s_2 <= 0)
    {f_1 = f_1 - 1e3} else
    {f_1 = f_1 + log((1-s_1-s_2)*(1-ED)^n[j]+s_2)}
  } else  {
    if (1 - (1-s_1-s_2)*(1-ED)^n[j] - s_2 <= 0)
    {f_2 = f_2 - 1e3} else
    {f_0 = f_0 + log(1 - (1-s_1-s_2)*(1-ED)^n[j] - s_2)}
  }
  }
  
  res = - (f_1 + f_0)

  return(res)
}

#Calculate hatED 
hatED = optim(par=0.5, fn=fobject, method="Brent", lower=0.001, upper=0.999)
hatED = hatED$par
q = rep(0,J)
for (j in 1:J) q[j]=(1-hatED)^(n[j]-1)

qq = rep(q,n)


#-------------------------------------------------------------------------------
#bandwidth selection
#-------------------------------------------------------------------------------

ZZ_qs = (qq*(1-s_1-s_2))^(-1) * (ZZ-s_2)

#Local constant estimator of m at a single point x
fcons = function(x,X,ZZ_qs,h)
{    
  a0 = dnorm((X-x)/h)

  S0 = sum(a0)
  T0 = sum(ZZ_qs*a0)

  res = T0/S0
  if(res<0) {
    res = 0
  } else if (res>1) {
    res = 1
  }
  
  return(res)
}


#Local constant estimator of m at a vector x
fconsvec =function(x,X,ZZ_qs,h)
{    
  sapply(x,fcons,X,ZZ_qs,h)
}


#Crossvalidation bandwidth with local constant estimator of m
CVcons = function(X,ZZ_qs,h)
{    
  index1 = cumsum(n1)+1
  index2 = cumsum(n)

  index1 = rep(index1[1:J],n)
  index2 = rep(index2,n)

  qantX = quantile(X,c(0.1,0.9))
  W = dunif(X,qantX[1],qantX[2]) * (qantX[2]-qantX[1])

  RS = 0
  for(i in 1:N)
  {   out = index1[i]:index2[i]

      if(i %in% which(W>0))
      {   RS = RS + (ZZ_qs[i] - fcons(X[i],X[-out],ZZ_qs[-out],h))^2
      }
  }
  
  return(RS)
}

CVobject = function(h){CVcons(X=X,ZZ_qs=ZZ_qs,h)}

h1 = optim(par=sd(X)*length(X)^(-0.2), fn=CVobject, method="Brent", lower=0.001, upper=max(x_eva)-min(x_eva), control=list(maxit=10))
h1 = h1$par

#Define the local constant estimator of m_{PILOT}(x)
mPIL0 = function(x){fconsvec(x,X=X,ZZ_qs=ZZ_qs,h=h1)}

#Calculate psi
qantX = quantile(X,c(0.1,0.9))

binsize = (qantX[2]-qantX[1])/100
mPIL0_X = t(as.matrix(mPIL0(seq(qantX[1],qantX[2],binsize))))
q_col = as.matrix(q)
ones_col = as.matrix(rep(1,J))
ones_row = t(as.matrix(rep(1,length(mPIL0_X))))

V = ((q_col*(1-s_1-s_2))^(-1)) %*% ((1-2*s_2)*mPIL0_X) + ((s_2-s_2^2)*(q_col*(1-s_1-s_2))^(-2)) %*% ones_row - ones_col %*% (mPIL0_X^2)

psi = (rowSums(abs(V)*binsize)+1e-4)^(-1)

psi = rep(psi,n)

#Global cubic estimator of m^{(2)}(x) at a point x
fcub = function(x,X,ZZ_qs)
{    
  S0 = N
  S1 = sum(X)
  S2 = sum(X^2)
  S3 = sum(X^3)
  S4 = sum(X^4)
  S5 = sum(X^5)
  S6 = sum(X^6) 

  T0 = sum(ZZ_qs)
  T1 = sum(ZZ_qs*X)
  T2 = sum(ZZ_qs*X^2)
  T3 = sum(ZZ_qs*X^3)

  a3 = (S3^3*T0 + S2^2*S5*T0 + S2^2*S4*T1 - S0*S4^2*T1 - S0*S2*S5*T2 + S1*(S4^2*T0 - S2*S5*T1 - S2*S4*T2) - S2^3*T3 + S0*S2*S4*T3 - S3^2*(S2*T1 + S1*T2 + S0*T3) + S3*(-2*S2*S4*T0 - S1*S5*T0 + S1*S4*T1 + S0*S5*T1 + S2^2*T2 + S0*S4*T2 + 2*S1*S2*T3) + S1^2*(S5*T2 - S4*T3))/
       (S3^4 + S2^2*S4^2 - S0*S4^3 + S1^2*S5^2 - S2^3*S6 - S1^2*S4*S6 - S3^2*(3*S2*S4 + 2*S1*S5 + S0*S6) + S2*(-2*S1*S4*S5 - S0*S5^2 + S0*S4*S6) + 2*S3*((S2^2 + S0*S4)*S5 + S1*(S4^2 + S2*S6)))
  a2 = (-S1*S4*S5*T0 - S2^2*S6*T0 + S3^3*T1 + S0*S4*S5*T1 - S0*S4^2*T2 - S1^2*S6*T2 + S1^2*S5*T3 - S3^2*(S4*T0 + S2*T2 + S1*T3) + S3*(S2*S5*T0 + S1*S6*T0 - S2*S4*T1 - S1*S5*T1 - S0*S6*T1 + 2*S1*S4*T2 + S2^2*T3 + S0*S4*T3) + S2*(S4^2*T0 + S1*S6*T1 + S0*S6*T2 - S1*S4*T3 - S0*S5*T3))/
       (S3^4 + S2^2*S4^2 - S0*S4^3 + S1^2*S5^2 - S2^3*S6 - S1^2*S4*S6 - S3^2*(3*S2*S4 + 2*S1*S5 + S0*S6) + S2*(-2*S1*S4*S5 - S0*S5^2 + S0*S4*S6) + 2*S3*((S2^2 + S0*S4)*S5 + S1*(S4^2 + S2*S6)))

  return(6*a3*x+2*a2)
}


#Global cubic estimator of m^{(2)}(x) at a vector x
fcubvec =function(x,X,ZZ_qs)
{    
  sapply(x,fcub,X,ZZ_qs)
}


#Calculate estimator of Theta^{(2)}
X_we = X[X>=qantX[1] & X<=qantX[2]]
Theta2 = sum((fcubvec(X_we,X,ZZ_qs))^2)/N


#Calculate the plug-in estimator of h_{PI}
hPI = (1/(2*sqrt(pi)*Theta2*sum(psi)))^(1/5)


#-------------------------------------------------------------------------------
#Local linear estimation
#-------------------------------------------------------------------------------
#Local linear estimator of m at a single point x
flin = function(x,X,ZZ_qs,psi,h)
{    
  a0 = dnorm((X-x)/h)
  a1 = a0*(X-x)/h
  a2 = a1*(X-x)/h

  S0 = sum(psi*a0)
  S1 = sum(psi*a1)
  S2 = sum(psi*a2)
  T0 = sum(psi*ZZ_qs*a0)
  T1 = sum(psi*ZZ_qs*a1)

  return((T0*S2-T1*S1)/(S0*S2-S1^2))
}


#Calculate the estimator of p(x) at a point x
hatp = function(x)
{    
  if(1-flin(x,X,ZZ_qs,psi,hPI)<0)
  {res=0} else if(1-flin(x,X,ZZ_qs,psi,hPI)>1)
  {res=1} else {res=1-flin(x,X,ZZ_qs,psi,hPI)}
  
  return(res)
}


#Calculate the estimator of p(x) at a vector x
hatpvec_naive = function(x)
{    
  sapply(x,hatp)
}

