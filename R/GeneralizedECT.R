rm(list=ls())
##Economic Design (ED) and Economic Statistical Design (ESD) of
##""Shewhart-Type"" Control Charts Under "Single"" Assignable Cause
##Under Different Process Failure Mechanisms (PFM) and Sampling Schemes (SS)
##Speciically, for
##Xbar (with Normal or Non-normal Distributions) and T2 Control Charts
##(The Lines alpha and beta Could be Modified for any Other Control Statistics.)
##Under Weibull and Chen PFMs
##(The Lines ETa, F[1]=p[1], h[j], and R[j] could be modified for Other PFMs.)
##Under Uniform and Non-uniform Sampling Schemes with Constant or Variable Integrated Hazard
##(The Line h[j] Could be modified for any Other SSs.)

##The Main Output Quantities of This Code is the Optimal Design Parameters
##(Sample Size n, Sampling intervals based on h1, Control Limit Coefficient L)
#n=x1;  h1=x2;  L=x3
##However, the Quantities Corresponding to each Sampling Interval Can be also obtained;
##If the Code does not Converge before m=m0, It will Stop at m0.
###############

##Control Statistics Distribution

INorm<-readline(prompt="Enter 1 for the Normal Xbar Chart, Otherwise Enter 0:")
INorm<-as.numeric(INorm)
if (INorm!=1 & INorm!=0) {print("Error: Enter 1 or 0")}

INNorm<-readline(prompt="Enter 1 or the Non-normal (Burr) Xbar Chart, Otherwise Enter 0:")
INNorm<-as.numeric(INNorm)
if (INNorm!=1 & INNorm!=0){print("Error: Enter 1 or 0")}

IT2<-readline(prompt="Enter 1 for T2 Chart, Otherwise Enter 0:")
IT2<-as.numeric(IT2)
if (IT2!=1 & IT2!=0){print("Error: Enter 1 or 0")}

if (sum(INorm,INNorm,IT2)!=1){print("Error: You Should choose one and only one of the above models. For Other Control Statistics, Modify the Lines alpha and beta.")}

##Shift Value/Shift Vector for Out-of-Control Status in Univariate/Multivariate Quality Data
##From Banerjee and Rahim (1988) and Yang and Rahim (2005)
delta<-0.5;deltaM=c(1,1.5)
##Burr Parameters
cb=5;kb=6;mb=0.65513;sb=0.16103
##Also, Parameters of the Multivariate Quality Characteristic Distributions
DF=2 #Chi-square Degrees o Freedom
sigmaM=matrix(c(2,1,1,2.5),nrow=2) #Known Covariance Matrix
NCP=as.numeric(t(deltaM)%*%solve(sigmaM)%*%deltaM) #Chi-square Non-central Parameter
###############

##Shock Model/Process Failure Mechanism (PFM) and its Parameters
##with its Expected Value ETa

IW<-readline(prompt="Enter 1 for the Weibull Shock Model, Otherwise Enter 0:")
IW<-as.numeric(IW)
if (IW!=1 & IW!=0){print("Error: Enter 1 or 0")}

ICH<-readline(prompt="Enter 1 for the Chen Shock Model, Otherwise Enter 0:")
ICH<-as.numeric(ICH)
if (ICH!=1 & ICH!=0){print("Error: Enter 1 or 0")}

if (sum(IW,ICH)!=1){print("Error: You Should choose one and only one of the above models. For other PFMs, Modify the Lines ETa, F[1]=p[1], h[j], and R[j].")}

lambda=0.05;nu=2
integral = function(t)
  {t*(lambda*nu*t^(nu-1)*exp(lambda*(1-exp(t^nu))+t^nu))}
ETaCH=integrate(integral,lower=0,upper=Inf)$value
ETaW=(1/lambda)^(1/nu)*gamma(1+(1/nu))
ETa=IW*(ETaW)+ICH*(ETaCH)
###############

##Sampling Scheme

IU<-readline(prompt="Enter 1 for the Uniform SS, Otherwise Enter 0:")
IU<-as.numeric(IU)
if (IU!=1 & IU!=0){print("Error: Enter 1 or 0")}

INUC<-readline(prompt="Enter 1 for Non-uniform SS with Constant Integrated Hazard, Otherwise Enter 0:")
INUC<-as.numeric(INUC)
if (INUC!=1 & INUC!=0){print("Error: Enter 1 or 0")}

INUD1<-readline(prompt="Enter 1 for the Non-uniform SS with (h[j]=[r^(j-1)]*h[1];0<r<1), Otherwise Enter 0:")
INUD1<-as.numeric(INUD1)
if (INUD1!=1 & INUD1!=0){print("Error: Enter 1 or 0")}

INUD2<-readline(prompt="Enter 1 for the Non-uniform SS with (h[j]=[1/(j^r)]*h[1];r>0), Otherwise Enter 0:")
INUD2<-as.numeric(INUD2)
if (INUD2!=1 & INUD2!=0){print("Error: Enter 1 or 0")}

INUD3<-readline(prompt="Enter 1 for the Non-uniform SS with (h[j]=[1/(1+(r/j))^(j-1))]*h[1];r>0), Otherwise Enter 0:")
INUD3<-as.numeric(INUD3)
if (INUD3!=1 & INUD3!=0){print("Error: Enter 1 or 0")}

INUD4<-readline(prompt="Enter 1 for the Non-uniform SS with (h[j]=[1/(1+r*ln(j)]*h[1];r>0), Otherwise Enter 0:")
INUD4<-as.numeric(INUD4)
if (INUD4!=1 & INUD4!=0){print("Error: Enter 1 or 0")}

if (sum(IU,INUC,INUD1,INUD2,INUD3,INUD4)!=1){print("Error: You Should choose one and only one of the above schemes. For other SSs, Modify the Line h[j].")}

##The Decreasing SS Parameter r
r1=0.97;r2=0.5;r3=1.8;r4=1.5
###############

##Model Parameters for Cost and Time
##From Banerjee and Rahim (1988) and Yang and Rahim (2005)

Z0<-0.25;Z1<-1;a<-20;b<-4.22;D0<-50;D1<-950;Y<-500;W<-1100
###############

##Convergence Criteria
TOL=10^-6;m0=5000

###############

##The Nain Function to Compute Expected Cost per Unit Time EC/ET

ECTG <- function(x)
{
  x1<-x[1] #x1 is the sample size n
  
  x2<-x[2] #x2 is the first sampling interval h1
  
  x3<-x[3] #x3 is the control limit coefficient L
  
  alpha<-INorm*(2*pnorm(-x3,0,1))+INNorm*(1+(1/((1+(mb+x3*sb)^cb)^kb))-(1/((1+(mb-x3*sb)^cb)^kb)))+IT2*(1- pchisq(x3, df=DF))
  beta<-INorm*((pnorm(sqrt(x1)*delta+x3,0,1)-pnorm(sqrt(x1)*delta-x3,0,1)))+INNorm*((1/((1+(mb-x3*sb-sb*delta*sqrt(x1))^cb)^kb))-(1/(1+(mb+x3*sb-sb*delta*sqrt(x1))^cb)^kb))+IT2*(pchisq(x3, df=DF, ncp =x1*NCP , lower.tail = TRUE))
  
  h=c();t=c();F=c();p=c();R=c();hR=c();hB=c();tF=c()
  h[1]=t[1]=hR[1]=hB[1]=x2
  F[1]=p[1]=IW*(1-exp(-lambda*(h[1]^nu)))+ICH*(1-exp(lambda*(1-exp((h[1]^nu)))))
  R[1]=1-F[1];tF[1]=t[1]*F[1]
  j=2
  S10=S11=S12=1;S0=S01=S02=0;m=0
##m0 is the number of repetition for while command for the cases the code does not stop till m0
  while(((S10-S0>TOL)|(S11-S01>TOL)|(S12-S02>TOL))&(m<m0)){
    h[j]=IU*h[1]+INUC*IW*((j^(1/nu)-(j-1)^(1/nu))*h[1])+INUC*ICH*((log(1-j*(1-exp(h[1]^nu))))^(1/nu)-(log(1-(j-1)*(1-exp(h[1]^nu))))^(1/nu))+INUD1*(r1^(j-1)*h[1])+INUD2*((1/(j^r2))*h[1])+INUD3*(1/((1+(r3/j))^(j-1))*h[1])+INUD4*((1/(1+r4*log(j)))*h[1])
    t[j]=sum(h[1:j])
    R[j]=IW*(exp(-lambda*(t[j]^nu)))+ICH*(exp(lambda*(1-exp((t[j]^nu)))))
    S0=sum(R[1:j-1])
    S10=sum(R[1:j])
    hR[j]=h[j]*R[j-1]
    S01=sum(hR[1:j-1])
    S11=sum(hR[1:j])
    hB[j]=h[j]*(beta^(j-1))
    F[j]=1-R[j]
    tF[j]=t[j]*(F[j]-F[j-1])
    S02=sum(tF[1:j-1])
    S12=sum(tF[1:j])
    p[j]= (F[j]-F[j-1])/(1-F[j-1])
    
    if(h[j] == "NaN"|t[j] == "NaN"|R[j] == "NaN"|hR[j] == "NaN"|hB[j] == "NaN"|F[j] == "NaN"|tF[j] == "NaN"){S0=1; S10=0; S01=1; S11=0; S02=1; S12=0; S03=1; S13=0; print("Error: There Is NaN")}
    else {S0=sum(R[1:j-1]); S10=sum(R[1:j]); S01=sum(hR[1:j-1]); S11=sum(hR[1:j]); S02=sum(hB[1:j-1]); S12=sum(hB[1:j]); S03=sum(tF[1:j-1]); S13=sum(tF[1:j])}
    
    j=j+1
    m=j-1
  }

  ##The First Series common in Both Expressions
  firstboth=sum(hR[1:m])
  
  ##The Average Number of Sampling When the Process Is in Control
  S=sum(R[1:m])
  
  ##The Third Series common in Both Expressions
  sebothj1=sum(hB[2:m])
  seboth=F[1]*sebothj1
  for(j in 2:(m-1)){
    for(i in (j+1):m){
      seboth = seboth + (beta^(i-j))*(F[j]-F[j-1])*h[i]
    }
  }
  
  ##The Third Series of the Second Expression
  thirdsecond=sum(tF[1:m])
  
  ET=firstboth+alpha*Z0*S+seboth+Z1
  EC=D0*firstboth+alpha*Y*S+(D0-D1)*ETa+(D1-D0)*thirdsecond+D1*seboth+(a+b*x1)*(1/(1-beta)+S)+W
  EC/ET
}

## load alabama package
hin.hs100<-function(x)
{
  f<-numeric(5)
  #f=rep(NA,1)
##Xbar Normal alpha and beta  
  f[1]<-1-(2*pnorm(-x[3],0,1))
  f[2]<-1-(pnorm(sqrt(x[1])*delta+x[3],0,1)-pnorm(sqrt(x[1])*delta-x[3],0,1))
##Xbar Non-Normal alpha and beta  
  #f[1]<-1-(1+(1/((1+(mb+x[3]*sb)^cb)^kb))-(1/((1+(mb-x[3]*sb)^cb)^kb)))
  #f[2]<-1-((1/((1+(mb-x[3]*sb-sb*delta*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta*sqrt(x[1]))^cb)^kb))
##T2 Hotelling alpha and beta  
  #f[1]<-1-(1- pchisq(x[3], df=2))
  #f[2]<-1-(pchisq(x[3], df=DF, ncp=x[1]*NCP, lower.tail=TRUE))
  
  return(f)
}

##Choose one of the following functions for optimization
##xbar starting points
o=auglag(par=c(20,3.2,1.26),ECTG,hin=hin.hs100, control.outer = list(method="nlminb",ilack.max=10))
##T2 starting points
#o=auglag(par=c(20,5,10.2),ECTG,hin=hin.hs100, control.outer = list(method="nlminb",ilack.max=10))

##To use function optim, the line h[j] is recommended to be splitted.
##Xbar starting points
#o=optim(par=c(30,1.5,2),fn=ECTG,method = "L-BFGS-B",lower=c(2,.1,.1),upper = c(100,60,3))
##T2 starting points
#o=optim(par=c(20,5,10.2),fn=ECTG,method = "L-BFGS-B",lower=c(1,.5,.5),upper = c(100,9,20))

o
x1=o$par[1];x2=o$par[2];x3=o$par[3]
x1
x2
x3
x=c(x1,x2,x3)

##Now, run the Code Inside the Function ECTG to
##Compute the Following Quantities

EC
ET
EC/ET
m
ETa
1-beta
alpha
S
ARL1=1/alpha
ARL1
ARL2=1/(1-beta)
ARL2

h[1:m]
t[1:m]
F[1:m]
R[1:m]
p[1:m]

hR[1:m]
hB[1:m]
tF[1:m]
S10-S0
S11-S01
S12-S02
S13-S03
