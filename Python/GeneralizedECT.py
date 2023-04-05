#!/usr/bin/env python
# coding: utf-8

# In[182]:


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

import numpy as np
from scipy.special import gamma
from scipy.integrate import quad
from scipy.stats import norm, ncx2, chi2
from scipy.optimize import minimize


# In[31]:


##Control Statistics Distribution

# Prompt user for input
INorm = int(input("Enter 1 for the Normal Xbar Chart, Otherwise Enter 0: "))
if INorm not in [0, 1]:
    print("Error: Enter 1 or 0")


# In[7]:


INNorm = int(input("Enter 1 or the Non-normal (Burr) Xbar Chart, Otherwise Enter 0: "))
if INNorm not in [0, 1]:
    print("Error: Enter 1 or 0")


# In[5]:


IT2 = int(input("Enter 1 for T2 Chart, Otherwise Enter 0: "))
if IT2 not in [0, 1]:
    print("Error: Enter 1 or 0")


# In[8]:


# Check if user input is valid
if sum([INorm, INNorm, IT2]) != 1:
    print("Error: You Should choose one and only one of the above models. For Other Control Statistics, Modify the Lines alpha and beta.")


# In[14]:


# Shift Value/Shift Vector for Out-of-Control Status in Univariate/Multivariate Quality Data
# From Banerjee and Rahim (1988) and Yang and Rahim (2005)
delta = 0.5
deltaM = np.array([1, 1.5])

# Burr Parameters
cb = 5
kb = 6
mb = 0.65513
sb = 0.16103

# Parameters of the Multivariate Quality Characteristic Distributions
DF = 2  # Chi-square Degrees o Freedom
sigmaM = np.array([[2, 1], [1, 2.5]])  # Known Covariance Matrix
NCP = float(np.transpose(deltaM) @ np.linalg.inv(sigmaM) @ deltaM)  # Chi-square Non-central Parameter


# In[15]:


# Shock Model/Process Failure Mechanism (PFM) and its Parameters
# with its Expected Value ETa
IW = int(input("Enter 1 for the Weibull Shock Model, Otherwise Enter 0: "))
if IW not in [0, 1]:
    print("Error: Enter 1 or 0")


# In[16]:


ICH = int(input("Enter 1 for the Chen Shock Model, Otherwise Enter 0: "))
if ICH not in [0, 1]:
    print("Error: Enter 1 or 0")


# In[17]:


# Check if user input is valid
if sum([IW, ICH]) != 1:
    print("Error: You Should choose one and only one of the above models. For other PFMs, Modify the Lines ETa, F[1]=p[1], h[j], and R[j].")


# In[18]:


lambda_= 0.05
nu = 2

def integral(t):
    return t * (lambda_ * nu * t ** (nu - 1) * np.exp(lambda_ * (1 - np.exp(t ** nu)) + t ** nu))

ETaCH = quad(integral, 0, np.inf)[0]
ETaW = (1 / lambda_) ** (1 / nu) * gamma(1 + (1 / nu))
ETa = IW * ETaW + ICH * ETaCH


# In[ ]:





# In[21]:


# Sampling Scheme

IU = input("Enter 1 for the Uniform SS, Otherwise Enter 0:")
IU = int(IU)
if IU != 1 and IU != 0:
    print("Error: Enter 1 or 0")


# In[22]:


INUC = input("Enter 1 for Non-uniform SS with Constant Integrated Hazard, Otherwise Enter 0:")
INUC = int(INUC)
if INUC != 1 and INUC != 0:
    print("Error: Enter 1 or 0")


# In[23]:


INUD1 = input("Enter 1 for the Non-uniform SS with (h[j]=[r^(j-1)]*h[1];0<r<1), Otherwise Enter 0:")
INUD1 = int(INUD1)
if INUD1 != 1 and INUD1 != 0:
    print("Error: Enter 1 or 0")


# In[24]:


INUD2 = input("Enter 1 for the Non-uniform SS with (h[j]=[1/(j^r)]*h[1];r>0), Otherwise Enter 0:")
INUD2 = int(INUD2)
if INUD2 != 1 and INUD2 != 0:
    print("Error: Enter 1 or 0")


# In[25]:


INUD3 = input("Enter 1 for the Non-uniform SS with (h[j]=[1/(1+(r/j))^(j-1))]*h[1];r>0), Otherwise Enter 0:")
INUD3 = int(INUD3)
if INUD3 != 1 and INUD3 != 0:
    print("Error: Enter 1 or 0")


# In[26]:


INUD4 = input("Enter 1 for the Non-uniform SS with (h[j]=[1/(1+r*ln(j)]*h[1];r>0), Otherwise Enter 0:")
INUD4 = int(INUD4)
if INUD4 != 1 and INUD4 != 0:
    print("Error: Enter 1 or 0")


# In[27]:


if sum([IU, INUC, INUD1, INUD2, INUD3, INUD4]) != 1:
    print("Error: You Should choose one and only one of the above schemes. For other SSs, Modify the Line h[j].")


# In[28]:


# The Decreasing SS Parameter r
r1 = 0.97
r2 = 0.5
r3 = 1.8
r4 = 1.5

# Model Parameters for Cost and Time
# From Banerjee and Rahim (1988) and Yang and Rahim (2005)
Z0 = 0.25
Z1 = 1
a = 20
b = 4.22
D0 = 50
D1 = 950
Y = 500
W = 1100

# Convergence Criteria
TOL = 10**(-6)
m0 = 5000


# In[187]:


def ECTG(x):
    x1, x2, x3 = x
#    x1 = x[0] #x1 is the sample size n
#    x2 = x[1] #x2 is the first sampling interval h1 (in python h[0])
#    x3 = x[2] #x3 is the control limit coefficient L
    
    alpha = INorm*(2*norm.cdf(-x3, 0, 1))+INNorm*(1+(1/((1+(mb+x3*sb)**cb)**kb))-(1/((1+(mb-x3*sb)**cb)**kb)))+IT2*(1- chi2.cdf(x3, df=DF))
    beta = INorm*((norm.cdf(np.sqrt(x1)*delta+x3, 0, 1)-norm.cdf(np.sqrt(x1)*delta-x3, 0, 1)))+INNorm*((1/((1+(mb-x3*sb-sb*delta*np.sqrt(x1))**cb)**kb))-(1/(1+(mb+x3*sb-sb*delta*np.sqrt(x1))**cb)**kb))+IT2*(ncx2.cdf(x3, df=DF, nc =x1*NCP))
    
    h,t,F,p,R,hR,hB,tF = [], [], [], [], [], [], [], []
    h.append(x2)
    t.append(x2)
    hR.append(x2)
    hB.append(x2)
    F.append(IW*(1-np.exp(-lambda_*(h[0]**nu)))+ICH*(1-np.exp(lambda_*(1-np.exp((h[0]**nu))))))
    p.append(F[0])
    R.append(1-F[0])
    tF.append(t[0]*F[0])
    j = 2
    S10, S11, S12 = 1, 1, 1
    S0, S01, S02 = 0, 0, 0
    m = 0
    
    while (((S10 - S0 > TOL) or (S11 - S01 > TOL) or (S12 - S02 > TOL)) and (m < m0)):
      h.append (IU * h[0] +
         INUC * IW * ((j ** (1/nu) - (j-1) ** (1/nu)) * h[0]) +
         INUC * ICH * ((np.log(1 - j * (1 - np.exp(h[0] ** nu)))) ** (1/nu) - (np.log(1 - (j-1) * (1 - np.exp(h[0] ** nu)))) ** (1/nu)) +
         INUD1 * (r1 ** (j-1) * h[0]) +
         INUD2 * ((1/(j ** r2)) * h[0]) +
         INUD3 * (1 / ((1 + (r3 / j)) ** (j-1)) * h[0]) +
         INUD4 * ((1 / (1 + r4 * np.log(j))) * h[0]))
      t.append(sum(h))
#      S0 = sum(R)
      R.append(IW * (np.exp(-lambda_ * (t[j-1] ** nu))) + ICH * (np.exp(lambda_ * (1 - np.exp((t[j-1] ** nu))))))
#      S10 = sum(R)
#      S01 = sum(hR)
      hR.append(h[j-1] * R[j-2])
#      S11 = sum(hR)
      hB.append(h[j-1] * (beta ** (j-1)))
      F.append(1 - R[j-1])
#      S02 = sum(tF)
      tF.append(t[j-1] * (F[j-1] - F[j-2]))
#      S12 = sum(tF)
      p.append((F[j-1] - F[j-2]) / (1 - F[j-2]))
    
      if (np.isnan(h[j-1]) or np.isnan(t[j-1]) or np.isnan(R[j-1]) or np.isnan(hR[j-1]) or np.isnan(hB[j-1]) or np.isnan(F[j-1]) or np.isnan(tF[j-1])):
        S0 = 1
        S10 = 0
        S01 = 1
        S11 = 0
        S02 = 1
        S12 = 0
        S03 = 1
        S13 = 0
        print("Error: There is NaN")
      else:
        S0 = sum(R[0:j-1])
        S10 = sum(R[0:j])
        S01 = sum(hR[0:j-1])
        S11 = sum(hR[0:j])
        S02 = sum(hB[0:j-1])
        S12 = sum(hB[0:j])
        S03 = sum(tF[0:j-1])
        S13 = sum(tF[0:j])
    
      j = j + 1
      m = j - 1
    
    ##The First Series common in Both Expressions
    firstboth = sum(hR[0:m])

    ##The Average Number of Sampling When the Process Is in Control
    S = sum(R[0:m])

    ##The Third Series common in Both Expressions
    sebothj1 = sum(hB[1:m])
    seboth = F[0] * sebothj1
    for j in range(2, m):
      for i in range(j+1, m+1):
        seboth += (beta**(i-j)) * (F[j-1]-F[j-2]) * h[i-1]

    ##The Third Series of the Second Expression
    thirdsecond = sum(tF[0:m])

    ET = firstboth + alpha * Z0 * S + seboth + Z1
    EC = D0 * firstboth + alpha * Y * S + (D0-D1) * ETa + (D1-D0) * thirdsecond + D1 * seboth + (a+b*x1) * (1/(1-beta) + S) + W
    return EC/ET


# In[189]:



def constraint0(x):
    x1, x2, x3 = x
    return 1 - (2 * norm.cdf(-x[2], 0, 1))

def constraint1(x):
    x1, x2, x3 = x
    return 1 - (norm.cdf(np.sqrt(x[0]) * delta + x[2], 0, 1) - norm.cdf(np.sqrt(x[0]) * delta - x[2], 0, 1))

    
    ## Xbar Normal alpha and beta  
    #f[0] = 1 - (2 * norm.cdf(-x3, 0, 1))
    #f[1] = 1 - (norm.cdf(np.sqrt(x1) * delta + x3, 0, 1) - norm.cdf(np.sqrt(x1) * delta - x3, 0, 1))
    
    ## Xbar Non-Normal alpha and beta  
    #f[0] = 1 - (1 + (1 / ((1 + (mb + x[2] * sb) ** cb) ** kb)) - (1 / ((1 + (mb - x[2] * sb) ** cb) ** kb)))
    #f[1] = 1 - ((1 / ((1 + (mb - x[2] * sb - sb * delta * np.sqrt(x[0])) ** cb) ** kb)) - (1 / (1 + (mb + x[2] * sb - sb * delta * np.sqrt(x[0])) ** cb) ** kb))
    
    ## T2 Hotelling alpha and beta  
    #f[0] = 1 - (1 - chi2.cdf(x[2], df=2))
    #f[1] = 1 - (ncx2.cdf(x3, df=DF, nc =x1*NCP, loc=0, scale=1) if x[2] >= 0 else 1 - chi2.cdf(-x[2], df=DF, nc=x[0]*NCP, loc=0, scale=1))
                                    
                                    
cons = ({'type': 'ineq', 'fun': constraint0},
        {'type': 'ineq', 'fun': constraint1},
        {'type': 'ineq', 'fun': lambda x: x[0]},  # x1 >= 0
        {'type': 'ineq', 'fun': lambda x: x[1]},  # x2 >= 0
        {'type': 'ineq', 'fun': lambda x: x[2]},  # x3 >= 0
        )

# Set the initial parameter values
par = [20, 3.2, 1.26]


# Run the optimization
result = minimize(ECTG, par, method='L-BFGS-B', constraints=cons)

# Print the optimized parameter values
print(result)


# In[169]:





# In[156]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




