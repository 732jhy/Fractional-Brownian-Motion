# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 00:06:57 2020

@author: Justin Yu

Implementation of Fractional Brownian Motion, Hosking's method.
"""
import numpy as np


def hosking(T, N, H):
    '''
    Generates sample paths of fractional Brownian Motion using the Davies Harte method
    
    args:
        T:      length of time (in years)
        N:      number of time steps within timeframe
        H:      Hurst parameter
    '''
    gamma = lambda k,H: 0.5*(np.abs(k-1)**(2*H) - 2*np.abs(k)**(2*H) + np.abs(k+1)**(2*H))  
    
    X = [np.random.standard_normal()]
    mu = [gamma(1,H)*X[0]]
    sigsq = [1 - (gamma(1,H)**2)]
    tau = [gamma(1,H)**2]
    
    d = np.array([gamma(1,H)])
    
    for n in range(1, N):
        
        F = np.rot90(np.identity(n+1))
        c = np.array([gamma(k+1,H) for k in range(0,n+1)])
                
        # sigma(n+1)**2
        s = sigsq[n-1] - ((gamma(n+1,H) - tau[n-1])**2)/sigsq[n-1]
        
        # d(n+1)
        phi = (gamma(n+1,H) - tau[n-1])/sigsq[n-1]
        d = d - phi*d[::-1]
        d = np.append(d, phi)        
        
        # mu(n+1) and tau(n+1)
        Xn1 = mu[n-1] + sigsq[n-1]*np.random.standard_normal()
        
        X.append(Xn1)
        sigsq.append(s)
        mu.append(d @ X[::-1])
        tau.append(c @ F @ d)
    
    fBm = np.cumsum(X)*(N**(-H))    
    return (T**H)*fBm




a = hosking(1,200, 0.2);#plt.plot(np.cumsum((a))*np.sqrt(1/200),lw=2)
plt.plot(a,lw=2)
plt.plot(fBm(1,200,0.2),color='green',lw=2)

a = hosking(1,200, 0.5);plt.plot(np.cumsum((a))*np.sqrt(1/200),lw=2)
plt.plot(fBm(1,200,0.5),color='green',lw=2)

a = hosking(1,200, 0.8);plt.plot(np.cumsum((a))*np.sqrt(1/200),lw=2)
plt.plot(fBm(1,200,0.8),color='green',lw=2)




np.random.normal(1.06,-1.304)


gamma = lambda k,H: 0.5*(np.abs(k-1)**(2*H) - 2*np.abs(k)**(2*H) + np.abs(k+1)**(2*H))  
gamma(1,0.2)

n = 3

F = np.rot90(np.identity(n+1))
c = np.array([gamma(k,0.2) for k in range(0,n+1)])

len(F)
len(c)

(F@c)@c



# Testing matrix functions in relation to uppercase gamma (cov matrix)
G = np.zeros((10,10));      H=0.35
for i in range(10):
    for j in range(10):
        G[i,j] = gamma(i-j,H)

c = np.array([gamma(k+1,H) for k in range(10)])


c @ np.linalg.inv(G)

np.linalg.inv(G) @ c


X = np.array([np.random.standard_normal()])
mu = [gamma(1,H)*X[0]]
sigsq = [1 - (gamma(1,H)**2)]
tau = [gamma(1,H)**2]

d = np.array([gamma(1,H)])

for n in range(1, 10):
    
    F = np.rot90(np.identity(n+1))
    c = np.array([gamma(k+1,H) for k in range(0,n+1)])
            
    # sigma(n+1)**2
    s = sigsq[n-1] - ((gamma(n+1,H) - tau[n-1])**2)/sigsq[n-1]
    
    # d(n+1)
    phi = (gamma(n+1,H) - tau[n-1])/sigsq[n-1]
    d = d - phi*d[::-1]
    d = np.append(d, phi)        
    
    # mu(n+1) and tau(n+1)
    Xn1 = mu[n-1] + sigsq[n-1]*np.random.standard_normal()
    
    X.append(Xn1)
    sigsq.append(s)
    mu.append(d @ X)
    tau.append(c @ F @ d)


d 






























