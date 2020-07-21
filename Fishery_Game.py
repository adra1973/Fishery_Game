#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 17:53:12 2019

@author: alandanielroblesaguilar
"""

# Program for modeling a stochastic and dynamical game among two enterprises 
# competing for fishery resources in a certain region.

from random import uniform

from scipy import *
from numpy import *

from scipy import linalg
import numpy as np

R=zeros((11,4,4))
S=zeros((11,4,4))

Q=zeros((8,11))
Q1=zeros((5,11))

T=zeros((11,11,4,4))
U=zeros((11,11,4,4))
V=zeros((11,11,4,4))

V0=zeros((11,1))
V1=zeros((11,1))
W0=zeros((11,1))
W1=zeros((11,1))

H=zeros((11,4,4))
I=zeros((11,4,4))

P=zeros(11)
P1=zeros(11)
P2=zeros(11)

A=zeros(4)
B=zeros(4)


# A[] and B[] are vectores of the available actions for player 1 and 2 respectively
# This actions corresponds to the quantity of fish to be captured by each player.
# The fact that the actions are the same for both players faciliotates the calculations
# and simplify the model.

A[0]=0.0
A[1]=0.1
A[2]=0.2
A[3]=0.3

B[0]=0.0
B[1]=0.1
B[2]=0.2
B[3]=0.3

# p is the value of probability for the Binomial Distribution.

p=0.4


# P[] is a vector of probabilities from the Binimial Distribution                   

for l in range(11):
    P[l]=(math.factorial(10)/(math.factorial(10-l)* \
         math.factorial(l)))*(p**l)*((1-p)**(10-l))    


# P1[] is the vector of probabilities from de Empiric Distribution.

z=0
while z<1000: # "z<" corresponds to the number of iterations for the Empiric Distribution.
    u = uniform(0,1)
    c = float(p)/float((1-p))
    pr = ((1-p)**10) 
    F = pr
    for l in range(11):
        if u<F:
            P2[l]+=1
            break
        pr=pr*c*float(10-l)/float(l+1) 
        F=F+pr               
    for m in range(11):
        P1[m]=float(P2[m])/float((z+1))           
    z+=1



# The function payoffs(A,B) calculate the payoffs matrices for both players.
# Considers as arguments the vector A and B of actions for each player.
# Inckudes a vecrtor S[] of the 11 states of quantity of fish available.

    
def payoffs(A,B):
    U=zeros((11,4,4))
    S=zeros(11)

    
    S[0]=1.0
    S[1]=0.9
    S[2]=0.8
    S[3]=0.7
    S[4]=0.6
    S[5]=0.5
    S[6]=0.4
    S[7]=0.3
    S[8]=0.2
    S[9]=0.1   
    S[10]=0.0
    

    for i in range(11):        
        for j in range(4):            
            for k in range(4):
                if B[k]+A[j]<=S[i]:
                    U[i,j,k]=(A[j])**0.5
                else:
                    U[i,j,k]=(S[i]*A[j]/(B[k]+A[j]))**0.5
                    
                 
                       
    return U


# The function tran10(A,B,P) calculates the Transition Matrices.
# Considers as arguments the vector A and B of actions for each player,
# P or P1 correspond to the probability vectors for Binomial Distribution "P"
# or Empiric Distribution "P1"
# and inckudes a vecrtor S[] of the 11 states of quantity of fish available and
# E[] the vector of discrepances factors, considered for the quantity of fish in each state.

def tran10(A,B,P):   
    
    S=zeros(11)
    E=zeros(11)
    
    
    F=zeros((11,11,4,4))
    G=zeros((11,11,4,4))
    H=zeros((11,11))
    I=zeros((11,11))
    

    S[0]=1.0
    S[1]=0.9
    S[2]=0.8
    S[3]=0.7
    S[4]=0.6
    S[5]=0.5
    S[6]=0.4
    S[7]=0.3
    S[8]=0.2
    S[9]=0.1   
    S[10]=0.0   
        
    E[0]=0.5
    E[1]=0.6
    E[2]=0.7
    E[3]=0.8
    E[4]=0.9
    E[5]=1.0
    E[6]=1.1
    E[7]=1.2
    E[8]=1.3
    E[9]=1.4
    E[10]=1.5

    alfa=0.6
    r1=0.0    
    
    
    
    for l in range(4):
        for m in range(4):
            
    
            for i in range(11):
                for j in range(11):
                    r1=E[i]*max(S[j]-A[l]-B[m],0)**alfa
                    F[j,i,l,m]=10  
                    
                    for k in range(10):
                        if r1> (S[10-k-1]+S[10-k])/2.0: 
                            F[j,i,l,m]=10-k-1  
                            
            
            for i in range(11):
                for j in range(11):
                    for k in range(11):
                        if F[i,k,l,m]==j:
                            G[i,j,l,m]=G[i,j,l,m]+P[k]
                        
                         
                
    
    return G

T=tran10(A,B,P)   # or T=tran10(A,B,P1)        

R=payoffs(A,B) 
S=payoffs(B,A)

beta1=0.75  # beta1 is the discount factor


#  The function sig10(A,n,u,v,w,x,y,z) corresponds to the calculation of the value
# of a game and the strategies of both players as a matricial multiplication.
# A is the game matrix, n is the current state, u,v,w are strategies for player 1; x,y,z 
# are strategies for player 2.

def sig10(A,n,u,v,w,x,y,z):    
    u1=x*(u*A[n,0,0]+v*A[n,1,0]+w*A[n,2,0]+(1-u-v-w)*A[n,3,0])+ y*(u*A[n,0,1]+v*A[n,1,1]+w*A[n,2,1]+(1-u-v-w)*A[n,3,1]) + \
      z*(u*A[n,0,2]+v*A[n,1,2]+w*A[n,2,2]+(1-u-v-w)*A[n,3,2])+ (1-x-y-z)*(u*A[n,0,3]+v*A[n,1,3]+w*A[n,2,3]+(1-u-v-w)*A[n,3,3])
    return u1

# The function val10(A,B,n,u,v,w,x,y,z) corresponds to the McKelvey formula utilized
# in calculation of Nash Equilibrium.
# A and B are payoffs matrices for both player 1 and 2 respectively, n is the current state, u,v,w are strategies for player 1; x,y,z 
# are strategies for player 2.
    
def val10(A,B,n,u,v,w,x,y,z):
    w1=(max(0,sig10(A,n,1,0,0,x,y,z)-sig10(A,n,u,v,w,x,y,z)))**2+(max(0,sig10(A,n,0,1,0,x,y,z)-sig10(A,n,u,v,w,x,y,z)))**2+ \
       (max(0,sig10(A,n,0,0,1,x,y,z)-sig10(A,n,u,v,w,x,y,z)))**2+(max(0,sig10(A,n,0,0,0,x,y,z)-sig10(A,n,u,v,w,x,y,z)))**2+ \
       (max(0,sig10(B,n,u,v,w,1,0,0)-sig10(B,n,u,v,w,x,y,z)))**2+(max(0,sig10(B,n,u,v,w,0,1,0)-sig10(B,n,u,v,w,x,y,z)))**2+ \
       (max(0,sig10(B,n,u,v,w,0,0,1)-sig10(B,n,u,v,w,x,y,z)))**2+(max(0,sig10(B,n,u,v,w,0,0,0)-sig10(B,n,u,v,w,x,y,z)))**2

    return w1



# The function  eq10(A,B,n) helps to calculate the minimum values for the McKelvey 
# function, wich are Nash Equilibriums.
# The arguments are the payoffs matrices A and B for both players, n is the current state.
       
def eq10(A,B,n):
    i1=0.0
    i2=1.0
    j1=0.0
    j2=1.0
    k1=0.0
    k2=1.0
    l1=0.0
    l2=1.0
    m1=0.0
    m2=1.0
    h1=0.0
    h2=1.0

    g=0
    t=0
    e=1.0
    h3=1.0

    a=0.0
    b=0.0
    c=0.0
    d=0.0
    e=0.0
    f=0.0
    
    a1=0.0
    b1=0.0
    c1=0.0
    d1=0.0
    e1=0.0
    f1=0.0

    while g<12: # "g<" corresponds to the precision for the calculation of Nash Equilibrium in decimal figures.
        g+=1
        s=1000.0
        i=i1
        while i <= i2: 
            j=j1
            while j<=min(j2,1.0-i):
                k=k1
                while k<=min(k2,1.0-i-j):

                        
                    t+=1                                                
                    fx=val10(A,B,n,i,j,k,i,j,k)
                    if fx<s:
                        s=fx
                        a1=i
                        b1=j
                        c1=k
                        d1=i
                        e1=j
                        f1=k

                    k+=10**(-g)
                j+=10**(-g) 
            i+=10**(-g)  
        
        

        l2=min(d1+10**(-g),1.0)
        k1=max(c1-10**(-g),0.0)
        k2=min(c1+10**(-g),1.0)
        j1=max(b1-10**(-g),0.0)
        j2=min(b1+10**(-g),1.0)    
        i1=max(a1-10**(-g),0.0)
        i2=min(a1+10**(-g),1.0)
    
    
        e2=((a-a1)**2+(b-b1)**2+(c-c1)**2+(d-d1)**2+(e-e1)**2+(f-f1)**2)**0.5
        f2=abs(val10(A,B,n,a,b,c,d,e,f)-val10(A,B,n,a1,b1,c1,d1,e1,f1))
        h3=val10(A,B,n,a1,b1,c1,d1,e1,f1)
        a=a1
        b=b1
        c=c1
        d=d1
        e=e1
        f=f1
        
    Q[0,n]=a1
    Q[1,n]=b1
    Q[2,n]=c1
    Q[3,n]=d1
    Q[4,n]=e1
    Q[5,n]=f1
    Q[6,n]=sig10(A,n,a1,b1,c1,d1,e1,f1)
    Q[7,n]=sig10(B,n,a1,b1,c1,d1,e1,f1)

    Q1[0,n]=a1
    Q1[1,n]=b1
    Q1[2,n]=c1
    Q1[3,n]=1-a1-b1-c1
    Q1[4,n]=Q[6,n]

    return t


# The output file correspond to a text file with rows for each stage of calculation, 
# the strategies for player 1, the value of the game for each state.
    
file= open("Fishery_Game.txt","w")    
    


z=0
f3=1.0

V0[0,0]=0.0
V0[1,0]=0.0
V0[2,0]=0.0
V0[3,0]=0.0
V0[4,0]=0.0
V0[5,0]=0.0
V0[6,0]=0.0
V0[7,0]=0.0
V0[8,0]=0.0
V0[9,0]=0.0
V0[10,0]=0.0

W0[0,0]=0.0
W0[1,0]=0.0
W0[2,0]=0.0
W0[3,0]=0.0
W0[4,0]=0.0
W0[5,0]=0.0
W0[6,0]=0.0
W0[7,0]=0.0
W0[8,0]=0.0
W0[9,0]=0.0
W0[10,0]=0.0


while z<20: # Here in "z<" we have the number of stages in the determination of Nash Equilibriums
   
    z+=1
    print "Stage", z

    for j in range(11):
        for x in range(4):
            for y in range(4):
                H[j,x,y]=R[j,x,y]+beta1*(V0[0,0]*T[j,0,x,y]+V0[1,0]*T[j,1,x,y]+ \
                 V0[2,0]*T[j,2,x,y]+V0[3,0]*T[j,3,x,y]+V0[4,0]*T[j,4,x,y] + \
                 V0[5,0]*T[j,5,x,y]+V0[6,0]*T[j,6,x,y]+V0[7,0]*T[j,7,x,y] + \
                 V0[8,0]*T[j,8,x,y]+V0[9,0]*T[j,9,x,y]+V0[10,0]*T[j,10,x,y] )
                
                I[j,x,y]=S[j,x,y]+beta1*(W0[0,0]*T[j,0,x,y]+W0[1,0]*T[j,1,x,y]+ \
                 W0[2,0]*T[j,2,x,y]+W0[3,0]*T[j,3,x,y]+W0[4,0]*T[j,4,x,y] + \
                 W0[5,0]*T[j,5,x,y]+W0[6,0]*T[j,6,x,y]+W0[7,0]*T[j,7,x,y] + \
                 W0[8,0]*T[j,8,x,y]+W0[9,0]*T[j,9,x,y]+W0[10,0]*T[j,10,x,y] )
                 
                

        t1=eq10(H,I,j)

    V1[0,0]=Q[6,0]
    V1[1,0]=Q[6,1]
    V1[2,0]=Q[6,2]
    V1[3,0]=Q[6,3]
    V1[4,0]=Q[6,4]
    V1[5,0]=Q[6,5]
    V1[6,0]=Q[6,6]
    V1[7,0]=Q[6,7]
    V1[8,0]=Q[6,8]
    V1[9,0]=Q[6,9] 
    V1[10,0]=Q[6,10] 
    
    W1[0,0]=Q[7,0]
    W1[1,0]=Q[7,1]
    W1[2,0]=Q[7,2]
    W1[3,0]=Q[7,3]
    W1[4,0]=Q[7,4]
    W1[5,0]=Q[7,5]
    W1[6,0]=Q[7,6]
    W1[7,0]=Q[7,7]
    W1[8,0]=Q[7,8]
    W1[9,0]=Q[7,9]
    W1[10,0]=Q[7,10] 

    f3=max(abs(V0[0,0]-V1[0,0]),abs(V0[1,0]-V1[1,0]),abs(V0[2,0]-V1[2,0]), \
          abs(V0[3,0]-V1[3,0]),abs(V0[4,0]-V1[4,0]),abs(V0[5,0]-V1[5,0]), \
          abs(V0[6,0]-V1[6,0]),abs(V0[7,0]-V1[7,0]),abs(V0[8,0]-V1[8,0]), \
          abs(V0[9,0]-V1[9,0]),abs(V0[10,0]-V1[10,0]))
    
    V0[0,0]=V1[0,0]
    V0[1,0]=V1[1,0]
    V0[2,0]=V1[2,0]
    V0[3,0]=V1[3,0]
    V0[4,0]=V1[4,0]
    V0[5,0]=V1[5,0]
    V0[6,0]=V1[6,0]
    V0[7,0]=V1[7,0]
    V0[8,0]=V1[8,0]
    V0[9,0]=V1[9,0]
    V0[10,0]=V1[10,0]
    
    W0[0,0]=W1[0,0]
    W0[1,0]=W1[1,0]
    W0[2,0]=W1[2,0]
    W0[3,0]=W1[3,0]
    W0[4,0]=W1[4,0]
    W0[5,0]=W1[5,0]
    W0[6,0]=W1[6,0]
    W0[7,0]=W1[7,0]
    W0[8,0]=W1[8,0]
    W0[9,0]=W1[9,0]
    W0[10,0]=W1[10,0]


    file.write('% s' %z+ ' , ')

    file.write('% s' %Q1[0,0]+ ' , ')
    file.write('% s' %Q1[1,0]+ ' , ')
    file.write('% s' %Q1[2,0]+ ' , ')
    file.write('% s' %Q1[3,0]+ ' , ')
    file.write('% s' %Q1[4,0]+ ' , ')

    file.write('% s' %Q1[0,1]+ ' , ')
    file.write('% s' %Q1[1,1]+ ' , ')
    file.write('% s' %Q1[2,1]+ ' , ')
    file.write('% s' %Q1[3,1]+ ' , ')
    file.write('% s' %Q1[4,1]+ ' , ')

    file.write('% s' %Q1[0,2]+ ' , ')
    file.write('% s' %Q1[1,2]+ ' , ')
    file.write('% s' %Q1[2,2]+ ' , ')
    file.write('% s' %Q1[3,2]+ ' , ')
    file.write('% s' %Q1[4,2]+ ' , ')

    file.write('% s' %Q1[0,3]+ ' , ')
    file.write('% s' %Q1[1,3]+ ' , ')
    file.write('% s' %Q1[2,3]+ ' , ')
    file.write('% s' %Q1[3,3]+ ' , ')
    file.write('% s' %Q1[4,3]+ ' , ')

    file.write('% s' %Q1[0,4]+ ' , ')
    file.write('% s' %Q1[1,4]+ ' , ')
    file.write('% s' %Q1[2,4]+ ' , ')
    file.write('% s' %Q1[3,4]+ ' , ')
    file.write('% s' %Q1[4,4]+ ' , ')

    file.write('% s' %Q1[0,5]+ ' , ')
    file.write('% s' %Q1[1,5]+ ' , ')
    file.write('% s' %Q1[2,5]+ ' , ')
    file.write('% s' %Q1[3,5]+ ' , ')
    file.write('% s' %Q1[4,5]+ ' , ')
    
    file.write('% s' %Q1[0,6]+ ' , ')
    file.write('% s' %Q1[1,6]+ ' , ')
    file.write('% s' %Q1[2,6]+ ' , ')
    file.write('% s' %Q1[3,6]+ ' , ')
    file.write('% s' %Q1[4,6]+ ' , ')
    
    file.write('% s' %Q1[0,7]+ ' , ')
    file.write('% s' %Q1[1,7]+ ' , ')
    file.write('% s' %Q1[2,7]+ ' , ')
    file.write('% s' %Q1[3,7]+ ' , ')
    file.write('% s' %Q1[4,7]+ ' , ')    
    
    file.write('% s' %Q1[0,8]+ ' , ')
    file.write('% s' %Q1[1,8]+ ' , ')
    file.write('% s' %Q1[2,8]+ ' , ')
    file.write('% s' %Q1[3,8]+ ' , ')
    file.write('% s' %Q1[4,8]+ ' , ')
    
    file.write('% s' %Q1[0,9]+ ' , ')
    file.write('% s' %Q1[1,9]+ ' , ')
    file.write('% s' %Q1[2,9]+ ' , ')
    file.write('% s' %Q1[3,9]+ ' , ')
    file.write('% s' %Q1[4,9]+ ' , ')
    
    file.write('% s' %Q1[0,10]+ ' , ')
    file.write('% s' %Q1[1,10]+ ' , ')
    file.write('% s' %Q1[2,10]+ ' , ')
    file.write('% s' %Q1[3,10]+ ' , ')
    file.write('% s' %Q1[4,10]+ ' , ')    


    file.write('\n')

file.close()