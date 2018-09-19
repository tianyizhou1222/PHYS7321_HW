# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 16:50:22 2018
@author: Tom K
"""


def forwardiff(f,a,b,N):
    h = (b-a)/N
    g=[]
    for k in range(1,N):
        slop = (f(a+k*h)-f(a+(k-1)*h))/h
        g.append(slop)
    return g
    
def backwardiff(f,a,b,N):
    h = (b-a)/N
    g=[]
    for k in range(0,N-1):
        slop = (f(a+(k+1)*h)-f(a+k*h))/h
        g.append(slop)
    return g

def centraldiff(f,a,b,N):
    h = (b-a)/N
    g=[]
    for k in range(0,N):
        slop = (f(a+(k+1/2)*h)-f(a+(k-1/2)*h))/h
        g.append(slop)
    return g

def hodiff(f,a,b,N):
    h = (b-a)/N
    g=[]
    for k in range(0,N):
        slop = (-2/3*f(a+(k+1)*h)+2/3*f(a+(k-1)*h)-1/12*f(a+(k+2)*h)+1/12*f(a+(k-2)*h))/h
        g.append(slop)
    return g

def central2diff(f,a,b,N):
    h = (b-a)/N
    g=[]
    for k in range(0,N):
        slop = (f(a+k*h)- 2*f(a) +f(a+k*h))/h**2
        g.append(slop)
    return g