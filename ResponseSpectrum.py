#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 16:30:41 2020

@author: deepakpant
"""

# import packages
import numpy as np
import matplotlib.pyplot as plt

# load input files and get data
X = np.loadtxt(fname = "infile.txt") #accleration file in two column format
T = np.loadtxt(fname = "period.txt") #load period for which the spectrum is to be computed
time = X[:,0]
accln = X[:,1]

# define damping ratio
damp = 0.05

# get number of points and delata t for time step
npts = len(time)
dt = time[1]-time[0]

# plot acceleration to check
plt.figure()
plt.plot(time,accln)

# function to calculate response sepctrum using Newmark's method
def RSNewmark(Ts,damp,accln,npts,dt,time):
    '''
    function to calculate response sepctrum using Newmark's emthod
    Ts: time
    damp: daming ratio
    accln: accleration vector
    npts: total number of points
    dt: time step
    time: time vector
    '''
    pi = 3.14159265359
    m = 1
    p = -m*accln
    omega = 2*pi/Ts
    c = 2*damp*m*omega
    k = omega**2*m
    khat = k+2*c/dt+4*m/(dt**2)
    a = 4*m/dt+2*c
    b = 2*m
    
    timenew = time[0:npts-1]
    
    acc = np.zeros(npts)
    vel = np.zeros(npts)
    dis = np.zeros(npts)
    p[0] = 0.0
    deltaphat = np.zeros(npts)
    deltadisp = np.zeros(npts)
    deltavel = np.zeros(npts)
    deltaacc = np.zeros(npts)
    
    for n, item in enumerate(timenew):
        deltaphat[n]=(p[n+1]-p[n])+a*vel[n]+b*acc[n];
        deltadisp[n]=deltaphat[n]/khat;
        deltavel[n]=2*deltadisp[n]/dt-2*vel[n];
        deltaacc[n]=4*(deltadisp[n]-dt*vel[n])/dt**2-2*acc[n];
    
        dis[n+1]=dis[n]+deltadisp[n]
        vel[n+1]=vel[n]+deltavel[n]
        acc[n+1]=acc[n]+deltaacc[n]
    
    resp = max(abs(dis))*981       
    
    return resp

# initialize variables    
resp = np.zeros(len(T))
pacc = np.zeros(len(T))
pvel = np.zeros(len(T))

# loop through each time step and compute response spectrum
for n, item in enumerate(T):
    resp[n] = RSNewmark(T[n],damp,accln,npts,dt,time)
    pi = 3.14159265359
    omega = 2*pi/T[n]   
    pacc[n]=resp[n]*omega**2/981
    pvel[n]=resp[n]*omega

# plot response spectra
plt.figure()
plt.plot(T,resp)
plt.figure()
plt.plot(T,pacc)
plt.figure()
plt.plot(T,pvel)


