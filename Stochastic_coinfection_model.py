# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 21:08:15 2018

@author: lubna
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import pylab as pl
import seaborn as sns; sns.set(color_codes=True)

# Model parameters
CF=10**0
b1,k1,d1,p1,c1,b2,k2,d2,p2,c2 = 3.2e-5*1e3/CF,4.6,5.2,4.6e-2*CF,5.2,3.2e-5/CF,4.6,5.2,4.6e-2*CF,5.2
#T0,E10, I10, V10, E20, I20, V20=4e8,0,0,7.5e-2*CF,0,0,7.5e-2*CF
T0,E10, I10, V10, E20, I20, V20=4e8,0,0,1,0,0,1
tau=0.001 #dt
Max=10 #maximum days of infection
INPUT = np.array((T0,E10, I10, V10, E20, I20, V20))

def stoc_eqs_tauleap(INP): 
	X = INP
	Rate=np.zeros((10))##propensity function [b1*T*V1 k1*E1 d1*I1 p1*I1 c1*V1 b2*T*V2 k2*E2 d2*I2 p2*I2 c2*V2]
	Transitions=np.zeros((10,7))##stoichiometric matrix, each row of which is a transition vector
	Rate[0] = b1*X[0]*X[3]; Transitions[0,:]=([-1, +1, 0, 0, 0, 0, 0]);
	Rate[1] = k1*X[1];  Transitions[1,:]=([0, -1, +1, 0, 0, 0, 0]);
	Rate[2] = d1*X[2];  Transitions[2,:]=([0, 0, -1, 0, 0, 0, 0]);
	Rate[3] = p1*X[2];  Transitions[3,:]=([0, 0, 0, +1, 0, 0, 0]);
	Rate[4] = c1*X[3];  Transitions[4,:]=([0, 0, 0, -1, 0, 0, 0]);
	Rate[5] = b2*X[0]*X[6]; Transitions[5,:]=([-1, 0, 0, 0, +1, 0, 0]);
	Rate[6] = k2*X[4];  Transitions[6,:]=([0, 0, 0, 0, -1, +1, 0]);
	Rate[7] = d2*X[5];  Transitions[7,:]=([0, 0, 0, 0, 0, -1, 0]);
	Rate[8] = p2*X[5];  Transitions[8,:]=([0, 0, 0, 0, 0, 0, +1]);
	Rate[9] = c2*X[6];  Transitions[9,:]=([0, 0, 0, 0, 0, 0, -1]);
	for i in range(10):
		leap=np.random.poisson(Rate[i]*tau);
		## To avoid negative numbers
		Use=min([leap, X[pl.find(Transitions[i,:]<0)]]);
		X=X+Transitions[i,:]*Use;
	return X

def Stoch_Iteration(INPUT):
    tt=0#initiate state vector to
    T=[0] #initiate state vector for all the cells
    E1=[0]
    I1=[0]
    V1=[0]
    E2=[0]
    I2=[0]
    V2=[0]
    for tt in time:
        res=stoc_eqs_tauleap(INPUT)
        T.append(INPUT[0])
        E1.append(INPUT[1])
        I1.append(INPUT[2])
        V1.append(INPUT[3])
        E2.append(INPUT[4])
        I2.append(INPUT[5])
        V2.append(INPUT[6])
        INPUT=res
    return [T, E1, I1, V1, E2, I2, V2]

plt.figure(figsize=(12, 9))
jpeak1=[]
jpeak2=[]
Vpeak1=[]
Vpeak2=[]
ipeak1=[]
ipeak2=[]
tpeak1=[]
tpeak2=[]
tc=[]
tcoin=[]
Extinction_co=[]
n_simulations=1000
for j in range (1,n_simulations+1):
    print('simulation=',j)
    time=np.arange(0.0, Max, tau)
    [T, E1, I1, V1, E2, I2, V2]=Stoch_Iteration(INPUT)

    t=np.array(time)
    
    tT=np.array(T)[1:,]
    
    tV1=np.array(V1)[1:,]
    tV2=np.array(V2)[1:,]
    
    tE1=np.array(E1)[1:,]
    tE2=np.array(E2)[1:,]
    
    tI1=np.array(I1)[1:,]
    tI2=np.array(I2)[1:,]
    #Coinfection duration
#    for i in range(len(t)):
#        if tV1[i]>1.0:
#            if tV2[i]>1.0:
#                tc.append(t[i])
#            else:
#                pass
#        else:
#            pass
#    if len(tc) == 0:
#        pass
#    else:
#        tcoin.append(max(tc)-min(tc))
    #print(tcoin)
    #Times of virus peak for both viruses
#    if max(tV1)>max(tV2):
#        jpeak1.append(j)
#    else:
#        jpeak2.append(j)
#    ipeak1.append(np.argmax(tV1))
#    ipeak2.append(np.argmax(tV2))
#    tpeak1 = [t[k] for k in ipeak1]
#    tpeak2 = [t[k] for k in ipeak2]
#    Vpeak1.append(max(V1))
#    Vpeak2.append(max(V2))
    #calculation of Extinction probability
    if max(tV1)>1.0:
        if max(tV2)>1.0:
           Extinction_co.append(j)
        else:
            pass
    else:
        pass
    #plt.figure(figsize=(12, 9))
    plt.semilogy(t,tV1,'r',t,tV2,'b',linewidth=2)
    plt.ylabel('Viral titer, [V]', fontsize=35)
    plt.xlabel('Time post infection (day)', fontsize=35)
    plt.tick_params(axis='both', which='major', labelsize=30)
    plt.tick_params(axis='both', which='minor', labelsize=25)
    plt.legend(('Virus 1','Virus 2',),fontsize=30)
#plt.text(5,10**6, r' $\alpha$ = $10^4$', fontsize=30)
#plt.text(5,10**6, r' No conversion', fontsize=30)
plt.show()
####################################################################### 
#print("Mean of tcoin=", np.mean(tcoin))
#print("Std of tcoin=",np.std(tcoin))
#
#print("Mean of tpeak1=", np.mean(tpeak1))
#print("Std of tpeak1=",np.std(tpeak1))
#
#print("Mean of tpeak2=", np.mean(tpeak2))
#print("Std of tpeak2=",np.std(tpeak2))
###################################################################
#Min, Maxx = 1e0, 1e7
##Min, Maxx = 1e0, 1e10
#
#plt.figure(figsize=(12, 9))  
#plt.hist(Vpeak1, bins=10**np.linspace(np.log10(Min), np.log10(Maxx), 10), color="indianred")
#plt.gca().set_xscale("log")
#plt.xlabel(r'Peak value of $Virus_{1}$, [V]', fontsize=35)
#plt.ylabel('No. of occurrences', fontsize=35)
##plt.title('Distribution for Virus 1',fontsize=20)
#plt.tick_params(axis='both', which='major', labelsize=30)
#plt.tick_params(axis='both', which='minor', labelsize=25)
#plt.show()
#
#
#plt.figure(figsize=(12, 9))  
#plt.hist(Vpeak2, bins=10**np.linspace(np.log10(Min), np.log10(Maxx), 10), color="salmon")
#plt.gca().set_xscale("log")
#plt.xlabel(r'Peak value of $Virus_{2}$, [V]', fontsize=35)
#plt.ylabel('No. of occurrences', fontsize=35)
##plt.title('Distribution for Virus 2',fontsize=20)
#plt.tick_params(axis='both', which='major', labelsize=30)
#plt.tick_params(axis='both', which='minor', labelsize=25)
#plt.show()
#
#plt.figure(figsize=(12, 9))  
#plt.hist(tcoin, bins=np.linspace(5, 7, 10), color="forestgreen")
#plt.xlabel('Duration of coinfection (day)', fontsize=35)
#plt.ylabel('No. of occurrences', fontsize=35)
##plt.title('Distribution for Virus 1',fontsize=20)
#plt.tick_params(axis='both', which='major', labelsize=30)
#plt.tick_params(axis='both', which='minor', labelsize=25)
#plt.show()
#
#plt.figure(figsize=(12, 9))  
#plt.hist(tpeak1, color="royalblue")
#plt.xlabel('Time of peak viral load of $Virus_{1}$ (day)', fontsize=35)
#plt.ylabel('No. of occurrences', fontsize=35)
##plt.title('Distribution for Virus 1',fontsize=20)
#plt.tick_params(axis='both', which='major', labelsize=30)
#plt.tick_params(axis='both', which='minor', labelsize=25)
#plt.show()
#
#
#plt.figure(figsize=(12, 9))  
#plt.hist(tpeak2, color="skyblue")
#plt.xlabel('Time of peak viral load of $Virus_{2}$ (day)', fontsize=35)
#plt.ylabel('No. of occurrences', fontsize=35)
##plt.title('Distribution for Virus 2',fontsize=20)
#plt.tick_params(axis='both', which='major', labelsize=30)
#plt.tick_params(axis='both', which='minor', labelsize=25)
#plt.show()
#
#print('##################################')
#plt.figure(figsize=(12, 9))
#plt.hist(jpeak1)
#plt.xlabel('Peak viral load of $Virus_{1}$', fontsize=35)
#plt.ylabel('No. of occurrences', fontsize=35)
##plt.title('Distribution for Virus 1',fontsize=20)
#plt.tick_params(axis='both', which='major', labelsize=30)
#plt.tick_params(axis='both', which='minor', labelsize=25)
#plt.show()
#print('##################################')
#plt.figure(figsize=(12, 9))
#plt.hist(jpeak2)
#plt.xlabel('Peak viral load of $Virus_{2}$', fontsize=35)
#plt.ylabel('No. of occurrences', fontsize=35)
##plt.title('Distribution for Virus 2',fontsize=20)
#plt.tick_params(axis='both', which='major', labelsize=30)
#plt.tick_params(axis='both', which='minor', labelsize=25)
#plt.show()
#################################################################
print('------Extinction probability --------')
a=len(Extinction_co)
b=a/n_simulations 
print('No of times both Virus have peaks above detection threshold=', a)
print('Probability of disease outbreak =', b) 
print('Probability of disease extinction=', 1-b)
#
#################################################################

#data1=np.column_stack((tcoin))
#np.savetxt("data_stoc_tcoin_fixed_V0.dat",data1)
#data2=np.column_stack((tpeak1))
#np.savetxt("data_stoc_tpeak_V1_fixed_V0.dat",data2)
#data3=np.column_stack((tpeak2))
#np.savetxt("data_stoc_tpeak_V2_fixed_V0.dat",data3)
#data2=np.column_stack((Vpeak1))
#np.savetxt("data_stoc_Vpeak1_fixed_V0.dat",data2)
#data3=np.column_stack((Vpeak2))
#np.savetxt("data_stoc_Vpeak2_fixed_V0.dat",data3)
#data4=np.column_stack((jpeak1))
#np.savetxt("data_stoc_jpeak_V1_fixed_V0.dat",data4)
#data5=np.column_stack((jpeak2))
#np.savetxt("data_stoc_jpeak_V2_fixed_V0.dat",data5)







