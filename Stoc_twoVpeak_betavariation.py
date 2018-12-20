import random
import matplotlib.pyplot as plt
import numpy as np

b2 = 3.2e-5
k1 = 4.6
k2 = 4.6
p1 = 4.6e-2
p2 = 4.6e-2
d1 = 5.2
d2 = 5.2
c1= 5.2 
c2= 5.2
T=[0]
E1=[0]
I1=[0]
V1=[0]
E2=[0]
I2=[0]
V2=[0]
tt=[0]
#Initial states of the model
T[0]=4.0e8#1.0#
E1[0]=0.0
I1[0]=0.0
V1[0]=7.5e-2#1.0#
E2[0]=0.0
I2[0]=0.0
V2[0]=7.5e-2#1.0#
tt[0]=0.0
Vpeak1=[]
Vpeak2=[]
ipeak1=[]
ipeak2=[]
tpeak1=[]
tpeak2=[]
tc=[]
tcoin=[]
imax=1000
br=np.logspace(-2,0.3,2)#b1 = 3.2e-5 
for ii in range (len(br)):
    b1=br[ii]
    for i in range (0,imax):#function goes from start to (end-1), inclusive
        dt=random.uniform(1e-6,1e-2)
        rem1 = np.fmod((b1*T[i]*V1[i]*dt),1.0)
        rem12 = np.fmod((b2*T[i]*V2[i]*dt),1.0)
        rem2 = np.fmod((k1*E1[i]*dt),1.0)
        rem22 = np.fmod((k2*E2[i]*dt),1.0)
        rem3 = np.fmod((d1*I1[i]*dt),1.0)
        rem32 = np.fmod((d2*I2[i]*dt),1.0)
        rem4 = np.fmod((p1*I1[i]*dt),1.0)
        rem42 = np.fmod((p2*I2[i]*dt),1.0)
        rem5 = np.fmod((c1*V1[i]*dt),1.0)
        rem52 = np.fmod((c2*V2[i]*dt),1.0)
        T_out1 = np.floor(b1*T[i]*V1[i]*dt)+random.uniform(0,rem1)
        T_out2 = np.floor(b2*T[i]*V2[i]*dt)+random.uniform(0,rem12)
        E_out1 = np.floor(k1*E1[i]*dt)+random.uniform(0,rem2)
        E_out2 = np.floor(k2*E2[i]*dt)+random.uniform(0,rem22)
        I_out1 = np.floor(d1*I1[i]*dt)+random.uniform(0,rem3)
        I_out2 = np.floor(d2*I2[i]*dt)+random.uniform(0,rem32)
        V_in1 = np.floor(p1*I1[i]*dt)+random.uniform(0,rem4)
        V_in2 = np.floor(p2*I2[i]*dt)+random.uniform(0,rem42)
        V_out1 = np.floor(c1*V1[i]*dt)+random.uniform(0,rem5) 
        V_out2 = np.floor(c2*V2[i]*dt)+random.uniform(0,rem52) 
        T.append(-T_out1-T_out2+T[i])
        E1.append(T_out1-E_out1+E1[i])
        E2.append(T_out2-E_out2+E2[i])
        I1.append(E_out1-I_out1+I1[i])
        I2.append(E_out2-I_out2+I2[i])
        V1.append(V_in1-V_out1+V1[i])
        V2.append(V_in2-V_out2+V2[i])
        tt.append(tt[i]+dt)
        if V1[i]>1.0:
            if V2[i]>1.0:
		      tc.append(tt[i])
            else:
                pass
        else:
		   pass
    #plt.semilogy(tt,V1,'r',tt,V2,'b')
    #plt.show(True)
    Vpeak1.append(max(V1))
    Vpeak2.append(max(V2))
    #Coinfection duration:
    tcoin.append(max(tc)-min(tc))
    #Time of peak virus load:
    ipeak2.append(np.argmax(V2))
    ipeak1.append(np.argmax(V1))
    tpeak1 = [tt[k] for k in ipeak1]
    tpeak2 = [tt[k] for k in ipeak2]
    
plt.figure()
plt.semilogy(tt,V1,'r',tt,V2,'b',linewidth=3)
#plt.ylim(1e-2,1e7)
#plt.xlim(0,10)
plt.ylabel('Viral titer', fontsize=20)
plt.xlabel('Time post infection (day)', fontsize=20)
plt.legend(('Virus 1', 'Virus 2'),'upper right',fontsize=10)
#plt.title('Two virus stochastic model',fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='minor', labelsize=10)
plt.show()

######END#####OF##### j for loop #######        
print '###############Variation in Vpeak with Variation in Growth Rates###################'        
plt.figure()
plt.loglog(br,Vpeak1,'r',br,Vpeak2,'b',linewidth=3)
#plt.ylim(1e-2,1e7)
#plt.xlim(0,10)
plt.ylabel('Viral peak load', fontsize=20)
plt.xlabel('Growth rate (infection rate)', fontsize=20)
plt.legend(('Virus 1', 'Virus 2'),'upper right',fontsize=18)
#plt.title('Two virus stochastic model',fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='minor', labelsize=10)
plt.show()

print '#######Variation in tpeak with Variation in Growth Rates########' 
print 'Time of peak1=', tpeak1
print 'Time of peak2=', tpeak2     
plt.figure()
plt.semilogx(br,tpeak1,'r',br,tpeak2,'b',linewidth=3)
#plt.ylim(1e-2,1e7)
#plt.xlim(0,10)
plt.ylabel('Time of peak', fontsize=20)
plt.xlabel('Growth rate (infection rate)', fontsize=20)
plt.legend(('Virus 1', 'Virus 2'),'upper right',fontsize=10)
#plt.title('Two virus stochastic model',fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='minor', labelsize=10)
plt.show()

print '#######Variation in tcoin with Variation in Growth Rates#######'
print 'Duration of coinfection=', tcoin     
plt.figure()
plt.semilogx(br,tcoin,'k',linewidth=3)
#plt.ylim(1e-2,1e7)
#plt.xlim(0,10)
plt.ylabel('coinfection duration', fontsize=20)
plt.xlabel('Growth rate (infection rate)', fontsize=20)
#plt.legend(('Virus 1', 'Virus 2'),'upper right',fontsize=10)
#plt.title('Two virus stochastic model',fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='minor', labelsize=10)
plt.show()


