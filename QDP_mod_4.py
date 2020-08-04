#!/usr/bin/env python
# coding: utf-8

# In[8]:


import matplotlib.pyplot as plt
import numpy as np
from math import exp
import sys
import time


# In[9]:


def Poten(x):
    C_0=C_2**2/(4*C_4)#term to make the potential positieve always
    V=K*(C_4*(x**4)-C_2*(x**2)+C_0)
    return V
def F(x,E):
    v=Poten(x)
    return 2*m*(v-E)/h**2


# In[10]:


def Numrov(x,q0,q1,psi1,f1,dx,E):
    q2=(dx**2)*f1*psi1+2*q1-q0
    f2=F(x,E)
    psi2=q2/(1-f2*dx**2)
    return q2,f2,psi2


# In[11]:


def initials(E=3/2,Xmin=-5,Xmax=5,psi_0=10**(-30),psi_1=10**(-29),div=10**5):
    '''
    Xmin,Xmax=minimum and maximum of the range
    div denotes the number of divisions for X
    '''
    X=np.linspace(Xmin,Xmax,div)
    dx=X[2]-X[1]
    f_0=F(X[0],E)
    f_1=F(X[1],E)
    q_0=psi_0*(1-dx**2*f_0/12)
    q_1=psi_1*(1-dx**2*f_1/12)
    psi=[psi_0,psi_1]
    f=[f_0,f_1]
    q=[q_0,q_1]
    return X,psi,f,q,dx


# In[12]:


def run_eq(X,q,f,psi,dx,E):
    #print(eps)
    for i in range(len(X)-2):
# q2,f2,psi2=Numrov(X[i+1],q[-2],q[-1],psi[-1],f[-1],dx,eps)
        x=X[i+1]
        f1=f[-1]
        psi1=psi[-1]
        q1=q[-1]
        q0=q[-2]
        q2,f2,psi2=Numrov(x,q0,q1,psi1,f1,dx,E)
        q.append(q2)
        f.append(f2)
        psi.append(psi2)
    psi_n=Normalize(X,psi)
    return X,psi_n


# In[13]:


def run_mult(range_eps,Ra):
    data=[]
    for eps in range_eps:
        X,psi,f,q,dx=initials(eps,Ra[0],Ra[1])
        X,n_psi=run_eq(X,q,f,psi,dx,eps)
        data.append([X,n_psi])
    return data


# In[14]:


def Eg_optim(E_init,num_prec,peak=[-1,1],Ra=[-10,10],div=.1):
    '''
    Energy optimizer function
    this function recieves an enegry assumption and optimizes this to an energy value which 
    will satisfy the boundary condition best
    this code follows similar logic the previously defined eigen finder
    
    E_init:-> initial assumption of E
    num_prec:-> the precision we require for of the eigen energy
    {num_pec denotes the number of decimal places till the eigen value is correct}
    higer the num_prec more accurate the eigen energy is 
    higer the bum_prec more time it takes to compute
    '''
    Dyn_E=E_init#We consider 2 values of E near to this Dynamic E
    dE=div#stores the range in which we check for an optimum(the values are:-> E-dE,E,E+dE)
    Dummy_data=run_mult([Dyn_E],Ra)#X is same for all
    x=Dummy_data[0][0]
    x_0=peak[0]
    x_l=peak[1]
    imin=np.where((x>x_0) & (x<(x_0+.1)))[0][0]# getting the range of indices near the peak region
    imax=np.where((x>x_l)&(x<(x_l+.1)))[0][0]#here our range is the peak range given
    opt_count=0#denotes the number of precision accuired, initialized to be zero
    print(f'Input value for optimization:{E_init}')
    #A loop is run until the needed optimization is achived 
    i=0
    while opt_count<=num_prec:
        '''
        Method: We run the different values of E and accquire an array called explo
        (explotion factor),the explo array keeps the ratio of last value of psi calculated 
        to the maxium found in the expected maximum range. the values in explo will be v
        ery large when the energy is not eigen.when a local minimum is aqcuired in Explo 
        that range of E is magnified to find a more precise value of Eigen Energy.
        '''
        E_range=[Dyn_E-dE,Dyn_E,Dyn_E+dE]#Our epsilon range 3 psi values for which we try the plotting
        q_s=run_mult(E_range,Ra)#getting the plot of quantum states for different energy value
        Explo=[]#array to store the fractional difference
        for j in [0,1,2]:#since we already know there are only 3 plots we calculated
            Y=q_s[j][1]#psi values are diffrent for different epsilons
            expec_max=0# The maximum near the expected reagion
            for y in Y[imin:imax]: #finding the maximum in our expected peak region
                if abs(y)>expec_max:
                    expec_max=abs(y)
            exp_f=abs(Y[-1]/expec_max)# the explotion factor for each energy
            Explo.append(exp_f)
        if Explo[2]<Explo[0]:# comparing the explotion factor between energies
            if Explo[2]<Explo[1]:# if the exp_f is lesser for the higher epsilon
                #print('\t 2<1\t2<0')
                Dyn_E=Dyn_E+dE #then shifts epsilon range to higer energy range
            else:
                #print('\t 1<2<0')
                dE=dE/10.0#if its least for the current energy we need to zoom into the energy range
                opt_count+=1#this zooming in counts as one optimization
        else:
            if Explo[0]<Explo[1]:
                #print('\t 0<2\t0<1')
                Dyn_E=Dyn_E-dE#shiftin our range towards the more accurate eigen energy
            else:
                #print('\t 1<0<2')
                dE=dE/10.0
                opt_count+=1#this zooming in counts as one optimization
        sys.stdout.flush()
        sys.stdout.write("\r{0}".format(f"\tEnergy Value after {1+i} optimizations:->{Dyn_E}"))
        i+=1
    if Dyn_E<0:
        print('Optimization failed')
    else:
        print(f'\nValue  after optimization =>{Dyn_E}')
        return Dyn_E


# In[66]:


def Explo_ary(E_array,peak_ind,Ra=[-10,10]):#function which will return the explotion array
    '''
    This function computes the explotion factor i,e the ration between last point in psi
    to the maximum near the expected region
    '''
    optim_eng=[]#stores the energies which are at closer to eigen enrgies
    imin=peak_ind[0]#the index marking the begining of the region where we ecpect a max
    imax=peak_ind[1]#the index marking the end of the region where we expect a maximum
    q_s=run_mult(E_array,Ra)#getting the plot of quantum states for different energy value
    Explo=[]#array to store the fractional difference
    lg=len(E_array)
    for j in range(lg):#since we already know there are only 3 plots we calculated
            Y=q_s[j][1]#psi values are diffrent for different epsilons
            expec_max=0# The maximum near the expected reagion
            for y in Y[imin:imax]: #finding the maximum in our expected peak region
                if abs(y)>expec_max:
                    expec_max=abs(y)
            exp_f=abs(Y[-1]/expec_max)# the explotion factor for each energy
            Explo.append(exp_f)
    for i in range(lg-2):
        if Explo[i+1]<Explo[i]:
            if Explo[i+1]<Explo[i+2]:
                optim_eng.append(E_array[i+1])
    return optim_eng


# In[98]:


def Energy_finder(E_est,num_prec=5,peak=[-1,1],Ra=[-10,10],div=.1):
    '''
    Returns the maximum possible accurate value of eigen energy based on the 
    prececuion given
    '''
    E_eigen=[]
    E_array=np.linspace(E_est-div,E_est+div,15)# an array with .01 differences
    dE=div#stores the range in which we check for an optimum(the values are:-> E-dE,E,E+dE)
    Dummy_data=run_mult([E_array[0]],Ra)#X is same for all
    x=Dummy_data[0][0]
    x_0=peak[0]
    x_l=peak[1]
    imin=np.where((x>x_0) & (x<(x_0+div)))[0][0]# getting the range of indices near the peak region
    imax=np.where((x>x_l)&(x<(x_l+div)))[0][0]#here our range is the peak range given
    Eg_ar=[]
    E_est=round(E_est,2)
    for i in range(num_prec):
        #print(E_array[0],'\t',E_est,'\t',E_array[-1],'\t',dE)
        Eg_ar=Explo_ary(E_array,[imin,imax],Ra)
        dE=dE/10
        if len(Eg_ar)>1:
            for e in Eg_ar:
                E_eigen.append(Eg_optim(e,(num_prec-i-1),peak,Ra,dE))
                #E=np.copy(E_eigen)
            return E_eigen#E.flatten()
        else:
            E_est=Eg_ar[0]
            E_array=np.linspace(E_est-3/2*dE,E_est+3/2*dE,20)
    return Eg_ar


# In[99]:


def Normalize(x,y,norml_Val=1):#function to normalize the function
    A=0
    norm_y=[]
    for i in range(len(x)-1):
        a=abs((x[i+1]-x[i])*(y[i]+y[i+1])/2)
        A=A+a
    for i in range(len(x)):
        norm_y.append(y[i]/A)
    return norm_y


# In[100]:


def Energy_loc(E_range,peak_r,X_range=[-10,10],div=.1):
    E_r=np.arange(E_range[0],E_range[1],div)
    d=run_mult(E_r,X_range)
    X=d[0][0]
    x_0=peak_r[0]
    x_l=peak_r[1]
    Eigen_E=[]#array to store eigen energies
    n=len(E_r)#number of energy values we are considering
    imin=np.where((X>x_0) & (X<(x_0+.1)))[0][0]# getting the range of indices near the peak region
    imax=np.where((X>x_l)&(X<(x_l+.1)))[0][0]#here our range is the peak range given
    Eng=Explo_ary(E_r,[imin,imax])
    ''' Explo=[]#arry to store the ratio of peak at expected range to the value of psi at x=10
    for i in range(n):
        Y=d[i][1]
        max_psi=0# The maximum ner origin
        for y in Y[imin:imax]: #finding the maximum near origin(our range(-2,2))
            if abs(y)>max_psi:
                max_psi=abs(y)
        exp_f=abs(Y[-1]/max_psi)# the maximum/maximum at peak range
        Explo.append(exp_f)
    
     
        "By looking at the trend of data in grad we get the approximate location
        of the minium which indicates a eigen enrgy state"
    
    for i in range(n-2):
        if Explo[i]>Explo[i+1]:
            if Explo[i+2]>Explo[i+1]:
    
                "The above 2 conditions indicate a minimum for explosion factor
                for the corresponding Energy value. This indicates the presence 
                of an eigen state here. So optimization of energy is done here"
                
                print(E_r[i+1])
                '''
    for E in Eng:         
        Eigen_E.append(Energy_finder(E,5,peak_r,X_range))
    return Eigen_E


# In[ ]:





# In[80]:


#E_1=Energy_matrix[1]
#E_0=Energy_matrix[0]
'''h=6.626070e-34
m=9.10938356e-31
w=5e-30'''
h=1
m=1
w=1
n=0.0
C_4=4
C_2=40
K=1/10
exp_max_r=[-2.5,-1.5]# the range between which we are more likely yo find a probability peak
#E=(n+1/2)*h*w
#E=.00001#Eigen_finder2(0.01,15,[-15,15],.0001)


# In[21]:


#Energy_matrix


# In[15]:


#d=run_mult([E_1,E_0],[-4,4])
#plt.plot(d[0][0],d[0][1],d[1][0],d[1][1])


# In[72]:





# In[102]:


b=time.time()
E=Energy_loc([0,10],exp_max_r)
TT=time.time()-b
print(TT)


# In[103]:


Eigen_Eg=[]
for e in E:
    for erg in e:
        Eigen_Eg.append(erg)


# In[87]:


Eigen_Eg


# In[114]:


EGR=Eigen_Eg[:]
tm=time.time()
D=run_mult(EGR,[-4,4])
tme=time.time()-tm
print(tme)
i=0
for d in D:
    plt.plot(d[0],d[1]+EGR[i],label=EGR[i])
    i+=1
#plt.legend()
plt.show()





