#!/usr/bin/env python
# coding: utf-8

# In[37]:


import matplotlib.pyplot as plt
import numpy as np
from math import exp
import sys
#get_ipython().run_line_magic('matplotlib', 'qt')


# In[2]:


def Poten(x):
    C_0=C_2**2/(4*C_4)#term to make the potential positieve always
    V=K*(C_4*(x**4)-C_2*(x**2)+C_0)
    return V
def F(x,E):
    v=Poten(x)
    return 2*m*(v-E)/h**2


# In[3]:


def Numrov(x,q0,q1,psi1,f1,dx,E):
    q2=(dx**2)*f1*psi1+2*q1-q0
    f2=F(x,E)
    psi2=q2/(1-f2*dx**2)
    return q2,f2,psi2


# In[4]:


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


# In[5]:


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


# In[6]:


def run_mult(range_eps,Ra):
    data=[]
    for eps in range_eps:
        X,psi,f,q,dx=initials(eps,Ra[0],Ra[1])
        X,n_psi=run_eq(X,q,f,psi,dx,eps)
        data.append([X,n_psi])
    return data


# In[40]:


def Eigen_finder2(eps_init,num_steps,Ra=[-5,5],div=.1):
    '''
    starts with an eigen values
    in all epsilon corresponding to non iegen energies the psi goes to infinty after
    the origin.Then find the  fractional difference between the higest point near origin and last value of psi
    Then this function tries to minimize this fractional difference by changing the epsilon.
    i,e it tries to find a psi which has least value near the end of our range(this our boundary condition)
    '''
    initial=eps_init
    Div=div
    print(f'Input value for optimization:{initial}')
    for i in range(num_steps):# an arbitary number of steps
        #print(f'step {i+1}: \n \t Initial guess for this step:{eps_init}')
        E_range=[eps_init-div,eps_init,eps_init+div]#Our epsilon range 3 psi values for which we try the plotting
        d=run_mult(E_range,Ra)#getting the values
        X=d[0][0]#X is same for all
        imin=np.where((X>-.1) & (X<-.09))[0][0]# getting the range of indices near origin
        imax=np.where((X>.1)&(X<.11))[0][0]#here our range is (-.1,.1)
        Grad=[]#array to store the fractional difference
        '''
        Method:
        1) find the fractional differencees for 3 epsilon values
        2)fins the minimum amoung the fractional differences
        3) I .Shifts the epsilon ranges toward the epsilon which gave minimum difference
          II .If the difference is less for the current epsilon the epsilon range
          is futher divided into more finer intervales
        '''
        for j in range(len(E_range)):
            Y=d[j][1]#psi values are diffrent for different epsilons
            max_psi=0# The maximum ner origin
            for y in Y[imin:imax]: #finding the maximum near origin(our range(-2,2))
                if abs(y)>max_psi:
                    max_psi=abs(y)
            g=abs(Y[-1]/max_psi)# the maximum
            Grad.append(g)
        if Grad[2]<Grad[0]:# comparing the fractional difference between different epsilons
            if Grad[2]<Grad[1]:# if the fractional difference is lesser for the higher epsilon
                #print('\t 2<1\t2<0')
                eps_init=eps_init+div#then shifts epsilon range to higer epsilon
            else:
                #print('\t 1<2<0')
                div=div/2# if its lest for current epsilon the range is further divided into finer intervals
        else:
            if Grad[0]<Grad[1]:
                #print('\t 0<2\t0<1')
                eps_init=eps_init-div
            else:
                #print('\t 1<0<2')
                div=div/2
        #print(f'\t value after optimization  ->{eps_init}')
        sys.stdout.flush()
        sys.stdout.write("\r{0}".format(f"\toptimized after step {i+1}->{eps_init}"))
    if eps_init<0:
        print('Optimization failed')
    else:
        E_final=eps_init
        print(f'Value  after optimization =>{E_final}')
        return E_final


# In[8]:


def Normalize(x,y,norml_Val=1):
    A=0
    norm_y=[]
    for i in range(len(x)-1):
        a=(x[i+1]-x[i])*(y[i]+y[i+1])/2
        A=A+a
    for i in range(len(x)):
        norm_y.append(y[i]/A)
    return norm_y


# In[9]:


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
#E=(n+1/2)*h*w
#E=.00001#Eigen_finder2(0.01,15,[-15,15],.0001)


# In[41]:


E=Eigen_finder2(1,30,[-10,10],.1)


# In[29]:


#for E in np.linspace(6,10,10):
X,psi,f,q,dx=initials(E,-4,4,1e-50,1e-50,10**5)


# In[32]:


X,norm_psi=run_eq(X,q,f,psi,dx,E)


# In[33]:


E_psi=np.copy(norm_psi)
for i in range(len(X)):
    E_psi[i]=norm_psi[i]+E
plt.plot(X,E_psi,label=E)
plt.legend()
plt.show()
print(E)


# In[34]:


plt.plot(X,Poten(X))


# In[ ]:




