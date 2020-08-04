import matplotlib.pyplot as plt
import numpy as np
from modsim import *
def nxt_point_RK4(i,param):#using RK4
    p=param.phi
    dp=param.dphi
    ddp=param.ddphi
    U=param.U
    u=U[i]#our U cordinate
    du=U[2]-U[1]#a small unit in which we evaluate  the equations
    np=p[i]+dp[i]*du#
    ndp=dp[i]+du*ddp[i]
    eps=param.eps
    '''
    k1 is rate of slope change at the begining of interval
    dk1 is change in slope with k1 as rate od change of slope
    k2 is rate the rate of change of slope if dk1 is the change in slope
    in one interval du
    '''
    k1=(u**2-eps)*(p[i])#k1,k2,k3,k4 at specific positions
    dk1=k1*du#change in dpsi due to k1+dpsi
    k2=((u+du/2)**2)*(p[i]+dk1*du/2)
    '''The right-end bracket denotes the new psi value at
        mid point of interval'''
    dk2=k2*du
    k3=((u+du/2)**2)*(p[i]+dk2*du/2)
    dk3=k3*du
    k4=((u+du)**2-eps)*(p[i]+dk3*du)
    nddp=(k1+2*k2+2*k3+k4)/6
    p.append(np)#*(2.718)**(-u**2))
    dp.append(ndp)
    ddp.append(nddp)
def nxt_point(i,param):
    p=param.phi
    dp=param.dphi
    ddp=param.ddphi
    U=param.U
    u=U[i]
    du=U[2]-U[1]
    np=p[i]+dp[i]*du
    ndp=dp[i]+du*ddp[i]
    eps=param.eps
    nddp=(u**2-eps)*p[i]
    p.append(np)#*(2.718)**(-u**2))
    dp.append(ndp)
    ddp.append(nddp)
def make_initials(phi_0,dphi_0):
    phi=[]
    dphi=[]
    ddphi=[]
    ddphi_0=(X[0]**2-eps)*phi_0
    phi.append(phi_0)
    dphi.append(dphi_0)
    ddphi.append(ddphi_0)
    return phi,dphi,ddphi
def run_eq(parameter):
    phi_0,dphi_0=parameter.phi_0,parameter.dphi_0
    U=parameter.U
    phi,dphi,ddphi=make_initials(phi_0,dphi_0)
    param_n=Params(phi=phi,dphi=dphi,ddphi=ddphi,
                   U=U,eps=parameter.eps)
    d=SweepSeries(index=U)
    max=0
    for i in range(len(U)-1):
        nxt_point(i,param_n)
        #print(U[i],phi[i],dphi[i],ddphi[i])
        d[U[i]]=phi[i]
        if(phi[i]>max):
            max=phi[i]
    return d,phi/max,dphi,ddphi
def Run_multiple(params,eps_range):
    for i in eps_range:
        params.eps=i
        phi,dphi,ddphi=run_eq(params)
        plot(phi)
h=6.63*10**(-34)
w=10**(-3)
m=10**(-31)
alpha=(h/(m*w))**(1/2)
X=np.linspace(-20,20,10**4)#here X is not the absolute cordinate it is x/alpha
U=X/alpha
du=U[2]-U[1]
n=0
E=h*w*(1/2+n)
eps=E/(h*w/2)
#eps=2*n+1#eps = E/(h*w/2)
phi_0=1*10**(-30)
dphi_0=1*10**(-3)#first defferential of phi
parameters=Params(eps=eps,phi_0=phi_0,dphi_0=dphi_0,
                  U=X,du=du)
#rang=np.arange(1,6,1)
#Run_multiple(parameters,rang)
'''for i in range(len(X)-1):
    nxt_point(i,X,phi,dphi,ddphi)'''
parameters.eps=1
d,phi,dphi,ddphi=run_eq(parameters)
plot(X,phi)
plt.show()
#Init=State(phi=phi_0,dphi=dphi_0,ddphi=ddphi_0,x=X[0])
#Data=SweepFrame(columns=Init.index)
#Data.row[X[0]]=Init
#for x in X:
#    try:
#        Data.row[x+dx]=nxt_state(Data.row[x])
#    except:
#        plot(Data.phi)
