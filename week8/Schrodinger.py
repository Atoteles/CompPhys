import numpy as np
from pylab import *

m=1
hbar=1
def Potential(delta_x,a):#创建势函数数组
    x = np.arange(0,a,delta_x)
    V=np.zeros(size(x,0))
    for i in range(0,size(x,0)):
        if x[i]>a*0.375 and x[i]<a*0.625:
            V[i]=-10
    return V

def CalBeta(delta_t,delta_x):#计算表达式里的一个常数
    return -1j*delta_t*hbar/(4*delta_x**2*m)

def InitialCondition(x,a,k):#创建初始条件
    return np.sqrt(1/np.pi)*exp(-np.power(x-40,2)/2+1j*k*x)

def CalCrankMatrix(delta_t,delta_x,Nx):#计算CrankNicolson方法中的矩阵
    V=Potential(delta_x,delta_x*Nx)
    beta=CalBeta(delta_t,delta_x)
    A=zeros((Nx,Nx),dtype=complex)
    A[0,0]=1-2*beta-1j*delta_t*V[0]/(2*hbar)
    A[0,1]=beta
    A[Nx-1,Nx-1]=1-2*beta-1j*delta_t*V[Nx-1]/(2*hbar)
    A[Nx-1,Nx-2]=beta
    for i in range(1,Nx-1):
        A[i,i]=1-2*beta-1j*delta_t*V[i]/(2*hbar)
        A[i,i-1]=beta
        A[i,i+1]=beta
    return A

def CrankNicolson(psi_0,T,delta_t,delta_x):
    Nt=int(floor(T/delta_t))
    Nx=size(psi_0,0)

    psi=zeros((Nx,Nt),dtype=complex)
    psi[:,0]=psi_0
    psi[0,0]=0
    psi[Nx-1,0]=0

    A1=CalCrankMatrix(delta_t,delta_x,Nx)
    A2=CalCrankMatrix(-delta_t,delta_x,Nx)
    A=matmul(inv(A1),A2)
    for n in range(1,Nt):
        psi[:,n]=matmul(A,psi[:,n-1])#时间演化
        psi[0,n]=0#设置边条
        psi[Nx-1,n]=0
        pause(0.1)
    return psi

def CalExplicitMatrix(delta_t,delta_x,Nx):#计算StableExplicitScheme中的矩阵
    V=Potential(delta_x,delta_x*Nx)
    beta=CalBeta(delta_t,delta_x)
    A=zeros((Nx,Nx),dtype=complex)
    A[0,0]=8*beta-2j*delta_t*V[0]/hbar
    A[0,1]=-4*beta
    A[Nx-1,Nx-1]=8*beta-2j*delta_t*V[Nx-1]/hbar
    A[Nx-1,Nx-2]=-4*beta
    for i in range(1,Nx-1):
        A[i,i]=8*beta-2j*delta_t*V[i]/hbar
        A[i,i-1]=-4*beta
        A[i,i+1]=-4*beta
    return A

def StableExplicit(psi_0,T,delta_t,delta_x):
    Nt=int(floor(T/delta_t))
    Nx=size(psi_0,0)

    psi=zeros((Nx,Nt),dtype=complex)
    psi[:,0]=psi_0
    psi[0,0]=0
    psi[Nx-1,0]=0

    A0=CalExplicitMatrix(delta_t,delta_x,Nx)
    A1=CalCrankMatrix(delta_t,delta_x,Nx)
    A2=CalCrankMatrix(-delta_t,delta_x,Nx)
    A=matmul(inv(A1),A2)
    psi[:,1]=matmul(A,psi[:,0])#先用CrankNicolson方法计算第一步
    psi[0,1]=0
    psi[Nx-1,1]=0

    for n in range(2,Nt):
        psi[:,n]=psi[:,n-2]+matmul(A0,psi[:,n-1])#时间演化
        psi[0,n]=0#设置边条
        psi[Nx-1,n]=0
    return psi

a=200
T=20
k=5
alpha=0.49#稳定条件
delta_x=0.3
delta_t=alpha*delta_x**2
print('delta_t=',delta_t)

t = np.arange(0,T,delta_t)
x = np.arange(0,a,delta_x)
psi_0=InitialCondition(x,a,k)
psi=CrankNicolson(psi_0,T,delta_t,delta_x)#Crank方法
p=abs(multiply(psi,conjugate(psi)))
psi1=StableExplicit(psi_0,T,delta_t,delta_x)#explicit方法
p1=abs(multiply(psi1,conjugate(psi1)))

figure
subplot(211)
cset = plt.imshow(p,aspect='auto') 
plt.colorbar(cset)
subplot(212)
cset = plt.imshow(p1,aspect='auto') 
plt.colorbar(cset)
plt.show()