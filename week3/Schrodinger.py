import numpy as np
from pylab import *
import scipy.linalg as SL

m=2
hbar=1
SampleNum=1025

def ConsS2(nu_p,s_p,nu_k,s_k):#由于所选基不正交形成的广义本征值问题中S矩阵的矩阵元计算
    return sqrt(nu_p*nu_k/(pi*(nu_k+nu_p)))*exp(-((s_k-s_p)**2)*nu_k*nu_p/(nu_k+nu_p))

def ConsH2(nu_p,s_p,nu_k,s_k):#势能V(x)=x^2的哈密顿量矩阵元计算
    return (exp(-((s_k-s_p)**2)*nu_k*nu_p/(nu_k+nu_p))*sqrt(nu_k*nu_p)*(m*nu_p*(1+2*s_p**2*nu_p)+nu_k*(m+4*m*s_k*s_p*nu_p+2*hbar*nu_p**2)+2*nu_k**2*(m*s_k**2+hbar*nu_p-2*hbar*(s_k-s_p)**2*nu_p**2)))/(2*m*sqrt(pi*(nu_k+nu_p))*(nu_k+nu_p)**2)

def ConsH4(nu_p,s_p,nu_k,s_k):#势能V(x)=x^4-x^2的哈密顿量矩阵元计算
    P1=4*nu_k**4*(m*s_k**2*(s_k**2-1)+hbar*nu_p-2*hbar*(s_k-s_p)**2*nu_p**2)
    P2=2*nu_k**3*(m*(-1+6*s_k**2)+4*m*s_k*(-s_k+(-1+2*s_k**2)*s_p)*nu_p+6*hbar*nu_p**2-8*hbar*(s_k-s_p)**2*nu_p**3)
    P3=nu_k**2*(3*m+2*nu_p*(-3*m+6*m*s_k*(s_k+2*s_p)+2*m*(-4*s_k*s_p-s_p**2+s_k**2*(-1+6*s_p**2))*nu_p+6*hbar*nu_p**2-4*hbar*(s_k-s_p)**2*nu_p**3))
    P4=m*nu_p**2*(3+2*nu_p*(-1+2*s_p**2*(3+(-1+s_p**2)*nu_p)))
    P5=2*nu_p*nu_k*(3*m+nu_p*(-3*m+2*m*s_p**2*(3-2*nu_p)+2*hbar*nu_p**2+4*m*s_k*s_p*(3+(-1+2*s_p**2)*nu_p)))
    return exp(-(s_k-s_p)**2*nu_k*nu_p/(nu_k+nu_p))*sqrt(nu_k*nu_p)*(P1+P2+P3+P4+P5)/(4*m*sqrt(pi)*np.power(nu_k+nu_p,4.5))

def GaussBasis(nu_,s_,x):#建立一个高斯函数的基
    return sqrt(nu_/pi)*exp(-nu_*(x-s_)**2)

def Basis(nu_,s_,x,N):#将N个基拼成矩阵
    B=zeros((N,SampleNum))
    for i in range(0,N):
        B[i,:]=GaussBasis(nu_[i],s_[i],x)
    return B
    
N=101#取N个基
step=0.15#每个基的均值s间隔step
nu_num=6#各个基有统一的nu

s=linspace(-(N-1)/2,(N-1)/2,N,endpoint=True)*step#基的均值的序列
nu=ones(N)*nu_num#基的nu的序列

xl=-(N-1)/2*step
xh=(N-1)/2*step
x=linspace(xl,xh,SampleNum)

#对于势能V(x)=x^2
S=zeros((N,N))
H=zeros((N,N))
for p in range(0,N):
    for k in range(p,N):#构造S和N矩阵
        S[p,k]=ConsS2(nu[p],s[p],nu[k],s[k])
        H[p,k]=ConsH2(nu[p],s[p],nu[k],s[k])
        S[k,p]=S[p,k]
        H[k,p]=H[p,k]

(E,C)=SL.eigh(H,S)#解广义本征值问题

print(E[0:3])
scatter(range(1,size(E)+1),E)#绘制本征值分布
xlabel("n")
ylabel("E")
show()

psi0=matmul(C[:,0]/sqrt(linalg.norm(C[:,0])),Basis(nu,s,x,N))#计算并绘制前三个本征态的波函数
psi1=matmul(C[:,1]/sqrt(linalg.norm(C[:,1])),Basis(nu,s,x,N))
psi2=matmul(C[:,2]/sqrt(linalg.norm(C[:,2])),Basis(nu,s,x,N))
xlabel("x")
ylabel("y")
plot(x,psi0,x,psi1,x,psi2)
legend(["E=0.5","E=1.0","E=1.5"])
show()

#对于势能V(x)=x^4-x^2
S=zeros((N,N))
H=zeros((N,N))
for p in range(0,N):
    for k in range(p,N):#构造S和N矩阵
        S[p,k]=ConsS2(nu[p],s[p],nu[k],s[k])
        H[p,k]=ConsH4(nu[p],s[p],nu[k],s[k])
        S[k,p]=S[p,k]
        H[k,p]=H[p,k]

(E,C)=SL.eigh(H,S)#解广义本征值问题

print(E[0:3])
scatter(range(1,size(E)+1),E)#绘制本征值分布
xlabel("n")
ylabel("E")
show()

psi0=matmul(C[:,0]/sqrt(linalg.norm(C[:,0])),Basis(nu,s,x,N))#计算并绘制前三个本征态的波函数
psi1=matmul(C[:,1]/sqrt(linalg.norm(C[:,1])),Basis(nu,s,x,N))
psi2=matmul(C[:,2]/sqrt(linalg.norm(C[:,2])),Basis(nu,s,x,N))
xlabel("x")
ylabel("y")
plot(x,psi0,x,psi1,x,psi2)
legend(["E=0.147","E=0.872","E=2.128"])
show()