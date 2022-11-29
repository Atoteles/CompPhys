import numpy as np
from pylab import *
import math

hbar=1.0546e-34
m=9.11e-31
charge=1.6e-19
a=5.29e-11
Ex2=hbar**2/(2*m)/(charge*a**2)#h^2/2m ev*1^2
x2E=1/Ex2
V0=charge/(4*pi*8.85e-12)/a

def k2(A,B,r,Num):#Numerov算法中的k^2
    if Num==1:
        return -A-x2E*(-1/r*V0)-B/r**2
    elif Num==2:
        Z=3
        rl=0.4
        c1=-14.0093922
        c2=9.5099073
        c3=-1.7532723
        c4=0.0834586
        V=-Z/r*math.erf(r/(sqrt(2)*rl))+exp(-0.5*(r/rl)**2)*(c1+c2*(r/rl)**2+c3*(r/rl)**4+c4*(r/rl)**6)
        return -A-x2E*(V*V0)-B/r**2

def Numerov(u0,rmax,N,A,B,Num):
    u=np.zeros((2,N))
    r=linspace(1e-6,rmax,N)
    delta_r=r[1]-r[0]
    u[0,:]=r
    u[1,N-2]=u0[0]#远端两点的边界条件
    u[1,N-1]=u0[1]
    for i in range(N-2,0,-1):
        u[1,i-1]=(2*(1-5/12*delta_r**2*k2(A,B,r[i],Num))*u[1,i]-(1+1/12*delta_r**2*k2(A,B,r[i+1],Num))*u[1,i+1])/(1+1/12*delta_r**2*k2(A,B,r[i-1],Num))
    return u

def Bisection(Amin,Amax,B,epsilon,u0,rmax,N,Num):#输入时需要注意u_0(Amin)<0,u_0(Amax)>0
    error=999
    while error>epsilon:#二分法求解u_0(A)=0的根
        Amid=0.5*(Amax+Amin)
        result=Numerov(u0,rmax,N,Amid,B,Num)
        if result[1,0]>0:
            Amax=Amid
        elif result[1,0]<0:
            Amin=Amid
        else:
            return Amid*Ex2
        error=np.abs(Amax-Amin)
    return Amid*Ex2


rmax=20
N=10000
epsilon=1e-6

Num=1#第一种势场
#1s
figure
subplot(131)
u0=array([1e-9,0])
AR=arange(12*x2E,14*x2E,0.01)
B=0
Res=zeros(size(AR,0))
for i in range(0,size(AR,0)):
    result=Numerov(u0,rmax,N,AR[i],B,Num)
    plot(result[0,:],result[1,:],'bo-',markersize='2')
    Res[i]=result[1,0]
result=Numerov(u0,rmax,N,13.56324*x2E,B,Num)
plot(result[0,:],result[1,:],'o-',markersize='2',color='red')#绘制打靶法示意图
xlabel('x')
ylabel('$\\psi_{1s}$')
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
print('V_1 E_1s:',Bisection(12*x2E,14*x2E,0,epsilon,u0,rmax,N,Num))#二分法求根
#2s
subplot(132)
u0=array([-1e-5,0])
AR=arange(2*x2E,4*x2E,0.01)
B=0
Res=zeros(size(AR,0))
for i in range(0,size(AR,0)):
    result=Numerov(u0,rmax,N,AR[i],B,Num)
    plot(result[0,:],result[1,:],'bo-',markersize='2')
    Res[i]=result[1,0]
result=Numerov(u0,rmax,N,3.3906*x2E,B,Num)
plot(result[0,:],result[1,:],'o-',markersize='2',color='red')#绘制打靶法示意图
xlabel('x')
ylabel('$\\psi_{2s}$')
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
print('V_1 E_2s:',Bisection(2*x2E,4*x2E,B,epsilon,u0,rmax,N,Num))#二分法求根
#2p
subplot(133)
u0=array([1e-8,0])
AR=arange(2*x2E,4*x2E,0.01)
Res=zeros(size(AR,0))
B=2
for i in range(0,size(AR,0)):
    result=Numerov(u0,rmax,N,AR[i],B,Num)
    plot(result[0,:],result[1,:],'bo-',markersize='2')
    Res[i]=result[1,0]
result=Numerov(u0,rmax,N,3.3907*x2E,B,Num)
plot(result[0,:],result[1,:],'o-',markersize='2',color='red')#绘制打靶法示意图
xlabel('x')
ylabel('$\\psi_{2p}$')
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()
print('V_1 E_2p:',Bisection(4*x2E,2*x2E,B,epsilon,u0,rmax,N,Num))#二分法求根
# n=1
u0=array([1e-5,0])
AR=arange(1*x2E,14*x2E,0.01)
B=0
Res=zeros(size(AR,0))
for i in range(0,size(AR,0)):
    result=Numerov(u0,rmax,N,AR[i],B,Num)
    Res[i]=result[1,0]
figure
plot(AR*Ex2,Res)
scatter([1.3544,3.3907,13.56324],[0,0,0],color='red')#绘制E-u_0关系图观察零点分布
xlabel('E')
ylabel('u(r=0)')
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()

Num=2
# 1s
figure
subplot(131)
u0=array([1e-15,0])
AR=arange(120*x2E,124*x2E,0.01)
B=0
Res=zeros(size(AR,0))
for i in range(0,size(AR,0)):
    result=Numerov(u0,rmax,N,AR[i],B,Num)
    plot(result[0,:],result[1,:],'bo-',markersize='2')
    Res[i]=result[1,0]
result=Numerov(u0,rmax,N,120.973*x2E,B,Num)
plot(result[0,:],result[1,:],'o-',markersize='2',color='red')#绘制打靶法示意图
xlabel('x')
ylabel('$\\psi_{1s}$')
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
print('V_2 E_1s:',Bisection(120*x2E,124*x2E,0,epsilon,u0,rmax,N,Num))#二分法求根
# 2s
subplot(132)
u0=array([-1e-5,0])
AR=arange(28*x2E,31*x2E,0.01)
B=0
Res=zeros(size(AR,0))
for i in range(0,size(AR,0)):
    result=Numerov(u0,rmax,N,AR[i],B,Num)
    plot(result[0,:],result[1,:],'bo-',markersize='2')
    Res[i]=result[1,0]
result=Numerov(u0,rmax,N,30.27*x2E,B,Num)
plot(result[0,:],result[1,:],'o-',markersize='2',color='red')#绘制打靶法示意图
xlabel('x')
ylabel('$\\psi_{2s}$')
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
print('V_2 E_2s:',Bisection(28*x2E,31*x2E,B,epsilon,u0,rmax,N,Num))#二分法求根
# 2p
subplot(133)
u0=array([1e-8,0])
AR=arange(28*x2E,31*x2E,0.01)
Res=zeros(size(AR,0))
B=2
for i in range(0,size(AR,0)):
    result=Numerov(u0,rmax,N,AR[i],B,Num)
    plot(result[0,:],result[1,:],'bo-',markersize='2')
    Res[i]=result[1,0]
result=Numerov(u0,rmax,N,30.44*x2E,B,Num)
plot(result[0,:],result[1,:],'o-',markersize='2',color='red')#绘制打靶法示意图
xlabel('x')
ylabel('$\\psi_{2p}$')
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()
print('V_2 E_2p:',Bisection(31*x2E,28*x2E,B,epsilon,u0,rmax,N,Num))#二分法求根
n=1
u0=array([1e-5,0])
AR=arange(10*x2E,125*x2E,0.1)
B=0
Res=zeros(size(AR,0))
for i in range(0,size(AR,0)):
    result=Numerov(u0,rmax,N,AR[i],B,Num)
    Res[i]=result[1,0]
figure
plot(AR*Ex2,Res)
scatter([Bisection(10*x2E,15*x2E,B,epsilon,u0,rmax,N,Num),30.27,120.973],[0,0,0],color='red')#绘制E-u_0关系图观察零点分布
xlabel('E')
ylabel('u(r=0)')
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()
