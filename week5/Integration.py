import numpy as np
from pylab import *

Z=14
Ze=Z*sqrt(Z)/(9*sqrt(3))
rh=40
rl=0
r0=5e-4
th=log(rh/r0+1)

def R3s(r):
    rho=2*Z*r/3
    return abs((6-6*rho+rho**2)*exp(-0.5*rho)*Ze)**2*r**2#待积函数

def SimpsonEquaSpace(n):
    N=2*n+1#保证取点个数为奇数
    r=linspace(rl,rh,N,endpoint=True)
    I=0
    for i in range(0,n):
        I=I+R3s(r[2*i])+4*R3s(r[2*i+1])+R3s(r[2*i+2])#代入Simpson公式
    return (rh-rl)*I/(3*N)

def SimpsonNonEquaSpace(n):
    N=2*n+1
    t=linspace(rl,th,N)
    r=r0*(exp(t)-1)#换元
    I=0
    for i in range(0,n):
        I=I+R3s(r[2*i])*exp(t[2*i])+4*R3s(r[2*i+1])*exp(t[2*i+1])+R3s(r[2*i+2])*exp(t[2*i+2])#换元积分时需要乘dr/dt
    return (th-rl)*I*r0/(3*N)

num=1000
N=(arange(num)+1)*2+1
Equa=zeros(num)
NonEqua=zeros(num)
for n in range(1,num+1):
    Equa[n-1]=SimpsonEquaSpace(n)
    NonEqua[n-1]=SimpsonNonEquaSpace(n)
print("equal:",Equa[num-1],"Nonequal:",NonEqua[num-1],"when N=",2*num+1)

xlabel("r")
ylabel("R*r^2")
r=linspace(rl,rh,1000)
plot(r,R3s(r))#绘制f(x)图像
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()

xlabel("t")
ylabel("R*r^2*(dr/dt)")#绘制g(t)图像
t=linspace(rl,th,1000)
plot(t,R3s(r0*(exp(t)-1))*r0*exp(t))
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()

xlabel("N")
ylabel("I")
plot(N,Equa,N,NonEqua)#绘制两种取点方式的积分结果随着取点个数增加的收敛情况
legend(["Equal","Nonequal"])
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()
