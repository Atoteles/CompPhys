import numpy as np
from pylab import *

def F(m,t):
    return tanh(m/t)-m

def Fm(m,t):#导函数
    return 1/t*(1-tanh(m/t)**2)-1

def Newton(mi,t,epsilon):
    while np.abs(F(mi,t)/Fm(mi,t))>epsilon:#两次间差异大于允许值时继续循环
        mi=mi-F(mi,t)/Fm(mi,t)
    else:
        return mi

k=0
num=5000
m=[0]*num
for i in linspace(0.01,1.5,num):#在t=0.01到t=1.5间取5000个点
    m[k]=Newton(1,i,1e-10)
    k=k+1

plot(linspace(0.01,1.5,num),m)#绘制m-t图
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
xlabel("t")
ylabel("m")
show()