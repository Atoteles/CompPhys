import numpy as np
from pylab import *

def F(x):
    return x**3-5*x+3

def dF(x):#导函数
    return 3*x**2-5

def Bisection(xl,xr,epsilon):#二分法
    x0=(xl+xr)*0.5
    while np.abs(xr-xl)>epsilon:
        if F(x0)*F(xl)>0:
            xl=x0
        elif F(x0)*F(xr)>0:
            xr=x0
        else:
            return x0
        x0=(xl+xr)*0.5
    else:
        return x0

def Newton(xi,epsilon):#牛顿法
    while np.abs(F(xi)/dF(xi))>epsilon:
        xi=xi-F(xi)/dF(xi)
    else:
        return xi

def Hybrid(xl,xr,epsilon):#混合法
    x0=(xl+xr)*0.5
    while np.abs(xr-xl)>epsilon:
        if dF(x0)==0:#若导数为0，则用二分法
            if F(x0)*F(xl)>0:
                xl=x0
            elif F(x0)*F(xr)>0:
                xr=x0
            else:
                return x0
            x0=(xl+xr)*0.5
        else:#导数不为0用牛顿法
            x0=x0-F(x0)/dF(x0)
            if (x0>xl)and(x0<xr):#新解在范围内则更新范围
                # print("in")
                if F(x0)*F(xl)>0:
                    xl=x0
                elif F(x0)*F(xr)>0:
                    xr=x0
                else:
                    return x0
            else:#新解不在范围内用二分法
                x0=(xl+xr)*0.5
                if F(x0)*F(xl)>0:
                    xl=x0
                elif F(x0)*F(xr)>0:
                    xr=x0
                else:
                    return x0
                x0=(xl+xr)*0.5
    else:
        return x0

X=np.linspace(0,5,100)
plt.plot(X,F(X))
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
plt.show()#绘制f(x)图像

#Ques1:第一个正根在(0,1)，第二个正根在(1,2)
x1_bis=Bisection(0,1,1e-5)
x2_bis=Bisection(1,2,1e-5)
print("The result of bisection are:\tx1=",round(x1_bis,4),",x2=",round(x2_bis,4))
#Ques2:
x1_g=Newton(x1_bis,1e-15)
x2_g=Newton(x2_bis,1e-15)
print("The result of NR method are:\tx1=",round(x1_g,14),",x2=",round(x2_g,14))
#Ques3
x1_h=Hybrid(0,1,1e-15)
x2_h=Hybrid(1,2,1e-15)
print("The result of hybrid method are:\tx1=",round(x1_h,14),",x2=",round(x2_h,14))