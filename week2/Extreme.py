import numpy as np
from pylab import *

def F(x,y):
    return sin(x+y)+cos(x+2*y)
def dF(x,y):#导函数
    return [cos(x+y)-sin(x+2*y),cos(x+y)-2*sin(x+2*y)]

def Gradient(x0,y0,step,epsilon):
    lg=sqrt((dF(x0,y0)[0])**2+(dF(x0,y0)[1])**2)#计算梯度的模长
    while (lg>epsilon):#若梯度的模长大于允许值则继续循环
        x0=x0-step*dF(x0,y0)[0]/lg
        y0=y0-step*dF(x0,y0)[1]/lg#用梯度更新点坐标
        lg=sqrt((dF(x0,y0)[0])**2+(dF(x0,y0)[1])**2)
    else:
        return [x0,y0]

N = 100
x = np.linspace(-5*pi, 5*pi, N)
y = np.linspace(-5*pi, 5*pi, N)
X, Y = np.meshgrid(x, y)
Z = F(X,Y)
fig, ax = plt.subplots()
cs = ax.contourf(X, Y, Z, cmap=plt.get_cmap('Spectral'))#画出g(x,y)图像
cbar = fig.colorbar(cs)
plt.show()

for i in linspace(-2*pi,2*pi,3):#改变initial guess的值 求解不同的极小值
    for j in linspace(-2*pi,2*pi,3):
        [x,y]=Gradient(i,j,1e-3,1e-3)
        [x,y]=Gradient(x,y,1e-5,1e-5)
        print("For initial guess(",round(i,4),round(j,4),")",", One local minima is: (",round(x,4),round(y,4),")")
