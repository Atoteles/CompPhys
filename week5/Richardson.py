import numpy as np
from pylab import *

def Richardson(x,h,N):
    D=zeros((N+1,N+1))
    for n in range(0,N+1):
        D[n,0]=0.5*(sin(x+h)-sin(x-h))/h#中心差分，计算表格第一列
        h=h*0.5#对表格的每一行改变h
    for m in range(1,N+1):
        p4=4**m
        for n in range(m,N+1):
            D[n,m]=D[n,m-1]+(D[n,m-1]-D[n-1,m-1])/(p4-1)
            p4=p4*4
    return D

x=pi/3
h=1
N=10
D=Richardson(x,h,N)
d=zeros(N+1)
for i in range(0,N+1):
    d[i]=D[i,i]
    if i!=0:
        if abs(d[i]-d[i-1])<1e-6:#比较表格对角线上的相邻两个数值之差是否小于精度要求
            break
print("N=",i,"\nD=",D[0:i+1,0:i+1])

xlabel("N")
ylabel("D[N,N]")#绘制导数值随N的收敛情况
plot(range(0,i+1),d[0:i+1],'bo-')
show()
