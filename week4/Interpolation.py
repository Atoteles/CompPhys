import numpy as np
from pylab import *

def f(x):
    return 1/(1+25*x**2)

def Newton(n,Y,step):
    DY=zeros((n,n))
    DY[0,:]=Y
    for i in range(0,n-1):
        DY[i+1,0:n-i-1]=diff(DY[i,0:n-i])/(step*(i+1))
    print("coefficients of Newton methods are:",DY[:,0])
    return DY[:,0]

def NewtonF(B,X0,x):
    n=size(B,0)
    res=0
    for i in range(0,n):
        res=res+B[n-i-1]
        if i!=n-1:
            res=res*(x-X0[n-2-i])
    return res

def GetMat(n,Y,step):
    A=zeros((n,n))
    b=zeros(n)
    A[0,0]=1
    A[n-1,n-1]=1
    b[0]=0
    b[n-1]=0
    for i in range(1,n-1):
        A[i,i-1]=step
        A[i,i]=4*step
        A[i,i+1]=step
        b[i]=6/step*(Y[i+1]+Y[i-1]-2*Y[i])
    return A,b

def Thomas(A,b):
    n=size(b)
    Ep=zeros(n)
    E=zeros(n)
    F=zeros(n)
    Fp=zeros(n)
    G=zeros(n)
    for k in range(0,n):
        F[k]=A[k,k]
        if k!=0:
            E[k]=A[k,k-1]
        if k!=n-1:
            G[k]=A[k,k+1]
    Fp[0]=F[0]
    for k in range(1,n):
        Ep[k]=E[k]/Fp[k-1]
        Fp[k]=F[k]-Ep[k]*G[k-1]

    Y=zeros(n)
    Y[0]=b[0]
    for k in range(1,n):
        Y[k]=b[k]-Ep[k]*Y[k-1]
    X=zeros(n)
    X[n-1]=Y[n-1]/Fp[n-1]
    for k in range(n-2,-1,-1):
        X[k]=(Y[k]-G[k]*X[k+1])/Fp[k]
    print("f''(x_i)=",X)
    return X

def CubicSplines(Xpp,X,Y,x):
    step=X[1]-X[0]
    ipre=floor((x-X[0])/step)
    i=ipre.astype('int16')
    return Xpp[i]/(6*step)*(X[i+1]-x)**3+Xpp[i+1]/(6*step)*(x-X[i])**3+(Y[i]/step-Xpp[i]*step/6)*(X[i+1]-x)+(Y[i+1]/step-Xpp[i+1]*step/6)*(x-X[i])


#------------cos------------
x_draw=linspace(0, pi, num=500,endpoint=False)
X=linspace(0, pi, num=10, endpoint=True)
Y=cos(X)
##-----------Newton---------
y_newton=NewtonF(Newton(10,Y,pi/9),X,x_draw)
##-----------Cubic----------
(A,b)=GetMat(10,Y,pi/9)
y_cubic=CubicSplines(Thomas(A,b),X,Y,x_draw)
scatter(X,Y)
y_theo=cos(x_draw)
xlabel("x")
ylabel("y")
plot(x_draw,y_newton,x_draw,y_cubic,x_draw,y_theo)
legend(["SamplePoints","Newton","CubicSplines","Exact"])
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()

xlabel("x")
ylabel("error")
plot(x_draw,y_newton-y_theo,x_draw,y_cubic-y_theo)
legend(["Newton","CubicSplines"])
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()
#------------F------------
x_draw=linspace(-1, 1, num=500,endpoint=False)
X=linspace(-1, 1, num=10, endpoint=True)
Y=f(X)
##-----------Newton---------
y_newton=NewtonF(Newton(10,Y,2/9),X,x_draw)
##-----------Cubic----------
(A,b)=GetMat(10,Y,2/9)
y_cubic=CubicSplines(Thomas(A,b),X,Y,x_draw)
scatter(X,Y)
y_theo=f(x_draw)
xlabel("x")
ylabel("y")
plot(x_draw,y_newton,x_draw,y_cubic,x_draw,y_theo)
legend(["SamplePoints","Newton","CubicSplines","Exact"])
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()

xlabel("x")
ylabel("error")
plot(x_draw,y_newton-y_theo,x_draw,y_cubic-y_theo)
legend(["Newton","CubicSplines"])
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()