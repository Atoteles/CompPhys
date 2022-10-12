import numpy as np
from pylab import *

def LinearCo(X,Y):
    n=size(X)
    c0=sum(X)
    c1=sum(multiply(X,X))
    c2=sum(Y)
    c3=sum(multiply(X,Y))
    return (c0*c3-c1*c2)/(c0**2-(n+1)*c1),(c0*c2-(n+1)*c3)/(c0**2-(n+1)*c1)

def LinearFunc(A,x):
    print("Linear:y=",A[0],"+",A[1],"x")
    return A[0]+A[1]*x

def ParabolicCo(X,Y):
    n=size(X)
    X2=multiply(X,X)
    A=zeros((3,3))
    A[0,0]=n
    A[1,1]=sum(multiply(X,X))
    A[2,2]=sum(multiply(X2,X2))
    A[0,1]=sum(X)
    A[1,0]=A[0,1]
    A[0,2]=sum(X2)
    A[2,0]=A[0,2]
    A[1,2]=sum(multiply(X,X2))
    A[2,1]=A[1,2]
    b=zeros(3)
    b[0]=sum(Y)
    b[1]=sum(multiply(X,Y))
    b[2]=sum(multiply(X2,Y))

    (Q,R)=np.linalg.qr(A)

    return np.matmul(np.matmul(linalg.inv(R),transpose(Q)),b)

def ParabolicFunc(B,x):
    print("Parabolic:y=",B[0],"+",B[1],"x+",B[2],"x^2")
    return B[0]+B[1]*x+B[2]*multiply(x,x)

x_draw=linspace(1,9, num=500,endpoint=False)
X=linspace(1,9,9,endpoint=True)
T=[14.6,18.5,36.6,30.8,59.2,60.1,62.2,79.4,99.9]
y_linear=LinearFunc(LinearCo(X,T),x_draw)
y_para=ParabolicFunc(ParabolicCo(X,T),x_draw)
scatter(X,T)
xlabel("x")
ylabel("T")
plot(x_draw,y_linear,x_draw,y_para)
legend(["SamplePoints","Linear","Parabolic"])
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()

y_linear=LinearFunc(LinearCo(X,T),X)
y_para=ParabolicFunc(ParabolicCo(X,T),X)
xlabel("x")
ylabel("error")
scatter(X,y_linear-T)
scatter(X,y_para-T)
legend(["Linear","Parabolic"])
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()
