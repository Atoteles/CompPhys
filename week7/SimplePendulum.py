import numpy as np
from pylab import *

g=9.8
l=1

def y_dot(y):#返回导数向量
    return np.array([y[1],-g/l*sin(y[0])])

def Energy(y):#计算总能量
    return 0.5*np.power(y[1]*l,2)+g*l*(1-np.cos(y[0]))

def Euler(delta_t,y0,total_time):
    N=int(np.floor(total_time/delta_t))
    y=np.zeros((2,N))
    t=np.zeros(N)
    t[0]=0
    y[:,0]=y0

    for i in range(0,N-1):
        y[:,i+1]=y[:,i]+y_dot(y[:,i])*delta_t#forward
        t[i+1]=(i+1)*delta_t
    
    result=np.zeros((3,N))
    result[0,:]=t
    result[1,:]=y[0,:]
    result[2,:]=Energy(y)
    return result

def MidPoint(delta_t,y0,total_time):
    N=int(np.floor(total_time/delta_t))
    y=np.zeros((2,N))
    t=np.zeros(N)
    t[0]=0
    y[:,0]=y0

    for i in range(0,N-1):
        delta_y=y_dot(y[:,i])*delta_t#计算delta_y
        y[:,i+1]=y[:,i]+y_dot(y[:,i]+0.5*delta_y)*delta_t
        t[i+1]=(i+1)*delta_t

    result=np.zeros((3,N))
    result[0,:]=t
    result[1,:]=y[0,:]
    result[2,:]=Energy(y)
    return result

def RK4(delta_t,y0,total_time):
    N=int(np.floor(total_time/delta_t))
    y=np.zeros((2,N))
    t=np.zeros(N)
    t[0]=0
    y[:,0]=y0

    for i in range(0,N-1):
        k1=y_dot(y[:,i])
        k2=y_dot(y[:,i]+0.5*k1*delta_t)
        k3=y_dot(y[:,i]+0.5*k2*delta_t)
        k4=y_dot(y[:,i]+k3*delta_t)
        y[:,i+1]=y[:,i]+1/6*(k1+2*k2+2*k3+k4)*delta_t
        t[i+1]=(i+1)*delta_t

    result=np.zeros((3,N))
    result[0,:]=t
    result[1,:]=y[0,:]
    result[2,:]=Energy(y)
    return result

def EulerTrapezoidal(delta_t,y0,total_time,epsilon):
    N=int(np.floor(total_time/delta_t))
    y=np.zeros((2,N))
    t=np.zeros(N)
    t[0]=0
    y[:,0]=y0

    for i in range(0,N-1):
        dy0=y_dot(y[:,i])
        y_former=y[:,i]+delta_t*dy0#predictor
        e=999
        while e>epsilon:#corrector
            y_later=y[:,i]+0.5*delta_t*(y_dot(y_former)+dy0)
            e=np.max(abs(y_former-y_later))
            y_former=y_later
        y[:,i+1]=y_later
        t[i+1]=(i+1)*delta_t

    result=np.zeros((3,N))
    result[0,:]=t
    result[1,:]=y[0,:]
    result[2,:]=Energy(y)
    return result

delta_t=0.001
total_time=10.001
epsilon=1e-6

#大角度
y0=np.array([np.pi/4,0])
result1=Euler(delta_t,y0,total_time)
result2=MidPoint(delta_t,y0,total_time)
result3=RK4(delta_t,y0,total_time)
result4=EulerTrapezoidal(delta_t,y0,total_time,epsilon)

figure
subplot(121)
plot(result1[0,:],result1[1,:])
plot(result2[0,:],result2[1,:])
plot(result3[0,:],result3[1,:])
plot(result3[0,:],result4[1,:])
legend(['Euler','MidPoint','RK4','EulerTrapezoidal'])
xlabel('t')
ylabel('$\\theta$')

subplot(122)
plot(result1[0,:],result1[2,:])
plot(result2[0,:],result2[2,:])
plot(result3[0,:],result3[2,:])
plot(result3[0,:],result4[2,:])
legend(['Euler','MidPoint','RK4','EulerTrapezoidal'])
xlabel('t')
ylabel('E')
show()

#小角度
y0=np.array([np.pi/100,0])
result1=Euler(delta_t,y0,total_time)
result2=MidPoint(delta_t,y0,total_time)
result3=RK4(delta_t,y0,total_time)
result4=EulerTrapezoidal(delta_t,y0,total_time,epsilon)

L=total_time
N=size(result1[1,:],0)
print(N)
K=linspace(2*pi/L*(-(N-1)/2),2*pi/L*((N-1)/2),N)#建立k空间坐标

fourier1=1/N*fftshift(fft(result1[1,:]))#对数据做FFT
F21=multiply(abs(fourier1),abs(fourier1))
fourier2=1/N*fftshift(fft(result2[1,:]))
F22=multiply(abs(fourier2),abs(fourier2))
fourier3=1/N*fftshift(fft(result3[1,:]))
F23=multiply(abs(fourier3),abs(fourier3))
fourier4=1/N*fftshift(fft(result4[1,:]))
F24=multiply(abs(fourier4),abs(fourier4))

print(K[F21.argsort()[-2:]])

show()
figure
subplot(121)
plot(result1[0,:],result1[1,:])
plot(result2[0,:],result2[1,:])
plot(result3[0,:],result3[1,:])
plot(result3[0,:],result4[1,:])
legend(['Euler','MidPoint','RK4','EulerTrapezoidal'])
xlabel('t')
ylabel('$\\theta$')

subplot(122)
xlabel("k")
ylabel("$|C_k|^2$")
plot(K[5001-20:5001+21],F21[5001-20:5001+21])#绘制傅里叶谱
plot(K[5001-20:5001+21],F22[5001-20:5001+21])
plot(K[5001-20:5001+21],F23[5001-20:5001+21])
plot(K[5001-20:5001+21],F24[5001-20:5001+21])
legend(['Euler','MidPoint','RK4','EulerTrapezoidal'])
show()