import numpy as np
from numpy.fft import fft,fftshift
from pylab import *
import scipy.linalg as SL


hbar = 1.0545e-34
m = 9.11e-31
a = 1e-9
b = 0.9e-9
V0 = 2
E0=2*hbar**2*pi**2/(m*a**2)/1.6e-19

N1=50
N = 2*N1+1

X = linspace(0, a, N, endpoint=False)#采样点位置
V = zeros(N)
for i in range(0, N):
    if X[i]<b:#为势能函数赋值
        V[i] = -V0

Vq = fftshift(fft(V))#对势能函数做FFT

H = zeros((N, N), dtype=complex)
for p in range(-N1, N1+1):
    for q in range(-N1, N1+1):
        if (p-q >= -N1) and (p-q <= N1):
            Vtemp = Vq[p-q+N1]/N#哈密顿量矩阵元的势能部分
        else:
            Vtemp = 0
        if p == q:
            T = E0*q**2#哈密顿量矩阵元的动能部分
        else:
            T=0
        H[p+N1, q+N1] = T+Vtemp

E = SL.eigvalsh(H)#求解本征值问题
print(E[0:3])#输出三个最低的本征值

xlabel("x(m)")
ylabel("V(x)(eV)")
plot(X,V)#绘制势能函数
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()

xlabel("q")
ylabel("Vq")
plot(range(-N1,N1+1),abs(Vq))#绘制势能函数的谱
ax = gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
show()