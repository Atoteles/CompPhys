import numpy as np
from pylab import *
import scipy.linalg as SL

data = np.loadtxt('./sunspots.txt')#加载数据
fourier=fftshift(fft(data[:,1]))#对数据做FFT
F2=multiply(abs(fourier),abs(fourier))

L=data[-1,0]
N=size(data,0)
K=linspace(2*pi/L*(-(N-1)/2),2*pi/L*((N-1)/2),N)#建立k空间坐标

print(K[F2.argsort()[-3:]])#返回F2数组由到大排序后的各元素原来的下标，取最后三个，赋予k空间坐标数组K，得到峰值位置

subplot(211)
xlabel("t")
ylabel("N")
plot(data[:,0],data[:,1])#绘制原始数据
subplot(212)
xlabel("k")
ylabel("|C_k|^2")
plot(K,F2)#绘制傅里叶谱
show()