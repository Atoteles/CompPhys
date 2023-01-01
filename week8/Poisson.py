import numpy as np
from pylab import *

def GenerateField(Lx,Ly,h,BC):#设置边条
    u=zeros((int(floor(Lx/h)),int(floor(Ly/h))))
    Nx=size(u,0)
    Ny=size(u,1)
    #左右下上
    u[0,:]=BC[0]
    u[Nx-1,:]=BC[1]
    u[:,0]=BC[2]
    u[:,Ny-1]=BC[3]
    return u

def GaussSeidel(rho,h,u,Max_iter_Num,epsilon,omega):
    Nx=size(u,0)
    Ny=size(u,1)    
    iter=0
    u[1:Nx-1,1:Nx-1]=random((Nx-2,Nx-2))
    rho_h=h**2*rho
    norm_f=999
    norm_l=0
    while (iter<Max_iter_Num) and (abs(norm_l-norm_f)>epsilon):#超过最大迭代次数或收敛
        norm_f=norm_l
        for i in range(1,Nx-1):
            for j in range(1,Ny-1):
                u[i,j]=u[i,j]+omega/4*(-rho_h+u[i-1,j]+u[i+1,j]-4*u[i,j]+u[i,j-1]+u[i,j+1])
        norm_l=sum(sum(multiply(u[1:Nx-1,1:Nx-1],u[1:Nx-1,1:Nx-1])))#计算矩阵的模
        iter+=1
        print(iter,abs(norm_l-norm_f))
    return u

Max_iter_Num=1e3
epsilon=1e-6
h=0.005
#Q1
BC1=array([0,0,0,1])
rho1=0
Lx1=1
Ly1=1.5

u1=GenerateField(Lx1,Ly1,h,BC1)
omega1=2/(1+pi/((size(u1,0)+size(u1,1))/2))
u1=GaussSeidel(rho1,h,u1,Max_iter_Num,epsilon,omega1)

#Q2
BC2=array([0,0,0,0])
rho2=1
Lx2=1
Ly2=1
u2=GenerateField(Lx2,Ly2,h,BC2)
omega2=2/(1+pi/((size(u2,0)+size(u2,1))/2))
u2=GaussSeidel(rho2,h,u2,Max_iter_Num,epsilon,omega2)

figure
subplot(121)
x = np.arange(0,Lx1,h)
y = np.arange(0,Ly1,h)
cset = plt.imshow(u1)
yticks([0,int(size(u1,0)/2),size(u1,0)],[x[0],x[int(size(u1,0)/2)-1],x[size(u1,0)-1]])
xticks([0,int(size(u1,1)/2),size(u1,1)],[y[0],y[int(size(u1,1)/2)-1],y[size(u1,1)-1]])
plt.colorbar(cset)

subplot(122)
x = np.arange(0,Lx2,h)
y = np.arange(0,Ly2,h)
cset = plt.imshow(u2) 
yticks([0,int(size(u2,0)/2),size(u2,0)],[x[0],x[int(size(u2,0)/2)-1],x[size(u2,0)-1]])
xticks([0,int(size(u2,1)/2),size(u2,1)],[y[0],y[int(size(u2,1)/2)-1],y[size(u2,1)-1]])
plt.colorbar(cset)
plt.show()
imshow()