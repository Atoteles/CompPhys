import matplotlib.pyplot as plt
import numpy as np
import math
import time

J=1

def GenerateVec():#随机生成三维单位矢量
    theta=np.random.random(1)*2*np.pi
    u=np.random.random(1)
    phi=np.arccos(1-2*u)

    x=np.sin(phi)*np.cos(theta)
    y=np.sin(phi)*np.sin(theta)
    z=np.cos(phi)
    return np.array([x,y,z])

def Initialize(nx,ny,nz):#初始化晶格
    spin=np.zeros((nx,ny,nz,3))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                spin[i,j,k,:]=np.squeeze(GenerateVec())
    return spin

def RandomChoose(nx,ny,nz):#随机选取格点进行改变
    return np.array([np.random.randint(0,nx),np.random.randint(0,ny),np.random.randint(0,nz)])

def ForwardIndex(Coord,nx,ny,nz):#周期性边条下标越界处理
    if Coord[0]==nx:
        Coord[0]=0
    elif Coord[0]==-1:
        Coord[0]=nx-1

    if Coord[1]==ny:
        Coord[1]=0
    elif Coord[1]==-1:
        Coord[1]=ny-1

    if Coord[2]==nz:
        Coord[2]=0
    elif Coord[2]==-1:
        Coord[2]=nz-1
    return Coord

def Acceptance(Coord,J,spin,beta,Vec_New,nx,ny,nz):
    list=np.zeros((12,3))
    list[0,:]=ForwardIndex(np.array([Coord[0]+1,Coord[1],Coord[2]]),nx,ny,nz)#分类寻找12个配位格点
    list[1,:]=ForwardIndex(np.array([Coord[0]-1,Coord[1],Coord[2]]),nx,ny,nz)
    list[2,:]=ForwardIndex(np.array([Coord[0],Coord[1]+1,Coord[2]]),nx,ny,nz)
    list[3,:]=ForwardIndex(np.array([Coord[0],Coord[1]-1,Coord[2]]),nx,ny,nz)

    list[6,:]=ForwardIndex(np.array([Coord[0],Coord[1],Coord[2]+1]),nx,ny,nz)
    list[7,:]=ForwardIndex(np.array([Coord[0],Coord[1],Coord[2]-1]),nx,ny,nz)
    if np.mod(Coord[1],2)==0:
        list[4,:]=ForwardIndex(np.array([Coord[0]-1,Coord[1]-1,Coord[2]]),nx,ny,nz)
        list[5,:]=ForwardIndex(np.array([Coord[0]-1,Coord[1]+1,Coord[2]]),nx,ny,nz)

        if np.mod(Coord[2],3)==0:
            list[8,:]=ForwardIndex(np.array([Coord[0],Coord[1]-1,Coord[2]-1]),nx,ny,nz)
            list[9,:]=ForwardIndex(np.array([Coord[0]-1,Coord[1]-1,Coord[2]-1]),nx,ny,nz)
            list[10,:]=ForwardIndex(np.array([Coord[0]-1,Coord[1]-1,Coord[2]+1]),nx,ny,nz)
            list[11,:]=ForwardIndex(np.array([Coord[0]-1,Coord[1],Coord[2]+1]),nx,ny,nz)
        elif np.mod(Coord[2],3)==1:
            list[8,:]=ForwardIndex(np.array([Coord[0],Coord[1]-1,Coord[2]-1]),nx,ny,nz)
            list[9,:]=ForwardIndex(np.array([Coord[0]+1,Coord[1],Coord[2]-1]),nx,ny,nz)
            list[10,:]=ForwardIndex(np.array([Coord[0]+1,Coord[1],Coord[2]+1]),nx,ny,nz)
            list[11,:]=ForwardIndex(np.array([Coord[0],Coord[1]+1,Coord[2]+1]),nx,ny,nz)
        elif np.mod(Coord[2],3)==2:
            list[8,:]=ForwardIndex(np.array([Coord[0]-1,Coord[1],Coord[2]-1]),nx,ny,nz)
            list[9,:]=ForwardIndex(np.array([Coord[0]-1,Coord[1]+1,Coord[2]-1]),nx,ny,nz)
            list[10,:]=ForwardIndex(np.array([Coord[0]-1,Coord[1]+1,Coord[2]+1]),nx,ny,nz)
            list[11,:]=ForwardIndex(np.array([Coord[0],Coord[1]+1,Coord[2]+1]),nx,ny,nz)

    elif np.mod(Coord[1],2)==1:
        list[4,:]=ForwardIndex(np.array([Coord[0]+1,Coord[1]-1,Coord[2]]),nx,ny,nz)
        list[5,:]=ForwardIndex(np.array([Coord[0]+1,Coord[1]+1,Coord[2]]),nx,ny,nz)

        if np.mod(Coord[2],3)==0:
            list[8,:]=ForwardIndex(np.array([Coord[0]+1,Coord[1]-1,Coord[2]-1]),nx,ny,nz)
            list[9,:]=ForwardIndex(np.array([Coord[0],Coord[1]-1,Coord[2]-1]),nx,ny,nz)
            list[10,:]=ForwardIndex(np.array([Coord[0],Coord[1]-1,Coord[2]+1]),nx,ny,nz)
            list[11,:]=ForwardIndex(np.array([Coord[0]-1,Coord[1],Coord[2]+1]),nx,ny,nz)
        elif np.mod(Coord[2],3)==1:
            list[8,:]=ForwardIndex(np.array([Coord[0]+1,Coord[1]-1,Coord[2]-1]),nx,ny,nz)
            list[9,:]=ForwardIndex(np.array([Coord[0]+1,Coord[1],Coord[2]-1]),nx,ny,nz)
            list[10,:]=ForwardIndex(np.array([Coord[0]+1,Coord[1],Coord[2]+1]),nx,ny,nz)
            list[11,:]=ForwardIndex(np.array([Coord[0]+1,Coord[1]+1,Coord[2]+1]),nx,ny,nz)
        elif np.mod(Coord[2],3)==2:
            list[8,:]=ForwardIndex(np.array([Coord[0]-1,Coord[1],Coord[2]-1]),nx,ny,nz)
            list[9,:]=ForwardIndex(np.array([Coord[0],Coord[1]+1,Coord[2]-1]),nx,ny,nz)
            list[10,:]=ForwardIndex(np.array([Coord[0],Coord[1]+1,Coord[2]+1]),nx,ny,nz)
            list[11,:]=ForwardIndex(np.array([Coord[0]+1,Coord[1]+1,Coord[2]+1]),nx,ny,nz)
    
    list=list.astype(int)

    dE=0
    for i in range(12):#计算改变一个自旋的能量变化
        dE+=-J*(np.dot(spin[list[i,0],list[i,1],list[i,2],:],np.squeeze(Vec_New))-np.dot(spin[list[i,0],list[i,1],list[i,2],:],spin[Coord[0],Coord[1],Coord[2],:]))
    # print(dE)
    return np.min([np.exp(-dE*beta),1.])

def Moveforward(Coord,spin,beta,nx,ny,nz):
    p=np.random.rand(1)#生成随机数
    Vec_New=GenerateVec()#生成一个新自旋
    acpt=Acceptance(Coord,J,spin,beta,Vec_New,nx,ny,nz)#计算接受概率
    if p<acpt:
        spin[Coord[0],Coord[1],Coord[2],:]=np.squeeze(Vec_New)#若接受则改变自旋
    return spin

def CalMag(spin):
    x=np.sum(np.sum(np.sum(spin[:,:,:,0])))
    y=np.sum(np.sum(np.sum(spin[:,:,:,1])))
    z=np.sum(np.sum(np.sum(spin[:,:,:,2])))
    return np.array([x,y,z])#计算晶格总磁化强度

bin_contain=48
step=6000
abondon=1000

b_num=20
T_l=0.02
T_step=0.5
T_h=T_l+T_step*b_num
T=np.linspace(T_l,T_h,b_num)
b=1/T

nx=4
ny=4
nz=3

M=np.zeros((b_num,3))

for k in range(0,b_num):
    start=time.perf_counter()
    Mlist=np.zeros((step,4))
    spin=Initialize(nx,ny,nz)
    for i in range(-abondon,step):#前abandon个采样舍去
        for j in range(0,bin_contain):#每bin_contain次改变作为一次采样
            Coord=RandomChoose(nx,ny,nz)#随机选取一个格点
            spin=Moveforward(Coord,spin,b[k],nx,ny,nz)
        if i>=0:
            Mlist[i,:]=CalMag(spin)
    M[k,:]=np.average(Mlist,0)
    end=time.perf_counter()
    print('time cost: {}s,beta_num={}'.format((end-start),k))
plt.figure
plt.plot(T,np.abs(M[:,0]),'-o')
plt.plot(T,np.abs(M[:,1]),'-o')
plt.plot(T,np.abs(M[:,2]),'-o')
plt.plot(T,np.sqrt(M[:,0]**2+M[:,1]**2+M[:,2]**2),'-o')
plt.legend(['$M_x$','$M_y$','$M_z$','$|M|$'])
plt.xlabel('$T$')
plt.ylabel('$M$')
plt.show()
np.savetxt(fname=".\M.csv", X=M, fmt="%f",delimiter=",")