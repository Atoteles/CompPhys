import matplotlib.pyplot as plt
import numpy as np
import math

def CalDist(x):#计算投点与原点的距离来判断是否在超球体内部
    if np.linalg.norm(x)>=1:
        return False
    else:
        return True

def CalVolume(dime,Point_Num):
    is_interior=0#落在超球体内部的点数
    for i in range(0,Point_Num):
        x=np.random.random(dime)
        if CalDist(x):
            is_interior+=1
    volume=is_interior/Point_Num*np.power(2,dime)#计算体积
    return volume

def TheoVal(dime):#d维超球体理论体积
    return np.power(np.pi,dime/2)/math.gamma(dime/2+1)

Point_Num_Range=np.arange(1000,100000,5000)
Simu_Num=np.size(Point_Num_Range,0)

MCResults=np.zeros((4,Simu_Num))
for i in [2,3,4,5]:#对d=2,3,4,5计算体积
    for j in range(Simu_Num):
        print('(i,j)=(',i,',',j,')')
        MCResults[i-2,j]=CalVolume(i,Point_Num_Range[j])

Theo2=TheoVal(2)
Theo3=TheoVal(3)
Theo4=TheoVal(4)
Theo5=TheoVal(5)

print('2D:Theoretical:',Theo2,',Numerical:',MCResults[0,-1],', Relative Error:',np.abs(Theo2-MCResults[0,-1])/Theo2)
print('3D:Theoretical:',Theo3,',Numerical:',MCResults[1,-1],', Relative Error:',np.abs(Theo3-MCResults[1,-1])/Theo3)
print('4D:Theoretical:',Theo4,',Numerical:',MCResults[2,-1],', Relative Error:',np.abs(Theo4-MCResults[2,-1])/Theo4)
print('5D:Theoretical:',Theo5,',Numerical:',MCResults[3,-1],', Relative Error:',np.abs(Theo5-MCResults[3,-1])/Theo5)

plt.figure
plt.subplot(221)
plt.plot(Point_Num_Range,MCResults[0,:],'-o')
plt.axhline(Theo2)
plt.xlabel('point number')
plt.ylabel('MC result of 2D sphere')


plt.subplot(222)
plt.plot(Point_Num_Range,MCResults[1,:],'-o')
plt.axhline(Theo3)
plt.xlabel('point number')
plt.ylabel('MC result of 3D sphere')

plt.subplot(223)
plt.plot(Point_Num_Range,MCResults[2,:],'-o')
plt.axhline(Theo4)
plt.xlabel('point number')
plt.ylabel('MC result of 4D sphere')

plt.subplot(224)
plt.plot(Point_Num_Range,MCResults[3,:],'-o')
plt.axhline(Theo5)
plt.xlabel('point number')
plt.ylabel('MC result of 5D sphere')
plt.show()