import numpy as np
from pylab import *

def Pivote(A,col_index):#col_index要从0开始，目的是把第col_index列在第col_index行之下的最大值换到第col_index行
    row_num=size(A,0)
    col_num=size(A,1)
    largest_num=abs(A[col_index,col_index])
    largest_index=col_index
    for i in range(col_index,row_num):#寻找第col_index列在第col_index行之下的最大值
        if abs(A[i,col_index])>largest_num:
            largest_num=abs(A[i,col_index])
            largest_index=i
    if largest_num==0:#若最大值是0，则不能pivote，直接返回原矩阵
        return 0,A
    elif largest_index==col_index:#若最大值即为第col_index行的值，则将该数值化为1
        A[col_index,col_index:col_num]=A[col_index,col_index:col_num]/A[col_index,col_index]
        return 1,A
    else:#交换两行并将最大值化为1
        largest_row=array(A[largest_index,col_index:col_num])
        A[largest_index,col_index:col_num]=A[col_index,col_index:col_num]
        A[col_index,col_index:col_num]=largest_row/largest_row[0]
        return 1,A

def ElemTrans(A,col_index):#初等变换
    row_num=size(A,0)
    for i in range(0,row_num):
        if i==col_index:
            continue
        else:
            A[i,:]=A[i,:]-A[col_index,:]*A[i,col_index]#第col_index行第col_index列已经为1
    return A

A=array([[2.,8.,4.,2.],[2.,5.,1.,5.],[4.,10.,-1.,1.]])
print(A,"\n")
row_num=size(A,0)
col_num=size(A,1)
for i in range(0,min(row_num,col_num)):
    (status,A)=Pivote(A,i)
    if status==0:#若该列全为0，则不进行初等变换
        continue
    else:
        A=ElemTrans(A,i)
print(A)