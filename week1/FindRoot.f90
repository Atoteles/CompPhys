Program FindRoot
!find root for ax^2+bx+c=0
implicit none

Character Conti
REAL a,b,c,Judge,x3,x4
Complex x1,x2,temp

Do
    Print*,"Input a,b and c of ax^2+bx+c=0"
    Read*,a,b,c
    Judge=b**2-4*a*c    !定义判别式，abc为用户键盘输入
    x3=(-b+Sqrt(Judge))/(2*a)   !代入求根公式
    x4=(-b-Sqrt(Judge))/(2*a)

    If(Judge>0) Then    !根据判别式判断根的情况
        Write(*,'("The equation has two different real roots:x1=",F0.5,",x2=",F0.5)')x3,x4
    Else If(Judge==0) Then
        Write(*,'("The equation has two identical real roots:x=",F0.5)')x3
    Else
        temp=Complex(Judge,0.0)   !把判别式数据类型转换为复数，便于下面开方运算
        x1=(-b+Sqrt(temp))/(2*a)
        x2=(-b-Sqrt(temp))/(2*a)
        Write(*,'("The equation has two different complex&
        roots:x1=",F0.5,SP,F0.5,"i,x2=",F0.5,SP,F0.5,"i")')real(x1),aimag(x1),real(x2),aimag(x2)
    End If

    1 Print*,"Do you want to continue?(y/n)"
    Read*,Conti   !用户输入conti变量决定是否继续运行
    If(Conti=='n') Then
        Exit
    Else If(Conti=='y') Then
        Cycle
    Else 
        Print*,"Invalid input!"
        Goto 1
    End if
End Do
End Program FindRoot

