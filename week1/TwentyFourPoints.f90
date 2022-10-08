Program TwentyFourPoints
!To solve the problem of twenty-four points
character Mode,conti
integer i,j,k,m,n,p,c
real poker1(4),poker2(3),poker3(2),poker4(1)
character(len=100) segment1,segment2,segmentmid1,segmentmid2,segmentfinal
do
    1 print*,"Do you want to input by yourself or generate randomly?(I/R)"
    read*,Mode
    If(Mode=='I') then
        print*,"Please input four integers in the range of [1,13]:"
        read*,poker1(1),poker1(2),poker1(3),poker1(4)
    else if(Mode=='R') then
        do i=1,4
            call random_number(temp)
            poker1(i)=1+int(temp*13)    !把[0,1)的实数随机数转换为1-13的整数随机数
        end do
        print*,poker1(1),poker1(2),poker1(3),poker1(4)
    else
        print*,"Invalid Input!"
        goto 1
    end if

    do i=1,4    !第一轮运算：数组中有4个数字
        do j=i+1,4
            do k=1,6
                write(segment1,*)poker1(i)
                write(segment2,*)poker1(j)  !生成两个操作数对应的字符串
                call operate(4,poker1,poker2,k,i,j,segment1,segment2,segmentmid1)   !运算与字符串拼接并生成新数组
                do m=1,3    !第二轮运算：数组中有3个数字
                    do n=m+1,3
                        do p=1,6
                            if(m==1) then
                                segment1=trim(segmentmid1)  !若这一轮选中的第一个操作数poker2(m)是上一轮的运算结果而非初始牌面中的数，则其字符串应为上一轮的拼接结果
                            else
                                write(segment1,*)poker2(m)
                            end if
                            write(segment2,*)poker2(n)
                            call operate(3,poker2,poker3,p,m,n,segment1,segment2,segmentmid2)
                            do c=1,6    !第三轮运算：数组中有2个数字
                                segment1=trim(segmentmid2)
                                if(m==1) then
                                    write(segment2,*)poker3(2)  !第一个操作数对应的字符串是上一轮的拼接结果
                                else
                                    segment2=trim(segmentmid1)
                                end if
                                call operate(2,poker3,poker4,c,1,2,segment1,segment2,segmentfinal)
                                if(abs(poker4(1)-24)<0.0001) then   !判断并输出
                                    Print*,trim(trim(segmentfinal)//"="),poker4(1)
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

    print*,"Do you want to continue?(y/n)"
    read*,conti
    if(conti=='y') then
        cycle
    else if(conti=='n') then
        exit
    else
        print*,"Invalid Input!"
    end if
end do
contains
    subroutine operate(num,pokerf,pokerl,oper,i,j,segment1,segment2,segmentmid)    !进行运算和字符串拼接，从pokerf数组生成pokerl数组。ij是被挑选的两个操作数，oper是运算符号，segment1,2是两个操作数对应的字符串，segmentmid是字符串拼接结果
        integer num,i,j,m,n,oper
        real pokerf(num),pokerl(num-1),step
        character(len=*) segment1,segment2,segmentmid

        select case(oper)   !减法和除法不满足交换律，于是两个数的运算有6种情况
            case(1)
                step=pokerf(i)+pokerf(j)    !运算
                segmentmid="("//trim(segment1)//"+"//trim(segment2)//")"    !拼接字符串
            case(2)
                step=pokerf(i)-pokerf(j)
                segmentmid="("//trim(segment1)//"-"//trim(segment2)//")"
            case(3)
                step=pokerf(i)*pokerf(j)
                segmentmid="("//trim(segment1)//"*"//trim(segment2)//")"
            case(4)
                if(pokerf(j)/=0) then   !解决了分母为0的情况。
                    step=pokerf(i)/pokerf(j)
                    segmentmid="("//trim(segment1)//"/"//trim(segment2)//")"
                end if
            case(5)
                step=pokerf(j)-pokerf(i)
                segmentmid="("//trim(segment2)//"-"//trim(segment1)//")"        
            case(6)
                if(pokerf(i)/=0) then
                    step=pokerf(j)/pokerf(i)
                    segmentmid="("//trim(segment2)//"/"//trim(segment1)//")"
                end if
        end select
        
        !把运算结果和pokerf中没有用到过的数字存到pokerl中  
        pokerl(1)=step
        do m=1,num
            if((m==i).or.(m==j)) then
                cycle
            else
                pokerl(2)=pokerf(m)
                exit
            end if
        end do
        if(num==4) then
            do n=1,num
                if((n==i).or.(n==j).or.(n==m)) then
                    cycle
                else
                    pokerl(3)=pokerf(n)
                    exit
                end if
            end do
        end if
    end subroutine
end program TwentyFourPoints
