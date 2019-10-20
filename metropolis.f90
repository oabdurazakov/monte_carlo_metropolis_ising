program ising

  implicit none 
  integer::i,neq,j,l,m,lmax
  integer, parameter::n=20
  integer, dimension(n)::s
  real(kind=8)::x,t,e,e2,mg,mg2,e_e,m_e


call set_seed


  do i=1,n

     call random_number(x)
     if(x-0.5>0)then 
        s(i)=1
     else
        s(i)=-1
     end if
     !write(*,*)s(i)
  end do



lmax=50000
t=0.0
do m=1,50
e=0.0
e2=0.0
mg2=0.0
mg=0.0
do l=1,lmax
do j=1,10000
call mc(s)
!write(*,*)j, mag(s),energy(s)
end do

mg=mg+mag(s)
mg2=mg2+(mag(s))**2
e=e+energy(s)
e2=e2+(energy(s))**2
end do
e2=e2/lmax
e=e/lmax
mg=mg/lmax
mg2=mg2/lmax
e_e=-tanh(1.0/t)
t=t+0.1
write(*,*)t,e/n,e_e,mg/n,sqrt((e2-e*e)/(lmax-1))/n,sqrt((mg2-mg*mg)/(lmax-1))/n
end do



contains

function energy(s)
integer::u1,u,sum
integer::energy
integer, dimension(n)::s
sum=0
do i=1,n
   if(i==n)then
   u1=s(i)*s(1)
   else
   u=s(i)*s(i+1)
   end if
   sum=sum+u
end do
energy=-(sum+u1)

end function energy

function mag(s)
integer::mag,rho
integer, dimension(n)::s
rho=0
do i=1,n
rho=rho+s(i)
end do
mag=rho
end function mag


subroutine mc(s)

  implicit none
  integer, dimension(n)::s
  real(kind=8)::r,y
  integer::delE
  call random_number(r)
  i=n*r+1

if(i==n) then
   delE=s(i)*(s(1)+s(i-1))
elseif (i==1)then
   delE=s(i)*(s(n)+s(2))
else
   delE=s(i)*(s(i+1)+s(i-1))
end if
call random_number(y)
if(delE<=0) then
   s(i)=-s(i)
elseif(exp(-2*delE/t)>=y)then
      s(i)=-s(i)
   else 
      s(i)=s(i)
   end if
   
end subroutine mc		

subroutine set_seed
    integer,allocatable :: seed(:)
    integer :: n_seed,clock,i
    real(kind=8) :: x
    
    ! get the number of seed integers (varies with system)
    !
    ! note use of 'keyword argument' size
    call random_seed(size=n_seed)
    
    ! allocate seed array
    allocate(seed(n_seed))
    
    ! initialize the seed from the time
    call system_clock(count=clock)

    do i=1,n_seed
       seed(i)=clock/i  ! could do other things here
    enddo

    ! actually initialize the generator
    call random_seed(put=seed)

    deallocate(seed)
  end subroutine 


end program ising
