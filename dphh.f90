!This function computes the angle between 2 vectors in degrees. remember to set L for pbc $
real*8 function angle(vec1,vec2,L) result(res)
implicit none
real*8 :: vec1(3),vec2(3)
real*8 :: xa,L
  vec1 = vec1 - L* nint(vec1/L)
  vec2 = vec2 - L* nint(vec2/L)
xa = dot_product(vec1,vec2)/sqrt(dot_product(vec1,vec1))
xa = xa/sqrt(dot_product(vec2,vec2))
if(( xa >= 1.0 ).or.( xa <= -1.0 ))xa=int(xa)
res = acos(xa)*(180.0/3.1415926536)
end function angle
!This function computes the distance between vectors remember to set L for pbc
real*8 function dis(vec1,vec2,L) result(res)
implicit none
real*8 :: vec1(3),vec2(3),vesc(3)
real*8 :: xa,L
  vesc= (vec1-vec2) - L* nint((vec1-vec2)/L)
res = sqrt(dot_product(vesc,vesc))
end function dis

! This function gives us the length of a vector
real*8 function xnorm(vec1,L) result(res)
implicit none
real*8 :: vec1(3)
real*8 :: x, L 
vec1 = vec1 - L* nint(vec1/L)
x = sqrt(dot_product(vec1,vec1))
res = x
end function xnorm


subroutine pbc(x,L)
implicit none
real*8, intent(inout), dimension(3) :: x
real*8 :: L
x = x - L*nint(x/L)
end subroutine pbc



! this subroutine strips spaces from a string
subroutine StripSpaces(string)
    character(len=*) :: string
    integer :: stringLen 
    integer :: last, actual
    stringLen = len (string)
    last = 1
    actual = 1
    do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        endif
    end do
end subroutine

!This subroutine sorts the the vectors y and z from largest to smallest based on the magnitude of the elements of y
subroutine indescendingorder(y,z,n)
implicit none
integer :: n
real*8, intent(inout), dimension(n) :: y,z
integer :: i, j
real*8 :: k, l
do i = 1,n
    do j = i + 1, n
        if ( y(i) < y(j) ) then
            k    = y(i)
            l    = z(i)
            y(i) = y(j)
            z(i) = z(j)
            y(j) = k
            z(j) = l 
        end if
    end do
end do
end subroutine indescendingorder


!This program encodes the connectivity matrix of a molecule up to the second shell in a string
Program weightedstring
implicit none
!!!!!!!!!!!!!!!!!input number of molecules in to kg!!!!!!!!
 Integer, parameter :: kg=1019
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Integer :: i, j=1,k=1,s,m=1,ii,nn,cnum=30
 real*8 ::  oh1(3)=0,oh2(3)=0,L
 real*8 ,dimension(kg*3) :: dp=0,hh=0
 character*6 :: ast
 integer     :: atnum,ig,jg,as,ios,ay, numframes, kgg=kg
 character*4 :: molnum
 character*4 :: attype
 character*3 :: molname
 character*1 :: chain
 character(7) :: bstring
 character(55) :: output
 real*8 :: x(3*kg,3), r1, r2, ang1, ang2, angw
 real*8, external :: xnorm, angle,dis

!The input of this file is a pdb file with one frame it can be modified to work if you have more than one file
!number of frames
 numframes=500
 L = 31.197

!wmicrostates.dat will contains the weights of the microstate in the string following the order of the microstates
	output = 'prd_space.pdb'
	call stripspaces(output)
	open(101,file=output,status='old',action='read')
	open(50,file= "dp_space.dat" , action='write')
	open(60,file= "hh_space.dat" , action='write')


 do ii=1 ,numframes
!Read the file and extract the adjacency matrix
	read(101,*)
	read(101,*)
	read(101,*)
	read(101,*)
	read(101,*)
	Do ig=1,kgg*3
	   read(101, *) ast,atnum,attype,molname,molnum,x(ig,:),r1,r2,chain 
	   
	end do
	read(101,*)
	read(101,*)
!Writing the dipole angles	
 dp=0
 hh=0
	do ig=1,kgg
 	  oh1 = x( 3*(ig-1)+1,:) - x(3*(ig-1)+2,:)
	  oh2 = x( 3*(ig-1)+1,:) - x(3*(ig-1)+3,:)        
	  call pbc(oh1,L)
	  call pbc(oh2,L)	  
	  dp(3*(ig-1)+1:3*(ig-1)+3) = oh1 + oh2
	  hh(3*(ig-1)+1:3*(ig-1)+3) = oh1 - oh2
	end do
	write(50,'(3057(F7.3,1X))')dp(:)
	write(60,'(3057(F7.3,1X))')hh(:)
 end do
	


end Program weightedstring

