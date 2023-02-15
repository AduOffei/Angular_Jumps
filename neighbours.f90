!This function computes the angle between 2 vectors in degrees. remember to set L for pbc $
real*8 function angle(vec1,vec2,L) result(res)
implicit none
real*8 :: vec1(3),vec2(3)
real*8 :: xa, L
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
real*8 :: xa, L 
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
        end if
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
Integer, parameter :: kg = 1019
Integer :: i, j = 1 , k = 1, s,m=1,ii,nn,cnum=30, defects(kg*6)=0,kk,ll
real*8 :: ang_ct,dist_ct
 character*6 :: ast
integer     :: atnum,ig,jg,as,ios,ay, numframes, adein(6)= 0 ,adeout(6)=0
 character*4 :: molnum
 character*4 :: attype
 character*3 :: molname
 character*1 :: chain
 character(7) :: bstring
 character(55) :: output
 real*8 :: x(3*kg,3), r1, r2, ang1, ang2, L,kgg=kg!, angw
 real*8, external :: xnorm, angle, dis
 integer :: quickmatrix(kg,kg),angw,indices(kg),si,indices2(kg)=0

!The input of this file is a pdb file with one frame it can be modified to work if you have more than one file
!the output is microstates.dat which contains a string encoding connectivity
!the input  is wmicrostates.dat which contains the weights of the connectivity
!The output atomstring.dat contains the list of atoms
!and mergedmicrostates2.dat contains all the atoms in the list p.s. I know this is redundant

!number of frames
numframes= 500
!angular cutoff
ang_ct = 30
!distance cutoff
dist_ct =3.5
L = 31.197




	open(unit = 50, file = 'neigh_space.dat', action = 'write')

	output = 'prd_space.pdb'
	call stripspaces(output)
	open(101,file=output,status='old',action='read')

	!allocate(quickmatrix(1019,1019))


do ii=1 ,numframes
	quickmatrix=0
	!Here I read the file and extract the adjacency matrix
	read(101,*)
	read(101,*)
	read(101,*)
	read(101,*)
	read(101,*)
	Do ig=1,kg*3
	   read(101, *) ast,atnum,attype,molname,molnum,x(ig,:),r1,r2,chain
	end do
	read(101,*)
	read(101,*)
	! here I compute the adjacency matrix
	! the cutoff for the distance is 3.5 and is binary
	! Adjacency matrix computed here
	do ig=1,kg-1
	 do jg=1+ig,kg

if(dis(x(3*(ig-1)+1,:),x(3*(jg-1)+1,:),L) <= dist_ct )then
	       quickmatrix(ig,jg)  = 1
	       quickmatrix(jg,ig)  = 1

	    
	!this cutoff is multiplied by the hydrogen bond angle which is non binary
	!but has a cutoff of 30
	! first I do this for the right out matrix ig to jg
	! i take the smallest of the hydrogen bond angles
	!this defined as O1 O2 to O1 O1h1_x x can be 1 or 2 
	   ang1=angle(x(3*(ig-1)+1,:)-x(3*(jg-1)+1,:),x(3*(ig-1)+1,:)-x(3*(ig-1)+2,:),L)
	   ang2=angle(x(3*(ig-1)+1,:)-x(3*(jg-1)+1,:),x(3*(ig-1)+1,:)-x(3*(ig-1)+3,:),L)
	   ang1=min(ang1,ang2)
!hear i define the angular weights
	   if(ang1 >= ang_ct)then
	       angw = 0
	   else
	       angw = 1 
	   end if

	   quickmatrix( ig, jg) = quickmatrix( ig, jg)*angw
	! first I do this for the right out matrix jg to ig
	! i take the smallest of the hydrogen bond angles
	!this defined as O1 O2 to O1 O1h1_x x can be 1 or 2 
	   ang1= angle(x(3*(jg-1)+1,:)-x(3*(ig-1)+1,:),x(3*(jg-1)+1,:)-x(3*(jg-1)+2,:),L)
	   ang2= angle(x(3*(jg-1)+1,:)-x(3*(ig-1)+1,:),x(3*(jg-1)+1,:)-x(3*(jg-1)+3,:),L)
	   ang1=min(ang1,ang2)
	   if(ang1>=ang_ct)then
	       angw = 0
	   else
	       angw = 1 !1.0 -(1.0/ang_ct)*(ang1-ang_ct)
	   end if
	       quickmatrix(jg,ig)= quickmatrix(jg,ig)*angw
end if
	 end do
	end do


	
!! the tag along stuff
do i=1,kg

        kk = 0
        ll = 0
        adein=0
        adeout =0

indices=0
indices2=0

indices = Merge( 0, [ ( si, si = 1, kg ) ], quickmatrix(i,:) == 0 )
adeout = Pack( indices, indices /= 0 )
indices2 = Merge( 0, [ ( si, si = 1, kg ) ], quickmatrix(:,i) == 0 )
adein = Pack( indices2, indices2 /= 0 )
           defects(1+(i-1)*6:3+(i-1)*6)= adeout(1:3)
           defects(4+(i-1)*6:6+(i-1)*6)= adein(1:3)
 end do
 write(50,*)defects(:)
 defects=0
end do

end Program weightedstring


!cat mergedmicrostates2.dat|awk 'NF>1'|awk 'NR%4==0'|awk '{for (i=1;i<=64;i++) g[i,NR]=$i;for (i=1;i<=64;i++) w[i,NR]=$(i+64)}END{for (j=1;j<=(NR-1);j++) for (i=j+1;i<=NR;i++) {d=0;for (k=1;k<=64;k++) if (g[k,i]==g[k,j]) {d=d+(w[k,i]-w[k,j])**2}else{d=d+1};print i,j,d}}'>distances_matrix_noweighted_every4
