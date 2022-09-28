cc=======================================================
	module basemod
	character*99 :: dataname                           ! file name containing ab initio points
	character*99 :: system                             ! 
	integer :: ntot                                    ! total number of ab initio points
	integer :: noders                                  ! number of radial grid
	integer :: nodea1,nodea2,nodeb1,nodeb2,nodec1      ! number of angle grid

	integer :: nG1                                     ! number of radius descriptor
	integer :: nG2                                     ! number of angular descriptor

	integer :: ja1max,ja2max,ja12max,k2,ma1max                ! angular quantum number

	integer,parameter :: nmax=1000000
	real*8 :: pi
	real*8,allocatable :: G1(:)                        ! stores radial descptor
	real*8,allocatable :: G2(:)                        ! stores angular descptor 
	real*8,allocatable :: fmatcg(:,:,:,:)              !Clebsch-Gordan coefficient
	
	real*8,allocatable :: a1val(:,:),a2val(:,:,:),a3val(:,:)
	integer :: ierr
	end module basemod


!======================================
!---
!======================================
       program main
		use basemod
       implicit none
!--------------------------------------
	
       open(5,file='input',status='old')
       open(7,file='ouput',status='unknown')
	   open(3,file='train.dat',status='unknown')
	   pi=dacos(-1.0d0)
	   read(5,*)system
	   read(5,*)dataname
	   read(5,*)
	   read(5,*)ja1max
	   read(5,*)
	   read(5,*)ja2max
	   if(system.eq.'v21') then
			call readv21
			print *,'readv21'
	   else if(system.eq. 'v31') then
			call readv31
			print *, 'readv31'
	   else if(system.eq. 'v22') then
			call readv22
			print *,'readv22'
	   else if(system.eq. 'v23') then
			call readv23
			print *, 'readv23'
	   else
			stop 'no system found in this program'
	   endif


       write(*,*)'program finished'
       stop
       end program main
       


	
!======================================================================
!  *   Program       legdmat                                                 
!  *   Function      generate legender data from function DM     
!======================================================================
	subroutine legdmat(a1)
	use basemod
	implicit none                                                                   
	integer j1,ma1
	real*8 beta,DM,a1
    	do 3100 j1=0,ja1max
    	do 3100 ma1=-ma1max,ma1max
        	beta=a1
        	a1val(j1,ma1)=DM(j1,0,ma1,beta) 
3100  continue 


	return 
	end
!======================================================================
!  *   Program       legdmat                                                 
!  *   Function      generate legender data from function DM     
!======================================================================
	subroutine legdmat1(a2)
	use basemod
	implicit none                                                                   
	integer j2,ma1
	real*8 beta,DM,a2
    	do 3100 j2=0,ja2max
    	do 3100 ma1=-ma1max,ma1max
        	beta=a2
        	a3val(j2,ma1)=DM(j2,0,ma1,beta) 
3100  continue 


	return 
	end
!======================================================================
!  *   Program       wigmat                                                 
!  *   Function      generate small wigner data from function DM     
!======================================================================
	subroutine wigmat(a2)
		use basemod
		implicit none                                                                   
		integer j2,ma1,k
		real*8 beta,DM,a2

		do 3000 j2=0,ja2max
		do 3000 ma1=-ma1max,ma1max
		do 3000 k=-ma1max,ma1max
			  beta=a2
			  a2val(k,j2,ma1)=DM(j2,k,ma1,beta)
3000  continue 
	
	!-----------------test--------------------------------------------- 	
		return 
		end
	
c  ===========================================================================
c      coeffcg       cgvdw
c  ===========================================================================
       subroutine    cgvdw
       use basemod
c      -----------------------------------------------------------------
       implicit      none
       real*8        FUNCTCG
       integer       j12,ja1,ja2,ma1,ma2
       
c      -----------------------------------------------------------------
     
c      -----------------------------------------------------------------
       allocate(fmatcg(-ma1max:ma1max,0:ja1max,
     &                       0:ja2max,0:ja12max),
     & stat=ierr)
c      -----------------------------------------------------------------
       fmatcg=0.0d0

          do 2000 ja1=0,ja1max
          do 2000 ja2=0,ja2max
          do 2000 j12=0,ja12max
          do 2000 ma1=-ja1,ja1
             ma2=-ma1
             fmatcg(ma1,ja1,ja2,j12)=FUNCTCG(ja1,ma1,ja2,ma2,j12)
 2000     continue
 
 		return
 		end
	
	
	
	
	
	
	
	
	
	
	
	
	
	
