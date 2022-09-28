!=====================================
!       generate input data
!=====================================
	
	subroutine readv31
	use basemod
	implicit none
	real*8 rs,a1,b1,ei
	real*8 minE
!-----------------------------------------------------------
	
!----------------------------------------------------------
	
!----------------------------------------------------------
	ntot=0
	minE=200
!-------------------------------------------------------------
	open(2,file=trim(dataname),status='old',action='read')

	do while (.true.)
        read(2,*,end=901) rs,a1,b1,ei
        
	if(ei.le. minE*10.d0) then 
		write(3,'(8f16.8)')rs,a1,b1,ei
        	a1=a1*pi/180.d0
        	b1=b1*pi/180.d0
        	call descriptorV31(rs,a1,b1,ei)
        	ntot=ntot+1
        endif
        enddo
901     close(2)
	close(3)
	close(7)
	print *,'ntrain',ntot
!---------------------------------------------------------------
!---------------------------------------------------------------------

	return 
	end


!=========================================================
!   generate the descriptor from ra and ang 
!=========================================================       
	subroutine descriptorV31(rs,a1,b1,ei)
	use basemod
	implicit none
!---------------------------------------------------------
	real*8 rs,a1,b1,temp,ei,vtemp
	integer ig1,ig2,ja1,ja2,j12,ma1,j12max,j12min
!-----------------------------------------------------------   
	nG1=5
	nG2=24
        ma1max=min(ja1max,ja2max)
	
!-----------------------------------------------------------
	allocate(G1(nG1),G2(nG2),stat=ierr)
	allocate(a1val(ja1max+1,-ma1max:ma1max),stat=ierr)
!--------------------------------------------------------------
	do ig1=1,nG1
	    G1(ig1)=exp(-rs/(0.95d0*ig1))
	enddo
!--------------------------------------------------------------------
	call legdmat(a1)
	ig2=1
!--------------------------------------------------------------------
	do 2000 ja1=1,ja1max
        do 2000 k2=0,ja1,2
            temp=0.0d0
            if(k2.eq. 0.0d0) then
                temp=a1val(ja1,k2)
            else
                temp=a1val(ja1,k2)*dcos(k2*b1)
            endif
	     print *,ig2
            G2(ig2)=temp
            ig2=ig2+1
            
	    nG2=ig2
2000   continue

        write(7,'(635E16.8)')G1,G2,ei

!	deallocate(G1,G2,a1val,a2val,stat=ierr)
	return 
	end

