!=====================================
!       generate input data
!=====================================
	
	subroutine readv23
	use basemod
	implicit none
	real*8 rs,a1,a2,b1,b2,ei
	real*8 minE
!-----------------------------------------------------------
	ja12max=ja1max+ja2max
!----------------------------------------------------------
	call cgvdw
!----------------------------------------------------------
	ntot=0
	minE=200
!-------------------------------------------------------------
	open(2,file=trim(dataname),status='old',action='read')

	do while (.true.)
        read(2,*,end=901) rs,a2,b2,a1,b1,ei
        
	if(ei.le. minE*10.d0) then 
		write(3,'(8f16.8)')rs,a1,b1,a2,b2,ei
        	a1=a1*pi/180.d0
        	a2=a2*pi/180.d0
        	b1=b1*pi/180.d0
        	b2=b2*pi/180.d0
        	call descriptorV23(rs,a1,a2,b1,b2,ei)
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
	subroutine descriptorV23(rs,a1,a2,b1,b2,ei)
	use basemod
	implicit none
!---------------------------------------------------------
	real*8 rs,a1,a2,b1,b2,temp,ei,vtemp
	integer ja12min
	integer ig1,ig2,ja1,ja2,j12,ma1,j12max,j12min
!-----------------------------------------------------------   
	nG1=3
	nG2=45
	ma1max=min(ja1max,ja2max)
	ja12max=ja1max+ja2max
!-----------------------------------------------------------
	allocate(G1(nG1),G2(nG2),stat=ierr)
	allocate(a1val(ja1max+1,-ma1max:ma1max),
     & a2val(-ma1max:ma1max,ja2max+1,-ma1max:ma1max),stat=ierr)
!--------------------------------------------------------------
	do ig1=1,nG1
	    G1(ig1)=exp(-rs/(0.95d0*ig1))
	enddo
!--------------------------------------------------------------------
	call legdmat(a1)
	call wigmat(a2)
	ig2=1
!--------------------------------------------------------------------
	do 2000 ja1=0,ja1max,2
	do 2000 ja2=0,ja2max
	    j12max=min(abs(ja1+ja2),ja12max)
	    j12min=max(abs(ja1-ja2),1)
!		print *,j12max,j12min
	    do 2000 j12=j12min,j12max
	    do 2000 k2=0,ja2,2
	    	temp=0.0d0
	    do 1800 ma1=-min(ja1,ja2),min(ja1,ja2)
	    	vtemp=dcos(ma1*b1+k2*b2)
	    	temp=temp
     &          +fmatcg(ma1,ja1,ja2,j12)*(a1val(ja1,ma1)
     &          *a2val(k2,ja2,-ma1)*vtemp)
1800    continue

        G2(ig2)=temp

        ig2=ig2+1
		nG2=ig2
2000   continue
        write(7,'(635E16.8)')G1,G2,ei

!	deallocate(G1,G2,a1val,a2val,stat=ierr)
	return 
	end

	
	

