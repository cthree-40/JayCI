!=====================================================================
!>acthv
! Subroutine to peform Hv=c
!---------------------------------------------------------------------
subroutine acthv(In_Vector,MOints1,MOints2,M1Len,M2Len,pString,    &
     pStep,pLocate,qString,qStep,qLocate,xRefList,ciDim,pDets,    &
     pDLen,qDets,qDLen,aDets,bDets,aElec,bElec,Orbitals,Diagonals,&
     nFrozen,nDocc,nActive,Out_Vector)
  use action_util,  only: hv_alpha,hv_beta
#ifdef DEBUGGING
  use construct
  use addressing
#endif
  IMPLICIT NONE
  integer,intent(IN)      ::M1Len,M2Len,ciDim,pDLen,qDLen,aDets,&
       bDets,aElec,bElec,Orbitals,nFrozen,nDocc,nActive
  integer,dimension(pDLen),intent(IN) ::pStep,pLocate,pDets
  integer,dimension(qDLen),intent(IN) ::qStep,qLocate,qDets
  integer,dimension(ciDim),intent(IN)     :: pString,qString
  integer,dimension(ciDim),intent(IN) ::xRefList
  real*8,dimension(M1Len),intent(IN)  ::MOints1
  real*8,dimension(M2Len),intent(IN)  ::MOints2
  real*8,dimension(ciDim),intent(IN)  ::In_Vector,Diagonals
  real*8,dimension(ciDim),intent(OUT)    ::Out_Vector
  integer     :: i
#ifdef DEBUGGING
  real*8, dimension(ciDim) :: hcolumn, invec, outvec
  integer                  :: testcol
  integer, dimension(ciDim):: detExp
  integer, dimension(aElec) :: alphastr
  integer, dimension(bElec) :: betastr
  hcolumn = 0d0
  invec = 0d0
  outvec = 0d0
  testcol=1
  invec(testcol) = 1d0
  open(unit=3,file="det.list",status="old")
  read(unit=3,fmt="(1x,I10)") (detExp(i), i=1, ciDim )
  close(unit=3)
  do i=1, ciDim
     hcolumn(i) = ham_element(detExp(i),detExp(testcol), &
          MOints1, M1Len, MOints2, M2Len, aElec, bElec,  &
          Orbitals )
  end do
  ! perform Hv
  do i=1, ciDim
     outvec(i) = outvec(i)+Diagonals(i)*invec(i)
  end do
  call hv_alpha(invec,MOints1,M1Len,MOints2,M2Len,pString,      &  
       pStep,pLocate,pDets,qString,qStep,qLocate,qDets,ciDim,      &   
       pDLen,qDLen,aDets,bDets,aElec,bElec,Orbitals,     &             
       nFrozen,nDocc,nActive,outvec)
  CALL hv_beta(invec,MOints1,M1Len,MOints2,M2Len,pString,      & 
       pStep,pLocate,pDets,qString,qStep,qLocate,qDets,ciDim,      & 
       pDLen,qDLen,aDets,bDets,aElec,bElec,Orbitals,     &           
       nFrozen,nDocc,nActive,xRefList,outvec)
  ! compare the two
  do i=1, ciDim
     write(*, "(1x,f20.12,A,f20.12)") hcolumn(i), "  ", outvec(i)
  end do
  call GenOrbString(1, aElec, Orbitals, aDets, alphastr )
  call GenOrbString(10263, bElec, Orbitals, bDets, betastr  )
  write(*,"(1x,7i3)") alphastr
  write(*,"(1x,7i3)") betastr
  call GenOrbString(571, bElec, Orbitals, bDets, betastr )
  write(*,"(1x,7i3)") betastr
  stop
#endif
  ! Zero out Out_vector
  Out_Vector=0d0
  ! Diagonal contribution
  do i=1,ciDim
     Out_Vector(i)=Out_Vector(i)+Diagonals(i)*In_Vector(i)
  end do !i
  ! Alpha string contribution
  print *, "     Computing alpha contribution "
  CALL hv_alpha(In_Vector,MOints1,M1Len,MOints2,M2Len,pString,      &
       pStep,pLocate,pDets,qString,qStep,qLocate,qDets,ciDim,      &
       pDLen,qDLen,aDets,bDets,aElec,bElec,Orbitals,     &
       nFrozen,nDocc,nActive,Out_Vector)
  ! Beta string contribution
  print *, "     Computing beta contribution "
  CALL hv_beta(In_Vector,MOints1,M1Len,MOints2,M2Len,pString,      &
       pStep,pLocate,pDets,qString,qStep,qLocate,qDets,ciDim,      &
       pDLen,qDLen,aDets,bDets,aElec,bElec,Orbitals,     &
       nFrozen,nDocc,nActive,xRefList,Out_Vector)
      
  RETURN
end subroutine acthv
!=====================================================================
