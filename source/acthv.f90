!=====================================================================
!>acthv
! Subroutine to peform Hv=c
!---------------------------------------------------------------------
subroutine acthv(In_Vector,MOints1,MOints2,M1Len,M2Len,pString,    &
     pStep,pLocate,qString,qStep,qLocate,xRefList,ciDim,pDets,    &
     pDLen,qDets,qDLen,aDets,bDets,aElec,bElec,Orbitals,Diagonals,&
     nFrozen,nDocc,nActive,Out_Vector)
  use action_util,  only: hv_alpha,hv_beta
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
