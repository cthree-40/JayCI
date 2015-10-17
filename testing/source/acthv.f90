subroutine acthv (in_vec, moints1, moints2, m1len, m2len, pString,   &
     pStep, pLocate, qString, qStep, qLocate, xRefList, ciDim, pDets,&
     pDLen, qDets, qDLen, aDets, bDets, aElec, bElec, orbitals,      &
     diagonals, nFrozen, nDocc, nActive, out_Vec)
  !=====================================================================
  ! acthv
  ! -----
  ! Subroutine to peform Hv=c
  !
  ! Input:
  !  in_vec  = v
  !  moints1 = 1-e integrals
  !  moints2 = 2-e integrals
  !  m1len   = len(moints1)
  !  m2len   = len(moints2)
  !  pString = q strings of determinenats ordered by p string index
  !  pStep   = number of q strings for p string
  !  pLocate = location of ith p string
  !  qString =
  !  qStep   =
  !  qLocate =
  !  xRefList= pString(i) |-> qString(j)
  !  ciDim   = expansion size
  !  pDets   = p strings
  !  pDLen   = length of pDets
  !  qDets   = q strings
  !  qDLen   = length of qDets
  !  aDets   = num of alpha strings without expansion restrictions
  !  bDets   = num of beta  strings without expansion restrictions
  !  aElec   = alpha electrons
  !  bElec   = beta  electrons
  !  orbitals= number of orbitals
  !  diagonals = diagonal elements of H
  !  nfrozen = number of frozen core orbitals
  !  ndocc   = number of docc orbitals
  !  nactive = number of active space orbitals
  !
  ! Output:
  !  out_vec
  !---------------------------------------------------------------------
  use action_util,  only: hv_alpha,hv_beta
  implicit none

  ! .. INPUT arguments ..
  integer, intent(in) :: m1len, m2len, ciDim, pDLen, qDLen, aDets
  integer, intent(in) :: bDets, aElec, bElec, orbitals
  integer, intent(in) :: nfrozen, ndocc, nactive
  integer, dimension(pDLen), intent(in) :: pStep,pLocate,pDets
  integer, dimension(qDLen), intent(in) :: qStep,qLocate,qDets
  integer, dimension(ciDim), intent(in) :: pString,qString
  integer, dimension(ciDim), intent(in) :: xRefList
  real*8,  dimension(m1len), intent(in) :: moints1
  real*8,  dimension(m2len), intent(in) :: moints2
  real*8,  dimension(ciDim), intent(in) :: in_vec, diagonals
  
  ! .. OUTPUT arguments ..
  real*8,  dimension(ciDim), intent(out) :: out_vec

  ! .. LOCAL arguments ..
  integer          :: i
  double precision :: tstart, tfinal
  
  ! Zero out Out_vector
  out_vec = 0d0

  ! get time
  call cpu_time(tstart)

  ! Diagonal contribution
  do i=1,ciDim
     out_vec(i) = out_vec(i) + diagonals(i)*in_vec(i)
  end do !i

  ! Alpha string contribution
  CALL hv_alpha(in_vec,moints1,m1len,moints2,m2Len,pString,  &
       pStep,pLocate,pDets,qString,qStep,qLocate,qDets,ciDim,&
       pDLen,qDLen,aDets,bDets,aElec,bElec,orbitals,         &
       nfrozen,ndocc,nactive,out_vec)
  
  ! Beta string contribution
  CALL hv_beta(in_vec,moints1,m1len,moints2,m2len,pString,   &
       pStep,pLocate,pDets,qString,qStep,qLocate,qDets,ciDim,&
       pDLen,qDLen,aDets,bDets,aElec,bElec,orbitals,         &
       nfrozen,ndocc,nactive,xRefList,out_vec)

  call cpu_time(tfinal)
  write(*, "(1x, f6.3,' seconds to perform Hv=c.')") (tfinal - tstart)
  
  RETURN
end subroutine acthv
