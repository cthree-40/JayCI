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
#ifdef HAMSOLVE
  use construct
  use addressing
#endif
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
  integer     :: i

#ifdef HAMSOLVE
  ! .. HAMSOLVE arguments ..
  real*8, dimension(ciDim, ciDim) :: hamiltonian
  integer,dimension(ciDim)        :: detlist
  real*8, dimension(ciDim)        :: test_in, test_out
  character*1 :: uplo, jobz, rnge
  real*8      :: abstol, dlamch, vl, vu
  integer     :: lwork, liwork, il, iu, info, eigfound
  integer, dimension(:), allocatable :: isuppz
  real*8,  dimension(:), allocatable :: work, iwork, eig_val
  real*8,  dimension(:,:), allocatable :: eig_vec
  integer :: p, q
  integer, dimension(aelec) :: pstring1
  integer, dimension(belec) :: qstring1
  character*300 :: fmt9, fmt8
  
  test_in = 0d0
  test_in(1) = 1d0
  
  write (*, "(1x,'Explicitly constructing Hamiltonian for diagonalization')")
  
  open (unit = 4, file = "det.list", status = "old", action= "read")
  read (unit = 4, fmt = "(1x,i10)") (detlist(i), i = 1, ciDim)
  close(unit = 4)

  do i = 1, ciDim
          test_out(i) = ham_element(detlist(i), 1, moints1, m1len, moints2, &
            m2len, aelec, belec, orbitals)
  end do
  !call exp_construct(moints1, m1len, moints2, m2len, ciDim, &
  !  aelec, belec, orbitals, detlist, hamiltonian)

  ! perform test multiplication
  !call dgemv ('n', ciDim, ciDim, 1d0, hamiltonian, &
  !  ciDim, in_vec, 1, 0d0, test_out, 1)
  
  uplo = 'l'
  jobz = 'v'
  rnge = 'a'
  abstol = dlamch('Safe minimum')
  lwork = 30*ciDim + 10
  liwork= 15*ciDim
  info = 0
  !allocate(isuppz(2 * ciDim))
  !allocate(work(lwork))
  !allocate(iwork(liwork))
  !allocate(eig_vec(ciDim, ciDim))
  !allocate(eig_val(ciDim))

  !call dsyevr (jobz, rnge, uplo, ciDim, hamiltonian, ciDim, vl, vu, il, &
  !  iu, abstol, eigfound, eig_val, eig_vec, ciDim, isuppz, &
  !  work, lwork, iwork, liwork, info)

  !write (*, "(1x,'Eigenvalue = ',f15.8)") eig_val(1)

  !deallocate(isuppz, work, iwork, eig_vec, eig_val)
  
  !stop
#endif
  
  ! Zero out Out_vector
  out_vec = 0d0
  
  ! Diagonal contribution
  do i=1,ciDim
     out_vec(i) = out_vec(i) + diagonals(i)*test_in(i)
  end do !i

  ! Alpha string contribution
  print *, "     Computing alpha contribution "
  CALL hv_alpha(test_in,moints1,m1len,moints2,m2Len,pString,  &
       pStep,pLocate,pDets,qString,qStep,qLocate,qDets,ciDim,&
       pDLen,qDLen,aDets,bDets,aElec,bElec,orbitals,         &
       nfrozen,ndocc,nactive,out_vec)

  ! Beta string contribution
  print *, "     Computing beta contribution "
  CALL hv_beta(test_in,moints1,m1len,moints2,m2len,pString,   &
       pStep,pLocate,pDets,qString,qStep,qLocate,qDets,ciDim,&
       pDLen,qDLen,aDets,bDets,aElec,bElec,orbitals,         &
       nfrozen,ndocc,nactive,xRefList,out_vec)
      
#ifdef HAMSOLVE
  write (fmt9,'("(1x,i10,i10,",i0,"i3,i10,",i0,"i3,5x,2f15.8)")') aelec, belec
  write (fmt8,'("(1x,a10,a10,a",i0,",a10,a",i0,",5x,2a15)")') aelec*3, belec*3 
  ! print and compare two vectors
  write (*, "(1x,A)") "Comparing Hv to acthv()"
  write (*, "(1x,A)") ""
  write (*, fmt8)  "Det Index", "a index", "a string", &
    "b index", "b string", "Hv", "acthv()"
  write (*, "(1x,A,A)") &
    "--------------------------------------------------------------------",&
    "--------------------------------------------------------------------"
  do i = 1, ciDim
          call K2Indc(detlist(i), belec, orbitals, p, q)
          call GenOrbString_Alpha(p, aelec, orbitals, aDets, pDets(pDLen),&
            pstring1)
          call GenOrbString_Beta (q, belec, orbitals, bDets, qDets(qDLen),&
            qstring1)
          write (*, trim(fmt9)) detlist(i), p, pstring1, q, qstring1, &
            test_out(i), out_vec(i)
  end do
8 format(1x,2a10,a30,a10,a30,5x,2a15)
9 format(1x,i10,i10,5i3,i10,5i3,5x,f15.8,f15.8)
  stop
#endif
  RETURN
end subroutine acthv
