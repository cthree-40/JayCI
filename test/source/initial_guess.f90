!********************************************************************
!>INITIALGUESS
! This module contains subroutines for generating the initial
!  basis for the davidson algorithm
!====================================================================
MODULE InitialGuess
  IMPLICIT NONE
CONTAINS
!====================================================================
!>genVec_dgemm
! Generates ritz vectors using DGEMM
!     (In_Vec)_nxm (matrix)_mxo = (Out_Vec)_nxo
!--------------------------------------------------------------------
  subroutine genVec_dgemm( imatrix, n, m, o, In_Vec, Out_Vec )
    IMPLICIT NONE
    integer,intent(IN)      ::n, m, o
    real*8,dimension(n,m),intent(IN)    ::In_Vec
    real*8,dimension(m,o),intent(IN)    ::imatrix
    real*8,dimension(n,o),intent(OUT)    ::Out_Vec
    character(len=2)  :: transa, transb
    integer           :: lda,ldb,ldc,p,q,r
    real*8            :: a,b
    transa      = "No"      !Use transpose of A
    transb      = "No"      !Use transpose of B
    p=n
    q=o
    r=m
    a=1d0
    b=0d0
    lda=n
    ldb=m
    ldc=n
    CALL DGEMM(transa,transb,p,q,r,a,In_Vec,lda,imatrix,ldb,b, &
      Out_Vec,ldc)
    RETURN
  end subroutine genVec_dgemm
  !====================================================================
  !>preDiagSubBlock
  ! Diagonalizes an explicitly constructed subblock of the Hamiltonian
  !--------------------------------------------------------------------
  subroutine preDiagSubBlock(ciDim,MOints1,M1Len,MOints2,M2Len,     &
    Det_Exp,subBlockDim,initGuessDim,aElec,bElec,Orbitals,&
    initial_vecs)
    use construct,    only:ham_element
    use david_util,    only:diag_dsyevr
    IMPLICIT NONE
    integer,intent(IN)      ::ciDim,M1Len,M2Len,subBlockDim,    &
      initGuessDim,aElec,bElec,Orbitals
    integer,dimension(ciDim),intent(IN) ::Det_Exp
    real*8,dimension(M1Len),intent(IN)  ::MOints1
    real*8,dimension(M2Len),intent(IN)  ::MOints2
    real*8,dimension(ciDim,initGuessDim),intent(OUT)::initial_vecs
    real*8,dimension(:,:),ALLOCATABLE   ::ham_block,unit_vecs,tmp_vecs,&
      tmp_eigvecs
    real*8,dimension(:),ALLOCATABLE     ::tmp_eigvals
    integer     ::i,j
    !Allocate arrays
    ALLOCATE(ham_block(subBlockDim,subBlockDim))
    ALLOCATE(unit_vecs(ciDim,subBlockDim))
    ALLOCATE(tmp_eigvecs(ciDim,subBlockDim))
    ALLOCATE(tmp_eigvals(subBlockDim))
    !Construct subBlockDim x subBlockDim Hamiltonian
    do i=1,subBlockDim
            do j=1,subBlockDim
                    ham_block(j,i)=ham_element(Det_Exp(j),Det_Exp(i),&
                      MOints1,M1Len,MOints2,M2Len,aElec,bElec,   &
                      Orbitals)
            end do
    end do
    !Diagonalize ham_block using DSYEVR
    print *, "Diagonalizing ham_block"
    CALL diag_dsyevr(ham_block,subBlockDim,tmp_eigvecs,tmp_eigvals)
    
    ! PP Block for explicit hamiltonian diagonalization
#ifdef HAMILTONIAN
    write(*,"(A,F15.8)") "Eigenvalue 1 = ", tmp_eigvals(1)
#endif
    
    ALLOCATE(tmp_vecs(ciDim,initGuessDim))
    tmp_vecs(1:ciDim,1:initGuessDim)=tmp_eigvecs(1:ciDim,1:initGuessDim)
    DEALLOCATE(tmp_eigvecs,tmp_eigvals)
    
    !Generate unit vectors
    unit_vecs=0d0
    do i=1,subBlockDim
            unit_vecs(i,i)=1d0
    end do
    !Generate ritz vectors
    CALL genVec_dgemm(tmp_vecs,ciDim,subBlockDim,initGuessDim,unit_vecs,&
      initial_vecs)
    !Deallocate arrays
    DEALLOCATE(ham_block,unit_vecs,tmp_vecs)
    RETURN !Return initial_vecs
  end subroutine preDiagSubBlock
!====================================================================
END MODULE InitialGuess
