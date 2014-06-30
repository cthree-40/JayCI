!>DAVIDSON
! Davidson algorithm subroutine
!====================================================================
SUBROUTINE Davidson(initGuessDim,initGuess_vecs,Diagonals,MOints1,  &
            M1Len,MOints2,M2Len,ciDim,pString,pLocate,pStep,qString,&
            qLocate,qStep,xRefList,pDets,pDLen,qDets,qDLen,aDets,   &
            bDets,aElec,bElec,Orbitals,KryMin,KryMax,res_tol,Roots, &
            Max_Iter,nFrozen,nDocc,nActive,Eig_Values,Eig_Vectors)
      use david_util
      IMPLICIT NONE
      integer,intent(IN)      ::initGuessDim,M1Len,M2Len,ciDim,pDLen, &
            qDLen,aDets,bDets,aElec,bElec,Orbitals,KryMin,KryMax,Roots,&
            Max_Iter,nFrozen,nDocc,nActive
      real*8,intent(IN) :: res_tol
      integer,dimension(ciDim,2),intent(IN)     ::qString,pString
      integer,dimension(qDLen),intent(IN)       ::qDets,qLocate,qStep
      integer,dimension(pDLen),intent(IN)       ::pDets,pLocate,pStep
      integer,dimension(ciDim),intent(IN)       ::xRefList
      real*8,dimension(ciDim),intent(IN)        ::Diagonals
      real*8,dimension(M1Len),intent(IN)        ::MOints1
      real*8,dimension(M2Len),intent(IN)        ::MOints2
      real*8,dimension(ciDim,initGuessDim),intent(IN) ::initGuess_vecs
      real*8,dimension(Roots),intent(OUT)       ::Eig_Values
      real*8,dimension(ciDim,Roots),intent(OUT) ::Eig_Vectors
      integer     ::i,j,k
      integer     ::curr_dim,curr_root
      real*8      ::resid_norm,ddot
      real*8,parameter  ::max_ovrlap=1.0D-16
      real*8,dimension(:,:),ALLOCATABLE   ::basis_vecs,hv_vectors,&
            subsp_ham,Kry_eigvec
      real*8,dimension(:),ALLOCATABLE     ::Kry_eigval,Residual,  &
            new_Vector
      !Make curr_root=1
      curr_root=1
      !Allocate arrays
      ALLOCATE(basis_vecs(ciDim,KryMax))
      ALLOCATE(hv_vectors(ciDim,KryMax))
      basis_vecs=0d0
      hv_vectors=0d0
      print *, initGuessDim
      !Equate basis_vecs and initial guess vectors
      basis_vecs(1:ciDim,1:KryMin)=initGuess_vecs(1:ciDim,1:KryMin)
      !MAIN LOOP --------------------------------
      main_loop:  do i=1, int(Max_Iter/(KryMax-KryMin))
            print *, "Performing Hv on initial space..."
            !Perform Hv on initial space
            CALL David_Initial(basis_vecs,KryMin,ciDim,MOints1,M1Len,&
                  MOints2,M2Len,pString,pStep,pLocate,qString,qStep, &
                  qLocate,xRefList,pDets,qDets,pDLen,qDLen,aDets,bDets,&
                  aElec,bElec,Orbitals,Diagonals,nFrozen,nDocc,nActive,&
                  hv_vectors)
            !Construct subspace hamiltonian
            if(ALLOCATED(subsp_ham))DEALLOCATE(subsp_ham)
            ALLOCATE(subsp_ham(KryMin,KryMin))
            CALL build_subham(basis_vecs,hv_vectors,ciDim,KryMin,subsp_ham)
            !Diagonalize this matrix
            if(ALLOCATED(Kry_eigvec))DEALLOCATE(Kry_eigvec)
            if(ALLOCATED(Kry_eigval))DEALLOCATE(Kry_eigval)
            ALLOCATE(Kry_eigvec(KryMin,KryMin))
            ALLOCATE(Kry_eigval(KryMin))
            CALL diag_dsyevr(subsp_ham,KryMin,Kry_eigvec,Kry_eigval)
            !Deallocate subsp_ham, it has been destroyed by DSYEVR()
            DEALLOCATE(subsp_ham)
            !Set space dimension
            curr_dim=KryMin
            !Enter sub_loop
            sub_loop:   do j=1,KryMax+1
                  !From residual vector
                  if(ALLOCATED(Residual))DEALLOCATE(Residual)
                  ALLOCATE(Residual(ciDim))
                  CALL gen_ResidVec(basis_vecs,hv_vectors,ciDim,curr_dim,&
                        Kry_eigvec,Kry_eigval,curr_root,Residual)
                  !Compute norm of residual
                  resid_norm=SQRT(ddot(ciDim,Residual,1,Residual,1))
                  print *, "========================================="
                  print *, "Iteration ",(i-1)*(KryMax+1)+j
                  print *, "Solving for root # ", curr_root
                  print *, "Eigenvalue:        ",Kry_eigval(curr_root)
                  print *, "Residual norm:     ",resid_norm
                  print *, "Current dimension: ",curr_dim
                  !Test if ||r|| < res_tol 
                  if(resid_norm.LT.res_tol)then
                        Eig_Values(curr_root)=Kry_eigval(curr_root)
                        EXIT main_loop
                  end if
                  !Deallocate Kry_eigvec
                  DEALLOCATE(Kry_eigvec)
                  !Generate new vector
                  if(ALLOCATED(new_Vector))DEALLOCATE(new_Vector)
                  ALLOCATE(new_Vector(ciDim))
                  CALL gen_newVector(Residual,Diagonals,Kry_eigval(curr_root),&
                        ciDim,new_Vector)
                  !Deallocate Kry_eigval
                  DEALLOCATE(Kry_eigval)
                  !Orthogonalize this vector to the basis space
                  CALL orthog_newvec(basis_vecs,new_Vector,ciDim,curr_dim)
                  !Add new vector to the basis. *Note:  curr_dim=curr_dim+1
                  ! after subroutine call
                  CALL append_newvec(basis_vecs,curr_dim,ciDim,new_Vector)
                  !Perform Hv on new vector
                  CALL acthv(basis_vecs(1,curr_dim),MOints1,MOints2,M1Len,M2Len,    &
                        pString,pStep,pLocate,qString,qStep,qLocate,xRefList,       &
                        ciDim,pDets,pDLen,qDets,qDLen,aDets,bDets,aElec,bElec,      &
                        Orbitals,Diagonals,nFrozen,nDocc,nActive,hv_vectors(1,      &
                        curr_dim))
                  !Build vHv subspace matrix
                  ALLOCATE(subsp_ham(curr_dim,curr_dim))
                  CALL build_subham(basis_vecs,hv_vectors,ciDim,curr_dim,subsp_ham)
                  ALLOCATE(Kry_eigval(curr_dim))
                  ALLOCATE(Kry_eigvec(curr_dim,curr_dim))
                  !Diagonalize subspace matrix
                  CALL diag_dsyevr(subsp_ham,curr_dim,Kry_eigvec,Kry_eigval)
                  !Deallocate subsp_ham, it has been destroyed by dsyevr
                  DEALLOCATE(subsp_ham)
                  !Test size of subspace
                  if(curr_dim.EQ.KryMax)then
                        CALL space_trunc(basis_vecs,ciDim,KryMin,KryMax,curr_root, &
                              Kry_eigvec )
                        curr_dim=KryMin
                        !Zero out hv_vectors
                        hv_vectors=0d0
                        !Deallocate Kry_eigvec and Kry_eigval
                        DEALLOCATE(Kry_eigvec,Kry_eigval)
                        !Exit sub_loop
                        EXIT sub_loop     !Reinitialize space
                  end if
            end do sub_loop
      end do main_loop
      RETURN

END SUBROUTINE Davidson
