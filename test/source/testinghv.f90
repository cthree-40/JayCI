module testinghv
  implicit none
contains
  subroutine buildham( ciDim, MOints1, M1Len, MOints2, M2Len, &
       Det_Exp, aElec, bElec, Orbitals, hamiltonian )
    use construct,  only: ham_element
    use david_util, only: diag_dsyevr
    implicit none
    integer, intent(in) :: ciDim, M1Len, M2Len, aElec, bEle,  &
         Orbitals
    real*8,  dimension(M1Len), intent(in) :: MOints1
    real*8,  dimension(M2Len), intent(in) :: MOints2
    real*8,  dimension(ciDim, ciDim), intent(in out) :: hamiltonian
    integer, dimension(ciDim), intent(in) :: Det_Exp
    integer :: i, j
    ! Build hamiltonian
    hamiltonian = 0d0
    do i=1, ciDim
       do j=1, ciDim
          hamiltonian(j,i) = ham_element( Det_Exp(j), Det_Exp(i), &
               MOints1, M1Len, MOints2, M2Len, aElec, bElec,      &
               Orbitals )
       end do
    end do
    return
  end subroutine buildham
!==================================================================
  subroutine buildhamcol( ciDim, MOints1, M1Len, MOints2, M2Len, &          
       Det_Exp, aElec, bElec, Orbitals, column, hamcol )                    
    use construct,  only: ham_element                                    
    use david_util, only: diag_dsyevr                                    
    implicit none                                                        
    integer, intent(in) :: ciDim, M1Len, M2Len, aElec, bEle,  &          
         Orbitals, column                                                        
    real*8,  dimension(M1Len), intent(in) :: MOints1                     
    real*8,  dimension(M2Len), intent(in) :: MOints2                     
    real*8,  dimension(ciDim), intent(in out) :: hamcol      
    integer, dimension(ciDim), intent(in) :: Det_Exp                     
    integer :: i, j                                                      
    ! Build hamiltonian                                                  
    hamcol = 0d0                                                                     
    do j=1, ciDim                                                     
       hamcol(j) = ham_element( Det_Exp(j), Det_Exp(column), &      
            MOints1, M1Len, MOints2, M2Len, aElec, bElec,      &      
            Orbitals )                                                
    end do                                                               
    return                                                               
  end subroutine buildhamcol
!=====================================================================
  subroutine diagham( dim, hamiltonian, eigval )
    use david_util, only: diag_dsyevr
    implicit none
    integer, intent(in)  :: dim
    real*8,  dimension(dim,dim), intent(in out) :: hamiltonian
    real*8,  intent(out) :: eigval
    real*8,  dimension(:,:), allocatable :: eigvecs
    real*8,  dimension(:),   allocatable :: eigvals
    ! diagonalize hamiltonian
    allocate( eigvecs(dim,dim) )
    allocate( eigvals(dim) )
    call diag_dsyevr( hamiltonian, dim, eigvecs, eigvals )
    eigval = eigvals(1)
    deallocate(eigvecs,eigvals)
    return
  end subroutine diagham
!=====================================================================
  subroutine testhvcol( dim, ntest, hcol )
    use action_util
    implicit none
    integer, intent(in) :: dim, ntest
    real*8, dimension(dim), intent(in) :: hcol
    real*8, dimension(:), allocatable :: column, out
    integer :: i
    ! Build unit vector
    allocate( column( dim ) )
    column = 0d0
    column(ntest) = 1d0
    ! perform hv
    call acthv( column, MOints1, MOints2




