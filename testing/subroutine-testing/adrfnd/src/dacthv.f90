! This subroutine debugs acthv. Specifically it tests for:
!   - hermiciticity: Perform Hv on vectors vi and vj, where vi_a = 0
!                    for all a except a=i and vj_a = 0 for all a 
!                    except a=j.
!                    Prints out: elements, determinant indices, etc.
!====================================================================
subroutine dacthv( moints1, moints2, max1e, max2e, actstr, actstrlen, &
  aelec, belec, totels, orbitals, adets, bdets, dgls, nfrzn, ndocc,   &
  ncas )

  use detci2, only: k2indc

  implicit none

! ...input integer scalars...
  integer, intent(in) :: max1e, max2e, actstrlen, aelec, belec, totels,&
                         orbitals, adets, bdets, nfrzn, ndocc, ncas

! ...input real*8 arrays...
  real*8, dimension(max1e) :: moints1
  real*8, dimension(max2e) :: moints2
  real*8, dimension(actstrlen) :: dgls

! ...input integer arrays...
  integer, dimension(actstrlen) :: actstr

! ...loop integer scalars...
  integer :: i, j, k, l
  integer :: elm1, elm2
  integer :: p, q

! ...loop real*8 scalars...
  real*8 ::  rm1, rm2

! ...loop real*8 arrays...
  real*8, dimension(actstrlen, 2) :: rdmvectors, hrdmvectors
  real*8, dimension(actstrlen) :: hartree, honhartree
!====================================================================
  print *, "  Debugging acthv. To do this we generate a singular vector."
  print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  hartree=0d0
  hartree(1) = 1d0
  call acthv( hartree, moints1, moints2, max1e, max2e, actstr, actstrlen,&
              aelec, belec, totels, orbitals, adets, bdets, nfrzn, &
              ndocc, ncas, dgls, honhartree )

  print *, " Debugging h on the hartree fock GS "
  print *, " Vector..."
  do i=1, actstrlen
    print *, " Element ", i, honhartree(i)
  end do
  print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print *, "   Debugging acthv. To do this we will generate two random"
  print *, " integers between 1 and $actstrlen. We then perform H on vi"
  print *, " and check Hermicity."
  print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

! Generate vectors.
  call random_seed
  call random_number(rm1)
  call random_number(rm2)
  elm1 = 1 + int( (actstrlen - 1)*rm1 )
  elm2 = 1 + int( (actstrlen - 1)*rm2 )
  rdmvectors=0d0
  rdmvectors(elm1,1)=1d0
  rdmvectors(elm2,2)=1d0
  hrdmvectors=0d0

! Perform Hv.
  print *, " Calling acthv..."
  do i=1, 2
    call acthv( rdmvectors(1,i), moints1, moints2, max1e, max2e, actstr,&
                actstrlen, aelec, belec, totels, orbitals, adets, bdets,&
                nfrzn, ndocc, ncas,dgls, hrdmvectors(1,i))
  end do

! Check hermiticity.
  print *, " Checking hermiticity..."
  print *, " Vector 1's ", elm1,"th element is nonzero."
  print *, "   - This corresponds to determinant ", actstr(elm1)
  print *, "   - Determinant ", actstr(elm1)," is composed of..."
! Find index of p and q strings of determinant
  call k2indc( actstr(elm1), belec, orbitals, p, q )
  print *, "      Alpha string, ", p, " and Beta string, ", q
  print *, " Vector2's ", elm2,"th element is nonzero."
  print *, "   - This corresponds to determinant ", actstr(elm2)
  print *, "   - Determinant ", actstr(elm2)," is composed of..."
! Find index of p and q strings of determinant
  call k2indc( actstr(elm2), belec, orbitals, p, q )
  print *, "      Alpha string, ", p, " and Beta string, ", q
  print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print *, hrdmvectors(elm2, 1), " SHOULD EQUAL ", hrdmvectors(elm1, 2)
  print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"


! Now check double excitations
  rdmvectors=0d0
  hrdmvectors=0d0
  call random_number(rm1)
  call random_number(rm2)
  elm1 = 500 + int( (actstrlen - 500)*rm1 )
  elm2 = 500 + int( (actstrlen - 500)*rm2 )
  
  
! Perform Hv.
  print *, " Calling acthv..."
  do i=1, 2
    call acthv( rdmvectors(1,i), moints1, moints2, max1e, max2e, actstr,&
                actstrlen, aelec, belec, totels, orbitals, adets, bdets,&
                nfrzn, ndocc, ncas, dgls, hrdmvectors(1,i))
  end do

! Check hermiticity.
  print *, " Checking hermiticity..."
  print *, " Vector 1's ", elm1,"th element is nonzero."
  print *, "   - This corresponds to determinant ", actstr(elm1)
  print *, "   - Determinant ", actstr(elm1)," is composed of..."
! Find index of p and q strings of determinant
  call k2indc( actstr(elm1), belec, orbitals, p, q )
  print *, "      Alpha string, ", p, " and Beta string, ", q
  print *, " Vector2's ", elm2,"th element is nonzero."
  print *, "   - This corresponds to determinant ", actstr(elm2)
  print *, "   - Determinant ", actstr(elm2)," is composed of..."
! Find index of p and q strings of determinant
  call k2indc( actstr(elm2), belec, orbitals, p, q )
  print *, "      Alpha string, ", p, " and Beta string, ", q
  print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print *, hrdmvectors(elm2, 1), " SHOULD EQUAL ", hrdmvectors(elm1, 2)
  print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

  return

end subroutine
