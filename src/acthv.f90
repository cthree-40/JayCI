! Subroutine to compute the action of H on some vector V
! This subroutine follows closely:
!  Theor. Chem. Acc. (2001) 106:339-351
!====================================================================
! By Christopher L. Malbon
! Yarkony Group
! Dept. of Chemistry, The Johns Hopkins University
!
! Last Edit: 11-04-13
!====================================================================
!Input:
! vector1     = input vector                      real*8  array  1-d
! moints1     = 1-e integrals                     real*8  array  1-d
! moints2     = 2-e integrals                     real*8  array  1-d
! max1e       = length of 1-e integral array      integer scalar
! max2e       = length of 2-e integral array      integer scalar
! actdets     = list of active determinants       integer array  1-d
! actdetslen  = length of active determinants     integer scalar
! aelec       = alpha electrons                   integer scalar
! belec       = beta electrons                    integer scalar
! totels      = total electrons                   integer scalar
! orbitals    = number of MO's in system          integer scalar
! adets       = alpha determinants                integer scalar
! bdets       = beta determinants                 integer scalar
! nfrzn       = frozen orbitals                   integer scalar
! ndocc       = DOCC orbitals                     integer scalar
! ncas        = CAS orbitals                      integer scalar
! dgls        = < K | H | K >                     real*8  array  1-d
!Output: 
! vector2     = output vector                     real*8  array  1-d             
!====================================================================
subroutine acthv( vector1, moints1, moints2, max1e, max2e, actdets, &
  actdetslen, aelec,belec, totels, orbitals, adets, bdets, nfrzn,   &
  ndocc, ncas, dgls,vector2 )

  use detci2
  use detci5

  implicit none

! ...input integer scalars...
  integer, intent(in) :: max1e, max2e, actdetslen, aelec, belec,      &
                         totels, orbitals, adets, bdets, nfrzn, ndocc,&
                         ncas

! ...input integer arrays...
  integer, dimension(actdetslen) :: actdets

! ...input real*8 arrays...
  real*8, dimension(actdetslen) :: vector1
  real*8, dimension(max1e) :: moints1
  real*8, dimension(max2e) :: moints2
  real*8, dimension(actdetslen) :: dgls

! ...output real*8 arrays...
  real*8, dimension(actdetslen) :: vector2

! ...loop integer scalars...
  integer :: i, j, k, l, m, n, o
  integer :: eps1, xindx1, vecindx1
  integer :: eps2, xindx2, vecindx2
  integer :: eps3, xindx3, vecindx3
  integer :: singexlen

! ...loop integer arrays...
  integer, dimension(adets, aelec ) :: alphmat
  integer, dimension(bdets, belec ) :: betamat
  integer, dimension(aelec) :: astring1, astring2
  integer, dimension(belec) :: bstring1, bstring2
  integer, dimension(orbitals-aelec) :: pexits1
  integer, dimension(orbitals-belec) :: qexits1

! This array contains information from the single excitations
!  srepinfo( det index, parity, bound orbital, replacement orbital )
  integer, dimension(:,:), allocatable :: srepinfo


! ...loop real*8 scalars...
  real*8 :: int1e1, int2e1, int2e2, int3e1weps, int3e2
!--------------------------------------------------------------------

#ifdef DACTHV
! This will create a file determinants.list that will contain
!  determinant = alpha . beta
  open( unit=10, file='determinants.list', status='unknown', position='rewind' )
  write( unit=10, fmt=103 ) " Determinant  alpha string  beta string "
103 format(1x,A)
  do i=1, actdetslen
    call k2indc( i, belec, orbitals, j, k )
    write(unit=10, fmt=102) i, j, k 
  end do
102 format(1x,I12, I14, I14)
  close(unit=10)
#endif

! ...DIAGONAL CONTRIBUTION...
  do i=1, actdetslen
    vector2(i) = vector1(i)*dgls(i)
  end do

#ifdef DACTHV
  print *, " Below is vector2 after contribution from diagonal elements..."
  do i=1, actdetslen
    write(*,fmt=101) i, vector2(i)
  end do
101 format(1x,I5,F10.6)
#endif  


!
! ...LOOP OVER ALPHA STRINGS...
#ifdef DACTHV
  print *, " Beginning loop over alpha determinants..."
#endif

  do i=1, adets

! Generate orbital index string
    if ( i .eq. 1 ) then
      do j=1, aelec
        astring1(j) = j
      end do
    else
      call nstrfnd( astring2, aelec, orbitals, adets, astring1 )
    end if
#ifdef DACTHV
    print *, " Here is our alpha string..."
    do j=1, aelec
      print *, astring1(j)
    end do
#endif

! Generate list of possible excitations
    call possex1( astring1, orbitals, aelec, (orbitals-aelec), &
                  pexits1 )

#ifdef DACTHV
    print *, " Here are our possible excitations..."
    do j=1, orbitals-aelec
      print *, pexits1(j)
    end do
#endif

! Loop over all single replacements
! First, allocate srepinfo
    singexlen = aelec*(orbitals-aelec)
    if ( allocated(srepinfo) ) deallocate(srepinfo)
    allocate( srepinfo( singexlen, 4 ) )

    do j=1, aelec
      do k=1, orbitals-aelec
        call singrepinfo( astring1, aelec, pexits1(k), j, orbitals, &
                          eps1, xindx1 )

        srepinfo((j-1)*(orbitals-aelec) + k, 1) = astring1(j)       ! Orbital exciting from
        srepinfo((j-1)*(orbitals-aelec) + k, 2) = pexits1(k)        ! Orbital exciting to
        srepinfo((j-1)*(orbitals-aelec) + k, 3) = eps1              ! Parity
        srepinfo((j-1)*(orbitals-aelec) + k, 4) = xindx1            ! Index of new string

        int1e1 = moints1(ind2val(astring1(j),pexits1(k)))
#ifdef DACTHV
        open(unit=15,file='1eintegrals',status='unknown',position='append')
        write(unit=15,fmt=107) "Integral for ", astring1(j), pexits1(k),":"
        write(unit=15,fmt=109) int1e1
107 format(1x,A,I4,I4,A) 
108 format(1x,A,A,A,I4,I4,A)
109 format(1x,F10.6)
        close(unit=15)
#endif
        int2e1 = 0d0
        do l=1,aelec
          if ( l .ne. j ) then
            int2e1 = int2e1 + moints2( index2e2( astring1(l), astring1(l),       &
                       astring1(j), pexits1(k))) - moints2( index2e2( astring1(l),&
                       astring1(j), astring1(l), pexits1(k) ))
          end if
        end do

#ifdef DACTHV
        open(unit=16,file='2eintegrals',status='unknown',position='append')
        write(unit=16,fmt=108) "Integral for ", "n","n",astring1(j),pexits1(k),":"
        write(unit=16,fmt=109) int2e1
        close(unit=16)
#endif

! Loop over all beta strings
        do l=1, bdets

! Generate orbital index string
          if ( l .eq. 1 ) then
            do m=1, belec
              bstring1(m) = m
            end do
          else
            call nstrfnd( bstring2, belec, orbitals, bdets, bstring1 )
          end if
          int2e2 = 0d0
          do n=1, belec
            int2e2 =  int2e2 + moints2( index2e2( bstring1(n), bstring1(n),   &
                        astring1(j), pexits1(k) ) )
          end do

! Find index of determinant C[r(pjv*),q]
          vecindx1 = indxk( xindx1, l, belec, orbitals )


! Note: this is constant -> v(p,q) = SUM_r,s < p, q | H | r, s > C( r, s)
          vecindx2 = indxk( i, l, belec, orbitals )

#ifdef DACTHV
          open(unit=11, file='singleconts.alpha',status='unknown',position='append')
          write(unit=11, fmt=104) " Single excitation. C[r,s]=",vecindx1," V(p,q)=",vecindx2
          write(unit=11, fmt=109) eps1*(int1e1 + int2e1 + int2e2)*vector1(vecindx1)
104 format(1x,A,I4,A,I4)
          close(unit=11)
#endif

! Add the stored integrals and multiply by parity
          vector2(vecindx2) = vector2(vecindx2) + eps1*( int1e1 + int2e1 + int2e2 )*vector1(vecindx1)

! Close loop over beta strings
          bstring2 = bstring1
          bstring1 = 0
        end do


! Loop over additional single replacements ( l > j; m > k )
#ifdef DACTHV
        print *, " We are now looping over additional single excitations..."
#endif
        do l=j+1, aelec
          do m=k+1, orbitals-aelec
            call doublerepinfo( astring1, aelec, pexits1(k), j, pexits1(m), l,&
                                orbitals, eps2, xindx2 )
            int3e1weps = eps2*( moints2( index2e2(astring1(j), pexits1(k), astring1(l), &
                           pexits1(m) )) - moints2( index2e2( astring1(j), pexits1(m),   &
                           astring1(l), pexits1(k) ) ))


! Loop over all beta strings without generating orbital sets
            do n=1, bdets
              vecindx1 = indxk( xindx2, n, belec, orbitals )
              vecindx2 = indxk(i, n, belec, orbitals )
              vector2(vecindx2) = vector2(vecindx2) + int3e1weps*vector1(vecindx1)
            end do
          end do

! Close loop over additional single replacements
        end do
      end do

! Close loop over single excitations
    end do

#ifdef DACTHV
    print *, " We are now computing 2 single excitations, these should be nonzero..."
#endif
! Loop over all beta strings. This is for the contribution of two single excitations
    do j=1, bdets
      if ( j .eq. 1 ) then
        do k=1, belec
          bstring1(k) = k
        end do
      else
        call nstrfnd( bstring2, belec, orbitals, bdets, bstring1 )
      end if

! Generate list of possible excitations
      call possex1( bstring1, orbitals, belec, (orbitals-belec), qexits1 )


! Loop over single excitations in the beta strings
      do k=1, belec
        do l=1, orbitals-belec
          call singrepinfo( bstring1, belec, qexits1(l), k, orbitals, eps3, &
                             xindx1 )


! Loop over single alpha string replacements. This information is stored in srepinfo(:,:)
          do m=1, singexlen
            int3e2 = moints2( index2e2( srepinfo( m, 1), srepinfo(m, 2), bstring1(k), &
                                        qexits1(l) ) )

            vecindx1 = indxk( srepinfo( m,4 ), xindx1, belec, orbitals )

            vecindx2 = indxk( i, j, belec, orbitals )


#ifdef DACTHV
            open( unit=20,file='dblesingle.ab',status='unknown',position='append')
            write(unit=20,fmt=111) "Integral for ", srepinfo(m,1), srepinfo(m,2), bstring1(k), qexits1(l),":"
            write(unit=20,fmt=112) int3e2
111 format(1x,A,I3,I3,I3,I3,A)
112 format(1x,F10.6)
#endif
            vector2(vecindx2) = vector2(vecindx2) + srepinfo(m,3)*eps3*int3e2*vector1(vecindx1)
          
          end do  ! Close loop over single alpha replacements

        end do    

      end do      ! Close loop over excitations
      bstring2=bstring1
      bstring1=0

    end do        ! Close loop over beta determinants
    astring2 = astring1
    astring1 = 0 
  end do          ! Close loop over alpha strings
 



! ...LOOP OVER BETA STRINGS...
! This loop follows the same structure as the loop(s) above
  do i=1, bdets

! Generate orbital index string
    if ( i .eq. 1 ) then
      do j=1, belec
        bstring1(j) = j
      end do
    else
      call nstrfnd( bstring2, belec, orbitals, bdets, bstring1 )
    end if


! Generate list of possible excitations
    call possex1( bstring1, orbitals, belec, (orbitals-belec), &
                  qexits1 )


! Loop over all single replacements
    do j=1, belec
      do k=1, orbitals-belec
        call singrepinfo( bstring1, belec, qexits1(k), j, orbitals, &
                          eps1, xindx1 )

        int1e1 = moints1( ind2val( bstring1(j), qexits1(k) ))
        int2e1 = 0d0
        do l=1, belec
          if ( l .ne. j ) then
            int2e1 = int2e1 + moints2( index2e2( bstring1(l), bstring1(l),        &
                       bstring1(j), qexits1(k))) - moints2( index2e2( bstring1(l),&
                       bstring1(j), bstring1(l), qexits1(k) ))
          end if
        end do



! Loop over all alpha strings
        do l=1, adets

! Generate orbital index string
          if ( l .eq. 1 ) then
            do m=1, aelec
              astring1(m) = m
            end do
          else
            call nstrfnd( astring2, aelec, orbitals, adets, astring1 )
          end if
          int2e2 = 0d0
          do n=1, aelec
            int2e2 = int2e2 + moints2( index2e2( astring1(n), astring1(n), &
                      bstring1(j), qexits1(k)) )
          end do


! Find index of determinant C[r(pjv*),q]
          vecindx1 = indxk( l, xindx1, belec, orbitals )

 
!    
          vecindx2 = indxk( l, i, belec, orbitals )


! Add the stored integrals and multiply by parity
          vector2(vecindx2) = vector2(vecindx2) + eps1*( int1e1 + int2e1 + int2e2 )*vector1(vecindx1)

! Close the loop over alpha strings
          astring2 = astring1
          astring1 = 0
        end do  

!
! Loop over additional single replacements ( l > j; m > k )
        do l=j+1, belec
          do m=k+1, orbitals-belec
            call doublerepinfo( bstring1, belec, qexits1(k), j, pexits1(m), l, &
                                orbitals, eps2, xindx2 )
            int3e1weps = eps2*( moints2( index2e2( bstring1(j), qexits1(k), bstring1(l), &
                           qexits1(m) )) - moints2( index2e2( bstring1(j), qexits1(m),   &
                           bstring1(l), qexits1(k) ) ))


! Loop over all alpha strings without generating orbital sets
            do n=1, adets
              vecindx1 = indxk( n, xindx2, belec, orbitals )
              vecindx2 = indxk( n, i, belec, orbitals )
              vector2(vecindx2) = vector2(vecindx2) + int3e1weps*vector1(vecindx1)
            end do
          end do

! Close loop over aditional single replacements
        end do
      end do

! Close loop over single excitations
    end do
 
! Close loop over beta strings
   bstring2=bstring1
   bstring1=0
   end do

! End subroutine
   return

        
end subroutine
