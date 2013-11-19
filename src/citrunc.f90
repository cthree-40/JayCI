! This subroutine truncates the CI expansion.
! It prints out a list of alpha and beta determinants in the expansion.
!====================================================================
! By Christopher Malbon
! Yarkony Group
! Dept. of Chem., The Johns Hopkins University
!
! Last edit: 11-13-13
!====================================================================
!Input:
! adets        = number of alpha strings            integer scalar
! bdets        = number of beta strings             integer scalar
! aelec        = number of alpha electrons          integer scalar
! belec        = number of beta electrons           integer scalar
! orbitals     = number of MO's                     integer scalar
! nfrozen      = number of frozen core orbitals     integer scalar
! ndocc        = number of double occupied orbitals integer scalar
! nactive      = number of active orbitals          integer scalar
! xlevel       = excitation level out of docc       integer scalar
! cidimension  = dimension of ci expansion          integer sclar
!====================================================================
subroutine citrunc( adets, bdets, aelec, belec, orbitals, nfrozen, &
  ndocc, nactive, xlevel, modalphalen, modbetalen, cidimension )

  use detci2
  use truncation

  implicit none

! ...input integer scalars...
  integer, intent(in) :: adets, bdets, aelec, belec, orbitals, nfrozen,&
                         ndocc, nactive, xlevel

! ...OUTPUT integer scalars...
  integer :: cidimension
  integer :: modalphalen, modbetalen

! ...loop integer scalars...
  integer :: i, j, k, l
  integer :: test, testa, testb
  integer :: p, q
  integer :: index1

! ...loop integer arrays...
  integer, dimension( aelec ) :: astring1, astring2
  integer, dimension( belec ) :: bstring1, bstring2
  integer, dimension( adets ) :: alphastrings
  integer, dimension( bdets ) :: betastrings

  integer, dimension( adets, aelec ) :: alphmat 
  integer, dimension( bdets, belec ) :: betamat
 
  integer, dimension(:), allocatable   :: modalpha, modbeta

  integer, dimension(:), allocatable   :: pstrdetlist, crssreflist

  integer, dimension(:,:), allocatable :: strngpr1, strngpr2, &
                                          refdpairs
  integer, dimension(:), allocatable :: plocate, qlocate, pstep, qstep

  integer, dimension(:), allocatable :: fnldets
  integer :: fnldetslen, remdet
  
  integer :: step
!--------------------------------------------------------------------

! Generate alpha and beta string lists
  call genstrings( alphastrings, betastrings, adets, bdets )
! ...ENFORCE FROZEN ORBITAL RESTRICTIONS...
! alpha strings
  call enffrozen( alphastrings, nfrozen, aelec, orbitals, adets )
! beta strings
  call enffrozen( betastrings, nfrozen, belec, orbitals, bdets )
! ...ENFORCE DOCC RESTRICTIONS...
  call enfdocc( alphastrings, adets, aelec, orbitals, nfrozen, &
                ndocc, xlevel )
  call enfdocc( betastrings, bdets, belec, orbitals, nfrozen, &
                ndocc, xlevel )
! ...ENFORCE CAS RESTRICTIONS
  call enfactive( alphastrings, adets, orbitals, aelec, nfrozen,&
                  ndocc, nactive, xlevel )
  call enfactive( betastrings, bdets, orbitals, belec, nfrozen, &
                  ndocc, nactive, xlevel )
  
! Write out the alpha string indices
  open( unit=2, file='alpha.dets', status='new' )
  do i=1, adets
    if ( alphastrings(i) .ne. 0 ) then
      write( unit=2, fmt=9 ) alphastrings(i)
    end if
  end do
  close(unit=2)

! Write out the beta string indices
  open( unit=3, file='beta.dets', status='new' )
  do i=1, bdets
    if ( betastrings(i) .ne. 0 ) then
      write( unit=3, fmt=9 ) betastrings(i)
    end if
  end do
  close(unit=3)

9 format(1x,I10)


! Get lengths of modalpha and modbeta
  call modspstrlen( alphastrings, adets, modalphalen )
  call modspstrlen( betastrings,  bdets, modbetalen  )

! Allocate modalpha and modbeta
  if ( allocated(modalpha)) deallocate(modalpha)
  allocate(modalpha(modalphalen))
  if ( allocated(modbeta)) deallocate(modbeta)
  allocate(modbeta(modbetalen))


! Read in index values from beta.dets and alpha.dets
! Alpha strings
  open( unit=2, file='alpha.dets', status='old', position='rewind' )
  read(unit=2, fmt=9) ( modalpha(i), i=1, modalphalen )
  close(unit=2)
! Beta strings
  open ( unit=3, file='beta.dets', status='old', position='rewind' )
  read(unit=3, fmt=9) ( modbeta(i), i=1, modbetalen )
  close(unit=3)

  fnldetslen=modalphalen*modbetalen
  allocate( fnldets(fnldetslen) )
! Generate determinant list
  call gendetlist( modalpha, modbeta, modalphalen, modbetalen, belec,&
                   orbitals, fnldets )
  print *, "FNLDETS() BEFORE ANY TRUNCATION-------------------------"
  do i=1, fnldetslen
    print *, fnldets(i)
  end do
! ...ENFORCE DOCC RESTRICTIONS ON DETERMINANTS...
  remdet=0
  call enfdoccdet( fnldets, fnldetslen, aelec, adets, belec, bdets, &
                   orbitals, nfrozen, ndocc, xlevel, remdet )
  print *, "FNLDETS() AFTER DOCC TRUCATION--------------------------"
  do i=1, fnldetslen
    print *, fnldets(i)
  end do
  print *, "DETERMINANTS REMOVED :    ", remdet
! ...ENFORCE CAS RESTRICTIONS ON DETERMINANTS...
  call enfactivedet( fnldets, fnldetslen, aelec, adets, belec, bdets,&
                     orbitals, nfrozen, ndocc, nactive, xlevel, remdet )
  print *, "FNLDETS() AFTER ACTIVE TRUNCATION------------------------"
  do i=1, fnldetslen
    print *, fnldets(i)
  end do

! Write determinant list to file
  open( unit=4, file='det.list', status='new', position='rewind' )
  do i=1, fnldetslen
    if ( fnldets(i) .ne. 0 ) then
      write(unit=4,fmt=9) fnldets(i)
    end if
  end do
  close(unit=4)
! Deallocate fnldets
  deallocate( fnldets )

! Compute number of determinants in ci expansion
  cidimension = (modalphalen*modbetalen) - remdet
! Allocate fnldets, but with adjusted size
  allocate(fnldets(cidimension))
! Read in determinant list from file
  open( unit=4, file='det.list', status='old', position='rewind' )
  read( unit=4, fmt=9 ) ( fnldets(i), i=1, cidimension )
  close(unit=4)
! ...FORM STRING PAIRS...
! String pairs (p,q)

  allocate( strngpr1(cidimension,2))
  allocate( strngpr2(cidimension,2))

  allocate( pstep(modalphalen))
! Loop over alpha strings
  l=1
  do i=1, modalphalen
    step=0
! Loop over beta strings
    do j=1, modbetalen
    
      index1 = indxk( modalpha(i), modbeta(j), belec, orbitals ) 

! Test if index1 is in expansion fnldets()
      do k=1, cidimension
        if ( index1 .eq. fnldets(k) ) then
          print *, index1
          strngpr1(l,1) = i
          strngpr1(l,2) = j
          l=l+1
          step=step+1
          exit
        end if
      end do
    end do
    pstep(i) = step
  end do

  print *, "Printing step information..."
  do i=1, modalphalen
    print *, pstep(i)
  end do
  print *, "Printing string pairings..."
  do i=1, cidimension
    print *, strngpr1(i,1), strngpr1(i,2), indxk( strngpr1(i,1),strngpr1(i,2),belec,orbitals)
  end do

  allocate(qstep(modbetalen))
! Loop over beta strings
  l=1
  do i=1, modbetalen
    step=0
! Loop over alpha strings
    do j=1, modalphalen
    
      index1 = indxk( j, i, belec, orbitals )

! Test if index1 is in expansion fnldets()
      do k=1, cidimension
        if ( index1 .eq. fnldets(k) ) then
          strngpr2(l,1) = i
          strngpr2(l,2) = j
          l=l+1
          step=step+1
          exit
        end if
      end do
    end do
    qstep(i) = step
  end do

  print *, "cidimension: ", cidimension
! ...WRITE strings to respective files...
  open(unit=5,file='pstring.list',status='new')
  do i=1,cidimension
    write(unit=5,fmt=10) strngpr1(i,1), strngpr1(i,2)
  end do
  close(unit=5)
  open(unit=6,file='qstring.list',status='new')
  do i=1,cidimension
    write(unit=6,fmt=10) strngpr2(i,1), strngpr2(i,2)
  end do
  close(unit=6)
10 format(1x,I10,I10)


! ...FORM DETERMINANT CROSS REFERENCE LIST...
! the 'p' list is the cannonical ordering of determinants, so this list
!  will list K'(K(p))

! Generate determinant list of determinants in stringpr1 list
  allocate( pstrdetlist(cidimension) )

  do i=1, cidimension
    pstrdetlist(i) = indxk( strngpr1(i,1),strngpr1(i,2), belec, orbitals ) 
  end do

  

  deallocate(modalpha)
  deallocate(modbeta)
  return

end subroutine
