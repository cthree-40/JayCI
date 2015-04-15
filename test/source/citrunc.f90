! This subroutine truncates the CI expansion.
! It prints out a list of alpha and beta determinants in the expansion.
!====================================================================
! By Christopher Malbon
! Yarkony Group
! Dept. of Chem., The Johns Hopkins University
!
! Added Mapping arrays: 04-30-2014
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
  ndocc, nactive, xlevel, modalphalen, modbetalen, cidimension, unit1, &
  unit2, unit3, unit4, unit5, unit6, unit7, unit8, unit9, unit10 )!,      &
  !unit11, unit12 )

  use addressing
  use truncation

  implicit none

! ...input integer scalars...
  integer, intent(in) :: adets, bdets, aelec, belec, orbitals, nfrozen,&
                         ndocc, nactive, xlevel, unit1, unit2, unit3,  &
                         unit4, unit5, unit6, unit7, unit8, unit9,     &
                         unit10!, unit11, unit12

! ...OUTPUT integer scalars...
  integer :: cidimension
  integer :: modalphalen, modbetalen

! ...loop integer scalars...
  integer :: i, j, k, l, m
  integer :: test, testa, testb
  integer :: p, q
  integer :: index1

! ...loop integer arrays...
  integer, dimension( aelec ) :: astring1, astring2
  integer, dimension( belec ) :: bstring1, bstring2
  integer, dimension( adets ) :: alphastrings
  integer, dimension( bdets ) :: betastrings

!  integer, dimension( adets, aelec ) :: alphmat 
!  integer, dimension( bdets, belec ) :: betamat
 
  integer, dimension(:), allocatable   :: modalpha, modbeta

  integer, dimension(:), allocatable   :: pstrdetlist, crssreflist

  integer, dimension(:,:), allocatable :: strngpr1, strngpr2, &
                                          refdpairs
  integer, dimension(:), allocatable :: plocate, qlocate, pstep, qstep

  integer, dimension(:), allocatable :: fnldets
  
  integer, dimension(:), allocatable :: qcrossreflist, pcrossreflist
  
  integer, dimension(:), allocatable :: qdeterms, pdeterms

  integer :: fnldetslen, remdet
  
  integer :: step

!  integer, dimension(:), allocatable :: alphaMap, betaMap ! Mapping arrays
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
                ndocc, nactive, xlevel, 'A' )
  call enfdocc( betastrings, bdets, belec, orbitals, nfrozen, &
                ndocc, nactive, xlevel, 'B' )
! ...ENFORCE CAS RESTRICTIONS
  call enfactive( alphastrings, adets, orbitals, aelec, nfrozen,&
                  ndocc, nactive, xlevel, 'A' )
  call enfactive( betastrings, bdets, orbitals, belec, nfrozen, &
                  ndocc, nactive, xlevel, 'B' )

! Write out the alpha string indices
! Allocate mapping array
!  if (allocated(alphaMap)) deallocate(alphaMap)
!  allocate(alphaMap(adets))
!  alphaMap=0      ! Zero out array
!  m=1             ! Indexing integer for alphaMap 
  open( unit=unit3, file='alpha.dets', status='new' )
  do i=1, adets
    if ( alphastrings(i) .ne. 0 ) then
      write( unit=unit3, fmt=9 ) alphastrings(i)
!      alphaMap(i) = m   ! Location in alpha strings of ith array
!      m = m + 1
    end if
  end do
  close(unit=unit3)

! Write out mapping array
!  open( unit=unit11, file='alpha.Map', status='new' )
!  do i=1, adets
!    write( unit=unit11, fmt=9 ) alphaMap(i)
!  end do
!  close(unit=unit11)
!  deallocate(alphaMap)

! Write out the beta string indices
! Allocate mapping array
 ! if (allocated(betaMap)) deallocate(betaMap)
!  allocate(betaMap(bdets))
!  betaMap=0      ! Zero out array
!  m=1             ! Indexing integer for alphaMap 
  open( unit=unit2, file='beta.dets', status='new' )
  do i=1, bdets
    if ( betastrings(i) .ne. 0 ) then
      write( unit=unit2, fmt=9 ) betastrings(i)
!      betaMap(i) = m
!      m = m + 1
    end if
  end do
  close(unit=unit2)

! Write out mapping array
!  open( unit=unit12, file='beta.Map', status='new' )
!  do i=1, bdets
!    write( unit=unit12, fmt=9 ) betaMap(i)
!  end do
!  close(unit=unit12)
!  deallocate(betaMap)

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
  open( unit=unit3, file='alpha.dets', status='old', position='rewind' )
  read(unit=unit3, fmt=9) ( modalpha(i), i=1, modalphalen )
  close(unit=unit3)
! Beta strings
  open ( unit=unit2, file='beta.dets', status='old', position='rewind' )
  read(unit=unit2, fmt=9) ( modbeta(i), i=1, modbetalen )
  close(unit=unit2)

  fnldetslen=modalphalen*modbetalen
  allocate( fnldets(fnldetslen) )
! Generate determinant list
  call gendetlist( modalpha, modbeta, modalphalen, modbetalen, belec,&
                   orbitals, fnldets )
! ...ENFORCE DOCC RESTRICTIONS ON DETERMINANTS...
  remdet=0
  call enfdoccdet( fnldets, fnldetslen, aelec, adets, belec, bdets, &
                   orbitals, nfrozen, ndocc, nactive, xlevel, remdet )
! ...ENFORCE CAS RESTRICTIONS ON DETERMINANTS...
  call enfactivedet( fnldets, fnldetslen, aelec, adets, belec, bdets,&
                     orbitals, nfrozen, ndocc, nactive, xlevel, remdet )
! Write determinant list to file
  open( unit=unit1, file='det.list', status='new', position='rewind' )
  do i=1, fnldetslen
    if ( fnldets(i) .ne. 0 ) then
      write(unit=unit1,fmt=9) fnldets(i)
    end if
  end do
  close(unit=unit1)
! Deallocate fnldets
  deallocate( fnldets )

! Compute number of determinants in ci expansion
  cidimension = (modalphalen*modbetalen) - remdet
! Allocate fnldets, but with adjusted size
  allocate(fnldets(cidimension))
! Read in determinant list from file
  open( unit=unit1, file='det.list', status='old', position='rewind' )
  read( unit=unit1, fmt=9 ) ( fnldets(i), i=1, cidimension )
  close(unit=unit1)
! ...FORM STRING PAIRS...
! String pairs (p,q)

  allocate( strngpr1(cidimension,2))
  allocate( strngpr2(cidimension,2))

  allocate( pstep(modalphalen))
  allocate( qstep(modbetalen))
  allocate( plocate(modalphalen))
  allocate( qlocate(modbetalen))
  allocate( pdeterms(cidimension))
  allocate( qdeterms(cidimension))
! Generate alpha string pairs
!  call genstrpairs( modalpha, modalphalen, modbeta, modbetalen, &
!                    fnldets, cidimension, belec, orbitals, 'a', &
!                    strngpr1, cidimension, pstep, modalphalen,  &
!                    pdeterms, plocate )
!  call genstrpairs( modbeta, modbetalen, modalpha, modalphalen, &
!                    fnldets, cidimension, belec, orbitals, 'b', &
!                    strngpr2, cidimension, qstep, modbetalen,   &
!                    qdeterms, qlocate )
  print *, "Final dets=",cidimension
  call gen_alphastrpr( fnldets, cidimension, belec, orbitals, &
                       strngpr1, pstep, plocate, modalphalen )
  call gen_betastrpr(  strngpr1, cidimension, pstep, modalphalen, &
                       plocate, strngpr2, qstep, qlocate, modbeta, modbetalen )

! ...WRITE strings to respective files...
  open(unit=unit4,file='pstring.list',status='new')
  do i=1,cidimension
    write(unit=unit4,fmt=10) strngpr1(i,1), strngpr1(i,2)
  end do
  close(unit=unit4)
  open(unit=unit8,file='pstep.list',status='new')
  do i=1,modalphalen
    write(unit=unit8,fmt=9) pstep(i)
  end do
  close(unit=unit8)
  open(unit=unit9,file='plocate.list',status='new')
  do i=1,modalphalen
    write(unit=unit9,fmt=9) plocate(i)
  end do
  close(unit=unit9)
  open(unit=unit5,file='qstring.list',status='new')
  do i=1,cidimension
    write(unit=unit5,fmt=10) strngpr2(i,1), strngpr2(i,2)
  end do
  close(unit=unit5)
  open(unit=unit6,file='qstep.list',status='new')
  do i=1,modbetalen
    write(unit=unit6,fmt=9) qstep(i)
  end do
  close(unit=unit6)
  open(unit=unit7,file='qlocate.list',status='new')
  do i=1,modbetalen
    write(unit=unit7,fmt=9) qlocate(i)
  end do
  close(unit=unit7)
10 format(1x,I10,I10)


! ...FORM DETERMINANT CROSS REFERENCE LIST...
! the 'p' list is the cannonical ordering of determinants, so this list
!  will list K'(K(p))

! Generate determinant list of determinants in stringpr1 list
  allocate(qcrossreflist(cidimension))
!  allocate(pcrossreflist(cidimension))
!  call detcrossref( pdeterms, qdeterms, fnldets, cidimension, pcrossreflist, &
!                    qcrossreflist )
  call xreflist_gen( cidimension, fnldets, strngpr2, belec, orbitals, qcrossreflist )

! Write cross reference list to file
!  open(unit=9,file='pcross.ref',status='new')
!  do i=1, cidimension
!    write(unit=9,fmt=9) pcrossreflist(i)
!  end do
!  close(unit=9)
  open(unit=10,file='xref.list',status='new')
  do i=1, cidimension
    write(unit=10,fmt=9) qcrossreflist(i)
  end do
  close(unit=10)
! Deallocate all arrays
  deallocate(qstep,pstep,strngpr1,strngpr2,modalpha,&
             modbeta,qdeterms,pdeterms,qcrossreflist)
  return

end subroutine
