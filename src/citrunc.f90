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
                ndocc, nactive, xlevel )
  call enfdocc( betastrings, bdets, belec, orbitals, nfrozen, &
                ndocc, nactive, xlevel )
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
! ...ENFORCE DOCC RESTRICTIONS ON DETERMINANTS...
  remdet=0
  call enfdoccdet( fnldets, fnldetslen, aelec, adets, belec, bdets, &
                   orbitals, nfrozen, ndocc, nactive, xlevel, remdet )
! ...ENFORCE CAS RESTRICTIONS ON DETERMINANTS...
  call enfactivedet( fnldets, fnldetslen, aelec, adets, belec, bdets,&
                     orbitals, nfrozen, ndocc, nactive, xlevel, remdet )

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

  call gen_alphastrpr( fnldets, cidimension, belec, orbitals, &
                       strngpr1, pstep, plocate, modalphalen )
  call gen_betastrpr(  strngpr1, cidimension, pstep, modalphalen, &
                       plocate, strngpr2, qstep, qlocate, modbeta, modbetalen )

! ...WRITE strings to respective files...
  open(unit=5,file='pstring.list',status='new')
  do i=1,cidimension
    write(unit=5,fmt=10) strngpr1(i,1), strngpr1(i,2)
  end do
  close(unit=5)
  open(unit=7,file='pstep.list',status='new')
  do i=1,modalphalen
    write(unit=7,fmt=9) pstep(i)
  end do
  close(unit=7)
  open(unit=11,file='plocate.list',status='new')
  do i=1,modalphalen
    write(unit=11,fmt=9) plocate(i)
  end do
  close(unit=11)
  open(unit=6,file='qstring.list',status='new')
  do i=1,cidimension
    write(unit=6,fmt=10) strngpr2(i,1), strngpr2(i,2)
  end do
  close(unit=6)
  open(unit=8,file='qstep.list',status='new')
  do i=1,modbetalen
    write(unit=8,fmt=9) qstep(i)
  end do
  close(unit=8)
  open(unit=12,file='qlocate.list',status='new')
  do i=1,modbetalen
    write(unit=12,fmt=9) qlocate(i)
  end do
  close(unit=12)
10 format(1x,I10,I10)


! ...FORM DETERMINANT CROSS REFERENCE LIST...
! the 'p' list is the cannonical ordering of determinants, so this list
!  will list K'(K(p))

! Generate determinant list of determinants in stringpr1 list
!  allocate(qcrossreflist(cidimension))
!  allocate(pcrossreflist(cidimension))
!  call detcrossref( pdeterms, qdeterms, fnldets, cidimension, pcrossreflist, &
!                    qcrossreflist )

! Write cross reference list to file
!  open(unit=9,file='pcross.ref',status='new')
!  do i=1, cidimension
!    write(unit=9,fmt=9) pcrossreflist(i)
!  end do
!  close(unit=9)
!  open(unit=10,file='qcross.ref',status='new')
!  do i=1, cidimension
!    write(unit=10,fmt=9) qcrossreflist(i)
!  end do
!  close(unit=10)
! Deallocate all arrays
  deallocate(qstep,pstep,strngpr1,strngpr2,modalpha,&
             modbeta,qdeterms,pdeterms)
  return

end subroutine
