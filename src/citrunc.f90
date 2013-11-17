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

  integer, dimension(:,:), allocatable :: strngpr1, strngpr2, &
                                          refdpairs
  integer, dimension(:), allocatable :: plocate, qlocate, pstep, qstep

  integer, dimension(adets*bdets) :: fnldets
  integer :: fnldetslen, remdet
  
  integer :: step
!--------------------------------------------------------------------

! Generate alpha and beta string lists
  do i=1, adets
    alphastrings(i) = i
  end do
  do i=1, bdets
    betastrings(i) = i
  end do


! ...ENFORCE FROZEN ORBITAL RESTRICTIONS...

! Loop through alpha string list removing all determinants without an
!  alpha electron in the frozen core orbitals

  lpa: do i=1, adets

    if ( i .eq. 1 ) then
      do j=1, aelec
        astring1(j) = j
      end do
    else
      call nstrfnd( astring2, aelec, orbitals, adets, astring1 )
    end if
    
! Test if the first nfrozen orbitals are occupied
    test = 0
    lpb: do j= 1, nfrozen
      if ( astring1(j) .eq. j ) then
        cycle lpb
      else
        test = test + 1
      end if
      if ( test .ne. 0 ) then
        alphastrings(i) = 0
        exit lpb
      end if
    end do lpb
    astring2 = astring1
    astring1 = 0
  end do lpa
        

! Loop through beta string list removing all determinants without a
!  beta electron in the forzen core orbitals

  lpc: do i=1, bdets
    
    if (i .eq. 1) then
      do j=1, belec
        bstring1(j) = j
      end do
    else
      call nstrfnd( bstring2, belec, orbitals, bdets, bstring1 )
    end if

! Test if the first nfrozen orbitals are occupied
    test = 0
    lpd: do j=1, nfrozen
      if ( bstring1(j) .eq. j ) then
        cycle lpd
      else
        test = test + 1
      end if
      if ( test .ne. 0 ) then
        betastrings(i) = 0
        exit lpd
      end if
    end do lpd
    bstring2 = bstring1
    bstring1 = 0
  end do lpc



 
! ...ENFORCE DOCC RESTRICTIONS...

! Loop through alpha strings
  lpe: do i=1, adets
! Test if string is in expansion
    if ( alphastrings(i) .eq. 0 ) then
      cycle lpe
    end if
   
    if ( alphastrings(i) .eq. 1 ) then
      do j=1, aelec
        astring1(j) = j
      end do
    else
      call nstrfnd( astring2, aelec, orbitals, adets, astring1 )
    end if

    testa=0
! Count excitations in astring1
    do j=1, aelec
      if ( astring1(j) > nfrozen+ndocc ) then
        testa = testa + 1
      end if
    end do
    
    if ( testa > xlevel ) then
      alphastrings(i) = 0         ! String is not in expansion
    end if

! Close loop over alpha strings
    astring2 = astring1
    astring1 = 0
  end do lpe

! Loop through beta strings
  lpf: do i=1, bdets
! Test if string is in expansion
    if ( betastrings(i) .eq. 0 ) then
      cycle lpf
    end if

    if ( betastrings(i) .eq. 1 ) then
      do j=1, belec
        bstring1(j) = j
      end do
    else
      call nstrfnd( bstring2, belec, orbitals, bdets, bstring1 )
    end if


    testb=0
! Count excitations in bstring1
    do j=1, belec
      if ( bstring1(j) > nfrozen+ndocc ) then
        testb = testb + 1
      end if
    end do

    if ( testb > xlevel ) then
      betastrings(i) = 0             ! String is not in expansion
    end if

! Close loop over betastrings
    bstring2 = bstring1
    bstring1 = 0
  end do lpf


! ...ENFORCE CAS RESTRICTIONS
! Loop over alpha determinants
  lpg: do i=1, adets
! Test if string is in expansion
    if ( alphastrings(i) .eq. 0 ) then
      cycle lpg
    end if

    if ( alphastrings(i) .eq. 1 ) then
      do j=1, aelec
        astring1(j) = j
      end do
    else
      call nstrfnd( astring2, aelec, orbitals, aelec, astring1 )
    end if


    testa=0
! Count excitations in astring1
    do j=1, aelec
      if ( astring1(j) > nfrozen+ndocc+nactive ) then
        testa = testa + 1
      end if
    end do

    if ( testa > xlevel ) then
      alphastrings(i) = 0
    end if

! Close loop over alpha strings
    astring2 = astring1
    astring1 = 0
  end do lpg

! Loop over beta determinants
  lph: do i=1, bdets
! Test if string is in expansion
    if ( betastrings(i) .eq. 0 ) then
      cycle lph
    end if

    if ( betastrings(i) .eq. 1 ) then
      do j=1, belec
        bstring1(j) = j
      end do
    else
      call nstrfnd( bstring2, belec, orbitals, belec, bstring1 )
    end if


    testb=0
! Count excitations in bstring1
    do j=1, belec
      if ( bstring1(j) > nfrozen+ndocc+nactive ) then
        testb = testb + 1
      end if
    end do

    if ( testb > xlevel ) then
      betastrings(i) = 0
    end if

! Close loop over beta strings
    bstring2 = bstring1
    bstring1 = 0
  end do lph

  
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
! Alpha strings
  modalphalen = 0
  do i=1, adets
    if ( alphastrings(i) .ne. 0 ) then
      modalphalen = modalphalen + 1
    end if
  end do
  
  modbetalen = 0
  do i=1, bdets
    if ( betastrings(i) .ne. 0 ) then
      modbetalen = modbetalen + 1
    end if
  end do

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


  fnldets=0
  fnldetslen = modalphalen*modbetalen
! Generate determinant list
  do i=1, modalphalen
    do j=1, modbetalen
      k = (i-1)*modbetalen + j
      fnldets(k) = indxk(modalpha(i),modbeta(j),belec,orbitals)
    end do
  end do
  remdet=0
! ...ENFORCE DOCC RESTRICTIONS ON DETERMINANTS...
  call strfnd( adets, aelec, orbitals, adets, alphmat )
  call strfnd( bdets, belec, orbitals, adets, betamat )

! Loop over determinants
  do i=1, fnldetslen      
    call k2indc( fnldets(i), belec, orbitals, p, q )
! Find strings
    do j=1, aelec
      astring1(j) = alphmat(p, j)
    end do
    do j=1, belec
      bstring1(j) = betamat(q, j)
    end do

    test = 0
! Check alpha string    
    do j=1,aelec
      if ( astring1(j) > nfrozen+ndocc ) then
        test = test + 1
      end if
    end do
! Check beta string
    do j=1,belec
      if ( bstring1(j) > nfrozen+ndocc ) then
        test = test + 1
      end if
    end do

! If test > xlevel, throw determinant away
    if ( test > xlevel ) then
      fnldets(i) = 0
      remdet = remdet + 1
    end if
  end do

! ...ENFORCE CAS RESTRICTIONS ON DETERMINANTS...
! Loop over determinants
  do i=1, fnldetslen
! Test if determinant is in expansion
    if ( fnldets(i) .eq. 0 ) then
      cycle
    end if
    call k2indc( fnldets(i), belec, orbitals, p, q )
! Find strings
    do j=1, aelec
      astring1(j) = alphmat(p,j)
    end do
    do j=1, belec
      bstring1(j) = betamat(p,j)
    end do

    test = 0
! Check alpha string
    do j=1, aelec
      if ( astring1(j) > nfrozen+ndocc+nactive ) then
        test = test + 1
      end if
    end do
! Check beta string
    do j=1, belec
      if ( bstring1(j) > nfrozen+ndocc+nactive ) then
        test = test + 1
      end if
    end do

! If test > xlevel, throw determinant away
    if ( test > xlevel ) then
      fnldets(i) = 0
      remdet = remdet + 1 
    end if
  end do


! Write determinant list to file
  open( unit=4, file='det.list', status='new', position='rewind' )
  do i=1, fnldetslen
    if ( fnldets(i) .ne. 0 ) then
      write(unit=4,fmt=9) fnldets(i)
    end if
  end do
  close(unit=4)


! Compute number of determinants in ci expansion
  cidimension = (modalphalen*modbetalen) - remdet
   

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
    
      index1 = indxk( i, j, belec, orbitals ) 

! Test if index1 is in expansion fnldets()
      do k=1, cidimension
        if ( index1 .eq. fnldets(k) ) then
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
  do i=1, modalphalen
    print *, strngpr1(i,1), strngpr1(i,2)
  end do


! ...Form determinant-cross reference list...


  deallocate(modalpha)
  deallocate(modbeta)
  return

end subroutine
