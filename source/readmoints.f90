! readmoints
! ----------
! Read 1 and 2 electron integrals in SIFS structure
! Please consult SIFS documentation for further details.
! structure:
!  read header 1
!  read header 2
!  read 1-e integrals
!  read 2-e integrals
! Input:
!  type1n   = type of integrals to read
!  orbitals = number of orbitals in system
!  m1len    = length of moints1
!  m2len    = length of moints2
! Output:
!  moints1  = 1-e integrals
!  moints2  = 2-e integrals
!  en(1) =  nuc_rep  = nuclear repulsion energy
!  en(2) =  fcenergy = frozen core energy (array, size = 10)
!--------------------------------------------------------------------
subroutine readmoints (moints1, moints2, type1n, orbitals, &
  m1len, m2len, en)
  implicit none
  
  ! **input
  integer,      intent(in) :: m1len, m2len, type1n, orbitals
  
  ! **output
  real*8, dimension(m1len), intent(out) :: moints1
  real*8, dimension(m2len), intent(out) :: moints2
  real*8, dimension(2),     intent(out) :: en

  ! ** molecular integral file
  character*20 :: mofile
  
  ! **header 1 variables for SIFS file
  ! version = sifs version number
  ! ntitle  = number of character*80 titles
  ! nsym    = number of symmetry blocks of integrals
  ! nbft    = total number of basis functions
  ! ninfo   = number of info values to be read in
  ! nenrgy  = number of energy values to be read
  ! nmap    = number of integer mapping vectors too be read
  integer :: ntitle, nsym, nbft, ninfo, nenrgy, nmap
  integer*4 :: ntitle4, nsym4, nbft4, ninfo4, nenrgy4, nmap4
  integer :: version
  
  ! **header 2 variables for SIFS file
  ! title(1:ntitle)  = identifying titles
  ! nbpsy(1:nsym)    = number of basis functions in each symmetry block
  ! slabel(1:nsym)   = character labels for the symmetry blocks
  ! info(1:ninfo)    = additional parameters describing the integral file
  !                    info(1) = fsplit: allows separation of 1 and 2 int rec.s
  !                                     1 = same file
  !                                     2 = 2-e records are on aoints2
  !                    info(2) = l1rec: 1-e integral record length in working
  !                                     precison
  !                    info(3) = n1max: max number of 1-e integrals in a record
  !                    info(4) = l2rec: 2-e integral record length
  !                    info(5) = n2max: max number of 2-e integrals in a record
  ! bfnlab(1:nbft)   = chracter labels for the individual basis functions.
  ! ietype(1:nenrgy) = energy types
  ! energy(1:nenrgy) = core energy contributions accordint to the corresponding
  !                   ietype(*) entries
  ! imtype(1:nmap)   = map vector types
  ! map(1:bft,1:nmap)= integer mapping vectors
  character*80, dimension(:),   allocatable :: title
  character*4,  dimension(:),   allocatable :: slabel
  character*8,  dimension(:),   allocatable :: bfnlab
  integer,      dimension(:),   allocatable :: nbpsy, info, ietype, imtype
  integer,      dimension(:,:), allocatable :: map
  integer,      dimension(:),   allocatable :: mapin
  real*8,       dimension(:),   allocatable :: energy
  
  !**dword variables
  integer :: num, itypea, itypeb, last, ifmt, ibvtyp, lab1
  !**record variables
  real*8,  dimension(:),   allocatable :: vals, ibitv, buf
  integer, dimension(:),   allocatable :: symb
  integer, dimension(:,:), allocatable :: labels
  real*8,  dimension(:,:), allocatable :: s1h1
  real*8 :: fcore, score, hcore
  integer, parameter :: msame = 0, nmsame = 1, nomore = 2
  integer, parameter :: iretbv =-1
  
  integer :: nipv ! nipv = 2: 1-e integrals, nipv = 4: 2-e integrals
  integer :: ios, ierr, i, index, j, l1, l2, l3, l4
  integer :: aoints, n1int, n2int
  integer, parameter :: verin = 2

  print *, "Reading molecular integrals."
  mofile = "moints"
  
  aoints = 16
  ! open integral file
  call trnfln (1, mofile)
  open (unit = aoints, file = mofile, form = "unformatted", status = "old",&
    iostat = ios)
  if (ios .ne. 0) stop "*** Could not open integral file! ***"
  
  ! read header 1 information
  read (unit = aoints, iostat = ierr) version, ntitle, nsym, nbft, ninfo, &
    nenrgy, nmap
  if ((version - verin) .ne. 0) then
    write (*,"(1x,A)") "WARNING: MO's generated with version != 2!"
  end if
  if (ierr .ne. 0) stop "*** Error reading header 1! ***"
  
  n1int = Ind2Val(nbft, nbft)
  n2int = Index2E(nbft, nbft, nbft, nbft)
  
  ! allocate arrays
  allocate (title(ntitle), stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating title()! ***"
  allocate (nbpsy(nsym),   stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating nbpsy()! ***"
  allocate (slabel(nsym),  stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating slabel()! ***"
  allocate (info(ninfo),   stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating info()! ***"
  allocate (bfnlab(nbft),  stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating bfnlab()! ***"
  allocate (ietype(nenrgy),stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating ietype()! ***"
  allocate (energy(nenrgy),stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating energy()! ***"
  allocate (imtype(nmap),  stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating imtype()! ***"
  allocate (map(nbft,nmap),stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating map()! ***"
  ! read header 2 information
  call sifreadh2 (aoints, ntitle, nsym, nbft, ninfo, nenrgy, nmap, title, &
    nbpsy, slabel, info, bfnlab, ietype, energy, imtype, map, ierr )
  if (ierr .ne. 0) stop "*** Error reading header 2! ***"
  
  ! set nuclear repulsion energy
  en(1) = energy(1)
  
  ! allocate arrays
  allocate (buf(info(2)), stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating buf()!***"
  allocate (labels(2,info(3)), stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating labels()!***"
  allocate (vals(info(3)), stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating vals()!***"
  allocate (ibitv(info(3) + 63), stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating ibitv()!***"
  allocate (s1h1(n1int, 2), stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating s1h1()!***"
  allocate (symb(nbft), stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating symb()! ***"
  allocate (mapin(n1int), stat = ierr)
  if (ierr .ne. 0) stop "*** Error allocating mapin()!***"
  do i = 1, n1int
    mapin(i) = i
  end do
  
  ! read 1 electron record
  moints1 = 0d0
  nipv = 2
  en(2) = 0d0
  
  if (type1n .eq. 1) then ! we are after H1(*)
    
    call sifrsh(aoints, info, buf, vals, labels, nsym, &
      nbpsy, mapin, n1int, s1h1, score, hcore, symb, ierr)
    if (ierr .ne. 0) stop "*** Error reading 1-e integrals!***"
    
    moints1(1:n1int) = s1h1(1:n1int,2)
    en(2) = hcore
    
  else if (type1n .eq. 3) then ! we are after the dipole moments
    
    do while (last .ne. nomore)
      call sifrd1(aoints, info, nipv, iretbv, buf, num, last, &
        itypea, itypeb, ibvtyp, vals, labels, fcore, ibitv, ierr)
      if (ierr .ne. 0) stop "***Error reading integral record!***"
      
      if (itypea .eq. 1 .and. itypeb .eq. 0 &
        .and. last .eq. 1) then
        ! <x> dipole moment
        ! collect values
        do i = 1, num
          l1 = labels(1,i)
          l2 = labels(2,i)
          index = Ind2Val(l1,l2)
          moints1(index) = moints1(index) + vals(i)
        end do
      else if (itypea .eq. 1 .and. itypeb .eq. 1 &
        .and. last .eq. 1) then
        ! <y> dipole moment
        ! collect values
        do i = 1, num
          l1 = labels(1,i)
          l2 = labels(2,i)
          index = Ind2Val(l1,l2)
          index = index + m1len
          moints1(index) = moints1(index) + vals(i)
        end do
      else if (itypea .eq. 1 .and. itypeb .eq. 2 &
        .and. last .eq. 1) then
        ! <z> dipole moment
        do i = 1, num
          l1 = labels(1,i)
          l2 = labels(2,i)
          index = Ind2Val(l1,l2)
          index = index + m1len
          index = index + m1len
          moints1(index) = moints1(index) + vals(i)
        end do
      end if
      
    end do
  end if
  
  ! deallocate arrays
  deallocate (buf, s1h1, mapin, labels, vals, ibitv)
  
  if (type1n .eq. 3) then
    close(aoints)
    return
  end if
  
  ! allocate arrays
  allocate (buf(info(4)))
  allocate (labels(4,info(5)))
  allocate (vals(info(5)))
  allocate (ibitv(info(5) + 63))
  ! read 2 electron records
  call trnfln (1, mofile)
  call sifo2f (aoints, aoints, mofile, info, i, ierr) ! i is dummy variable
  if (ierr .ne. 0) stop "*** Error positioning moints2 file! ***"
  moints2 = 0d0
  last = msame
  nipv = 4
  do while (last .ne. nomore)
    call sifreadd2 (aoints, info, ninfo, nipv, iretbv, buf, info(4), &
      num, last, itypea, itypeb, ibvtyp, vals, info(5), labels, info(5),&
      ibitv, (info(5) + 63), ierr)
    if (ierr .ne. 0) stop "*** Error reading integral record! ***"
    
    do i = 1, num
      l1 = int(labels(1,i))
      l2 = int(labels(2,i))
      l3 = int(labels(3,i))
      l4 = int(labels(4,i))
      index = Index2E(l1,l2,l3,l4)
      moints2(index) = moints2(index) + vals(i)
    end do
  end do
  
  ! deallocate arrays
  deallocate(buf, ibitv, vals, labels)
  deallocate(title, nbpsy,slabel, info, bfnlab, ietype, energy, imtype, map)
  
  close(aoints)

  return
contains
  !====================================================================
  !>sifreadd2
  ! ---------
  ! Read and decode 2-e integral record
  subroutine sifreadd2 (aoints, info, ninfo, nipv, iretbv, buf, lenbuf, &
    num, last, itypea, itypeb, ibvtyp, vals, nvals, labels, nlabels,    &
    ibitv, lenbitv, ierr)
    implicit none

    integer, intent(in) :: aoints, nipv, iretbv, ninfo, lenbuf, nlabels
    integer, intent(in) :: nvals, lenbitv
    integer, dimension(ninfo),  intent(in)       :: info
    real*8,  dimension(lenbuf), intent(inout)    :: buf
    real*8,  dimension(nvals),  intent(inout)    :: vals
    integer, dimension(4,nlabels), intent(inout) :: labels
    real*8,  dimension(lenbitv),   intent(inout) :: ibitv
    integer, intent(inout) :: num, last, itypea, itypeb, ibvtyp, ierr

    ! check value of split
    if (info(1) .ne. 1) then
            ierr = 100
            print *, "2-e records are separate. Not implemented."
            return
    end if

    ! read buffer
    read (aoints) buf

    ! unpack 2-e buffer
    call sifd2 (info, nipv, iretbv, buf, num, last, itypea, itypeb, &
      ibvtyp, vals, labels, ibitv, ierr)
    if (ierr .ne. 0) stop "*** Error unpacking 2-e buffer! ***"
    
    return

  end subroutine sifreadd2
  !====================================================================
  !>sifreadh2
  ! ---------
  ! Read header 2 of SIF formatted file
  subroutine sifreadh2 (aoints, ntitle, nsym, nbft, ninfo, nenrgy, nmap, &
    title, nbpsy, slabel, info, bfnlab, ietype, energy, imtype, map, ierr)
    implicit none

    integer, intent(in) :: aoints, ntitle, nsym, nbft, ninfo, nenrgy
    integer, intent(in) :: nmap
    integer, intent(inout) :: ierr
    character*80, dimension(ntitle),    intent(inout) :: title
    integer,      dimension(nsym),      intent(inout) :: nbpsy
    character*4,  dimension(nsym),      intent(inout) :: slabel
    integer,      dimension(ninfo),     intent(inout) :: info
    character*8,  dimension(nbft),      intent(inout) :: bfnlab
    integer,      dimension(nenrgy),    intent(inout) :: ietype
    real*8,       dimension(nenrgy),    intent(inout) :: energy
    integer,      dimension(nmap),      intent(inout) :: imtype
    integer,      dimension(nbft,nmap), intent(inout) :: map

    integer :: mapdim
    integer, parameter :: wrnerr = 0, nfterr = 1, faterr = 2
    
    ! set mapdim
    mapdim = max (nmap, 1)

    if (nmap .eq. 0) then
            read (aoints, iostat = ierr) &
              title, nbpsy, slabel, info, bfnlab, ietype, energy
    else
            read (aoints, iostat = ierr) &
              title, nbpsy, slabel, info, bfnlab, ietype, energy &
              , imtype, map
    end if

    return
  end subroutine sifreadh2
  !====================================================================
  !>Ind2Val
  ! Computes index of element i,j packed in lower triangle array
  !--------------------------------------------------------------------
  integer function Ind2Val(i,j)
    IMPLICIT NONE
    integer,intent(IN)      ::i,j
    if(i.GE.j)then
       Ind2Val=(i-1)*i/2+j
    else
       Ind2Val=(j-1)*j/2+i
    end if
  end function Ind2Val
  !====================================================================
  !>Index2E
  ! Computes index of 2-e integrals
  !--------------------------------------------------------------------
  integer function Index2E(i,j,k,l)
    IMPLICIT NONE
    integer,intent(IN)      ::i,j,k,l
    Index2E=Ind2Val(Ind2Val(i,j),Ind2Val(k,l))
  end function Index2E
  !====================================================================
end subroutine readmoints
