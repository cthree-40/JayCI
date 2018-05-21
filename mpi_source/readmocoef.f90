subroutine readmocoef (c, clen)
  implicit none

  ! ..input..
  integer, intent(in) :: clen
  
  ! ..output..
  double precision, dimension(clen), intent(out) :: c

  ! ..local variables..
  integer, parameter :: len = 1046529
  character*255 :: mocoeffl
  integer       :: mofl
  integer       :: ios, i

  ! ..moread() variables..
  ! mxtitl = maximum number of title cards
  ! filerr = mo file error number
  ! syserr = i/o system error
  ! nsym   = number of irreps in point group
  ! flcod  = file code (what to read from file)
  ! nbfpsy = number of basis functions per mo of each irrep
  ! nmopsy = number of mos per irrep
  ! labels = symmetry labels for each irrep
  integer, parameter :: mxtitl = 20
  integer :: filerr, ntitle, nsym, flcod, syserr
  integer, dimension(len) :: nbfpsy, nmopsy
  character*80, dimension(mxtitl) :: titles
  character*80 :: afmt
  character*4 :: labels

  print *, "Reading molecular orbital coefficents."
  mocoeffl = "mocoef_mc.sp"
  mofl = get_flunit()

  ! open molecular orbital coefficient file
  open(file = mocoeffl, unit = mofl, status = "old", form = "formatted", &
          action = "read", iostat = ios)
  if (ios .ne. 0) stop "*** Could not open molecular coefficient file! ***"

  ! Call colib subroutine: moread to read in molecular orbital coefficients
  flcod = 10
  call moread(mofl, flcod, filerr, syserr, mxtitl, ntitle, titles, afmt, nsym,&
          nbfpsy, nmopsy, labels, len, c)
  write(*,*) 'Titles from mocoef file:'
  do i = 1, ntitle
          write(*,*) titles(i)
  end do
  write(*,900) nsym
  write(*,901) (nbfpsy(i),i=1,nsym)
  write(*,902) (nmopsy(i),i=1,nsym)
  
900 format(' number of irreducible representations:',i4)
901 format(' number of basis functions per irrep:',8i4)
902 format(' number of molecular orbitals per irrep:',8i4)

  write (*,*) 'reading mocoefficients, flcod = 20'
  flcod = 20
  call moread(mofl, flcod, filerr, syserr, mxtitl, ntitle, titles, afmt, nsym,&
          nbfpsy, nmopsy, labels, len, c)
  if (filerr .eq. 0 .and. syserr .eq. 0) then
          write (*,*) 'coefficients successfully read.'
  else
          write (*,*) 'Error occured!', filerr, syserr
  end if
  
  return

contains

  ! get_flunit: get available unit number for file
  function get_flunit() result(u)
    implicit none

    ! ..output..
    integer :: u

    ! ..local variables..
    integer :: i
    logical :: uexist, uopen
    u = 0
    do i = 15, 99999
            inquire(unit = i, exist = uexist, opened = uopen)
            if (uexist .and. .not. uopen) then
                    u = i
                    exit
            end if
    end do
    if (u .eq. 0) stop "get_flunit: failed to find available unit."
  end function get_flunit
  
  
end subroutine readmocoef
