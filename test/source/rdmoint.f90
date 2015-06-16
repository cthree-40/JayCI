module rdmoint
  ! module for reading in 1-e and 2-e integrals
  !----------------------------------------------------
  ! moint1 = 1-e integral array
  ! moint2 = 2-e integral array
  ! m1len  = len(moint1)
  ! m2len  = len(moint2)
  use integral
  implicit none
  
  real*8, dimension(:), allocatable :: moint1, moint2
  integer                           :: m1len,  m2len

contains
  !>rdheader1
  ! subroutine to read header 1 of sifs integral file
  subroutine rdheader1( mointsfl, ntitle, nsym, nbfcns, &
    ninfo, nenrgy, nmap, ierr )
    ! mointsfl = unit number of molecular integral file
    ! ntitle   = number of titles
    ! nsym     = number of symmetry blocks
    ! nbfcns   = total number of basis functions
    ! ninfo    = number of record-definition parameters
    ! nenrgy   = number of core energies.
    ! nmap     = number of optional map vectors
    !
    ! sifvrs   = routine library version number.
    ! ntitmx   = maximum number of titles allowed.
    ! ninchk   = minimum number of info(*) elemets.
    ! lrecmx   = maximum record length allowed.
    implicit none
    integer, intent(in) :: moints
    integer, intent(out):: ntitle, nsym, nbft, ninfo, &
      nenrgy, nmap, ierr
    integer             :: sifvrs, ntitmx, ninchk, lrecmx
    integer             :: ioerr

    ! read header 1 from file
    read( unit=moints, format="unformatted", position="rewind", &
      action="read", iostat=ioerr )                             &
      sifvrs, ntitle, nsym, nbfcns, nenrgy, nmap

    ! test for errors
    if ( ioerr .ne. 0 ) then
            ierr = ioerr
            return ! return with no errors
    else
            write(*,"(A,i3)") "Could not open file. unit= ", moints
            stop "Exiting..."
    end if

  end subroutine rdheader1
  !>rdheader2
  ! subroutine to read header 2 of SIFS integral file.
  subroutine rdheader2( moints, ntitle, nsym, nbfcns, ninfo, &
    nenrgy, nmap, title, nbfpsy, symlabel, info, bsfnlab,    &
    ietype, energy, imtype, map, ierr )
  end subroutine
end module rdmoint
