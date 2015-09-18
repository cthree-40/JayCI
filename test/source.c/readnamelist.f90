subroutine readnamelist(nmlist, nmlstr, err)
  ! readnamelist
  ! ------------
  ! Reads namelist.
  ! Input:
  !  nmlist = 1: general
  !           2:
  ! Output:
  !  nmlstr = namelist string
  implicit none

  ! .. input arguments ..
  integer, intent(in) :: nmlist

  ! .. output arguments ..
  character*300, dimension(30), intent(out) :: nmlstr
  integer, intent(out) :: err
  
  ! .. namelists ..
  integer :: gen_nml = 1
  
  ! .. &general arguments ..
  integer :: electrons, orbitals
  integer :: nfrozen, ndocc, nactive, nfrzvirt
  integer :: xlevel, printlvl

  namelist /general/ electrons, orbitals, nfrozen, ndocc, nactive, &
    nfrzvirt, xlevel, printlvl

  ! initialize error flag
  err = 0
  if (nmlist .eq. gen_nml) then
    electrons = 0
    orbitals  = 0
    nfrozen   = 0
    ndocc     = 0
    nfrzvirt  = 0
    xlevel    = 0
    printlvl  = 0
    open(file = "jayci.in", unit = 10, action = "read", status = "old", &
      iostat = err)
    if (err .ne. 0) return

    read(10, nml=general)

    ! write values to namelist array
    write(nmlstr(1),9) electrons
    write(nmlstr(2),9) orbitals
    write(nmlstr(3),9) nfrozen
    write(nmlstr(4),9) ndocc
    write(nmlstr(5),9) nactive
    write(nmlstr(6),9) xlevel
    write(nmlstr(7),9) nfrzvirt
    write(nmlstr(7),9) printlvl
    
    close(10)
    return
  else
    ! uknown namelist flag
    err = 99
    return
  endif
9 format(i10)
end subroutine readnamelist
