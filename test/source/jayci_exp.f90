program jayci_exp
  !=============================================================================
  ! jayci_exp
  ! ---------
  !
  ! Purpose: generate the wavefunction expansion for jayci.x
  !
  ! Written by: cthree-40
  ! Dept. of Chemistry, The Johns Hopkins University
  !=============================================================================
  use input_proc, only: alpha_beta
  use combinatorial
  use truncation
  implicit none

  ! -- IN/OUT FILES ------------------------------------------------------------
  ! ipfl_name  = input file name
  ! ipfl_unit  = input file unit
  ! outfl_name = output file name
  ! outfl_unit = output file unit
  ! expfl_name = expansion file name
  ! expfl_unit = expansion file unit
  character*20 :: ipfl_name = "jayci.in", outfl_name = "jayci.out"
  character*20 :: expfl_name = "expansion.def"
  integer      :: ipfl_unit = 10, outfl_unit = 11, expfl_unit = 12

  ! -- general NAMELIST ---------------------------------------------------------
  ! electrons = electrons in system
  ! orbitals  = orbitals in system
  ! nfrozen   = number of frozen core orbitals
  ! ndocc     = number of doubly-occupied orbitals
  ! nactive   = number of active orbitals
  ! xlevel    = excitation level of expansion
  ! printlvl  = print level
  ! nrzvirt   = number of frozen virtuals
  integer :: electrons, orbitals, nfrozen, ndocc, nactive, xlevel, nfrzvirt
  integer :: printlvl
  namelist /general/ electrons, orbitals, nfrozen, ndocc, nactive, xlevel, &
    printlvl, nfrzvirt

  ! .. LOCAL scalars ..
  ! aelec   = alpha electrons
  ! belec   = beta  electrons
  ! ierr    = error handling
  ! astrlen = number of alpha strings after truncation
  ! bstrlen = number of beta  strings after truncation
  ! dtrmlen = number of determinants  after truncation
  integer :: aelec, belec
  integer :: ierr
  integer :: astrlen, bstrlen, dtrmlen
  
  ! default &general namelist values
  electrons = 0
  orbitals  = 0
  nfrozen   = 0
  ndocc     = 0
  nactive   = 0
  nfrzvirt  = 0
  xlevel    = 2 ! This won't return error. Default CI program is SDCI.
  printlvl  = 0
  
  ! open output file
  open(unit = outfl_unit, file = outfl_name, status = "unknown", &
    position = "rewind", action = "write", iostat = ierr)
  if (ierr .ne. 0) stop "*** Cannot open output file! ***"
  
  ! write header
  call write_header(outfl_unit)
  
  ! read input file
  open(unit = ipfl_unit, file = ipfl_name, status = "old", action = "read", &
    iostat = ierr)
  if (ierr .ne. 0) stop "*** Cannot open input file! ***"
  read(unit = ipfl_unit, nml = general)
  close(ipfl_unit)
  ! check namelist input values
  if (electrons .eq. 0) then
          stop "*** electrons = 0! ***"
  else if (orbitals .eq. 0) then
          stop "*** orbitals = 0! ***"
  end if

  ! write input values
  call write_gennml(electrons, orbitals, nfrozen, ndocc, nactive, xlevel, &
    nfrzvirt, outfl_unit)
  ! set alpha/beta electron numbers
  call alpha_beta(electrons, aelec, belec)

  ! truncate ci space
  write(*, "(1x, A)") "Truncating CI space..."
  call citrunc1(aelec, belec, orbitals, nfrozen, ndocc, nactive, nfrzvirt, &
    xlevel, astrlen, bstrlen, dtrmlen, ierr)
  if (ierr .ne. 0) stop "*** Error during truncation! ***"

  write(*, "(1x, A, i15)") "Determinants =", dtrmlen
  write(*, "(1x, A, i15)") "Alpha strings=", astrlen
  write(*, "(1x, A, i15)") "Beta  strings=", bstrlen

  ! generate input for jayci.x
  write (*, "(1x, A)") "Generating input for jayci.x..."
  call gen_input(dtrmlen, astrlen, bstrlen, aelec, belec, orbitals, nfrozen, &
    nfrzvirt, ierr)
  
contains

  subroutine gen_input(num_det, num_astr, num_bstr, aelec, belec, orbs, &
    nfrzn, nfvrt, ierr)
    !===========================================================================
    ! gen_input
    ! ---------
    ! Purpose: Write input file for jayci.x: input.jayci
    !
    ! Input:
    !  num_det  = number of determinants  in expansion
    !  num_astr = number of alpha strings in expansion
    !  num_bstr = number of beta  strings in expansion
    !  aelec    = number of alpha electrons
    !  belec    = number of beta  electrons
    !  orbs     = number of orbitals
    !  nfrzn    = number of frozen orbitals
    !  nfvrt    = number of frozen virtual orbitals
    ! Output:
    !  ierr     = error handling
    !---------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: num_det, num_astr, num_bstr
    integer, intent(in) :: aelec, belec, orbs
    integer, intent(in) :: nfrzn, nfvrt

    ! .. OUTPUT arguments ..
    integer, intent(out) :: ierr

    ! .. LOCAL arrays ..
    ! flname = name of input file
    character*25 :: flname = "input.jayci"
    
    ! .. LOCAL scalars ..
    ! flunit   = unit of input file name
    ! ci_orbs  = orbs  - nfrzn
    ! ci_aelec = aelec - nfrzn
    ! ci_belec = belec - nfrzn
    integer :: flunit = 14
    integer :: ci_orbs, ci_aelec, ci_belec

    ci_orbs  = orbitals - nfrozen - nfvrt
    ci_aelec = aelec - nfrozen
    ci_belec = belec - nfrozen
    
    ! open input file
    open(file = flname, unit = flunit, status = "unknown", position = "rewind",&
      action = "write", iostat = ierr)
    if (ierr .ne. 0) stop "*** Error writing input for jayci! ***"

    ! write file
    write(flunit, 10) "&expansion"
    write(flunit, 11) "ci_orbitals   =", ci_orbs
    write(flunit, 11) "ci_electrons  =", (ci_aelec + ci_belec)
    write(flunit, 11) "astr_len   =", num_astr
    write(flunit, 11) "bstr_len   =", num_bstr
    write(flunit, 11) "dtrm_len   =", num_det
    write(flunit, 10) "/"

10  format(A)
11  format(1x,A,i15)
  end subroutine gen_input

  subroutine write_gennml(elec, orb, nfrzn, ndocc, nactv, xlvl, nfvrt, outfl_unit)
    !===========================================================================
    ! write_gennml
    ! ------------
    ! Purpose: Write &general namelist input values.
    !
    ! elec  = number of electrons
    ! orb   = number of orbitals
    ! nfrzn = number of frozen core orbitals
    ! ndocc = number of docc orbitals
    ! nactv = number of active orbitals
    ! xlvl  = excitation level
    !---------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: outfl_unit, elec, orb, nfrzn, ndocc, nactv, nfvrt
    integer, intent(in) :: xlvl

    ! write input values
    write(outfl_unit, 10) ""
    write(outfl_unit, 10) "&general input values:"
    write(outfl_unit, 11) "electrons =", elec
    write(outfl_unit, 11) "orbitals  =", orb
    write(outfl_unit, 11) "nfrozen   =", nfrzn
    write(outfl_unit, 11) "ndocc     =", ndocc
    write(outfl_unit, 11) "nactive   =", nactv
    write(outfl_unit, 11) "nfrzvirt  =", nfvrt
    write(outfl_unit, 11) "xlevel    =", xlvl
    write(outfl_unit, 10) ""

    return
10  format(1x,A)
11  format(2x,A,i5)
  end subroutine write_gennml
  subroutine write_header(outfl_unit)
    !===========================================================================
    ! write_header
    ! ------------
    ! Purpose: Write header for jayci_exp.
    !---------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: outfl_unit

    write(unit = outfl_unit, fmt = "(A)") " ------------------------------------"
    write(unit = outfl_unit, fmt = "(A)") " jayci_exp.x "
    write(unit = outfl_unit, fmt = "(A)") " ------------------------------------"

    return
  end subroutine write_header
end program jayci_exp
  
