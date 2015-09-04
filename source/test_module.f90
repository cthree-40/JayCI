module test_module
  ! Testing module for jayci.
  implicit none

contains
  !*
  !*
  subroutine test_cannon4()
    !===========================================================================
    ! test_cannon4()
    ! --------------
    ! Purpose: test cannon4() subroutine
    !---------------------------------------------------------------------------
    use string_util, only: cannon4
    implicit none

    ! .. cannon4 input arguments ..
    integer :: strlen, indx1

    ! .. cannon4 input/output arguments ..
    integer, dimension(5) :: string, order
    integer, dimension(10):: scratch
    
    ! .. cannon4 output arguments ..
    integer :: perm

    ! .. local arrays ..
    character*25 :: inflnm = "test.in"

    ! .. local scalars ..
    integer :: influn = 9
    integer :: subval

    ! .. test_cannon4 namelist ..
    namelist /testcannon4/ strlen, indx1, subval

    ! open input file
    open(file = inflnm, unit = influn, status = "old", action = "read")
    read(influn, nml = testcannon4)

    ! write namelist input
    write(*,10) strlen, indx1, subval
10  format('strlen =', i5, /, 'indx1 =', i3, /, &
      'subval =', i5)
    close(influn)
    
    if (strlen .ne. 5) stop "strlen .ne. 5 !!"
    ! make substitution
    string = (/0, 2, 4, 6, 8/)
    string(indx1) = subval
    write(*,"('old string =', 5i3)") string

    ! call cannon4
    order = (/1, 2, 3, 4, 5/)
    perm = 1
    call cannon4(string, strlen, indx1, scratch, order, perm)

    ! print final string
    write(*,"('new string 2=', 5i3)") string
    write(*,"('new order   =', 5i3)") order
    write(*,"('permutation =', i3)" ) perm
  end subroutine test_cannon4
  !*
  !*
  subroutine test_citrunc()
    !===========================================================================
    ! test_citrunc()
    ! --------------
    ! Purpose: test citrunc() subroutine.
    !---------------------------------------------------------------------------
    use truncation, only: citrunc1
    implicit none

    ! .. CITRUNC1 INPUT arguments ..
    integer :: aelec, belec, orbitals, nfrozen, ndocc, nactive, xlevel, nfrzv

    ! .. CITRUNC1 OUTPUT arguments ..
    integer :: astr_len, bstr_len, dtrm_len, ierr

    ! .. LOCAL arrays ..
    character*25 :: inflnm = "test.in"

    ! .. LOCAL scalars ..
    integer :: influn = 9

    ! .. TESTCITRUNC namelist ..
    namelist /testcitrunc/ aelec,   belec, orbitals, &
                           nfrozen, ndocc,  nactive, nfrzv, xlevel
    nfrzv = 0
    aelec = 0
    belec = 0
    orbitals = 0
    nfrozen = 0
    ndocc = 0
    nactive = 0
    xlevel = 0
    
    ! open input file and read namelist
    open(file = inflnm, unit = influn, action = "read", status = "old", &
      iostat = ierr)
    if (ierr .ne. 0) stop "*** Error opening input file: test.in! ***"
    read(unit = influn, nml = testcitrunc)
    close(influn)
    
    ! write namelist input to output stream
    write(*, 10) "aelec", aelec, "belec", belec, "orbitals", orbitals, &
      "nfrozen", nfrozen, "ndocc", ndocc, "nactive", nactive, "xlevel",&
      xlevel
10  format(7(a15, " =", i5, /))

    ! call citrunc
    call citrunc1(aelec, belec, orbitals, nfrozen, ndocc, nactive, &
      nfrzv, xlevel, astr_len, bstr_len, dtrm_len, ierr)

    ! print output arguments
    write(*, "(a)") ""
    write(*, 11) "Alpha strings", astr_len, "Beta strings", bstr_len, &
      "Determinants", dtrm_len
11  format(3(a15, " =", i15, /))
    
    return
  end subroutine test_citrunc
  !*
  !*
end module test_module
    
