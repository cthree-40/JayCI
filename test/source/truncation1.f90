module truncation
  !=============================================================================
  ! truncation
  ! ----------
  ! Purpose: Module containing subroutines to truncate CI space according to
  !          user-defined inputs.
  !
  ! Subroutines:
  !   citrunc1
  !   enf_actv_det
  !   enf_actv_str
  !   enf_docc_det
  !   enf_docc_str
  !   find_maxstr
  !   write_list2file
  ! ----------------------------------------------------------------------------
  ! Written by: cthree-40
  ! Dept. of Chemistry, The Johns Hopkins University
  !-----------------------------------------------------------------------------
  implicit none
contains
  subroutine citrunc1 (aelec, belec, orb, nfrzn, ndocc, nactv, xlvl, &
    astr_len, bstr_len, dtrm_len, ierr)
    ! citrunc1
    ! --------
    ! Purpose: driver for truncation of ci expansion
    !
    ! Input:
    !  aelec = alpha electrons
    !  belec = beta  electrons
    !  orb   = orbitals
    !  nfrzn = number of frozen core orbitals
    !  ndocc = number of docc orbitals
    !  nactv = number of active orbitals
    !  xlvl  = excitaiton level
    !
    ! Output:
    !  ierr       = error flag
    !  astr_len   = number of alpha strings after truncation
    !  bstr_len   = number of beta  strings after truncation
    !  dtrm_len   = number of determinants after truncation
    use combinatorial
    use addressing, only: strfind2, k2indc, indexk
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: aelec, belec, orb, nfrzn, ndocc, nactv, xlvl

    ! .. OUTPUT arguments ..
    integer, intent(out) :: ierr
    integer, intent(out) :: astr_len, bstr_len, dtrm_len

    ! .. LOCAL scalars ..
    ! ci_orbs    = (orbitals) - (frozen core)
    ! ci_aelec   = (alpha electrons) - (frozen core)
    ! ci_belec   = (beta  electrons) - (frozen core)
    ! axdocc     = alpha string docc excitations
    ! axactv     = alpha string active orbital exciations
    ! bxdocc     = beta  string docc excitations
    ! bxactv     = beta  string actv excitations
    ! deti       = determinant index
    ! astrings   = binom(ci_orbs, ci_aelec)
    ! bstrings   = binom(ci_orbs, ci_belec)
    ! tot_dets   = astrings * bstrings
    ! dtlst_unit = unit of determinant list file
    integer :: ci_orbs, ci_aelec, ci_belec
    integer :: axdocc, axactv, bxdocc, bxactv
    integer :: astrings, bstrings, tot_dets
    integer :: dtlst_unit = 12, strlst_unit = 13
    integer :: i, j
    integer :: deti
    
    ! .. Local arrays ..
    ! astring1    = alpha string
    ! astring2    = alpha string
    ! bstring1    = beta  string
    ! bstring2    = beta  string
    ! dtlst_flnm  = determinant list file name
    integer, dimension(:), allocatable :: astring1, astring2
    integer, dimension(:), allocatable :: bstring1, bstring2, bstring_start
    character*25 :: dtlst_flnm  = "det.list"
    character*25 :: strlst_flnm = "str.list"
    
    ! if frozen-core calculation, remove non-interacting electrons
    ci_orbs  = orb - nfrzn
    ci_aelec = aelec - nfrzn
    ci_belec = belec - nfrzn

    ! check input

    if (ndocc .gt. ci_belec) then
            write (*, "('Error: Unoccupied DOCC orbitals!')")
            ierr = -100
            return
    end if
    ! compute total, untruncated (excepting frozen-core truncation just
    !  performed) number of alpha/beta determinant strings
    astrings = binom(ci_orbs, ci_aelec)
    bstrings = binom(ci_orbs, ci_belec)
    ! compute total number of determinants
    tot_dets = astrings * bstrings
    write (*, "(1x,A,i15)") "Untruncated expansion size = ", tot_dets
    
    ! allocate arrays
    allocate(astring1(ci_aelec))
    allocate(astring2(ci_aelec))
    allocate(bstring1(ci_belec))
    allocate(bstring2(ci_belec))
    allocate(bstring_start(ci_belec))
    ! initialize (a/b)string1
    do i = 1, ci_aelec
            astring1(i) = i
    end do
    do i = 1, ci_belec
            bstring1(i) = i
    end do
    bstring_start = bstring1

    ! open expansion file
    open(file = dtlst_flnm, unit = dtlst_unit, status = "unknown", &
      position = "rewind", action = "write", iostat = ierr)
    if (ierr .ne. 0) stop "*** Cannot open expansion file for writing! ***"
    ! open string flie
    open(file = strlst_flnm, unit = strlst_unit, status = "unknown", &
      position = "rewind", action = "write", iostat = ierr)
    if (ierr .ne. 0) stop "*** Cannot open string list file for writing! ***"
    
    ! first strings are obviously included
    write(dtlst_unit, 10) 1, 1, 1
    write(strlst_unit,11) "BETA"
    write(strlst_unit,12) 1
    astr_len = 1
    bstr_len = 1
    dtrm_len = 1
    !**** FIRST LOOP P=1, Q = 1 .. bstrings
    ! loop over beta strings, pairing with first alpha string
    do i = 2, bstrings
            ! generate alpha string
            call strfind2(bstring1, ci_belec, ci_orbs, bstrings, bstring2)

            ! test docc restrictions on string
            if (ndocc .ne. 0) then
                    call enf_docc_str(bstring2, ci_belec, ndocc, bxdocc)
                    if (bxdocc .gt. xlvl) then
                            bstring1 = bstring2
                            cycle
                    end if
            else
                    bxdocc = 0
            end if

            ! test active orbital restrictions on string
            if (nactv .ne. 0) then
                    call enf_actv_str(bstring2, ci_belec, ndocc, nactv, bxactv)
                    if (bxactv .gt. xlvl) then
                            bstring1 = bstring2
                            cycle
                    end if
            else
                    bxactv = 0
            end if
            
            ! string is valid
            bstr_len = bstr_len + 1
            write(strlst_unit, 12) i
            ! compute determinant and write it to file
            deti = indexk(1, i, ci_orbs, ci_belec)
            write(dtlst_unit, 10) 1, i, deti
            dtrm_len = dtrm_len + 1
            bstring1 = bstring2
    end do

    ! **** SECOND LOOP P = 2 .. astrings, Q = 1.. bstrings
    ! loop over alpha strings
    write(strlst_unit, 11) "ALPHA"
    write(strlst_unit, 12) 1
    do i = 2, astrings
            ! generate alpha string
            call strfind2(astring1, ci_aelec, ci_orbs, astrings, astring2)

            ! test docc restrictions on string
            if (ndocc .ne. 0) then
                    call enf_docc_str(astring2, ci_aelec, ndocc, axdocc)
                    if (axdocc .gt. xlvl) then
                            astring1 = astring2
                            cycle
                    end if
            else
                    axdocc = 0
            end if

            ! test active orbital restrictions on string
            if (nactv .ne. 0) then
                    call enf_actv_str(astring2, ci_aelec, ndocc, nactv, axactv)
                    if (axactv .gt. xlvl) then
                            astring1 = astring2
                            cycle
                    end if
            else
                    axactv = 0
            end if
            ! if here, then p = i can be paired with q = 1 (no exciation)
            
            ! compute determanint
            astr_len = astr_len + 1
            write(strlst_unit, 12) i
            deti = indexk(i, 1, ci_orbs, ci_belec)
            write(dtlst_unit, 10) i, 1, deti
            dtrm_len = dtrm_len + 1

            bstring1 = bstring_start
            ! now we loop over bstrings 2 .. bstrings
            do j = 2, bstrings
                    ! generate beta string
                    call strfind2(bstring1, ci_belec, ci_orbs, bstrings, &
                      bstring2)
                    
                    ! test docc restrictions on string
                    if (ndocc .ne. 0) then
                            call enf_docc_str(bstring2, ci_belec, ndocc, bxdocc)
                            if (bxdocc .gt. xlvl) then
                                    bstring1 = bstring2
                                    cycle
                            end if
                    else
                            bxdocc = 0
                    end if

                    ! test active orbital restrictions on string
                    if (nactv .ne. 0) then
                            call enf_actv_str(bstring2, ci_belec, ndocc, nactv, bxactv)
                            if (bxactv .gt. xlvl) then
                                    bstring1 = bstring2
                                    cycle
                            end if
                    else
                            bxactv = 0
                    end if

                    ! if here than both strings are valid, so we must test
                    !  if the determinant they form is.

                    ! first, docc restriction test
                    if ((bxdocc + axdocc) .gt. xlvl) then
                            bstring1 = bstring2
                            cycle
                    end if

                    ! then active orbital test
                    if ((bxactv + axactv) .gt. xlvl) then
                            bstring1 = bstring2
                            cycle
                    end if

                    ! test together
                    if ((bxactv + axactv) + (bxdocc + axdocc) .gt. xlvl) then
                            bstring1 = bstring2
                            cycle
                    end if

                    ! if here, write determinant
                    deti = indexk(i, j, ci_orbs, ci_belec)
                    write(dtlst_unit, 10), i, j, deti
                    dtrm_len = dtrm_len + 1
                    bstring1 = bstring2
            end do
            
            astring1 = astring2
    end do

    deallocate(astring1, astring2, bstring1, bstring2, bstring_start)
    close(dtlst_unit)
    return
10  format(i15,i15,i15)
11  format(A5)
12  format(i15)
  end subroutine citrunc1

  subroutine enf_actv_str(string, elec, ndocc, nactv, test)
    !===========================================================================
    ! enf_docc_str
    ! ------------
    ! Purpose: enforce doubly-occupied orbital restrictions on alpha/beta
    !          strings
    !
    ! Input:
    !  string = alpha/beta string to test
    !  elec   = number of electrons
    !  ndocc  = number of docc orbitals
    !  nactv  = number of active orbitals
    !
    ! Output:
    !  test = 1 if string is valid
    !         0 if string is not valid
    !---------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: elec, ndocc, nactv
    integer, dimension(elec), intent(in) :: string

    ! .. OUTPUT arguments ..
    integer, intent(out) :: test

    ! .. LOCAL scalars ..
    integer :: i

    ! loop over string
    test = 0
    do i = 1, elec
            if (string(i) .gt. (ndocc + nactv)) then
                    test = test + 1
            end if
    end do
    
    return
  end subroutine enf_actv_str

  subroutine enf_docc_str(string, elec, ndocc, test)
    !===========================================================================
    ! enf_docc_str
    ! ------------
    ! Purpose: enforce doubly-occupied orbital restrictions on alpha/beta
    !          strings
    !
    ! Input:
    !  string = alpha/beta string to test
    !  elec   = number of electrons
    !  ndocc  = number of docc orbitals
    !  xlvl   = excitation level of expansion
    !
    ! Output:
    !  test   = number of excitations
    !---------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: elec, ndocc
    integer, dimension(elec), intent(in) :: string

    ! .. OUTPUT arguments ..
    integer, intent(out) :: test

    ! .. LOCAL scalars ..
    integer :: i

    ! loop over string
    test = 0
    do i = 1, ndocc
            if (string(i) .gt. ndocc) then
                    test = test + 1
            end if
    end do
    
    return
  end subroutine enf_docc_str

  subroutine find_maxstr(list, len, max)
    !===========================================================================
    ! find_maxstr
    ! -----------
    ! Purpose: find max index in string
    !
    ! Input:
    !  list = list of integers (1 or 0)
    !  len  = length of list
    !
    ! Output:
    !  max  = greatest index of list .ne. 0
    !---------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer,                 intent(in) :: len
    integer, dimension(len), intent(in) :: list

    ! .. OUTPUT arguments ..
    integer, intent(out) :: max

    ! .. LOCAL scalars ..
    integer :: i

    max = 1
    do i = len, 1, -1
            if (list(i) .ne. 0) then
                    max = i
                    return
            end if
    end do

    ! if here, error
    max = 0
    return
  end subroutine find_maxstr
  
  subroutine write_list2file(list, len, flname, flunit, num_nz, ierr)
    !===========================================================================
    ! write_list2file
    ! ---------------
    ! Purpose: write integer list to file.
    !          list is a 1-d integer array containing 1 or 0 for each element.
    !
    ! Input:
    !  list   = list of integers (1 or 0)
    !  len    = length of list
    !  flname = name of file to print to
    !  flunit = unit of file to print to
    !
    ! Output:
    !  num_nz = number of nonzero elements
    !  ierr   = error handling
    !---------------------------------------------------------------------------
    implicit none

    ! .. INPUT arguments ..
    integer,                 intent(in) :: flunit, len
    integer, dimension(len), intent(in) :: list
    character*25,            intent(in) :: flname

    ! .. OUTPUT arguments ..
    integer, intent(out) :: ierr, num_nz

    ! .. LOCAL scalars ..
    integer :: i

    ! open file
    open(file = trim(flname), unit = flunit, status = "unknown", &
      action = "write", position = "rewind", iostat = ierr)
    if (ierr .ne. 0) return

    ! write list
    num_nz = 0
    do i = 1, len
            if (list(i) .ne. 0) then
                    write(flunit, 10) i
                    num_nz = num_nz + 1
            end if
    end do

    ! close file
    close(flunit)

    return
10  format(i15)
  end subroutine write_list2file
end module truncation
