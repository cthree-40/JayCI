subroutine readnamelist(nmlist, nmlstr, err)
  ! readnamelist
  ! ------------
  ! Reads namelist.
  ! Input:
  !  nmlist = 1: general
  !           2: dalginfo
  !           3: dysonorb
  ! Output:
  !  nmlstr = namelist string
  implicit none

  ! .. input arguments ..
  integer, intent(in) :: nmlist
  
  ! .. output arguments ..
  character*300, dimension(30), intent(out) :: nmlstr
  integer, intent(out) :: err
  
  ! .. namelists ..
  integer :: gen_nml = 1, dai_nml = 2, dys_wf0 = 3, dys_wf1 = 4
  integer :: dys_orb = 5
  
  ! .. &general arguments ..
  integer :: electrons, orbitals
  integer :: nfrozen, ndocc, nactive, nfrzvirt
  integer :: xlevel, printlvl, printwvf

  ! .. &dalginfo arguments ..
  integer :: maxiter, krymin, krymax, nroots
  integer :: prediagr, refdim, buflen
  real*8  :: restol

  ! .. &wavefcn0 & &wavefcn1 unique arguments ..
  character*300 :: wvfcn_file0, wvfcn_file1
  integer :: nstates

  ! .. &dysonorbital arguments ..
  integer, dimension(5) :: states0, states1

  integer :: i, j
  
  namelist /general/ electrons, orbitals, nfrozen, ndocc, nactive, &
          nfrzvirt, xlevel, printlvl, printwvf
  namelist /dalginfo/ maxiter, krymin, krymax, nroots, prediagr, &
          restol, refdim, buflen
  namelist /wavefcn0/ electrons, orbitals, nfrozen, ndocc, nactive, &
          nfrzvirt, xlevel, nstates
  namelist /wavefcn1/ electrons, orbitals, nfrozen, ndocc, nactive, &
          nfrzvirt, xlevel, nstates
  namelist /dysonorbital/ states0, states1
  
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
          printwvf  = 0
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
          write(nmlstr(8),9) printlvl
          write(nmlstr(9),9) printwvf
          
          close(10)
          return
  else if (nmlist .eq. dai_nml) then
          ! dalginfo namelist
          maxiter   = 20
          krymin    =  3
          krymax    =  8
          nroots    =  1
          prediagr  =  1
          restol    = 1.0d-5
          refdim    =  3
          buflen    = 1000
          open(file = "jayci.in", unit = 10, action = "read", status = "old", &
                  iostat = err)
          if (err .ne. 0) return
          
          read(10, nml=dalginfo)
          
          write(nmlstr(1),9) maxiter
          write(nmlstr(2),9) krymin
          write(nmlstr(3),9) krymax
          write(nmlstr(4),9) nroots
          write(nmlstr(5),9) prediagr
          write(nmlstr(6),8) restol
          write(nmlstr(7),9) refdim
          write(nmlstr(8),9) buflen
          
          close(10)
          return
  else if (nmlist .eq. dys_wf0) then
          electrons = 0
          orbitals  = 0
          nfrozen   = 0
          ndocc     = 0
          nfrzvirt  = 0
          xlevel    = 0
          nstates   = 0
          open(file = "dycicalc.in", unit = 10, action = "read", status = "old", &
                  iostat = err)
          if (err .ne. 0) return

          read(10, nml=wavefcn0)

          ! write values to namelist array
          write(nmlstr(1),9) electrons
          write(nmlstr(2),9) orbitals
          write(nmlstr(3),9) nfrozen
          write(nmlstr(4),9) ndocc
          write(nmlstr(5),9) nactive
          write(nmlstr(6),9) xlevel
          write(nmlstr(7),9) nfrzvirt
          write(nmlstr(8),9) nstates
          
          close(10)
          return
  else if (nmlist .eq. dys_wf1) then
          electrons = 0
          orbitals  = 0
          nfrozen   = 0
          ndocc     = 0
          nfrzvirt  = 0
          xlevel    = 0
          nstates   = 0
          open(file = "dycicalc.in", unit = 10, action = "read", status = "old",&
                  iostat = err)
          if (err .ne. 0) return

          read(10, nml=wavefcn1)

          ! write values to namelist array
          write(nmlstr(1),9) electrons
          write(nmlstr(2),9) orbitals
          write(nmlstr(3),9) nfrozen
          write(nmlstr(4),9) ndocc
          write(nmlstr(5),9) nactive
          write(nmlstr(6),9) xlevel
          write(nmlstr(7),9) nfrzvirt
          write(nmlstr(8),9) nstates
          
          close(10)
          return

  else if (nmlist .eq. dys_orb) then
          states0 = (/ 1, 0, 0, 0, 0 /)
          states1 = (/ 1, 0, 0, 0, 0 /)
          open(file="dycicalc.in", unit=10, action="read", status="old", &
                  position="rewind", iostat=err)
          if (err .ne. 0) return

          read(10, nml=dysonorbital)

          ! write values to namelist array
          do i = 1, 5
                  j = i + 5
                  write(nmlstr(i),9) states0(i)
                  write(nmlstr(j),9) states1(i)
          end do

          close(10)

          return
  else
          ! uknown namelist flag
          err = 99
          return
  endif
8 format(f10.5)
9 format(i10)
end subroutine readnamelist
