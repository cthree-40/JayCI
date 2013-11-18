program test

  implicit none

  integer :: adets, bdets, aelec, belec, orbitals, nfrozen, &
             ndocc, nactive, xlevel, cidimension, adetlen, bdetlen

!-------------------------

   print *, "Starting..."

   adets    = 20
   bdets    = 20
   aelec    = 3
   belec    = 3
   orbitals = 6
   nfrozen  = 1
   ndocc    = 2
   nactive     = 2
   xlevel   = 1

   call citrunc( adets, bdets, aelec, belec, orbitals, nfrozen,&
                 ndocc, nactive, xlevel, adetlen, bdetlen, cidimension )
   
   print *, cidimension
end program
