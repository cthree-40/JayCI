program test
  
  implicit none
  integer :: adets, bdets, aelec, belec, orbitals, nfrozen, &
             ndocc, nactive, xlevel, cidimension, adetlen, bdetlen

!-------------------------

   print *, "Starting..."

   adets    = 15504
   bdets    = 15504
   aelec    = 5
   belec    = 5
   orbitals = 20
   nfrozen  = 1
   ndocc    = 2
   nactive  = 5
   xlevel   = 2

   call citrunc( adets, bdets, aelec, belec, orbitals, nfrozen,&
                 ndocc, nactive, xlevel, adetlen, bdetlen, cidimension )
   
end program
