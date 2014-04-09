module detci3new
  use detci1
  use detci2
  implicit none
! This module contains the subroutines required to truncate the CI
!  expansion based upon frozen core, docc and active space restrictions
! [FC][DOCC][CAS][VIRTS]
!   =[No excitations out of][ xlevel excitations into CAS and VIRTS]
!			[Any number of excitations within and xlevel excitations into VIRTS]
contains
  subroutine frozen(nfrzn, aelec, belec, totels, orbitals, adets, & 
                      bdets, alphmat, betamat, strlen, alphfzl, betafzl, str)
    implicit none
    integer,intent(in) :: nfrzn, aelec, belec, orbitals, totels, strlen, &
      alphfzl, betafzl, adets, bdets
    integer,dimension(adets,aelec) :: alphmat 
    integer,dimension(bdets,belec) :: betamat
    integer,dimension(strlen) :: str
!Input:
! nfrzn    = number of frozen orbitals
! aelec    = number of alpha electrons
! belec    = number of beta electrons
! orbitals = number of orbitals
! adets    = number of alpha determinant strings
! bdets    = number of beta determinant strings
! alphmat  = alpha string array

! betamat  = beta string array
! strlen   = length of str
!Output:

! str      = string of determinant indices satisifying the requirement
!            that the first nfrzn orbitals are doubly occupied
! alphfzl  = length of astr
! betafzl  = legnth of bstr
    integer,dimension(alphfzl):: astr
    integer,dimension(betafzl):: bstr
    integer :: i, k, p, q, iter, test, j
! # Generate the string of determinant addresses for alpha/beta strings
    iter=1
20  if (iter .le. alphfzl) then
      do i = iter, adets
        test = 0
        do j = 1, nfrzn
          if (alphmat(i,j) == j) then
            test = test + 1 
          end if
        end do
        if (test == nfrzn) then
          astr(iter) = i
          iter = iter + 1
          go to 20
        end if
      end do
    end if
    iter = 1
30  if (iter .le. betafzl) then
			do i = iter, bdets
				test = 0
				do j = 1, nfrzn
					if (betamat(i,j) == j) then
						test = test + 1
					end if
				end do
				if (test == nfrzn) then
					bstr(iter) = i
					iter = iter + 1
					go to 30
				end if
			end do
		end if
! # Generate list of determinants with correct occupations
		do i = 1, alphfzl
			do j = 1, betafzl
				str(j + betafzl*(i-1)) = indxk(astr(i),bstr(j),belec, orbitals)
			end do
		end do
		return
	end subroutine
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  subroutine docc(numdocc, xlevel, nfrzn, aelec, belec, totels, orbitals, &
										adets, bdets, alphmat, betamat, fzdlist, fzdlen, &
										dcclist, dcclen)
! # This subroutine truncates the determinant list further by requiring
!   the numdocc orbitals after the frozen orbitals be occupied by at most
!   2*numdocc electrons and at least 2*numdocc-xlevel electrons
!Input:
! numdocc  = number of docc ORBITALS
! xlevel   = excitation level
! aelec    = alpha electrons
! belec    = beta electrons
! totels   = aelec + belec
! orbitals = number of orbitals
! adets    = number of alpha determinant strings
! bdets    = number of beta determinant strings
! alphmat  = number of alpha determinants
! betamat  = number of beta determinants
! fzdlist  = list of determinants satisfying the frozen orb. requirement
! fzdlen   = length of fzdlist
!Output:
! dcclist  = list of determinants satisfying the frozen orbital requirement
! dcclen   = length of fzdlist
		implicit none
		integer,intent(in) :: numdocc, xlevel, aelec, belec, totels, orbitals, &
														adets, bdets, nfrzn, fzdlen
		integer,dimension(adets,aelec) :: alphmat
		integer,dimension(bdets,belec) :: betamat
		integer,dimension(fzdlen) :: fzdlist,dlist
		integer,dimension(:),allocatable :: dcclist
		integer,intent(out) :: dcclen
		integer :: i, p, q, nils, j, test
! # Loop through the determinant indices in fzdlist making indices 0 that do not
!   satisfy the DOCC requirement and counting how many are not made zero
!   dlist(fzdlen) is the list generated from ^

    print *, "Welcome to docc()"

    nils = 0
    dcclen = 0

    do i=1, fzdlen
			call k2indc(fzdlist(i),belec,orbitals,p,q)
			test = 0
			do j = nfrzn+1,nfrzn+numdocc
				if (alphmat(p,j) > nfrzn+numdocc) then ! Test alpha excitations
					test = test + 1
				end if
			end do
			do j = nfrzn+1,nfrzn+numdocc
				if (betamat(q,j) > nfrzn+numdocc) then ! Test beta excitations
					test = test + 1
				end if
			end do
			if (test > xlevel) then  ! Check if the number of excitations > allowed
			  dlist(i)=0 						! If so make index 0
			  nils = nils + 1
			else
				dlist(i) = fzdlist(i) ! If not do nothing
				dcclen = dcclen + 1
			end if
		end do
! # Allocate dcclist
		allocate(dcclist(dcclen))
		i=1
		do j = 1, dcclen
40 		if ( i .le. fzdlen ) then
				if (dlist(i) .ne. 0) then
					dcclist(j) = dlist(i)
					i = i + 1
					cycle
				else
					i = i + 1
					go to 40
				end if
			end if
		end do
		return
	end subroutine
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	subroutine active(numcas,numdocc,nfrzn,xlevel,aelec, belec, totels, &
											orbitals, adets, bdets, alphmat, betamat, &
											dcclist, dcclen, actlist, &
											actlen)
		implicit none
		integer,intent(in) :: numcas,numdocc,nfrzn,xlevel,aelec,belec,totels,&
														orbitals,adets,bdets,dcclen
		integer,dimension(adets,aelec)::alphmat
		integer,dimension(bdets,belec)::betamat
		integer,dimension(dcclen)::dcclist
		integer,intent(out)::actlen
		integer,dimension(:),allocatable::actlist
! # Note: dcclist can be either a list modified by frzn, or docc.
!         both numdocc and nfrzn are asked for for loop parameters
!Input:
! numcas   = number of active orbitals
! numdocc  = number of docc orbitals
! nfrzn    = number of frozen orbitals
! xlevel   = excitation level
! aelec    = alpha electrons
! belec    = beta electrons
! totels   = total electrons
! orbitals = number of orbitals
! adets    = alpha determinants
! bdets    = beta determinants
! alphmat  = alpha string array
! betamat  = beta string array
! dcclist  = list that has been truncated by either/both frzn and docc
! dcclen   = length of dcclist
!Output:
! actlist  = list of determinant indices truncated by active
! actlen   = length of actlist
! ---------------------------
! CAS allows for any excitations within it, but only xlevel excitations out
!  of it
! # Most/all of this program follows the structure of docc()
		integer,dimension(dcclen) :: alist
		integer :: i, p, q, nils, j, test
! # Loop through the detemrinant indices in dcclist making indices 0 that
!   do not satisfy the CAS requirement and counting how many are not made zero
!   alist(dcclen) is the list generated from ^
		nils = 0
		actlen = 0
		do i=1, dcclen
			call k2indc(dcclist(i),belec,orbitals,p,q)
			test = 0
			do j = nfrzn+numdocc+1,aelec
				if (alphmat(p,j) > nfrzn+numdocc+numcas) then
					test = test + 1
				end if
			end do
			do j = nfrzn+numdocc+1,belec
				if (betamat(q,j) > nfrzn+numdocc+numcas) then
					test=test+1
				end if
			end do
			if (test > xlevel) then
				alist(i) = 0
				nils = nils + 1
			else
				alist(i) = dcclist(i)
				actlen = actlen + 1
			end if
		end do
! # Allocate actlist
		allocate(actlist(actlen))
		i=1
		do j=1,actlen
40		if (i .le. dcclen ) then
				if (alist(i) .ne. 0 ) then
					actlist(j) = alist(i)
					i = i + 1
					cycle
				else
					i = i + 1
					go to 40
				end if
			end if
		end do
		return
	end subroutine
end module detci3new
		
					
