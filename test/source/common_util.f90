!********************************************************************
! This module contains combinatorial utilities used in JayCI
!====================================================================
MODULE combinatorial
  IMPLICIT NONE
  PUBLIC      :: binom
  integer,dimension(:,:),ALLOCATABLE  :: Binom_Data
  integer     :: Binom_Data_MaxN,Binom_Data_MaxK
CONTAINS
  !=====================================================================
  !>Factl(m,n) returns the product of m-->n
  !---------------------------------------------------------------------
  integer function factl(m,n)
    IMPLICIT NONE
    integer,intent(IN)                  ::m,n
    integer                             ::i
    if(n.LT.0)then
       Factl=0
    else if(n.EQ.0)then
       Factl=1
    else
       Factl=PRODUCT((/ (i,i=m,n)/))
    end if
  end function factl
  !====================================================================
  !>Binom(n,k) returns binomial coefficient of n and k
  !--------------------------------------------------------------------
  integer function binom(n,k)
    IMPLICIT NONE
    integer,intent(IN)                  ::n,k
    !integer,dimension(:,:),ALLOCATABLE  ::tmp_BD
    integer     ::i,j
    !Test if array has been initialized
    if(.NOT.ALLOCATED(Binom_Data))then !Initialization
       Binom_Data_MaxN=130
       Binom_Data_MaxK=10
       ALLOCATE(Binom_Data(0:Binom_Data_MaxN,0:Binom_Data_MaxK))
       do i=0,Binom_Data_MaxN
          do j=0, MIN(i,Binom_Data_MaxK)
             if(j.EQ.0 .OR. j.EQ.i)then
                Binom_Data(i,j)=1
             else
                Binom_Data(i,j)= Binom_Data(i-1,j-1) + &
                     Binom_Data(i-1,j)
             end if
          end do
       end do
    end if
    !if(n.LE.Binom_Data_MaxN .AND. k.LE.Binom_Data_MaxK)then
    binom=Binom_Data(n,k)
  end function binom
  !====================================================================
END MODULE combinatorial
!********************************************************************
! This module contains addressing utilities used in JayCI
!====================================================================
MODULE addressing
  IMPLICIT NONE
  PUBLIC      :: GenOrbString_Alpha,GenOrbString_Beta
  integer,dimension(:,:),ALLOCATABLE  :: alpha_strings,beta_strings
  integer,dimension(:),ALLOCATABLE    :: alpha_indices,beta_indices
  integer     ::str_peg_int_a, str_peg_int_b     !String peg interval
  
CONTAINS
  !====================================================================
  !>IndexK
  ! Computes the index of determinant |p,q>
  !--------------------------------------------------------------------
  integer function IndexK(p,q,Orbitals,bElec)
    use combinatorial, only: Binom
    IMPLICIT NONE
    integer,intent(IN)      ::p,q,Orbitals,bElec
    IndexK=(p-1)*Binom(Orbitals,bElec)+q
  end function IndexK
  !====================================================================
  !>K2Indc
  ! Returns alpha and beta string indices of determiant K
  !--------------------------------------------------------------------
  subroutine K2Indc(k,bElec,Orbitals,p,q)
    use combinatorial, only: Binom
    IMPLICIT NONE
    integer,intent(IN)      ::k,bElec,Orbitals
    integer,intent(OUT)     ::p,q
    p=ceiling(dble(k)/dble(Binom(Orbitals,bElec)))
    q=k-(p-1)*Binom(Orbitals,bElec)
    RETURN
  end subroutine K2Indc
  !====================================================================
  !>AdrFind
  ! Returns the address of an orbital index string
  !--------------------------------------------------------------------
  subroutine AdrFind(string,elecs,Orbitals,address)
    use combinatorial, only: Binom
    IMPLICIT NONE
    integer,intent(IN)      ::elecs,Orbitals
    integer,dimension(elecs),intent(IN) ::string
    integer,intent(OUT)     ::address
    integer     ::i,j
    address=1                  !Start at 1
    do i=2,elecs
       do j=string(i-1)+1,string(i)-1
          address=address+Binom(Orbitals-j,elecs-i)
       end do!j
    end do!i
    do j=1,string(1)-1
       address=address+Binom(Orbitals-j,elecs-1)
    end do
    RETURN            !Return address
  end subroutine AdrFind
  !===================================================================
  !>StrFind1
  ! Generates strings up to index given
  !-------------------------------------------------------------------
  subroutine StrFind1(address,spElecs,orbs,spDets,array)
    IMPLICIT NONE
    integer,intent(IN)      ::address,spElecs,orbs,spDets
    integer,dimension(spElecs,address),intent(OUT)::array
    integer     ::i,j
    integer,dimension(:),ALLOCATABLE::tmp
    !Initialize array(1,i)
    do i=1,spElecs
       array(i,1)=i
    end do
    do i=2,address
       CALL StrFind2(array(1,i-1),spElecs,orbs,spDets,array(1,i))
    end do
    RETURN      !Return array
  end subroutine StrFind1
  !====================================================================
  !>StrFind2
  ! Generates index+1 string from index string
  !--------------------------------------------------------------------
  subroutine StrFind2(string1,spElecs,orbs,spDets,string2)
    IMPLICIT NONE
    integer,intent(IN)      ::spElecs,orbs,spDets
    integer,dimension(spElecs),intent(IN)     ::string1
    integer,dimension(spElecs),intent(OUT)    ::string2
    integer::i,j
    !Loop over electrons
    
    do i=0,spElecs-1
       if(string1(spElecs-i).EQ.(orbs-i))then
          CYCLE !i->i+1
       else !Determine new string
          do j=1,(spElecs-i-1)
             string2(j)=string1(j)
          end do
          do j=(spElecs-i),spElecs
             string2(j)=string1((spElecs-i))+1+j-(spElecs-i)
          end do
       end if
       EXIT  !i
    end do !i
    RETURN
  end subroutine StrFind2
  !====================================================================
  !>StrFind3
  ! Generates index+n string from index string
  !--------------------------------------------------------------------
  subroutine StrFind3(string1,index1,spElecs,orbs,spDets,index2,string2)
    IMPLICIT NONE
    integer,intent(IN)      ::spElecs,orbs,spDets,index1,index2
    integer,dimension(spElecs),intent(IN)     ::string1
    integer,dimension(spElecs),intent(OUT)    ::string2
    integer,dimension(:,:),ALLOCATABLE    ::temp
    integer::i,j
    ALLOCATE(temp(spElecs,index2-index1+1))       !Allocate temporary array
    temp(1:spElecs,1)=string1(1:spElecs)
    do i=1,index2-index1
       CALL StrFind2(temp(1,i),spElecs,orbs,spDets,temp(1,i+1))
    end do
    string2=temp(1:spElecs,index2-index1+1)
    DEALLOCATE(temp)              !Deallocate temporary array
    RETURN
  end subroutine StrFind3
  
  !====================================================================
  !>GenOrbString
  ! Generates orbital string i
  !--------------------------------------------------------------------
  subroutine GenOrbString(str_index,str_elec,Orbitals,str_Dets, &
       str_string)
    IMPLICIT NONE
    integer,intent(IN)      ::str_index,str_elec,Orbitals, &
         str_Dets
    integer,dimension(str_elec),intent(OUT)   :: str_string
    integer     ::i
    integer,dimension(:),ALLOCATABLE    ::temp_string
    ! Allocate temporary string
    ALLOCATE(temp_string(str_elec))
    ! Initialize
    do i=1,str_elec
       str_string(i)=i
    end do
    do i=2,str_index
       temp_string=str_string
       CALL StrFind2(temp_string,str_elec,Orbitals,str_Dets,&
            str_string)
    end do
    DEALLOCATE(temp_string)
    RETURN !Return str_string
  end subroutine GenOrbString
  !====================================================================
  !>GenOrbString_Alpha
  ! Generates alpha strings and their respective indices
  ! This subroutine holds 100 strings in the alpha string expansion.
  ! To look up and index it starts from the closest, then builds the strings
  !  from there. Note: 100 strings of the UNTRUNCATED alpha string expansion.
  ! INPUT:
  !  str_index      = index of alpha string to find
  !  str_elec       = number of alpha electrons
  !  Orbitals       = number of orbitals
  !  str_Dets       = Binom(Orbitals,str_elec)
  !  pdets_MAX      = Maximum value of alpha string in truncated expansion
  !  str_string     = resultant string (str_index)
  !--------------------------------------------------------------------
  subroutine GenOrbString_Alpha(str_index,str_elec,Orbitals,str_Dets,&
       pdets_MAX,str_string)
    IMPLICIT NONE
    integer,intent(IN)      ::str_index,str_elec,Orbitals,str_Dets,&
         pdets_MAX
    integer,dimension(str_elec),intent(OUT)   ::str_string
    integer     ::i,j
    integer     ::peg_start       !peg to start at
    integer     ::a_str_size      !number of indices kept

    if(.NOT. ALLOCATED(alpha_strings))then    !Initialize

       if(pdets_MAX.LT.100)then
          a_str_size=pdets_MAX
       else
          a_str_size=100
       end if!pdets_MAX.LT.100
       ALLOCATE(alpha_strings(0:str_elec-1,0:a_str_size-1))
       ALLOCATE(alpha_indices(0:a_str_size-1))
       str_peg_int_a=INT(pdets_MAX/a_str_size) !Always rounds down
       do i=0,a_str_size-1
          alpha_indices(i)=(i)*str_peg_int_a+1      !Set up index array
       end do
       do i=0,str_elec-1
          alpha_strings(i,0)=i+1    !Initialize first string
       end do
       do i=1,a_str_size-1
          CALL StrFind3(alpha_strings(0,i-1),alpha_indices(i-1),&
               str_elec,Orbitals,str_Dets,alpha_indices(i),    &
               alpha_strings(0,i))
       end do

    end if!if alpha_strings has not been allocated
    peg_start=MAX(0,INT(str_index/str_peg_int_a)-2)
    ! this is a bounds check
    peg_start=min(peg_start,99)
    
    if(alpha_indices(peg_start).EQ.str_index)then
       str_string(1:str_elec)=alpha_strings(0:str_elec-1,peg_start)
    else
       CALL StrFind3(alpha_strings(0,peg_start),alpha_indices(peg_start),&
            str_elec,Orbitals,str_Dets,str_index,str_string)
    end if
    RETURN
  end subroutine GenOrbString_Alpha
  !====================================================================
  !>GenOrbString_Beta
  ! Generates beta strings and their respective indices
  !--------------------------------------------------------------------
  subroutine GenOrbString_Beta(str_index,str_elec,Orbitals,str_Dets,&
       qdets_MAX,str_string)
    IMPLICIT NONE
    integer,intent(IN)      ::str_index,str_elec,Orbitals,str_Dets,&
         qdets_MAX
    integer,dimension(str_elec),intent(OUT)   ::str_string
    integer     ::i,j
    integer     ::peg_start       !peg to start at
    integer     ::b_str_size      !number of indices kept 

    if(.not. allocated(beta_strings))then    !Initialize

            if(qdets_MAX.LT.100)then
                    b_str_size=qdets_MAX
            else
                    b_str_size=100
            end if

            ALLOCATE(beta_strings(0:str_elec-1,0:b_str_size-1))
            ALLOCATE(beta_indices(0:b_str_size-1))

            str_peg_int_b=INT(qdets_MAX/b_str_size) !Always rounds down
            do i=0,b_str_size-1
                    beta_indices(i)=(i)*str_peg_int_b+1      !Set up index array
            end do
            do i=0,str_elec-1
                    beta_strings(i,0)=i+1    !Initialize first string
            end do
            do i=1,b_str_size-1
                    CALL StrFind3(beta_strings(0,i-1),beta_indices(i-1),&
                      str_elec,Orbitals,str_Dets,beta_indices(i),    &
                      beta_strings(0,i))
            end do

    end if

    peg_start=MAX(0,INT(str_index/str_peg_int_b) - 2)

    ! this is a bounds check
    peg_start=MIN(99,peg_start)

    if(beta_indices(peg_start).EQ.str_index)then
       str_string(1:str_elec)=beta_strings(0:str_elec-1,peg_start)
    else
       CALL StrFind3(beta_strings(0,peg_start),beta_indices(peg_start),&
            str_elec,Orbitals,str_Dets,str_index,str_string)
    end if
    RETURN
  end subroutine GenOrbString_Beta
END MODULE addressing
!********************************************************************
! This module contains integral location utilities
! Also contains subroutines to read in molecular integrals
!====================================================================
MODULE integral
  IMPLICIT NONE
CONTAINS
  ! readmoints
  ! ----------
  ! Read 1 and 2 electron integrals in SIFS structure
  ! Please consult SIFS documentation for further details.
  ! structure:
  !  read header 1
  !  read header 2
  !  read 1-e integrals
  !  read 2-e integrals
  ! Input:
  !  type1n   = type of integrals to read
  !  orbitals = number of orbitals in system
  !  mofile   = name of molecular orbital file
  !  m1len    = length of moints1
  !  m2len    = length of moints2
  ! Output:
  !  moints1  = 1-e integrals
  !  moints2  = 2-e integrals
  !  nuc_rep  = nuclear repulsion energy
  !  fcenergy = frozen core energy (array, size = 10)
  !--------------------------------------------------------------------
  subroutine readmoints (moints1, moints2, type1n, orbitals, mofile, &
    m1len, m2len, nuc_rep, fcenergy)
    implicit none

    ! **input
    integer,      intent(in) :: m1len, m2len, type1n, orbitals
    character*20, intent(in) :: mofile

    ! **output
    real*8, dimension(m1len), intent(out) :: moints1
    real*8, dimension(m2len), intent(out) :: moints2
    real*8,                   intent(out) :: fcenergy
    real*8,                   intent(out) :: nuc_rep
    
    ! **header 1 variables for SIFS file
    ! version = sifs version number
    ! ntitle  = number of character*80 titles
    ! nsym    = number of symmetry blocks of integrals
    ! nbft    = total number of basis functions
    ! ninfo   = number of info values to be read in
    ! nenrgy  = number of energy values to be read
    ! nmap    = number of integer mapping vectors too be read
    integer :: ntitle, nsym, nbft, ninfo, nenrgy, nmap
    integer*4 :: ntitle4, nsym4, nbft4, ninfo4, nenrgy4, nmap4
    integer :: version
    
    ! **header 2 variables for SIFS file
    ! title(1:ntitle)  = identifying titles
    ! nbpsy(1:nsym)    = number of basis functions in each symmetry block
    ! slabel(1:nsym)   = character labels for the symmetry blocks
    ! info(1:ninfo)    = additional parameters describing the integral file
    !                    info(1) = fsplit: allows separation of 1 and 2 int rec.s
    !                                     1 = same file
    !                                     2 = 2-e records are on aoints2
    !                    info(2) = l1rec: 1-e integral record length in working
    !                                     precison
    !                    info(3) = n1max: max number of 1-e integrals in a record
    !                    info(4) = l2rec: 2-e integral record length
    !                    info(5) = n2max: max number of 2-e integrals in a record
    ! bfnlab(1:nbft)   = chracter labels for the individual basis functions.
    ! ietype(1:nenrgy) = energy types
    ! energy(1:nenrgy) = core energy contributions accordint to the corresponding
    !                   ietype(*) entries
    ! imtype(1:nmap)   = map vector types
    ! map(1:bft,1:nmap)= integer mapping vectors
    character*80, dimension(:),   allocatable :: title
    character*4,  dimension(:),   allocatable :: slabel
    character*8,  dimension(:),   allocatable :: bfnlab
    integer,      dimension(:),   allocatable :: nbpsy, info, ietype, imtype
    integer,      dimension(:,:), allocatable :: map
    integer,      dimension(:),   allocatable :: mapin
    real*8,       dimension(:),   allocatable :: energy

    !**dword variables
    integer :: num, itypea, itypeb, last, ifmt, ibvtyp, lab1
    !**record variables
    real*8,  dimension(:),   allocatable :: vals, ibitv, buf
    integer, dimension(:),   allocatable :: symb
    integer, dimension(:,:), allocatable :: labels
    real*8,  dimension(:,:), allocatable :: s1h1
    real*8 :: fcore, score, hcore
    integer, parameter :: msame = 0, nmsame = 1, nomore = 2
    integer, parameter :: iretbv =-1

    integer :: nipv ! nipv = 2: 1-e integrals, nipv = 4: 2-e integrals
    integer :: ios, ierr, i, index, j, l1, l2, l3, l4
    integer :: aoints, n1int, n2int
    integer, parameter :: verin = 2
    
    aoints = 16
    ! open integral file
    call trnfln (1, mofile)
    open (unit = aoints, file = mofile, form = "unformatted", status = "old",&
      iostat = ios)
    if (ios .ne. 0) stop "*** Could not open integral file! ***"

    ! read header 1 information
    read (unit = aoints, iostat = ierr) version, ntitle, nsym, nbft, ninfo, &
      nenrgy, nmap
    if ((version - verin) .ne. 0) then
            write (*,"(1x,A)") "WARNING: MO's generated with version != 2!"
    end if
    if (ierr .ne. 0) stop "*** Error reading header 1! ***"

    n1int = Ind2Val(nbft, nbft)
    n2int = Index2E(nbft, nbft, nbft, nbft)
    
    ! allocate arrays
    allocate (title(ntitle), stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating title()! ***"
    allocate (nbpsy(nsym),   stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating nbpsy()! ***"
    allocate (slabel(nsym),  stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating slabel()! ***"
    allocate (info(ninfo),   stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating info()! ***"
    allocate (bfnlab(nbft),  stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating bfnlab()! ***"
    allocate (ietype(nenrgy),stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating ietype()! ***"
    allocate (energy(nenrgy),stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating energy()! ***"
    allocate (imtype(nmap),  stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating imtype()! ***"
    allocate (map(nbft,nmap),stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating map()! ***"
    ! read header 2 information
    call sifreadh2 (aoints, ntitle, nsym, nbft, ninfo, nenrgy, nmap, title, &
      nbpsy, slabel, info, bfnlab, ietype, energy, imtype, map, ierr )
    if (ierr .ne. 0) stop "*** Error reading header 2! ***"

    ! set nuclear repulsion energy
    nuc_rep = energy(1)
    
    ! allocate arrays
    allocate (buf(info(2)), stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating buf()!***"
    allocate (labels(2,info(3)), stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating labels()!***"
    allocate (vals(info(3)), stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating vals()!***"
    allocate (ibitv(info(3) + 63), stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating ibitv()!***"
    allocate (s1h1(n1int, 2), stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating s1h1()!***"
    allocate (symb(nbft), stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating symb()! ***"
    allocate (mapin(n1int), stat = ierr)
    if (ierr .ne. 0) stop "*** Error allocating mapin()!***"
    do i = 1, n1int
            mapin(i) = i
    end do
    
    ! read 1 electron record
    moints1 = 0d0
    nipv = 2
    fcenergy = 0d0
    
    if (type1n .eq. 1) then ! we are after H1(*)

            call sifrsh(aoints, info, buf, vals, labels, nsym, &
              nbpsy, mapin, n1int, s1h1, score, hcore, symb, ierr)
            if (ierr .ne. 0) stop "*** Error reading 1-e integrals!***"
            
            moints1(1:n1int) = s1h1(1:n1int,2)
            fcenergy = hcore

    else if (type1n .eq. 3) then ! we are after the dipole moments

            do while (last .ne. nomore)
                    call sifrd1(aoints, info, nipv, iretbv, buf, num, last, &
                      itypea, itypeb, ibvtyp, vals, labels, fcore, ibitv, ierr)
                    if (ierr .ne. 0) stop "***Error reading integral record!***"

                    if (itypea .eq. 1 .and. itypeb .eq. 0 &
                      .and. last .eq. 1) then
                            ! <x> dipole moment
                            ! collect values
                            do i = 1, num
                                    l1 = labels(1,i)
                                    l2 = labels(2,i)
                                    index = Ind2Val(l1,l2)
                                    moints1(index) = moints1(index) + vals(i)
                            end do
                    else if (itypea .eq. 1 .and. itypeb .eq. 1 &
                      .and. last .eq. 1) then
                            ! <y> dipole moment
                            ! collect values
                            do i = 1, num
                                    l1 = labels(1,i)
                                    l2 = labels(2,i)
                                    index = Ind2Val(l1,l2)
                                    index = index + m1len
                                    moints1(index) = moints1(index) + vals(i)
                            end do
                    else if (itypea .eq. 1 .and. itypeb .eq. 2 &
                      .and. last .eq. 1) then
                            ! <z> dipole moment
                            do i = 1, num
                                    l1 = labels(1,i)
                                    l2 = labels(2,i)
                                    index = Ind2Val(l1,l2)
                                    index = index + m1len
                                    index = index + m1len
                                    moints1(index) = moints1(index) + vals(i)
                            end do
                    end if

            end do
    end if

    ! deallocate arrays
    deallocate (buf, s1h1, mapin, labels, vals, ibitv)

    if (type1n .eq. 3) then
            close(aoints)
            return
    end if

    ! allocate arrays
    allocate (buf(info(4)))
    allocate (labels(4,info(5)))
    allocate (vals(info(5)))
    allocate (ibitv(info(5) + 63))
    ! read 2 electron records
    call trnfln (1, mofile)
    call sifo2f (aoints, aoints, mofile, info, i, ierr) ! i is dummy variable
    if (ierr .ne. 0) stop "*** Error positioning moints2 file! ***"
    moints2 = 0d0
    last = msame
    nipv = 4
    do while (last .ne. nomore)
            call sifreadd2 (aoints, info, ninfo, nipv, iretbv, buf, info(4), &
              num, last, itypea, itypeb, ibvtyp, vals, info(5), labels, info(5),&
              ibitv, (info(5) + 63), ierr)
            if (ierr .ne. 0) stop "*** Error reading integral record! ***"
            
            do i = 1, num
                    l1 = int(labels(1,i))
                    l2 = int(labels(2,i))
                    l3 = int(labels(3,i))
                    l4 = int(labels(4,i))
                    index = Index2E(l1,l2,l3,l4)
                    moints2(index) = moints2(index) + vals(i)
            end do
    end do

    ! deallocate arrays
    deallocate(buf, ibitv, vals, labels)
    deallocate(title, nbpsy,slabel, info, bfnlab, ietype, energy, imtype, map)

    close(aoints)
    
    return
  end subroutine readmoints
  !====================================================================
  !>sifreadd2
  ! ---------
  ! Read and decode 2-e integral record
  subroutine sifreadd2 (aoints, info, ninfo, nipv, iretbv, buf, lenbuf, &
    num, last, itypea, itypeb, ibvtyp, vals, nvals, labels, nlabels,    &
    ibitv, lenbitv, ierr)
    implicit none

    integer, intent(in) :: aoints, nipv, iretbv, ninfo, lenbuf, nlabels
    integer, intent(in) :: nvals, lenbitv
    integer, dimension(ninfo),  intent(in)       :: info
    real*8,  dimension(lenbuf), intent(inout)    :: buf
    real*8,  dimension(nvals),  intent(inout)    :: vals
    integer, dimension(4,nlabels), intent(inout) :: labels
    real*8,  dimension(lenbitv),   intent(inout) :: ibitv
    integer, intent(inout) :: num, last, itypea, itypeb, ibvtyp, ierr

    ! check value of split
    if (info(1) .ne. 1) then
            ierr = 100
            print *, "2-e records are separate. Not implemented."
            return
    end if

    ! read buffer
    read (aoints) buf

    ! unpack 2-e buffer
    call sifd2 (info, nipv, iretbv, buf, num, last, itypea, itypeb, &
      ibvtyp, vals, labels, ibitv, ierr)
    if (ierr .ne. 0) stop "*** Error unpacking 2-e buffer! ***"
    
    return

  end subroutine sifreadd2
  !====================================================================
  !>sifreadh2
  ! ---------
  ! Read header 2 of SIF formatted file
  subroutine sifreadh2 (aoints, ntitle, nsym, nbft, ninfo, nenrgy, nmap, &
    title, nbpsy, slabel, info, bfnlab, ietype, energy, imtype, map, ierr)
    implicit none

    integer, intent(in) :: aoints, ntitle, nsym, nbft, ninfo, nenrgy
    integer, intent(in) :: nmap
    integer, intent(inout) :: ierr
    character*80, dimension(ntitle),    intent(inout) :: title
    integer,      dimension(nsym),      intent(inout) :: nbpsy
    character*4,  dimension(nsym),      intent(inout) :: slabel
    integer,      dimension(ninfo),     intent(inout) :: info
    character*8,  dimension(nbft),      intent(inout) :: bfnlab
    integer,      dimension(nenrgy),    intent(inout) :: ietype
    real*8,       dimension(nenrgy),    intent(inout) :: energy
    integer,      dimension(nmap),      intent(inout) :: imtype
    integer,      dimension(nbft,nmap), intent(inout) :: map

    integer :: mapdim
    integer, parameter :: wrnerr = 0, nfterr = 1, faterr = 2
    
    ! set mapdim
    mapdim = max (nmap, 1)

    if (nmap .eq. 0) then
            read (aoints, iostat = ierr) &
              title, nbpsy, slabel, info, bfnlab, ietype, energy
    else
            read (aoints, iostat = ierr) &
              title, nbpsy, slabel, info, bfnlab, ietype, energy &
              , imtype, map
    end if

    return
  end subroutine sifreadh2
  !====================================================================
  !>Ind2Val
  ! Computes index of element i,j packed in lower triangle array
  !--------------------------------------------------------------------
  integer function Ind2Val(i,j)
    IMPLICIT NONE
    integer,intent(IN)      ::i,j
    if(i.GE.j)then
       Ind2Val=(i-1)*i/2+j
    else
       Ind2Val=(j-1)*j/2+i
    end if
  end function Ind2Val
  !====================================================================
  !>Index2E
  ! Computes index of 2-e integrals
  !--------------------------------------------------------------------
  integer function Index2E(i,j,k,l)
    IMPLICIT NONE
    integer,intent(IN)      ::i,j,k,l
    Index2E=Ind2Val(Ind2Val(i,j),Ind2Val(k,l))
  end function Index2E
  !====================================================================
END MODULE integral
MODULE search_fcns
  ! This module contains functions for searching through lists
  IMPLICIT NONE
CONTAINS
  !=====================================================================
  !>int_search
  ! Searches for an integer match in a list and locates the matching element
  ! Returns 0if element is not found
  ! ** int_list is ORDERED
  !---------------------------------------------------------------------
  subroutine int_search( int_lookfor, int_list, list_length, int_found_el )
    IMPLICIT NONE
    integer,intent(IN)      :: int_lookfor, list_length
    integer,dimension(list_length),intent(IN) :: int_list
    integer,intent(OUT)     :: int_found_el
    integer                 :: start_pt      !Where to begin searches
    integer                 :: i,j
    int_found_el=0
    
    if(MOD(list_length,2).EQ.0)then     !If list length is even
       start_pt=(list_length)/2
    else
       start_pt=(list_length+1)/2
    end if
    if(int_list(start_pt).EQ.int_lookfor)then !Match found
       int_found_el=start_pt
       RETURN !Leave
    else if(int_list(start_pt).GT.int_lookfor)then !Search below
       do i=1,start_pt-1
          if(int_list(i).EQ.int_lookfor)then
             int_found_el=i
             RETURN      !Leave
          end if
       end do
       RETURN !Leave int_found_el=0
    else
       do i=start_pt+1,list_length
          if(int_list(i).EQ.int_lookfor)then
             int_found_el=i
             RETURN      !Leave
          end if
       end do
       RETURN !Leave int_found_el=0
    end if
  end subroutine int_search
  !=====================================================================
  !>int_search2
  ! Faster integer list searcher - Calls subroutine to truncate search space
  !  subroutine being called is str_search1
  !---------------------------------------------------------------------
  subroutine int_search2(int_lookfor,int_list,list_length,int_found_el)
    implicit none
    integer,intent(IN)      ::int_lookfor,list_length
    integer,dimension(list_length),intent(in) ::int_list
    integer,intent(inout)     ::int_found_el
    integer     ::lwr_bnd,upp_bnd !upper and lower bounds (for sub call)
    integer     ::len_substr      !length of substrings
    integer     ::i
    int_found_el = 0
    ! if list_length = 1, check if int_lookfor == int_length(1)
    if ( list_length == 1 ) then
       if ( int_lookfor == int_list(1) ) then
          int_found_el = 1
          return 
       else
          int_found_el = 0
          return
       end if
    end if

    len_substr=list_length        !Initialize at full string
    lwr_bnd=1
    upp_bnd=list_length
    int_found_el=0
    do while(len_substr.GT.3)
       call str_search1(int_lookfor,int_list(lwr_bnd:upp_bnd), &
            (upp_bnd-lwr_bnd+1),lwr_bnd,upp_bnd,len_substr,   &
            int_found_el) 
    end do
    do i=lwr_bnd,upp_bnd
       if (int_list(i).eq.int_lookfor)then
          int_found_el=i
          RETURN
       end if
    end do
    RETURN
  end subroutine int_search2
  subroutine str_search1(int_lookfor,int_list,list_length,lower,upper, &
       srch_length, matching_element)
    IMPLICIT NONE
    integer,intent(IN)      ::int_lookfor,list_length
    integer,dimension(list_length),intent(in) ::int_list
    integer,intent(in out)  ::lower,upper
    integer,intent(in out)  ::srch_length,matching_element
    integer     ::start_point
    !real*8,parameter  ::golden_ratio=1.61803398875
    !start_point=int(REAL(list_length)/golden_ratio)+1       !use golden ratio to search
    start_point=list_length*2/3+1
    !if(MOD(list_length,2).EQ.0)then
    !      start_point=list_length/2
    !else
    !      start_point=(list_length+1)/2
    !end if
    
    if(int_list(start_point).EQ.int_lookfor)then
       srch_length=0
       matching_element=start_point+lower-1
       RETURN      !Leave subroutine
    else if(int_list(start_point).GT.int_lookfor)then !Search below
       upper=start_point-1+lower
       srch_length=upper-lower+1
       RETURN      !Leave subroutine with bounds adjusted for a
       ! search below the current midpoint
    else if(int_list(start_point).LT.int_lookfor)then !Search above
       lower=start_point+lower-1
       srch_length=upper-lower+1
       RETURN      !Leave subroutine with bounds adjusted for a
       ! search above the current midpoint
    end if
    RETURN
  end subroutine str_search1
  !====================================================================
END MODULE search_fcns
!********************************************************************
! This module contains string manipulation utilities
!====================================================================
MODULE string_util
  IMPLICIT NONE      
CONTAINS
  !===================================================================
  ! cannon1swp
  ! This subroutine orders a string and returns its permutation index 
  !  for a single index swap
  !--------------------------------------------------------------------
  subroutine cannon1swp(string,length,indswp,permind)
    IMPLICIT NONE
    integer,intent(IN)                        ::length,indswp
    integer,dimension(length),intent(INOUT)   ::string
    integer,intent(OUT)                       ::permind
    integer,dimension(:),ALLOCATABLE          ::tmp1
    integer     :: i,j,test
    ! Set permind to 1
    permind=1
    if(length-indswp.EQ.0)then ! The last index was changed
       ALLOCATE(tmp1(length))
       tmp1=string
       do i=length,2,-1
          if(string(i).GT.string(i-1))then
             RETURN ! Finished
          else
             permind=permind*(-1) ! Swap columns
             string(i-1)=tmp1(i)
             string(i)=tmp1(i-1)
             tmp1=string
          end if !string(i).GT.string(i-1)
       end do !i
    else  ! indswp .NE. length
       ALLOCATE(tmp1(length))
       tmp1=string
       !IF correct position is down string
       if(string(indswp).GT.string(indswp+1))then
          do i=indswp+1,length
             if(string(i-1).GT.string(i))then
                permind=permind*(-1) !Swap columns
                string(i-1)=tmp1(i)
                string(i)=tmp1(i-1)
                tmp1=string
             else
                RETURN ! Finished
             end if !string(i-1).GT.string(i)
          end do !i
       else !Correct position is up string
          do i=indswp-1,1,-1
             if(string(i+1).LT.string(i))then
                permind=permind*(-1) !Swap columns
                string(i+1)=tmp1(i)
                string(i)=tmp1(i+1)
                tmp1=string !update
             else
                RETURN ! Finished
             end if !string(i+1).LT.string(i)
          end do !i
       end if !string(indswp).GT.string(indswp+1)
    end if !length-indswp=0
    DEALLOCATE(tmp1) !Deallocate arrays
    RETURN
  end subroutine cannon1swp
  !====================================================================
  !>possex1
  ! This subroutine returns possible excitations for a string
  ! INPUT:
  !  occstring = string of occupied orbitals for alpha/beta electrons
  !  orbs = number of MO's
  !  elecs = number of alpha/beta electrons
  !  pexitslen = orbs-elecs
  !  pexits = list of possible excitations
  !--------------------------------------------------------------------
  subroutine possex1(occstr,orbs,elecs,pexitslen,pexits)
    IMPLICIT NONE
    integer,intent(IN)      ::orbs,elecs,pexitslen
    integer,dimension(elecs),intent(IN)       ::occstr
    integer,dimension(pexitslen),intent(OUT)  ::pexits
    integer,dimension(:),ALLOCATABLE          ::totorbs
    integer :: i,j,e1,e2
    ! Construct totorbs
    ALLOCATE(totorbs(orbs))
    do i=1,orbs
       totorbs(i)=i
    end do
    ! Loop through occstr making totorbs(i)=0
    do i=1,elecs
       totorbs(occstr(i))=0
    end do!i
    ! Construct pexits
    j=1
    do i=1,orbs
       if(totorbs(i).NE.0)then
          pexits(j)=totorbs(i)
          j=j+1!j->j+1
       end if
    end do!i
    RETURN!Return pexits
  end subroutine possex1
  !====================================================================
  !>cannon2
  ! This suborutine returns a string to cannonical ordering and gives
  !  the character of the permutation for 2 index swaps
  ! INPUT:
  !  string = string to be ordered (IN: ordered,noswaps; OUT: ordered,swaps)
  !  length = length of string
  !  indx1  = first index of swap
  !  indx2  = second index of swap
  !  swp1 = new value of string(indx1)
  !  swp2 = new value of string(indx2)
  !  permind = parity
  !
  !--------------------------------------------------------------------
  subroutine cannon2(string,length,indx1,swp1,indx2,swp2,permind)
    use search_fcns,  only:int_search2
    IMPLICIT NONE
    integer, intent(IN) ::length,indx1,indx2,swp1,swp2
    integer, dimension(length), intent(INOUT) ::string
    integer, intent(OUT)                      ::permind
    integer,dimension(length)   :: order
    integer,dimension(2*length) :: scratch
    integer ::i,j,k,tperm
    !Initialize
    permind=1
    do i = 1, length
            order(i) = i
    end do
    ! make first substitution
    string(indx1) = swp1

    call cannon4(string, length, indx1, scratch, order, permind)
    if (permind .eq. 0) stop "ERROR!!! permind = 0"

    ! make second substitution
    string(order(indx2)) = swp2
    call cannon4(string, length, order(indx2), scratch, order, permind)
    if (permind .eq. 0) stop "ERROR!!! permind = 0"

    RETURN!Return ordered string and parity of permutation
  end subroutine cannon2
  !====================================================================
  !>cannon3
  ! Returns cannonocalized string and parity
  !--------------------------------------------------------------------
  subroutine cannon3(string,length,index_sub,permind)
    implicit none
    integer,intent(in)                        ::length,index_sub
    integer,intent(inout)                     ::permind
    integer,dimension(length),intent(inout)   ::string
    integer,dimension(:),allocatable          ::temp
    integer                                   ::i,j,test

    !Check if index_sub is length or 1
    if(index_sub.EQ.length)then !Element to change is last
       if(string(index_sub).GT.string(index_sub-1))then !Do nothing
          RETURN
       else if(string(index_sub).LT.string(1))then !Move last element to first position
          allocate(temp(length-1))
          temp(1:length-1)=string(1:length-1)
          string(1)=string(index_sub)
          string(2:length)=temp(1:length-1)
          deallocate(temp)
          permind=permind*(-1**(length-1))
          RETURN
       else
          do i=index_sub-2,1,-1
             if(string(i).LT.string(index_sub))then
                allocate(temp(length-i-1))
                temp(1:length-i-1)=string(i+1:index_sub-1)
                string(i+1)=string(index_sub)
                string(i+2:length)=temp(1:length-i-1)
                deallocate(temp)
                permind=permind*((-1)**(length-i-1))
                RETURN
             end if
          end do
       end if
    else if(index_sub.EQ.1)then!Element to change is first
       if(string(index_sub).LT.string(index_sub+1))then !Do nothing
          RETURN
       else if(string(index_sub).GT.string(length))then !Put first element last
          allocate(temp(length-1))
          temp(1:length-1)=string(2:length)
          string(length)=string(index_sub)
          string(1:length-1)=temp(1:length-1)
          deallocate(temp)
          permind=permind*((-1)**(length-1))
          RETURN
       else
          do i=3,length
             if(string(index_sub).LT.string(i))then
                allocate(temp(1:i-2))
                temp(1:i-2)=string(2:i-1)
                string(i-1)=string(index_sub)
                string(1:i-2)=temp(1:i-2)
                deallocate(temp)
                permind=permind*((-1)**(i-2))
                RETURN
             end if
          end do
       end if
    else
       if(string(index_sub).GT.string(index_sub+1))then!Go above
          call order2('u',string(index_sub),(length-index_sub+1),permind)
          RETURN
       else if(string(index_sub).LT.string(index_sub-1))then!Go below
          call order2('d',string(1),index_sub,permind)
          RETURN
       end if
    end if
    RETURN
  end subroutine cannon3
  !=====================================================================
  !>order2
  ! ------
  ! find location and swap first or last element in a list
  subroutine order2 (direction, string, length, permind)
    implicit none

    character*1, intent(in) :: direction
    integer,     intent(in) :: length
    integer,     intent(inout) :: permind
    integer, dimension(length), intent(inout) :: string
    integer, dimension(:), allocatable :: temp
    integer :: i

    allocate(temp(length))

    ! up or down
    if (direction .eq. 'd') then ! permuting string(length)
            temp(1:(length - 1)) = string(1:(length-1))
            ! loop over temp finding where to place string(length)
            do i = (length - 1), 1, -1
                    if (string(length) .gt. temp(i)) then
                            string(1:i) = temp(1:i)
                            string(i+1) = string(length)
                            string(i+2:length) = temp(i+1:(length-1))
                            deallocate(temp)
                            permind = permind * ((-1) ** (length - i - 1))
                            return
                    end if
            end do
            ! if here, string(length) is less than string(1)
            string(1) = string(length)
            string(2:length) = temp(1:(length - 1))
            deallocate(temp)
            permind = permind * ((-1) ** (length - 1))
    else if (direction .eq. 'u') then ! permute string(1)
            temp(1:(length - 1)) = string(2:length)
            ! loop over temp finding where to place string(1)
            do i = 1, (length - 1)
                    if (string(1) .lt. temp(i)) then
                            string(i) = string(1)
                            string(1:i-1) = temp(1:i-1)
                            string(i+1:length) = temp(i:(length-1))
                            deallocate(temp)
                            permind = permind * ((-1) ** (i-1))
                            return
                    end if
            end do
            ! if here, string(1) is greater than string(length)
            string(length) = string(1)
            string(1:(length - 1)) = temp(1:(length - 1))
            deallocate(temp)
            permind = permind * ((-1) ** (length - 1))
    else
            STOP "ERROR READING PARAMETER 1 in ORDER()"
    end if
    if (allocated(temp)) deallocate(temp)
    return
  end subroutine order2
  !=====================================================================
  !>order
  ! Subroutine to find location, and swap, first or last element in a list
  !---------------------------------------------------------------------
  subroutine order(direction,string,length,permind)
    implicit none
    character(len=1),intent(IN)   ::direction
    integer,intent(IN)            ::length
    integer,intent(inout)         ::permind
    integer,dimension(length),intent(in out)  :: string
    integer,dimension(:),allocatable    ::temp
    integer     ::i
    allocate(temp(length))
    !Up or down?
    if(direction.eq.'d')then !string(length) is element to permute
            if(string(length).LT.string(1))then !Move element to first place
!          allocate(temp(1:length-1))    !for string(1:length-1)
                    temp(1:length-1)=string(1:length-1)
                    string(1)=string(length)      !perform swap
                    string(2:length)=temp(1:length-1)
                    deallocate(temp)
                    permind=permind*((-1)**(length-1))
                    RETURN
            else
                    do i=length-1,2,-1
                            if(string(length).GT.string(i))then !Found spot
                                    !                allocate(temp(1:length-i))    !for string(i+1,length-1)
                                    print *, "length = ", length
                                    print *, "i      = ", i
                                    print *, "length -i = ", (length - i)
                                    print *, "i+1 = ", (i + 1)
                                    print *, "len(temp) = ", size(temp)
                                    print *, "len(string) = ", size(string)
                                    print *, "string = ", string
                                    temp(1:length-i)=string(i+1:length-1)
                                    string(i+1)=string(length)
                                    string(i+2:length)=temp(1:length-i)
                                    deallocate(temp)
                                    permind=permind*((-1)**(length-i-1))
                                    RETURN
                            end if
                    end do!i
            end if!string(length).lt.string(1)
    else if(direction.eq.'u')then !string(1) is element to permute
            if(string(1).GT.string(length))then !Move elemnent to last place
                    !          allocate(temp(1:length-1))    !for string(2:length)
                    temp(1:length-1)=string(2:length)
                    string(length)=string(1)
                    string(1:length-1)=temp(1:length-1)
                    deallocate(temp)
                    permind=permind*((-1)**(length-1))
                    RETURN
            else
                    do i=2,length           !move up string
                            if(string(i).GT.string(1))then      !Found spot
                                    !                allocate(temp(1:i-2))         !for string(2:i-1)
                                    temp(1:i-2)=string(2:i-1)
                                    string(i-1)=string(1)
                                    string(1:i-2)=temp(1:i-2)
                                    deallocate(temp)
                                    permind=permind*((-1)**(i-2)) !(i-2) is length of temp. aka num. of switches
                                    RETURN
                            end if
                    end do!i
            end if! string(1).gt.string(length)
    else
            STOP "Error reading parameter1 in order()"
    end if
    if (allocated(temp)) deallocate(temp)
    RETURN
  end subroutine order
  !*
  !*
  subroutine cannon4(string, strlen, indx1, scratch, order, perm)
    !===========================================================================
    ! cannon4
    ! -------
    ! Purpose: order a integer string with 1 substitution
    !
    ! Input:
    !  string  = unordered string
    !  indx1   = index of 1 substitutions
    !  scratch = scratch array
    
    ! Order:
    !  order = new ordering
    !  perm  = parity of permutations
    !---------------------------------------------------------------------------
    implicit none

    ! .. input arguments ..
    integer, intent(in) :: strlen, indx1

    ! .. input/output arguments ..
    integer, dimension(strlen),   intent(inout) :: string, order
    integer, dimension(2*strlen), intent(inout) :: scratch
    ! .. output arguments ..
    integer, intent(inout) :: perm

    ! .. local scalars ..
    integer :: i, optr, htest, ltest
    
    ! set order ptr
    optr = strlen + 1
    ! zero out scratch
    scratch = 0
    
    ! bounds checks
    htest = min(strlen, (indx1 + 1))
    ltest = max(1, indx1 - 1)
    
    ! check direction to travel
    if (string(indx1) .gt. string(htest)) then
            ! string(indx1) is in wrong position, must move up
            call orderstr(1, string(indx1), order(indx1), &
              (strlen - indx1 + 1), scratch(1), scratch(optr), perm)
    else if (string(indx1) .lt. string(ltest)) then
            ! string(indx1) is in wrong position, must move down
            call orderstr(-1, string(1), order(1), indx1, &
              scratch(1), scratch(optr), perm)
    end if
    
    return
    
  contains

    subroutine orderstr(dir, str, ord, l, scr, oscr, p1)
      implicit none

      integer, intent(in) :: dir ! direction
      integer, intent(in) :: l       ! length of str and ord
      integer, dimension(l), intent(inout) :: str, ord, scr, oscr
      integer, intent(inout) :: p1   ! permutational parity

      integer :: i, htest

      oscr = ord
      scr  = str
      ! if we are moving first element up string
      if (dir .gt. 0) then
              ! test if we need only one switch
              htest = min(l, 3)
              if (l .eq. 2 .or. str(1) .lt. str(htest)) then
                      str(2) = str(1)
                      str(1) = scr(2)
                      ord(2) = oscr(1)
                      ord(1) = oscr(2)
                      p1 = p1 * (-1)
                      return
              end if
              ! find position
              do i = 3, l

                      if (str(1) .lt. str(i)) then

                              p1 = p1 * ((-1) ** (i-2))
                              str(i-1) = scr(1)
                              ord(1)   = oscr(i-1)

                              str(1:i-2) = scr(2:i-1)
                              ord(2:i-1) = oscr(1:i-2)
                              
                              return

                      else if (str(1) .eq. str(i)) then

                              p1 = 0

                              return

                      end if
              end do

              ! if here we must place first point last

              p1 = p1 * ((-1) ** (l-1))

              str(l) = str(1)
              ord(1) = ord(l)

              str(1:l-1) = scr(2:l)
              ord(2:l)   = oscr(1:l-1)
              return
      else ! we are moving down list
              ! test if we need only one switch
              htest = max(1, l-2)
              if (l .eq. 2 .or. str(l) .gt. str(htest)) then
                      str(l-1)= str(l)
                      str(l)  = scr(l-1)
                      ord(l-1) = ord(l)
                      ord(l)   = oscr(l-1)
                      p1 = p1 * (-1)
                      return
              end if
              do i = l-2, 1, -1
                      if (str(l) .gt. str(i)) then
                              p1 = p1 * ((-1) ** (l - i + 1))

                              str(i+1)  = str(l)
                              str(i+2:l)= scr(i+1:l-1)

                              ord(l)  = ord(i+1)
                              ord(i+1:l-1)= oscr(i+2: l)
              
                              return
                      else if (str(l) .eq. str(i)) then
                              p1 = 0
                              return
                      end if
              end do
              ! if here we must place last point first
              p1 = p1 * ((-1) ** (l - 1))

              str(1)   = str(l)
              str(2:l) = scr(1:l-1)

              ord(l)    = ord(1)
              ord(1:l-1)= oscr(2:l)
              return
      end if
      ! if here return error
      p1 = 0
      return
    end subroutine orderstr
  end subroutine cannon4
  subroutine strdiff(str1, str2, length, diffs)
    !===========================================================================
    ! strdiff
    ! -------
    ! Purpose: compute number of differences between two integer strings
    !
    ! Input:
    !  str1 = string1
    !  str2 = string2
    !  length = dimension of strings 1 and 2
    !
    ! Output
    !  diffs = number of differences
    !---------------------------------------------------------------------------
    use search_fcns,  only: int_search2
    implicit none

    ! .. INPUT arguments ..
    integer, intent(in) :: length
    integer, dimension(length), intent(in) :: str1, str2

    ! .. INPUT/OUTPUT arguments ..
    integer, intent(inout) :: diffs

    ! .. LOCAL scalars ..
    integer :: i, found

    ! loop over elements of str1
    diffs = 0
    do i = 1, length
            call int_search2(str1(i), str2, length, found)
            if (found .eq. 0) then
                    diffs = diffs + 1
            end if
    end do
    return
  end subroutine strdiff
END MODULE string_util
!********************************************************************
MODULE input_proc
  IMPLICIT NONE
  integer     :: aElec, bElec
CONTAINS
  !====================================================================
  !>alpha_beta
  ! Generates number of alpha and beta electrons for a system
  !--------------------------------------------------------------------
  subroutine alpha_beta(Electrons,aElec,bElec)
    IMPLICIT NONE
    integer,intent(IN)      :: Electrons
    integer,intent(OUT)     :: aElec,bElec
    ! Check if Electrons is even
    if(MOD(Electrons,2)==0)then
       aElec=Electrons/2
       bElec=Electrons/2
    else
       aElec=(Electrons+1)/2
       bElec=(Electrons-1)/2
    end if
    RETURN
  end subroutine alpha_beta
  !====================================================================
END MODULE input_proc
!********************************************************************
MODULE orthogroutines
  IMPLICIT NONE
CONTAINS
  !====================================================================
  subroutine modgramschmidt( matrix, matdim, arrdim, lda )
    
    implicit none
    
    !   ...input integer scalars...
    integer, intent(in) :: matdim, arrdim, lda
    
    !   ...input/output real*8 array...
    real*8, dimension(lda, matdim), intent(inout) :: matrix

    !   ...loop integer scalars...
    integer :: i, j, k
    
    !   ...loop real*8 scalars...
    real*8 :: overlap, ddot, norm
    !   ...loop real*8 arrays...
    real*8, dimension( :, : ),ALLOCATABLE :: scratch
    real*8, dimension( : ), ALLOCATABLE :: vec_scratch
    !--------------------------------------------------------------------
    ALLOCATE(scratch(lda,matdim))
    ALLOCATE(vec_scratch(lda))
    !  Loop over vectors
    do i=1, matdim
       norm =  sqrt( ddot( lda, matrix(1,i), 1, matrix(1,i), 1 ) )
       do j=1, lda
          vec_scratch(j) = matrix(j,i) / norm
       end do
       do j= 1, i-1
          call orthogvector( matrix, j, matdim, lda, vec_scratch )
       end do
       do j=1, lda
          matrix(j,i) = vec_scratch(j)
       end do
    end do
    
    DEALLOCATE(scratch, vec_scratch)
    
    return
    
  end subroutine modgramschmidt
  
  
  !--------------------------------------------------------------------
  ! Orthogonalizes a vector to a space
  !==================================================================== 
  !Input:
  ! ( See above)
  ! vector = vector to orthonormalize          real*8 array 1-d
  !Output:
  ! vector = now orthogonormalize              real*8 array 1-d
  !====================================================================
  subroutine orthogvector( matrix, matdim, arrdim, lda, vector )
    
    implicit none
    
    !   ...input integer scalars...
    integer, intent(in) :: matdim, arrdim, lda
    
    !   ...input real*8 arrays...
    real*8, dimension(lda, arrdim),intent(in) :: matrix
    
    !   ...input/output real*8 arrays...
    real*8, dimension(lda),intent(in out) :: vector
    
    !   ...loop integer scalars...
    integer :: i, j, k
    
    !   ...loop real*8 scalars...
    real*8 :: norm, overlap, ddot
    
    !--------------------------------------------------------------------
    
    ! Loop over vectors
    do i=1, matdim
       overlap = ddot( lda, matrix(1,i), 1, vector, 1 )
       do j=1, lda
          vector(j) = vector(j) - overlap*matrix(j,i)
       end do
    end do
    
    ! Normalize
    norm = sqrt( ddot( lda, vector, 1, vector, 1 ) )
    do i=1, lda
       vector(i) = vector(i) / norm
    end do
    
    return
  end subroutine orthogvector
END MODULE orthogroutines
!*********************************************************************










            
            
            
