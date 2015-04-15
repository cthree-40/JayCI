!********************************************************************
! This module contains combinatorial utilities used in JayCI
!====================================================================
MODULE combinatorial
  IMPLICIT NONE
  PUBLIC      :: Binom
  integer,dimension(:,:),ALLOCATABLE  :: Binom_Data
  integer     :: Binom_Data_MaxN,Binom_Data_MaxK
CONTAINS
  !=====================================================================
  !>Factl(m,n) returns the product of m-->n
  !---------------------------------------------------------------------
  integer function Factl(m,n)
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
  end function Factl
  !====================================================================
  !>Binom(n,k) returns binomial coefficient of n and k
  !--------------------------------------------------------------------
  integer function Binom(n,k)
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
    Binom=Binom_Data(n,k)
    !else if( n.GT.Binom_Data_MaxN .AND. k.LE.Binom_Data_MaxK)then
    !      ALLOCATE(tmp_BD(0:Binom_Data_MaxN,0:Binom_Data_maxK))
    !      tmp_BD=Binom_Data
    !      DEALLOCATE(Binom_Data)
    !      ALLOCATE(Binom_Data(0:n,0:Binom_Data_MaxK))
    !      Binom_Data(0:Binom_Data_maxN,0:Binom_Data_MaxK)=tmp_BD
    !      DEALLOCATE(tmp_BD)
    !      do i=Binom_Data_MaxN+1,n
    !            do j=0, MIN(i,Binom_Data_MaxK)
    !                  if(j.EQ.0 .OR. j.EQ.i)then
    !                        Binom_Data(i,j)=1
    !                  else
    !                        Binom_Data(i,j)=Binom_Data(i-1,j-1) + &
    !                              Binom_Data(i-1,j)
    !                  end if
    !            end do
    !      end do
    !      Binom=Binom_Data(n,k)
    !      Binom_Data_MaxN=n
    !else if(k.GT.Binom_Data_MaxK .AND. n.LE.Binom_Data_MaxN)then
    !      ALLOCATE(tmp_BD(0:Binom_Data_MaxN,0:Binom_Data_MaxK))
    !      tmp_BD=Binom_Data
    !      DEALLOCATE(Binom_Data)
    !      ALLOCATE(Binom_Data(0:Binom_Data_MaxN,0:k))
    !      Binom_Data(0:Binom_Data_MaxN,0:Binom_Data_MaxK)=tmp_BD
    !      DEALLOCATE(tmp_BD) !Deallocate temporary array
    !      do i=0,Binom_Data_MaxN
    !            do j=Binom_Data_MaxK+1,MIN(i,k)
    !                  if(j.EQ.0 .OR. j.EQ.i)then
    !                        Binom_Data(i,j)=1
    !                  else
    !                        Binom_Data(i,j)=Binom_Data(i-1,j-1)+ &
    !                              Binom_Data(i-1,j)
    !                  end if
    !            end do
    !      end do
    !      Binom=Binom_Data(n,k)
    !      Binom_Data_MaxK=k
    !else if(n.GT.Binom_Data_MaxN .AND. k.GT.Binom_Data_MaxK)then
    !      !Allocate temporary array
    !      ALLOCATE(tmp_BD(0:Binom_Data_MaxN,0:Binom_Data_MaxK))
    !      tmp_BD=Binom_Data             !Hold previous values of Binom_Data in tmp array
    !      DEALLOCATE(Binom_Data)        !Deallocate Binom_Data
    !      ALLOCATE(Binom_Data(0:n,0:k)) !Reallocate Binom_Data
    !      Binom_Data(0:Binom_Data_MaxN,0:Binom_Data_MaxK)=tmp_BD !Transfer prev. values
    !      DEALLOCATE(tmp_BD)            !Deallocate temporary array
    !      !n first
    !      do i=Binom_Data_MaxN+1,n
    !            do j=0, MIN(i,Binom_Data_MaxK)
    !                  if(j.EQ.0 .OR. j.EQ.i)then
    !                        Binom_Data(i,j)=1
    !                  else
    !                        Binom_Data(i,j)=Binom_Data(i-1,j-1) + &
    !                              Binom_Data(i-1,j)
    !                  end if
    !            end do
    !      end do
    !      Binom_Data_MaxN=n
    !      !k second
    !      do i=0,n
    !            do j=Binom_Data_MaxK+1,MIN(i,k)
    !                  if(j.EQ.0 .OR. j.EQ.i)then
    !                        Binom_Data(i,j)=1
    !                  else
    !                        Binom_Data(i,j)=Binom_Data(i-1,j-1) + &
    !                              Binom_Data(i-1,j)
    !                  end if
    !            end do
    !      end do
    !      Binom_Data_MaxK=k
    !      Binom=Binom_Data(n,k)
    !end if
  end function Binom
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
       print *, "    Initializing alpha strings array "
       if(pdets_MAX.LT.100)then
          a_str_size=10
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
       print *, "    Finished initializing alpha strings array."
    end if!if alpha_strings has not been allocated
    peg_start=MAX(0,INT(str_index/str_peg_int_a)-1)
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
    if(.NOT. ALLOCATED(beta_strings))then    !Initialize
       print *, "    Initializing beta strings array "
       if(qdets_MAX.LT.100)then
          b_str_size=10
       else
          b_str_size=100
       end if!pdets_MAX.LT.100
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
       print *, "    Finished initializing beta strings array "
    end if!if alpha_strings has not been allocated
    peg_start=MAX(0,INT(str_index/str_peg_int_b)-1)
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
!====================================================================
MODULE integral
  IMPLICIT NONE
CONTAINS
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
    integer,intent(out)     ::int_found_el
    integer     ::lwr_bnd,upp_bnd !upper and lower bounds (for sub call)
    integer     ::len_substr      !length of substrings
    integer     ::i
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
  !====================================================================
  !>cannon1swp
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
    integer,dimension(:),allocatable          ::temp_list
    integer ::i,j,k,tperm
    !Initialize
    permind=1
    ! Make subsitutions
    if(swp1.gt.swp2)then    !Swap subsituting orbitals
       permind=permind*(-1)
       string(indx1)=swp2
       string(indx2)=swp1
    else
       string(indx1)=swp1
       string(indx2)=swp2
    end if
    ! Order second swap
    call cannon3(string,length,indx2,tperm)
    permind=permind*tperm
    ! Order first swap
    call cannon3(string,length,indx1,tperm)
    permind=permind*tperm
    
    RETURN!Return ordered string and parity of permutation
  end subroutine cannon2
  !====================================================================
  !>cannon3
  ! Returns cannonocalized string and parity
  !--------------------------------------------------------------------
  subroutine cannon3(string,length,index_sub,permind)
    implicit none
    integer,intent(in)                        ::length,index_sub
    integer,intent(out)                       ::permind
    integer,dimension(length),intent(in out)  ::string
    integer,dimension(:),allocatable          ::temp
    integer                                   ::i,j,test
    permind=1
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
          call order('u',string(index_sub),(length-index_sub+1),permind)
          RETURN
       else if(string(index_sub).LT.string(index_sub-1))then!Go below
          call order('d',string(1),index_sub,permind)
          RETURN
       end if
    end if
    RETURN
  end subroutine cannon3
  !=====================================================================
  !>order
  ! Subroutine to find location, and swap, first or last element in a list
  !---------------------------------------------------------------------
  subroutine order(direction,string,length,permind)
    implicit none
    character(len=1),intent(IN)   ::direction
    integer,intent(IN)            ::length
    integer,intent(in out)        ::permind
    integer,dimension(length),intent(in out)  :: string
    integer,dimension(:),allocatable    ::temp
    integer     ::i
    !Up or down?
    if(direction.eq.'d')then !string(length) is element to permute
       if(string(length).LT.string(1))then !Move element to first place
          allocate(temp(1:length-1))    !for string(1:length-1)
          temp(1:length-1)=string(1:length-1)
          string(1)=string(length)      !perform swap
          string(2:length)=temp(1:length-1)
          deallocate(temp)
          permind=permind*((-1)**(length-1))
          RETURN
       else
          do i=length-1,2,-1
             if(string(length).GT.string(i))then !Found spot
                allocate(temp(1:length-i))    !for string(i+1,length-1)
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
          allocate(temp(1:length-1))    !for string(2:length)
          temp(1:length-1)=string(2:length)
          string(length)=string(1)
          string(1:length-1)=temp(1:length-1)
          deallocate(temp)
          permind=permind*((-1)**(length-1))
          RETURN
       else
          do i=2,length           !move up string
             if(string(i).GT.string(1))then      !Found spot
                allocate(temp(1:i-2))         !for string(2:i-1)
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
    RETURN
  end subroutine order
  
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










            
            
            
