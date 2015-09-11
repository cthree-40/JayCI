!*****************************************************************
!>module action_util
! This module contains utilities for performing Hv=c
!=================================================================
module action_util
  implicit none
contains
!=================================================================
!>hv_diagonals
! Computes contribution of diagonal matrix elements for c
!-----------------------------------------------------------------
  subroutine hv_diagonals( inVector, diagonals, ciDimension, & 
    outVector )
    ! Input:
    !  ciDimension = dimension of CI expansion
    !  inVector    = v
    !  diagonals   = H[i,i], for i=1, ciDimension
    ! Output:
    !  outVector   = c
    implicit none
    integer,                  intent(in)     :: ciDimension
    real*8, dimension(ciDimension), intent(in)     :: inVector, diagonals
    real*8, dimension(ciDimension), intent(in out) :: outVector
    integer :: i 
    ! Loop over diagonals
    do i=1, ciDimension
       outVector(i) = outVector(i) + diagonals(i)*inVector(i)
    end do
    return
  end subroutine hv_diagonals
!==================================================================
!>hv_alpha
! Computes contribution of excitations in alpha strings to c
!------------------------------------------------------------------
  subroutine hv_alpha( InVector, MOints1, M1Len, MOints2, M2Len,  &
       pString, pStep, pLocate, pDets, qString, qStep, qLocate,   &
       qDets, ciDim, pDetsTrunc, qDetsTrunc, aDets, bDets, aElec, &
       bElec, Orbitals, nFrozen, nDOCC, nActive, OutVector)
    ! Input:
    !  InVector = v
    !  MOints1  = 1-e integrals
    !  MOints2  = 2-e integrals
    !  pString  = beta string component of determinants ordered by 
    !              alpha strings. Ex: pString(1) => | 1,1 > => 1
    !                                 pSTring(2) => | 1,2 > => 2 .. etc
    !  qString  = see pString def.
    !  pLocate  = location of first occurance of ith alpha string in
    !             pString list
    !  qLocate  = see qString def.
    !  pStep    = number of beta strings associated with ith alpha string
    !  qStep    = see pStep def.
    !  pDets    = alpha strings in expansion
    !  qDets    = see pDets def.
    !  ciDim    = dimension of CI expansion
    !  aDets    = untruncated number of alpha strings
    !  bDets    = see aDets def.
    !  Orbitals = number of MO's
    !  nFrozen  = number of frozen orbitals
    !  nDOCC    = number of docc orbitals
    !  nActive  = number of active orbitals
    ! Output:
    !  OutVector = c
    use addressing,   only:StrFind3,GenOrbString_Alpha,GenOrbString_Beta
    use string_util,  only:possex1                                      
    use integral,     only:Index2E                                      
    use search_fcns,  only:int_search2                                  
    implicit none                                                       
    integer, intent(in)  :: M1Len, M2Len, ciDim, pDetsTrunc, qDetsTrunc, &
         aDets, bDets, aElec, bElec, Orbitals, nFrozen, nDOCC, nActive                                   
    integer, dimension(ciDim),      intent(in) :: pString, qString                
    integer, dimension(pDetsTrunc), intent(in) :: pStep, pLocate, pDets       
    integer, dimension(qDetsTrunc), intent(in) :: qStep, qLocate, qDets       
    real*8,  dimension(M1Len),      intent(in) :: MOints1                         
    real*8,  dimension(M2Len),      intent(in) :: MOints2                         
    real*8,  dimension(ciDim),      intent(in) :: InVector                        
    real*8,  dimension(ciDim),  intent(in out) :: OutVector
    type sRep ! Single replacement information linked list
       type(sRep), pointer :: next, prev
       integer             :: exIndex       ! index of |p*,q>
       integer             :: prty          ! parity returning p* to cannonical order
       integer             :: orbital1      ! excitation FROM this orbital
       integer             :: orbital2      ! excitation TO this orbital
       integer             :: detLoc        ! Location of determinant |p*,q>
    end type sRep
    type (sRep), pointer               :: sRepll, currsRep ! linked list, current node
    type (sRep), pointer               :: oldNode
    integer, dimension(:), allocatable :: aString, bString, pExct1, &
         qExct1, tmpString, new_bstring, new_astring
    integer :: i, j, k, l, m, n, o, p
    integer :: eps1, eps2, oVIndx, iVIndx, cLoc, cLoc2, dxindx1, sxindx1, sxindx2
    real*8  :: int1e1, int2e1, int2e2, int1e3, int3e2
    integer :: DeAllocStat
    ! Allocate temporary orital string. This array is used to generate a new
    ! string from a previous one.
    allocate(tmpString(aElec))
    ! allocate alpha string and beta string
    allocate(aString(aElec))
    allocate(new_astring(aElec))
    allocate(bString(bElec))
    allocate(new_bstring(bElec))
    ! allocate array of possible excitations
    allocate(pExct1(Orbitals - aElec))
    allocate(qExct1(Orbitals - bElec))
    ! Loop over alpha strings in expansion
    astrings: do i=1, pDetsTrunc
       ! generate string pDets(i)
       if ( i .ne. 1 ) then
          call StrFind3( tmpString, pDets(i-1), aElec, Orbitals, aDets, &
               pDets(i), aString )
       else
          do j=1, aElec
             aString(j) = j
          end do
       end if
       ! allocate array of single replacements in alpha strings (linked list)   
       allocate(sRepll)                                                         
       nullify(sRepll%prev)                                                    
       nullify(sRepll%next)                                                    
       currsRep => sRepll                                                       
       ! generate possible excitations
       call possex1( aString, Orbitals, aElec, (Orbitals - aElec), &
            pExct1 )
       ! loop over single excitations
       loop_elec: do j = 1, aElec
          loop_Orbitals: do k = 1, Orbitals - aElec
             ! generate information for current single replacement
             call SRepInfo( aString, aElec, pExct1(k), j, Orbitals, &
                  eps1, sxindx1, new_astring )
             cLoc = 0
             ! Test if sxindx1 is in expansion and is greater than pDets(i)
             !if ( sxindx1 .gt. pDets(i) ) then
             !   cycle loop_Orbitals
             !end if
             ! Test if sxindx1 is greater than the largest element of pDets
             if ( sxindx1 .gt. pDets(pDetsTrunc) ) then
                cycle loop_Orbitals
             end if
             ! search for determinant in string
             call int_search2( sxindx1, pDets(1), (pDetsTrunc), cLoc )
             ! if determinant is not found, cycle
             if ( cLoc .eq. 0 ) then
                cycle loop_Orbitals
             end if
             ! cLoc gives location of matching index above pDets(i) in pDets
             !cLoc = i
             ! store single replacment information
             currsRep%exIndex  = sxindx1
             currsRep%prty   = eps1
             currsRep%orbital1 = aString(j)
             currsRep%orbital2 = pExct1(k)
             currsRep%detLoc   = cLoc
             allocate(currsRep%next) ! allocate next node
             currsRep%next%prev => currsRep
             nullify(currsRep%next%next)
             currsRep => currsRep%next ! assign currsRep to next node
             ! compute first single excitation term
             call eval_singlex1( aString, pExct1, aElec, Moints1, Moints2, &
                  M1Len, M2Len, Orbitals, j, k, int1e1, int2e1 )
             ! loop over q strings
             loop_qstring1: do m = pLocate(cLoc) + 1, pLocate(cLoc) + pStep(cLoc)
                ! find q' string index in q string list
                call int_search2( pString(m), pString(pLocate(i) + 1), pStep(i), &
                     oVIndx )
                ! if string was not found
                if ( oVIndx .eq. 0 ) then
                   cycle loop_qstring1
                end if
                oVIndx = oVIndx + pLocate(i)
                ! compute contribution
                call eval_singlex2('A', pString(m), bElec, Orbitals, bDets, &
                     aString(j), pExct1(k), Moints1, Moints2, M1Len, M2Len, &
                     qDets(qDetsTrunc), int2e2 )
                OutVector(oVIndx) = OutVector(oVIndx) + eps1*(int1e1 + int2e1 &
                     + int2e2)*InVector(m)
                !OutVector(m)      = OutVector(m)      + eps1*(int1e1 + int2e1 &
                !     + int2e2)*InVector(oVIndx)
             end do loop_qstring1
             ! loop over additional replacments in alpha strings
             loop_elec2: do m = 1, j - 1
                loop_Orbitals2: do n = k + 1, Orbitals-aElec
                   ! Note: m must always be less than j
                   ! get double replacment information
                   call DRepInfo( aString, aElec, pExct1(n), m, pExct1(k), j, &
                        Orbitals, eps2, dxindx1 )
                   !call DRepInfo( aString, aElec, pExct1(k), j, pExct1(n), m, &
                   !     Orbitals, eps2, dxindx1 )
                   cLoc2 = 0
                   ! Test if dxindx1 is in expansion and greater than pDets(i)
                   !if ( dxindx1 .gt. pDets(i) ) then
                   !   cycle loop_Orbitals2
                   !end if
                   ! Test if dxindx1 is greater than largest element of pDets
                   if ( dxindx1 .gt. pDets(pDetsTrunc) ) then
                      cycle loop_Orbitals2
                   end if
                   ! search for determinant in string
                   call int_search2( dxindx1, pDets(1), (pDetsTrunc), cLoc2 )
                   if ( cLoc2 .eq. 0 ) then
                      cycle loop_Orbitals2
                   end if
                   !cLoc2 =
                   ! Compute integral
                   int1e3 = eps2*( Moints2( Index2E(aString(m),pExct1(n),aString(j), &
                        pExct1(k)) ) - Moints2( Index2E(aString(m), pExct1(k),       &
                        aString(j), pExct1(n) ) ) )
                   ! Loop over q's
                   loop_qstring2: do o = pLocate(cLoc2) + 1, pLocate(cLoc2) + pStep(cLoc2)
                      ! find index in stinrg: pLocate(i)...pLocate(i+pStep(i))
                      p = 0
                      call int_search2( pString(o), pString(pLocate(i)+1), pStep(i), p )
                      if ( p .eq. 0 ) then
                         cycle loop_qstring2
                      end if
                      iVIndx = o
                      oVIndx = p+pLocate(i)
                      OutVector(oVIndx) = OutVector(oVIndx) + int1e3*InVector(iVIndx)
                    !  OutVector(iVIndx) = OutVector(iVIndx) + int1e3*InVector(oVIndx)
                   end do loop_qstring2
                end do loop_Orbitals2
             end do loop_elec2
          end do loop_Orbitals
       end do loop_elec
       ! Remove empty linked list element
       if (associated(currsRep%prev)) then ! We are at head node
          currsRep => currsRep%prev
       end if
       if ( associated(currsRep%next) ) then ! We are at end of list
          nullify(currsRep%next)
       end if

       ! loop over single replacments stored from previous loop
       currsRep => sRepll
       singEx: do while ( associated ( currsRep ) )
          ! loop over q-strings
          loop_qstring3: do k = pLocate(i) + 1, pLocate(i) + pStep(i)
             ! generate orbital string for q
             call GenOrbString_Beta( pString(k), bElec, Orbitals, &
                  bDets, qDets(qDetsTrunc), bString )
             ! generate list of possible excitations
             call possex1( bString, Orbitals, bElec, (Orbitals-bElec), &
                  qExct1 )
             ! loop over single excitations in q strings
             loop_elec3: do l = 1, bElec
                loop_Orbitals3: do m = 1, Orbitals - bElec
                   ! Generate single replacement information
                   call SRepInfo( bString, bElec, qExct1(m), l, Orbitals, &
                        eps2, sxindx2, new_bstring )
                   ! loop over q indices that correspond to the single excitation
                   !  r(pjv*) in q*>q
                   !if ( sxindx2 .gt. pString(k) ) then
                   !   cycle loop_Orbitals3
                   !end if
                   ! Check if sxindx2 is in expansion
                   if ( sxindx2 .gt. pString(pLocate(i)+pStep(i)) ) then
                      cycle loop_Orbitals3
                   end if
                   cLoc2=0
                   ! find determinant location in pString(sRepll%detLoc)...
                   !  ...pString(sRepll%detLoc + pStep(sRepll%detLoc))
                   call int_search2( sxindx2, pString(pLocate(currsRep%detLoc) + 1 ), &
                        pStep(currsRep%detLoc), cLoc2 )
                   if ( cLoc2 .eq. 0 ) then
                      cycle loop_Orbitals3
                   end if
                   ! Compute contribution
                   int3e2 = Moints2( Index2E( currsRep%orbital1, currsRep%orbital2, &
                     bString(l), qExct1(m) ) )
                   iVIndx = pLocate( currsRep%detLoc ) + cLoc2
                   oVIndx = k
                   OutVector(oVIndx) = OutVector(oVIndx) + int3e2 * currsRep%prty * eps2 * &
                        InVector(iVIndx)
                   !OutVector(iVIndx) = OutVector(iVIndx) + int3e2 * currsRep%prty * eps2 * &
                   !     InVector(oVIndx)
                end do loop_Orbitals3
             end do loop_elec3
          end do loop_qstring3
          if ( associated(currsRep%next) ) then
             currsRep => currsRep%next
          else
             exit singEx
          end if
       end do singEx ! while associated
       ! Move aString to tmpSTring for generation of next string at start of loop
       tmpString = aString
       ! deallocate linked list            
       currsRep => sRepll
       do while( associated(currsRep%next) )
          currsRep => currsRep%next
          sRepll => currsRep%prev
          deallocate(sRepll)
       end do
       sRepll => currsRep
       deallocate(sRepll)
    end do astrings

    ! deallocate arrays
    deallocate( aString, bString, pExct1, qExct1, tmpString, new_astring )
    return
  end subroutine hv_alpha
!====================================================================
!>hv_beta
! Computes beta string contributions to Hv=c
!--------------------------------------------------------------------
  subroutine hv_beta( InVector, MOints1, M1Len, MOints2, M2Len,  &
       pString, pStep, pLocate, pDets, qString, qStep, qLocate,  &
       qDets, ciDim, pDetsTrunc, qDetsTrunc, aDets, bDets,       &
       aElec, bElec, Orbitals, nFrozen, nDOCC, nActive, XRefList,&
       OutVector )
    ! Input:
    !  (see hv_alpha)
    !  XRefList = mapping of determinants to alpha-ordered list
    ! Output:
    !  (see hv_alpha)
    use addressing,   only: StrFind3, GenOrbString_Alpha, &
                                      GenOrbString_Beta
    use string_util,  only: possex1
    use integral,     only: Index2E
    use search_fcns,  only: int_search2
    implicit none
    integer, intent(in)  :: M1Len, M2Len, ciDim, pDetsTrunc, qDetsTrunc, &    
         aDets, bDets, aElec, bElec, Orbitals, nFrozen, nDOCC, nActive        
    integer, dimension(ciDim),      intent(in) :: pString, qString, XRefList            
    integer, dimension(pDetsTrunc), intent(in) :: pStep, pLocate, pDets       
    integer, dimension(qDetsTrunc), intent(in) :: qStep, qLocate, qDets       
    real*8,  dimension(M1Len),      intent(in) :: MOints1                     
    real*8,  dimension(M2Len),      intent(in) :: MOints2                     
    real*8,  dimension(ciDim),      intent(in) :: InVector                    
    real*8,  dimension(ciDim),  intent(in out) :: OutVector                   
    integer, dimension(:), allocatable :: aString, bString, pExct1, &  
     qExct1, tmpString, new_bString                                             
    integer :: i, j, k, l, m, n, o, p
    integer :: eps1, eps2, oVIndx, iVIndx, sxindx1, sxindx2, cLoc, cLoc2, dxindx1
    real*8  :: int1e1, int2e1, int1e2, int2e2, int1e3
    ! Allocate temporary orbital string. This array will be used to generate
    ! a new string from a previoius one
    allocate(tmpString(bElec))
    ! allocate alpha string and beta string
    allocate(aString(aElec))
    allocate(bString(bElec))
    allocate(new_bString(bElec))
    ! allocate array of possible excitations
    allocate(pExct1(Orbitals - aElec))
    allocate(qExct1(Orbitals - bElec))

    ! loop over beta strings
    bstrings: do i=1, qDetsTrunc
       ! generate string qDets(i)
       if ( i .ne. 1 ) then
          call StrFind3( tmpString, qDets(i-1), bElec, Orbitals, bDets, &
               qDets(i), bString )
       else
          do j=1, bElec
             bString(j) = j
          end do
       end if
       ! generate possible excitations
       call possex1( bString, Orbitals, bElec, (Orbitals - bElec), &
            qExct1 )
       ! loop over single excitations
       ! Allocate replacement string new_bstring
       loop_elec: do j = 1, bElec
          loop_Orbitals: do k = 1, Orbitals - bElec
             ! generate information for current replacement
             call SRepInfo( bString, bElec, qExct1(k), j, Orbitals, eps1, &
                  sxindx1, new_bstring )
             ! Test if sxindx1 is in expansion and is greater than qDets(i)
             if ( sxindx1 .gt. qDets(i) ) then
                cycle loop_Orbitals
             end if
             ! Test if sxindx1 is within expansion
             if ( sxindx1 .gt. qDets(qDetsTrunc) ) then
                cycle loop_Orbitals
             end if
             ! search for determinant in string
             cLoc = 0
             call int_search2( sxindx1, qDets(1), (qDetsTrunc), cLoc )
             ! if determinant is not found, cycle
             if ( cLoc .eq. 0 ) then
                cycle loop_Orbitals
             end if
             ! cLoc gives location of matching index above qDets(i) in qDets
             cLoc = cLoc + 0
             ! compute first single excitation term
             call eval_singlex1( bString, qExct1, bElec, Moints1, Moints2, &
                  M1Len, M2Len, Orbitals, j, k, int1e1, int2e1 )
             ! loop over p strings
             loop_pstring1: do m = qLocate(cLoc) + 1, qLocate(cLoc) + qStep(cLoc)
                     ! find p' string index in p string list
                call int_search2( qString(m), qString(qLocate(i) + 1), qStep(i), &
                     oVIndx )
                ! if string was not found
                if ( oVIndx .eq. 0 ) then
                   cycle loop_pstring1
                end if
                oVIndx = XRefList(oVIndx + qLocate(i))
                iVIndx = XRefList(m)
                ! compute contribution
                call eval_singlex2( 'B', qString(m), aElec, Orbitals, aDets, &
                     bString(j), qExct1(k), Moints1, Moints2, M1Len, M2Len,  &
                     pDets(pDetsTrunc), int2e2 )
                OutVector(oVIndx) = OutVector(oVIndx) + eps1*(int1e1 + int2e1 &
                     + int2e2)*InVector(iVIndx)                                     
                OutVector(iVIndx)      = OutVector(iVIndx) + eps1*(int1e1 + int2e1 &
                     + int2e2)*InVector(oVIndx)
             end do loop_pstring1
             ! loop over additional replacements in beta strings
             loop_elec2: do m = 1 , j - 1
                loop_Orbitals2: do n=k+1, Orbitals-bElec !
                   ! Note: m must always be less than j
                   ! get double replacement information
                   call DRepInfo( bString, bElec, qExct1(n), m, qExct1(k), j, &
                     Orbitals, eps2, dxindx1 )
!                   call DRepInfo( bString, bElec, qExct1(k), j, qExct1(n), m, &
!                        Orbitals, eps2, dxindx1 )

                   cLoc2 = 0
                   ! test if dxindx1 is in expansion and greater than qDets(i)
                   if ( dxindx1 .gt. qDets(i) ) then
                      cycle loop_Orbitals2
                   end if
                   ! test if dxindx1 is within expansion
                   if ( dxindx1 .gt. qDets(qDetsTrunc) ) then
                      cycle loop_Orbitals2
                   end if
                   ! search for determinant in string qDets
                   call int_search2( dxindx1, qDets(1), qDetsTrunc, cLoc2 )
                   if ( cLoc2 .eq. 0 ) then
                      cycle loop_Orbitals2
                   end if
                   cLoc2 = cLoc2
                   ! Compute integral : ( m, k, j, n ) - ( m, n, j, k )
                   int1e3 = eps2 * ( Moints2( Index2E(bString(m),qExct1(n),bString(j), &
                        qExct1(k)) ) - Moints2( Index2E(bString(m), qExct1(k),         &
                        bString(j), qExct1(n) ) ) )
                   ! loop over p's
                   loop_pstring2: do o = qLocate(cLoc2) + 1, qLocate(cLoc2) + qStep(cLoc2)
                      ! find index in string: qLocate(i) .... qLocate (i + qStep(i) )
                      p = 0
                      call int_search2( qString(o), qString(qLocate(i)+1), qStep(i), p )
                      if ( p .eq. 0 ) then
                         cycle loop_pstring2
                      end if
                      iVIndx = XRefList(o)
                      oVIndx = XRefList(p + qLocate(i))
                      OutVector(oVIndx) = OutVector(oVIndx) + int1e3*InVector(iVIndx)
                      OutVector(iVIndx) = OutVector(iVIndx) + int1e3*InVector(oVIndx)
                      
                   end do loop_pstring2
                end do loop_Orbitals2
             end do loop_elec2
          end do loop_Orbitals
       end do loop_elec
       tmpString = bString
    end do bstrings
    ! deallocate arrays
    deallocate( bString, aString, pExct1, qExct1,tmpString, new_bstring )
    return
  end subroutine hv_beta
!===========================================================================
!====================================================================      
!>SRepInfo                                                                 
! Returns single excitation information                                    
! INPUT:                                                                   
!  String   =string with no excitations                                    
!  SLen     =length of String                                              
!  NewOrb   =orbital exciting into                                         
!  IndxSwp  =electron involved in excitation                               
!  Orbitals =number of MO's                                                
!  Eps1     =parity of excitation                                          
!  Indx     =New string address                                            
!--------------------------------------------------------------------      
  subroutine SRepInfo(String,SLen,NewOrb,IndxSwp,Orbitals,Eps1,Indx, &
       newString )   
    use addressing, only: AdrFind                                  
    use string_util, only: cannon1swp                              
    IMPLICIT NONE                                                  
    integer,intent(IN)      ::SLen,NewOrb,IndxSwp,Orbitals         
    integer,dimension(SLen),intent(IN)    ::String 
    integer,intent(OUT)                   ::Eps1,Indx                            
    integer,dimension(SLen),intent(INOUT) ::newString                
    !Allocate newString                   
    newString=String                                               
    !Swap orbitals                                                 
    newString(IndxSwp)=NewOrb                                      
    !Return to cannonical ordering                                 
    CALL cannon1swp(newString,SLen,IndxSwp,Eps1)                   
    !Locate address                                                
    CALL AdrFind(newString,SLen,Orbitals,Indx)
    RETURN                                                         
  end subroutine SRepInfo
!======================================================================
!====================================================================           
!>DRepInfo                                                                      
! Returns double excitation information                                         
! INPUT: (see above)                                                            
!--------------------------------------------------------------------           
  subroutine DRepInfo(String,SLen,NewOrb1,IndxSwp1,NewOrb2,IndxSwp2,&       
       Orbitals,Eps1,Indx)                                           
    use addressing, only: AdrFind                                       
    use string_util, only: cannon2                                      
    IMPLICIT NONE                                                       
    integer,intent(IN)      ::SLen,NewOrb1,IndxSwp1,NewOrb2,IndxSwp2,&  
         Orbitals                                                      
    integer,dimension(SLen),intent(IN)  ::String                        
    integer,intent(OUT)     ::Eps1,Indx                                 
    integer,dimension(:),ALLOCATABLE    ::newString                     
    !Allocate newString                                                 
    ALLOCATE(newString(SLen))                                           
    newString=String
    !Orbitals are swapped in subroutine cannon2                         
    CALL cannon2(newString,SLen,IndxSwp1,NewOrb1,IndxSwp2,NewOrb2,&     
         Eps1)                                                         
    !Find address
    CALL AdrFind(newString,SLen,Orbitals,Indx)                          
    DEALLOCATE(newString)                                               
    RETURN                                                              
  end subroutine DRepInfo
!====================================================================               
!>diag_element1                                                                     
! Computes diagonal element < k | H | k >                                           
!--------------------------------------------------------------------               
  double precision function diag_element1(index_K,MOints1,M1Len, &              
       MOints2,M2Len,aElec,bElec,Orbitals)                               
    use addressing,         only: StrFind1, K2Indc,GenOrbString             
    use combinatorial,      only: Binom                                     
    use construct,          only: Diagonal_Mat                              
    IMPLICIT NONE                                                           
    integer,intent(IN)      ::index_K,M1Len,M2Len,aElec,bElec, &            
         Orbitals                                                          
    real*8,dimension(M1Len),intent(IN)  ::MOints1                           
    real*8,dimension(M2Len),intent(IN)  ::MOints2                           
    integer     ::total_aDets,total_bDets,alpha_index,beta_index            
    integer,dimension(:),ALLOCATABLE    ::alpha_string,beta_string          
    !Allocate arrays                                                        
    ALLOCATE(alpha_string(aElec))                                           
    ALLOCATE(beta_string(bElec))                                            
    !Compute total number of alpha and beta determinants                    
    total_aDets=Binom(Orbitals,aElec)                                       
    total_bDets=Binom(Orbitals,bElec)                                       
    !Find determinant string indices for index_K                            
    CALL K2Indc(index_K,bElec,Orbitals,alpha_index,beta_index)              
    !Generate orbital strings for alpha_index and beta_index                
    CALL GenOrbString(alpha_index,aElec,Orbitals,total_aDets,alpha_string)  
    CALL GenOrbString(beta_index,bElec,Orbitals,total_bDets,beta_string)    
    !Compute diagonal contribution                                          
    diag_element1=Diagonal_Mat(alpha_string,beta_string,aElec,bElec, &      
         MOints1,M1Len,MOints2,M2Len)                                      
    DEALLOCATE(alpha_string,beta_string)                                    
  end function diag_element1
!====================================================================                             
!>eval_singlex1                                                                                   
! Evaluates single excitation matrix elements                                                     
!--------------------------------------------------------------------                             
  subroutine eval_singlex1(sp_string,sp_excites,sp_elecs,MOints1, &                           
       MOints2,M1Len,M2Len,Orbitals,bound_elec,excited_elec,contr_1elec, &                   
       contr_2elec)                                                                          
    use integral,         only: Ind2Val, Index2E
    IMPLICIT NONE                                                                         
    integer,intent(IN)      ::sp_elecs,M1Len,M2Len,Orbitals,bound_elec,&                  
         excited_elec                                                                    
    integer,dimension(sp_elecs),intent(IN)    ::sp_string                                 
    integer,dimension(Orbitals-sp_elecs),intent(IN)  ::sp_excites                         
    real*8,dimension(M1Len),intent(IN)        ::MOints1                                   
    real*8,dimension(M2Len),intent(IN)        ::MOints2                                   
    real*8,intent(OUT)                        ::contr_1elec,contr_2elec                   
    integer     ::i,j,k                                                                   
!    integer,dimension(:),ALLOCATABLE    ::tmp_string                                      
!    ALLOCATE(tmp_string(sp_elecs))                                                        
    !One electron contribution                                                            
    contr_1elec=MOints1(Ind2Val(sp_string(bound_elec),sp_excites(excited_elec)))          
    !Two electron contribution                                                            
!    tmp_string=sp_string                                                                  
!    tmp_string(bound_elec)=sp_excites(excited_elec)                                       
    !Compute contribution                                                                 
    contr_2elec=0d0                                                                       
    do i=1, bound_elec-1                                                                  
!       contr_2elec=contr_2elec+MOints2(Index2E(sp_string(bound_elec),                   &
!            sp_excites(excited_elec),sp_string(i),sp_string(i))) - &          
!            MOints2(Index2E(sp_string(bound_elec),sp_string(i),sp_string(i),            &        
!            sp_excites(excited_elec)))                                                
       contr_2elec = contr_2elec + MOints2(Index2E(sp_string(i),sp_string(i),   &
            sp_string(bound_elec),sp_excites(excited_elec))) -                  &
            MOints2(Index2E(sp_string(i),sp_string(bound_elec),sp_string(i),    &
            sp_excites(excited_elec)))
    end do
    do i=bound_elec+1,sp_elecs                                                            
!       contr_2elec=contr_2elec+MOints2(Index2E(sp_string(bound_elec),                   &
!            sp_excites(excited_elec),sp_string(i),sp_string(i))) - &          
!            MOints2(Index2E(sp_string(bound_elec),sp_string(i),sp_string(i),            &        
!            sp_excites(excited_elec)))                                                
       contr_2elec = contr_2elec + MOints2(Index2E(sp_string(i),sp_string(i),   &
            sp_string(bound_elec),sp_excites(excited_elec))) -                  &
            MOints2(Index2E(sp_string(i),sp_string(bound_elec),sp_string(i),    &
            sp_excites(excited_elec)))
    end do
!    DEALLOCATE(tmp_string)                                                                
    RETURN !Return values                                                                 
  end subroutine eval_singlex1
!=====================================================================                              
!>eval_singlex2                                                                                     
! Computes contribution of third integral of single excitations                                     
!---------------------------------------------------------------------                              
  subroutine eval_singlex2( sp, sp2index, sp2elec, orbitals, sp2_dets, sp1orbital_bnd,&         
       sp1orbital_ex , moints1, moints2, moints1len, moints2len, sp2_MAX, cont_2elec)          
    !Input:                                                                                 
    ! sp             = beta/alpha                                                           
    ! sp2index       = beta/alpha string index                                              
    ! sp2elec        = beta/alpha electron number                                           
    ! orbitals       = # of MO's                                                            
    ! sp2_dets       = non-truncated beta/alpha determinants                                
    ! sp2_string     = beta/alpha string                                                    
    ! sp1orbital_bnd = orbital of electron being excited FROM                               
    ! sp1orbital_ex  = orbital electron is being excited TO                                 
    ! cont_2elec     = 2-e contribution                                                     
    use addressing,   only: GenOrbString_Alpha,GenOrbString_Beta                            
    use integral,      only: Index2E                                                        
    IMPLICIT NONE                                                                           
    character(len=1),intent(IN)   :: sp                                                     
    integer, intent(in) :: sp2index, sp2elec, orbitals, sp2_dets, sp1orbital_bnd, &         
         sp1orbital_ex, moints1len, moints2len,sp2_MAX                            
    real*8, dimension(moints1len), intent(in) :: moints1                                    
    real*8, dimension(moints2len), intent(in) :: moints2                                    
    real*8, intent(out) :: cont_2elec                                                       
    integer, dimension(:),ALLOCATABLE :: sp2_string                                         
    integer :: i                                                                            
    !----------------------------------------------------------------                       
    ! Generate orbital index string                                                         
    ALLOCATE(sp2_string(sp2elec))                                                           
    if(sp.EQ.'A')then                                                                       
       call GenOrbString_Beta( sp2index, sp2elec, orbitals, sp2_dets, sp2_MAX, &         
            sp2_string )                                                                
    else                                                                                    
       call GenOrbString_Alpha(  sp2index, sp2elec, orbitals, sp2_dets, sp2_MAX, &       
            sp2_string )                                                                
    end if
    ! Compute contribution                                                                
    cont_2elec = 0d0                                                                      
    do i=1, sp2elec                                                                       
       cont_2elec = cont_2elec + moints2( Index2E( sp2_string(i), sp2_string(i), &     
            sp1orbital_bnd, sp1orbital_ex ))                                    
    end do
    DEALLOCATE(sp2_string)                                                                
    RETURN !Return                                                                        
  end subroutine eval_singlex2
                                                                                                  
END MODULE action_util                                                                            
