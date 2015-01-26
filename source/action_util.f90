!********************************************************************
!>MODULE action_util
! This module contains utilities for performing Hv=c
!====================================================================
MODULE action_util
      IMPLICIT NONE
CONTAINS
!====================================================================
!>hv_diagonals
! Computes contribution of diagonal matrix elements to c
!--------------------------------------------------------------------
  subroutine hv_diagonals(InVector,Diagonals,ciDim,OutVector)
    IMPLICIT NONE
    integer,intent(IN)      ::ciDim
    real*8,dimension(ciDim),intent(IN)   ::InVector,Diagonals
    real*8,dimension(ciDim),intent(INOUT)::OutVector
    integer     ::i
    do i=1,ciDim
       OutVector(i)=OutVector(i)+Diagonals(i)*InVector(i)
    end do
    RETURN
  end subroutine hv_diagonals
!====================================================================
!>hv_alpha
! Computes contributions of alpha strings to c
! INPUT:
!  InVector = v
!  Moints1(M1Len) = 1-e integral array
!  Moints2(M2Len) = 2-e integeral array
!  pString(ciDim,2) = Determinants ordered by alpha strings ([p,q]: columns)
!  qString(ciDim,2) = Determinents ordered by beta strings  ([p,q]: rows)
!  pLocate(pDetstrunc) = Location of alpha strings in pString
!  pStep(pDetstrunc)   = Number of q-strings corresponding to current p-string
!  pDets(pDetstrunc)   = p-strings
!
!--------------------------------------------------------------------
  subroutine hv_alpha(InVector,MOints1,M1Len,MOints2,M2Len,pString, &
       pStep,pLocate,pDets,qString,qStep,qLocate,qDets,ciDim,      &
       pDetsTrunc,qDetsTrunc,aDets,bDets,aElec,bElec,Orbitals,     &
       nFrozen,nDOCC,nActive,OutVector)
    use addressing,   only:StrFind3,GenOrbString_Alpha,GenOrbString_Beta
    use string_util,  only:possex1
    use integral,     only:Index2E
    use search_fcns,  only:int_search2
    IMPLICIT NONE
    integer,intent(IN)      ::M1Len,M2Len,ciDim,pDetsTrunc,     &
         qDetsTrunc,aDets,bDets,aElec,bElec,Orbitals,  &
         nFrozen,nDOCC,nActive
    integer,dimension(ciDim),intent(IN)::pString,qString
    integer,dimension(pDetsTrunc),intent(IN)::pStep,pLocate,pDets
    integer,dimension(qDetsTrunc),intent(IN)::qStep,qLocate,qDets
    real*8,dimension(M1Len),intent(IN)::MOints1
    real*8,dimension(M2Len),intent(IN)::MOints2
    real*8,dimension(ciDim),intent(IN)::InVector
    real*8,dimension(ciDim),intent(INOUT)::OutVector
    integer,dimension(:,:),ALLOCATABLE::SingleRepInfo
    integer,dimension(:),ALLOCATABLE::aString,bString,pExct1,qExct1,&
         tmpString
    real*8::int1e1,int2e1,int2e2,int1e3,int2e3,int3e2
    integer ::i,j,k,l,m,n,o,p
    integer ::SRILen,cLoc,cLoc2
    integer ::eps1,eps2,sxindx1,dxindx1,sRIndx,iVIndx,oVIndx,sxindx2
    !Allocate temporary orbital string. This is used to generate a new
    ! string from a previous one
    ALLOCATE(tmpString(aElec))
    !tmpString=0
    !Allocate alpha string and beta string
    ALLOCATE(aString(aElec))
    ALLOCATE(bString(bElec))
    !Allocate array of possible excitations
    ALLOCATE(pExct1(Orbitals-aElec))
    ALLOCATE(qExct1(Orbitals-belec))
    !Allocate array of single replacements in alpha strings
    SRILen=(Orbitals-aElec)*(aElec-nFrozen)
    ALLOCATE(SingleRepInfo(5,SRILen)) !Index;
    !Loop over alpha strings in expansion
    do i=1,pDetsTrunc
       !Generate string pDets(i)
       if (i.NE.1)then
          CALL StrFind3(tmpString,pDets(i-1),aElec,Orbitals,aDets, &
               pDets(i),aString)
          
       else
          do j=1,aElec
             aString(j)=j
          end do
       end if
       !Generate possible excitations
       CALL possex1(aString,Orbitals,aElec,(Orbitals-aElec),pExct1)
       !Loop over single excitations
       sRIndx=1
       do j=nFrozen+1,aElec
          do k=1,Orbitals-aelec
             CALL SRepInfo(aString,aElec,pExct1(k),j,Orbitals, &
                  eps1,sxindx1)
             cLoc=0
             !Test if sxindx1 is in expansion and greater than pDets(i)
             if(sxindx1.LT.pDets(i))then
                CYCLE!k->k+1
             end if!sxindx1.LT.pDets(i)
             ! Subroutine to return element of matching value in pDets
             call int_search2(sxindx1,pDets(i+1),(pDetsTrunc-i),cLoc)
             if(cLoc.EQ.0)then
                CYCLE!k->k+1
             end if
             !do l=i+1,pDetsTrunc
             !  if(sxindx1.EQ.pDets(l))then
             !        cLoc=l
             !        EXIT!l
             !  end if!sxindx.EQ.pDets(l)
             !end do!l
             !if(cLoc.EQ.0)then
             !      CYCLE !k->k+1
             !end if!cLoc.EQ.0
             ! Store single replacement information
             SingleRepInfo(1,sRIndx)=sxindx1
             SingleRepInfo(2,sRIndx)=eps1
             SingleRepInfo(3,sRIndx)=aString(j)
             SingleRepInfo(4,sRIndx)=pExct1(k)
             SingleRepInfo(5,sRIndx)=cLoc
             sRIndx=sRIndx+1 !sRIndx->sRIndx+1
             CALL eval_singlex1(aString,pExct1,aElec,Moints1,Moints2,&
                  M1Len,M2Len,Orbitals,j,k,int1e1,int2e1)    
             !Loop over corresponding q strings
             !do m=pLocate(cLoc)+1,pLocate(cLoc)+pStep(cLoc)
             !  call int_search2(pString(m),pString(pLocate(i)+1),pStep(i),oVIndx)
             !  if(oVIndx.EQ.0)then
             !    cycle !m->m+1
             !  else
             !    !Compute contribution
             !    CALL eval_singlex2('A',pString(m),bElec,Orbitals,bDets,&
             !      aString(j),pExct1(k),Moints1,Moints2,M1Len,M2Len, &
             !      qDets(qDetsTrunc),int2e2)
             !    OutVector(oVIndx)=OutVector(oVIndx)+eps1*(int1e1+int2e1+&
             !      int2e2)*InVector(m)
             !    OutVector(m)=OutVector(m)+eps1*(int1e1+int2e1+int2e2)*InVector(oVIndx)
             !  end if !oVIndx.EQ.0
             !end do!
             
             
             !Loop over corresponding q strings
             do m=pLocate(cLoc)+1,pLocate(cLoc)+pStep(cLoc)
                do n=pLocate(i)+1, pLocate(i)+ pStep(i)
                   if(pString(m).EQ.pString(n))then
                      !Compute contribution
                      CALL eval_singlex2('A',pString(m),bElec,&
                           Orbitals,bDets,aString(j),pExct1(k),Moints1,  &
                           Moints2,M1Len,M2Len,qDets(qDetsTrunc),int2e2)
                      !Find index of determinant c[r(pjv*,q)]
                      iVIndx=m
                      oVIndx=n
                      OutVector(oVIndx)=OutVector(oVIndx)+eps1*(int1e1+int2e1&
                           +int2e2)*InVector(iVIndx)
                      OutVector(iVIndx)=OutVector(iVIndx)+eps1*(int1e1+int2e1&
                           +int2e2)*InVector(oVIndx)
                      EXIT !n
                   end if!pstring(plocate(cLoc)+m,2).EQ.pString(pLocate(cLoc)+m,2)
                end do !n
             end do !m
             !Loop over  additional replacements in alpha strings
             do m=j+1,aElec
                do n=k+1,orbitals-aElec
                   !m must always be greater than j
                   CALL DRepInfo(aString,aElec,pExct1(k),j,pExct1(n),m,Orbitals,eps2,dxindx1)
                   cLoc2=0
                   !Test if dxindx1 is in expansion and greater than pDets(i)
                   if(dxindx1.LT.pDets(i))then
                      CYCLE !n->n+1
                   end if!dxindx1.LT.pDets(i)
                   !Call sub to find location of dxindx1 in pDets
                   call int_search2(dxindx1,pDets(i+1),(pDetsTrunc-i),cLoc2)
                   !do o=i+1,pDetsTrunc
                   !  if(dxindx1.EQ.pDets(o))then
                   !        cLoc2=o
                   !        EXIT!o
                   !  end if!dxindx1.EQ.pDets(o)
                   !end do!o
                   if(cLoc2.EQ.0)then
                      CYCLE !n->n+1
                   end if!cLoc2.EQ.0
                   int1e3=eps2*(Moints2(Index2E(aString(j),pExct1(k),aString(m),pExct1(n))) &
                        - Moints2(Index2E(aString(j),pExct1(n),aString(m),pExct1(k))))
                   !Loop over q's
                   do o=pLocate(cLoc2)+1,pLocate(cLoc2)+pstep(cLoc2)
                      do p=pLocate(i)+1,pLocate(i)+pstep(i)
                         if(pString(o).EQ.pString(p))then
                            iVIndx=o
                            oVIndx=p
                            OutVector(oVIndx)=OutVector(oVIndx)+int1e3*InVector(iVIndx)
                            OutVector(iVIndx)=OutVector(iVIndx)+int1e3*InVector(oVIndx)
                            EXIT !p
                         end if
                      end do !p
                   end do !o
                end do !n
             end do !m
          end do !k
       end do !j
       ! Loop over single excitations stored from loop j
       do j=1,sRIndx-1
          ! Loop over q-strings: |p,q> and |p*,q*> must be in expansion
          do k=pLocate(i)+1,pLocate(i)+pstep(i)
             !Generate an orbital string for q
             CALL GenOrbString_Beta(pString(k),bElec,Orbitals,bDets, &
                  qDets(qDetsTrunc),bString)
             !Generate a list of possible excitations
             CALL possex1(bString,Orbitals,bElec,(Orbitals-bElec),qExct1)
             !Loop over single excitations
             do l=nFrozen+1,bElec
                do m=1,Orbitals-belec
                   CALL SRepInfo(bString,bElec,qExct1(m),l,Orbitals,eps2,sxindx2)
                   !Loop over q indices that correpsond to the single excitation r(pjv*) in p strings
                   ! q*>q
                   if(sxindx2.LT.pString(k))then
                      CYCLE !m->m+1
                   end if!sxindx2.LT.pString(pLocate(i)+k,2)
                   cLoc2=0
                   call int_search2(sxindx2,pString(pLocate(SingleRepInfo(5,j))+1), &
                        pStep(SingleRepInfo(5,j)),cLoc2)
                   !do n=1,pstep(SingleRepInfo(5,j))
                   !  if(sxindx2.EQ.pString(pLocate(SingleRepInfo(5,j))+n,2))then
                   !    cLoc2=n
                   !    EXIT !n
                   !  end if!sxindx2.EQ.pString(pLocate(SingleRepInfo(5,j)),2)
                   !end do !n
                   if(cLoc2.EQ.0)then
                      CYCLE !m->m+1
                   end if!cLoc2.EQ.0
                   int3e2=Moints2(Index2E(SingleRepInfo(3,j),SingleRepInfo(4,j),bString(l),qExct1(m)))
                   iVIndx=pLocate(SingleRepInfo(5,j))+cLoc2
                   oVIndx=k
                   OutVector(oVIndx)=OutVector(oVIndx)+int3e2*SingleRepInfo(2,j)*eps2*InVector(iVIndx)
                   OutVector(iVIndx)=OutVector(iVIndx)+int3e2*SingleRepInfo(2,j)*eps2*InVector(oVIndx)
                   EXIT !m
                end do!m
             end do!l
          end do!k
       end do !j
       tmpString=aString ! For next iteration
    end do !i
    !Deallocate arrays
    DEALLOCATE(aString,bString,pExct1,qExct1,SingleRepInfo)
    RETURN !Return OutVector
  end subroutine hv_alpha
  !====================================================================
  !>hv_beta
! Computes beta string contribution to Hv=c
!--------------------------------------------------------------------
      subroutine hv_beta(InVector,Moints1,M1Len,Moints2,M2Len,pString,  &
            pStep,pLocate,pDets,qString,qStep,qLocate,qDets,ciDim,      &
            pDetsTrunc,qDetsTrunc,aDets,bDets,aElec,bElec,Orbitals,     &
            nFrozen,nDOCC,nActive,XRefList,OutVector)
            use addressing,   only:GenOrbString_Alpha,GenOrbString_Beta,StrFind3
            use string_util,  only:possex1
            use integral,     only:Index2E
            use search_fcns,  only:int_search2
            IMPLICIT NONE
            integer,intent(IN)      ::M1Len,M2Len,ciDim,pDetsTrunc,qDetsTrunc, &
                  aDets,bDets,aElec,bElec,Orbitals,nFrozen,nDOCC,nActive
            integer,dimension(pDetsTrunc),intent(IN)  ::pDets,pLocate,pStep
            integer,dimension(qDetsTrunc),intent(IN)  ::qDets,qLocate,qStep
            integer,dimension(ciDim),intent(IN)     ::pString,qString
            integer,dimension(ciDim),intent(IN)       ::XRefList
            real*8,dimension(M1Len),intent(IN)        ::Moints1
            real*8,dimension(M2Len),intent(IN)        ::Moints2
            real*8,dimension(ciDim),intent(IN)        ::InVector
            real*8,dimension(ciDim),intent(INOUT)     ::OutVector
            integer,dimension(:),ALLOCATABLE::bString,aString,qExct1,pExct1,&
                  tmpString
            real*8::int1e1,int2e1,int2e2,int1e3,int3e2
            integer ::i,j,k,l,m,n,o,p
            integer ::cLoc,cLoc2
            integer ::eps1,eps2,sxindx1,dxindx1,iVIndx,oVIndx
            !Allocate alpha string and beta string
            ALLOCATE(aString(aElec))
            ALLOCATE(bString(bElec))
            ALLOCATE(tmpString(bElec))
            !Allocate possible excitation arrays
            ALLOCATE(pExct1(Orbitals-aElec))
            ALLOCATE(qExct1(Orbitals-bElec))
            !Loop over beta strings in expansion
            do i=1,qDetsTrunc
                  !Generate string qDets(i)
                  !CALL GenOrbString(qDets(i),bElec,Orbitals,bDets,bString)
                   !Generate string pDets(i)
                  if (i.NE.1)then
                        CALL StrFind3(tmpString,qDets(i-1),bElec,Orbitals,bDets, &
                        qDets(i),bString)
                      
                  else
                        do j=1,bElec
                              bString(j)=j
                        end do
                  end if
                  !Generate possible excitations
                  CALL possex1(bString,Orbitals,bElec,(Orbitals-bElec),qExct1)
                  !Loop over single excitations
                  do j=nFrozen+1,bElec
                    do k=1,Orbitals-bElec
                      CALL SRepInfo(bString,bElec,qExct1(k),j,Orbitals,eps1,sxindx1)
                      cLoc=0
                      !Test if sxindx1 is in expansion and greater than qDets(i)
                      if(sxindx1.LT.qDets(i))then
                            CYCLE !k->k+1
                      end if!sxindx1.LT.qDets(i)
                      call int_search2(sxindx1,qDets(i+1),(qDetsTrunc-i),cLoc)
                      !do l=i,qDetsTrunc
                      !  if(sxindx1.EQ.qDets(l))then
                      !        cLoc=l
                      !        EXIT!l
                      !  end if!sxindx1.EQ.qDets(l)
                      !end do!l
                      if(cLoc.EQ.0)then
                            CYCLE!k->k+1
                      end if!cLoc.EQ.0
                      !Continue...
                      CALL eval_singlex1(bString,qExct1,bElec,Moints1,Moints2,M1Len,&
                            M2Len,Orbitals,j,k,int1e1,int2e1)
                      !Loop over corresponding p strings
                      do m=qLocate(cLoc)+1,qLocate(cLoc)+qStep(cLoc)
                        do n=qLocate(i)+1,qLocate(i)+qStep(i)
                          if(qString(m).EQ.qString(n))then
                            !Compute contribution
                            CALL eval_singlex2('B',qString(m),bElec,Orbitals,bDets,&
                                  bString(j),qExct1(k),Moints1,Moints2,M1Len,M2Len, &
                                  pDets(pDetsTrunc),int2e2)
                            !Index of element
                            iVIndx=XRefList(m)
                            oVIndx=XRefList(n)
                            OutVector(oVIndx)=OutVector(oVIndx)+eps1*(int1e1+int2e1+int2e2)*&
                                  InVector(iVIndx)
                            OutVector(iVIndx)=OutVector(iVIndx)+eps1*(int1e1+int2e1+int2e2)*&
                                  InVector(oVIndx)
                            EXIT !n
                          end if!qString(qLocate(cLoc)+m,2).EQ.qString(qLocate(i)+n,2)
                        end do !n
                      end do !m
                      !Loop over additional replacements
                      do m=j+1,bElec
                        do n=k+1,Orbitals-bElec
                          ! m must always be greater than j
                          CALL DRepInfo(bString,bElec,qExct1(k),j,qExct1(n),m,Orbitals,&
                                eps2,dxindx1)
                          int3e2=eps2*(Moints2(Index2E(bString(j),qExct1(k),bString(m),&
                                qExct1(n))) - Moints2(Index2E(bString(j),qExct1(n),    &
                                bString(m),qExct1(k))))
                          !Test if q''>q and dxindx1 is in expansion
                          cLoc2=0
                          if(dxindx1.LT.qDets(i))then
                                CYCLE !n->n+1
                          end if!dxindx1.LT.qDets(i)
                          call int_search2(dxindx1,qDets(i+1),(qDetsTrunc-i),cLoc2)
                          !do o=i+1,qDetsTrunc
                          !  if(dxindx1.EQ.qDets(o))then
                          !        cLoc2=o
                          !        EXIT !o
                          !  end if!dxindx1.EQ.qDets(o)
                          !end do !o
                          if(cLoc2.EQ.0)then
                                CYCLE !n->n+1
                          end if!cLoc2.EQ.0
                          int1e3=eps2*(Moints2(Index2E(bString(j),qExct1(k),bString(m),qExct1(n)))&
                                - Moints2(Index2E(bString(j),qExct1(n),bString(m),qExct1(k))))
                          !Loop over p's
                          do o=qLocate(cLoc2)+1,qLocate(cLoc2)+qstep(cLoc2)
                            do p=qLocate(i)+1,qLocate(i)+qstep(i)
                              if(qString(o).EQ.qString(p))then
                                    iVIndx=XRefList(o)
                                    oVIndx=XRefList(p)
                                    OutVector(oVIndx)=OutVector(oVIndx)+int1e3*InVector(iVIndx)
                                    OutVector(iVIndx)=OutVector(iVIndx)+int1e3*InVector(oVIndx)
                                    EXIT !p
                              end if
                            end do !p
                          end do !o
                        end do !n
                      end do !m
                    end do !k
                  end do !j
                  tmpString=bString
            end do !i
            !Dellocate arrays
            DEALLOCATE(aString,bString,pExct1,qExct1)
            RETURN !Return OutVector
      end subroutine hv_beta
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
      subroutine SRepInfo(String,SLen,NewOrb,IndxSwp,Orbitals,Eps1,Indx)
            use addressing, only: AdrFind
            use string_util, only: cannon1swp
            IMPLICIT NONE
            integer,intent(IN)      ::SLen,NewOrb,IndxSwp,Orbitals
            integer,dimension(SLen),intent(IN)  ::String
            integer,intent(OUT)     ::Eps1,Indx
            integer,dimension(:),ALLOCATABLE    ::newString
            !Allocate newString
            ALLOCATE(newString(SLen))
            newString=String
            !Swap orbitals
            newString(IndxSwp)=NewOrb
            !Return to cannonical ordering
            CALL cannon1swp(newString,SLen,IndxSwp,Eps1)
            !Locate address
            CALL AdrFind(newString,SLen,Orbitals,Indx)
            DEALLOCATE(newString)
            RETURN
      end subroutine SRepInfo
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
            integer,dimension(:),ALLOCATABLE    ::tmp_string
            ALLOCATE(tmp_string(sp_elecs))
            !One electron contribution
            contr_1elec=MOints1(Ind2Val(sp_string(bound_elec),sp_excites(excited_elec)))
            !Two electron contribution
            tmp_string=sp_string
            tmp_string(bound_elec)=sp_excites(excited_elec)
            !Compute contribution
            contr_2elec=0d0
            do i=1, bound_elec-1
                  contr_2elec=contr_2elec+MOints2(Index2E(sp_string(i),sp_string(i), &
                        sp_string(bound_elec),sp_excites(excited_elec))) -             &
                        MOints2(Index2E(sp_string(i),sp_string(bound_elec),sp_string(i), &
                        sp_excites(excited_elec)))
            end do
            do i=bound_elec+1,sp_elecs
                  contr_2elec=contr_2elec+Moints2(Index2E(sp_string(i),sp_string(i),     &
                        sp_string(bound_elec),sp_excites(excited_elec))) -                &
                        MOints2(Index2E(sp_string(i),sp_string(bound_elec),sp_string(i), &
                        sp_excites(excited_elec)))
            end do
            DEALLOCATE(tmp_string)
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
