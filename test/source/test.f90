! PROGRAM FOR TESTINGS SUBROUTINES
PROGRAM test
      use addressing
      use combinatorial
      use integral
      use string_util
      use input_proc
      use action_util
      use search_fcns
      integer     ::i,j
      
      integer     ::test_factl_m, test_factl_n, test_factl_out, &
                  test_factl_test
      !-----------------------------
      integer     ::test_binom_test,test_binom_n,test_binom_k, &
                  test_binom_out
      !-----------------------------
      integer     ::test_k2indc_testp,test_k2indc_testq, &
                  test_k2indc_k,test_k2indc_orbitals, test_k2indc_belecs,&
                  test_k2indc_p,test_k2indc_q
      !-----------------------------
      integer     ::test_indexk_p,test_indexk_q,test_indexk_orbitals, &
                  test_indexk_belec,test_indexk_testindex, &
                  test_indexk_out
      !-----------------------------
      integer     ::test_adrfind_p,test_adrfind_elec,test_adrfind_orbitals,&
                  test_adrfind_testaddress,test_adrfind_out
      integer,dimension(:),allocatable    ::test_adrfind_string
      !-----------------------------
      integer     ::test_strfind1_addr,test_strfind1_elec,test_strfind1_orb,&
                  test_strfind1_spdets,test_strfind1_test
      integer,dimension(:,:),allocatable  ::test_strfind1_out_str, &
                                          test_strfind1_test_str
      !-----------------------------
      integer     ::test_strfind2_elec,test_strfind2_orb,test_strfind2_spdets, &
                  test_strfind2_test
      integer,dimension(:),allocatable    ::test_strfind2_str1,test_strfind2_str2,&
                                          test_strfind2_test_str2
      !-----------------------------
      integer     ::test_gorbstr_index,test_gorbstr_elec,test_gorbstr_orb, &
                  test_gorbstr_dets,test_gorgstr_test
      integer,dimension(:),allocatable    ::test_gorbstr_out,test_gorbstr_teststr
      !-----------------------------
      integer     ::test_ind2val_i,test_ind2val_j,test_ind2val_out, &
                  test_ind2val_test
      !-----------------------------
      integer     ::test_cannon1swp_len,test_cannon1swp_indswp,test_cannon1swp_permind,&
                  test_cannon1swp_test_permind
      integer,dimension(:),ALLOCATABLE    ::test_cannon1swp_string
      !-----------------------------
      integer     ::test_possex1_orbs,test_possex1_elec,test_possex1_pexitslen, &
                  test_possex1_test
      integer,dimension(:),ALLOCATABLE    ::test_possex1_occstr,test_possex1_test_pexits,&
                                          test_possex1_pexits
      !-----------------------------
      integer     ::test_cannon2_len,test_cannon2_indx1,test_cannon2_swp1, &
                  test_cannon2_indx2,test_cannon2_swp2,test_cannon2_permind,&
                  test_cannon2_test_permind
      integer,dimension(:),ALLOCATABLE    ::test_cannon2_string,test_cannon2_test_str
      !-----------------------------
      integer     ::test_alpha_beta_elec,test_alpha_beta_aelec,test_alpha_beta_belec, &
                  test_alpha_beta_testa,test_alpha_beta_testb
      !=============================
      integer     ::hv_d_ciDim,hv_d_test
      real*8,dimension(:),ALLOCATABLE     ::hv_d_D,hv_d_IV, &
                                          hv_d_OV
      real*8,dimension(:),ALLOCATABLE     ::hv_d_OV_test
      !-----------------------------
      integer     ::srI_sLen,srI_NewOrb,srI_IndxSwp,srI_Orbitals,srI_Eps,srI_Indx,&
                  srI_test_Eps,srI_test_Indx
      integer,dimension(:),ALLOCATABLE    ::srI_String
      !-----------------------------
      integer     ::drI_sLen,drI_no1,drI_no2,drI_ind1,drI_ind2,drI_Orbitals,drI_eps1,&
                  drI_indx,drI_test_indx,drI_test_eps1
      integer,dimension(:),ALLOCATABLE    ::drI_string
      !------------------------------
      integer     ::int_search_length, int_search_test, int_search_out, int_search_find
      integer,dimension(:),ALLOCATABLE    :: int_search_list
      !===============================================================
      print *, "COMMON_UTIL TESTING..."
      !==== sub Factl(m,n) ===========================================
      print *, "Testing Factl(m,n) function..... "
      test_factl_m=97
      test_factl_n=100
      test_factl_test=94109400
      test_factl_out=Factl(test_factl_m,test_factl_n)
      if(test_factl_out.NE.test_factl_test)then
            print *, "  Factl(m,n) test failed."
      else
            print *, "  Factl(m,n) test succeeded."
      end if
      !==== sub Binom(n,k) ===========================================
      print *, "Testing Binom(n,k) function..... "
      test_binom_test=17310309456440
      test_binom_n=100
      test_binom_k=10
      test_binom_out=Binom(test_binom_n,test_binom_k)
      if(test_binom_out.NE.test_binom_test)then
            print *, "  Binom() test failed."
      else
            print *, "  Binom() test succeeded."
      end if
      !test_binom_test=1175759825268299040
      !test_binom_n=110
      !test_binom_k=15
      !test_binom_out=Binom(test_binom_n,test_binom_k)
      !if(test_binom_out.NE.test_binom_test)then
      !      print *,"  Binom() test 2 failed."
      !else
      !      print *,"  Binom() test 2 succeeded."
      !end if
      !===============================================================
      !==== sub K2Indc(K,orbitals,belec,p,q) =========================
      print *, "Testing K2Indc(K,orbitals,belec,p,q) subroutine... "
      test_k2indc_testp=5
      test_k2indc_testq=10
      test_k2indc_orbitals=20
      test_k2indc_belecs=4
      test_k2indc_k = IndexK(test_k2indc_testp,test_k2indc_testq,test_k2indc_orbitals,&
                        test_k2indc_belecs)
      CALL K2Indc(test_k2indc_k,test_k2indc_belecs,test_k2indc_orbitals,test_k2indc_p,test_k2indc_q)
      if(test_k2indc_p.EQ.test_k2indc_testp .AND. test_k2indc_q.EQ.test_k2indc_testq)then
            print *, "  K2Indc() test succeeded."
      else if(test_k2indc_p.NE.test_k2indc_testp .AND. test_k2indc_q.EQ.test_k2indc_testq)then
            print *, "  K2Indc() test FAILED: test_k2indc_p .NE. test_k2indc_testp. "
      else if(test_k2indc_p.EQ.test_k2indc_testp .AND. test_k2indc_q.NE.test_k2indc_testq)then
            print *, "  K2Indc() test FAILED: test_k2indc_q .NE. test_k2indc_testq. "
      else if(test_k2indc_p.NE.test_k2indc_testp .AND. test_k2indc_q.NE.test_k2indc_testq)then
            print *, "  K2Indc() test FAILED: test_k2indc_q .NE. test_k2indc_testq .AND. test_k2indc_p .NE. test_k2indc_testp. "
      end if
      !===============================================================
      !==== sub IndexK(p,q,Orbitals,belec) ===========================
      print *, "Testing IndexK(p,q,Orbitals,belec) subroutine... "
      test_indexk_orbitals=20
      test_indexk_belec=5
      test_indexk_testindex=13
      CALL K2Indc(test_indexk_testindex,test_indexk_belec,test_indexk_orbitals, &
            test_indexk_p,test_indexk_q)
      test_indexk_out=IndexK(test_indexk_p,test_indexk_q,test_indexk_orbitals, &
            test_indexk_belec)
      if(test_indexk_out.NE.test_indexk_testindex)then
            print *, "  IndexK() test failed."
      else
            print *, "  Indexk() test succeeded."
      end if
      !==============================================================
      !==== sub AdrFind(string,elecs,Orbitals,address) ==============
      print *, "Testing AdrFind(string,elecs,Orbitals,address) subroutine... "
      test_adrfind_elec=3
      test_adrfind_orbitals=6
      ALLOCATE(test_adrfind_string(test_adrfind_elec))
      test_adrfind_string=(/2,3,5/)
      test_adrfind_testaddress=12
      CALL AdrFind(test_adrfind_string,test_adrfind_elec,test_adrfind_orbitals, &
            test_adrfind_out)
      if(test_adrfind_out.NE.test_adrfind_testaddress)then
            print *, "  AdrFind() test failed."
      else
            print *, "  AdrFind() test succeeded."
      end if
      DEALLOCATE(test_adrfind_string) !Deallocate array
      !============================================================== 
      !==== sub StrFind1(address,spElecs,orbs,spDets,array) =========       
      print *, "Testing StrFind1() subroutine..."
      test_strfind1_addr=20
      test_strfind1_elec=3
      test_strfind1_orb=6
      test_strfind1_spdets=20
      ALLOCATE(test_strfind1_test_str(test_strfind1_elec,test_strfind1_addr))
      test_strfind1_test_str(1:3,1) =(/ 1,2,3 /)
      test_strfind1_test_str(1:3,2) =(/ 1,2,4 /)
      test_strfind1_test_str(1:3,3) =(/ 1,2,5 /)
      test_strfind1_test_str(1:3,4) =(/ 1,2,6 /)
      test_strfind1_test_str(1:3,5) =(/ 1,3,4 /)
      test_strfind1_test_str(1:3,6) =(/ 1,3,5 /)
      test_strfind1_test_str(1:3,7) =(/ 1,3,6 /)
      test_strfind1_test_str(1:3,8) =(/ 1,4,5 /)
      test_strfind1_test_str(1:3,9) =(/ 1,4,6 /)
      test_strfind1_test_str(1:3,10)=(/ 1,5,6 /)
      test_strfind1_test_str(1:3,11)=(/ 2,3,4 /)
      test_strfind1_test_str(1:3,12)=(/ 2,3,5 /) 
      test_strfind1_test_str(1:3,13)=(/ 2,3,6 /)
      test_strfind1_test_str(1:3,14)=(/ 2,4,5 /)
      test_strfind1_test_str(1:3,15)=(/ 2,4,6 /)
      test_strfind1_test_str(1:3,16)=(/ 2,5,6 /)
      test_strfind1_test_str(1:3,17)=(/ 3,4,5 /)
      test_strfind1_test_str(1:3,18)=(/ 3,4,6 /)
      test_strfind1_test_str(1:3,19)=(/ 3,5,6 /)
      test_strfind1_test_str(1:3,20)=(/ 4,5,6 /)
      ALLOCATE(test_strfind1_out_str(test_strfind1_elec,test_strfind1_addr))
      CALL StrFind1(test_strfind1_spdets,test_strfind1_elec,test_strfind1_orb,&
            test_strfind1_spdets,test_strfind1_out_str)
      test_strfind1_test=0
      do i=1,test_strfind1_addr
            do j=1,test_strfind1_elec
                  if( &
                  test_strfind1_out_str(j,i).NE.test_strfind1_test_str(j,i) )then
                        test_strfind1_test=test_strfind1_test+1
                  end if
            end do
      end do
      if(test_strfind1_test.GT.0)then
            print *, "  StrFind1() test failed."
      else
            print *, "  StrFind1() test succeeded."
      end if
      DEALLOCATE(test_strfind1_test_str,test_strfind1_out_str)
      !===============================================================
      !==== StrFind2(string1,spElecs,orbs,spDets,string2) ============
      print *, "Testing StrFind2() subroutine..."
      test_strfind2_elec=3
      test_strfind2_orb=6
      test_strfind2_spdets=20
      ALLOCATE(test_strfind2_str1(test_strfind2_elec))
      ALLOCATE(test_strfind2_str2(test_strfind2_elec))
      ALLOCATE(test_strfind2_test_str2(test_strfind2_elec))
      test_strfind2_str1=(/ 1,5,6 /)
      test_strfind2_test_str2=(/ 2,3,4 /)
      CALL StrFind2(test_strfind2_str1,test_strfind2_elec,test_strfind2_orb, &
            test_strfind2_spdets, test_strfind2_str2)
      test_strfind2_test=0
      do i=1,test_strfind2_elec
            if(test_strfind2_str2(i).NE.test_strfind2_test_str2(i))then
                  test_strfind2_test=test_strfind2_test+1
            end if
      end do
      if( test_strfind2_test.NE.0)then
            print *, "  StrFind2() test failed."
      else
            print *, "  StrFind2() test succeeded."
      end if
      DEALLOCATE(test_strfind2_str1,test_strfind2_str2,test_strfind2_test_str2)
      !===============================================================
      !==== GenOrbString(str_index,str_elec,Orbitals,str_Dets,str_string)
      print *, "Testing GenOrbString() subroutine... "
      test_gorbstr_elec=3
      test_gorbstr_orb=6
      test_gorbstr_dets=20
      test_gorbstr_index=18
      ALLOCATE(test_gorbstr_teststr(test_gorbstr_elec))
      ALLOCATE(test_gorbstr_out(test_gorbstr_elec))
      test_gorbstr_teststr=(/ 3,4,6 /)
      CALL GenOrbString(test_gorbstr_index,test_gorbstr_elec,test_gorbstr_orb,&
                  test_gorbstr_dets,test_gorbstr_out)
      test_gorbstr_test=0
      do i=1,test_gorbstr_elec
            if(test_gorbstr_teststr(i).NE.test_gorbstr_out(i))then
                  test_gorbstr_test=test_gorbstr_test+1
            end if
      end do
      if( test_gorbstr_test.NE.0)then
            print *, "  GenOrbString() test failed."
      else
            print *, "  GenOrbString() test succeeded."
      end if
      DEALLOCATE(test_gorbstr_teststr,test_gorbstr_out)
      !===============================================================
      !==== GenOrbString_Alpha() =====================================
      print *, "Testing GenOrbString_Beta()"
      ALLOCATE(test_gorbstr_teststr(test_gorbstr_elec))
      ALLOCATE(test_gorbstr_out(test_gorbstr_elec))
      test_gorbstr_teststr=(/ 3,4,6 /)
      CALL GenOrbString_Beta(test_gorbstr_index,test_gorbstr_elec,test_gorbstr_orb,&
                  test_gorbstr_dets, 20,test_gorbstr_out)
       test_gorbstr_test=0
      do i=1,test_gorbstr_elec
            if(test_gorbstr_teststr(i).NE.test_gorbstr_out(i))then
                  test_gorbstr_test=test_gorbstr_test+1
            end if
      end do
      if( test_gorbstr_test.NE.0)then
            print *, "  GenOrbString_Alpha() test failed."
            print *, test_gorbstr_out
      else
            print *, "  GenOrbString_Alpha() test succeeded."
      end if
   
      DEALLOCATE(test_gorbstr_teststr,test_gorbstr_out)
      
      !===============================================================
      !==== Ind2Val(i,j) =============================================
      print *, "Testing Ind2Val(i,j) subroutine..."
      test_ind2val_i=13
      test_ind2val_j=2
      test_ind2val_test=test_ind2val_i*(test_ind2val_i-1)/2+test_ind2val_j
      if(Ind2Val(13,2).NE.test_ind2val_test .OR. &
            Ind2Val(2,13).NE.test_ind2val_test)then
            print *, "  Ind2Val() test failed."
      else
            print *, "  Ind2Val() test succeeded."
      end if
      !==============================================================
      !==== cannon1swp(string,length,indswp,permind =================
      print *, "Testing cannon1swp() subroutine..."
      test_cannon1swp_len=5
      test_cannon1swp_indswp=3
      test_cannon1swp_test_permind=-1
      ALLOCATE(test_cannon1swp_string(test_cannon1swp_len))
      test_cannon1swp_string=(/5,7,8,10,13/)
      test_cannon1swp_string(test_cannon1swp_indswp)=11
      CALL cannon1swp(test_cannon1swp_string,test_cannon1swp_len,test_cannon1swp_indswp,&
            test_cannon1swp_permind)
      if(test_cannon1swp_permind.NE.test_cannon1swp_test_permind)then
            print *, "  cannon1swp() test  failed."
      else
            print *, "  cannon1swp() test  succeeded."
      end if
      CALL cannon1swp(test_cannon1swp_string,test_cannon1swp_len,test_cannon1swp_indswp,&
            test_cannon1swp_permind)
      if(test_cannon1swp_permind.NE.1)then
            print *, "  cannon1swp() test2 failed."
      else
            print *, "  cannon1swp() test2 succeeded."
      end if
      DEALLOCATE(test_cannon1swp_string)
      !===============================================================
      !==== possex1(occstr,obs,elecs,pexitslen,pexits) ===============
      print *, "Testing possex1() subroutine..."
      test_possex1_elec=3
      test_possex1_orbs=6
      test_possex1_pexitslen=3
      ALLOCATE(test_possex1_occstr(test_possex1_elec))
      ALLOCATE(test_possex1_test_pexits(test_possex1_pexitslen))
      ALLOCATE(test_possex1_pexits(test_possex1_pexitslen))
      test_possex1_occstr=(/ 2,5,6 /)
      test_possex1_test_pexits=(/ 1,3,4 /)
      CALL possex1(test_possex1_occstr,test_possex1_orbs,test_possex1_elec,&
            test_possex1_pexitslen,test_possex1_pexits)
      test_possex1_test=0
      do i=1, test_possex1_pexitslen
            if(test_possex1_pexits(i).NE.test_possex1_test_pexits(i))then
                  test_possex1_test=test_possex1_test+1
            end if
      end do
      if(test_possex1_test.NE.0)then
            print *, "  possex1() test failed."
      else 
            print *, "  possex1() test succeeded."
      end if
      DEALLOCATE(test_possex1_occstr,test_possex1_test_pexits,test_possex1_pexits)
      !===============================================================
      !==== cannon2(string,length,indx1,swp1,indx2,swp2,permind) =====
      print *, "Testing cannon2() subroutine..."
      test_cannon2_len=10
      test_cannon2_indx1=3
      test_cannon2_indx2=8
      test_cannon2_swp1=12
      test_cannon2_swp2=11
      test_cannon2_test_permind=-1
      ALLOCATE(test_cannon2_string(test_cannon2_len))
      test_cannon2_string=(/1,2,3,4,5,6,7,8,9,10/)
      CALL cannon2(test_cannon2_string,test_cannon2_len,test_cannon2_indx1,&
            test_cannon2_swp1,test_cannon2_indx2,test_cannon2_swp2, &
            test_cannon2_permind)
      !Should be -1
      if(test_cannon2_permind.NE.test_cannon2_test_permind)then
            print *, "  cannon2() test failed."
      else
            print *, "  cannon2() test succeeded."
      end if
      DEALLOCATE(test_cannon2_string)
      !===============================================================
      !==== alpha_beta(Electrons,aElec,bElec)
      print *, "Testing alpha_beta() subroutine..."
      test_alpha_beta_elec=5
      test_alpha_beta_testa=3
      test_alpha_beta_testb=2
      CALL alpha_beta(test_alpha_beta_elec,test_alpha_beta_aelec, &
            test_alpha_beta_belec)
      if( test_alpha_beta_aelec.NE.test_alpha_beta_testa .OR. &
            test_alpha_beta_belec.NE.test_alpha_beta_testb )then
            print *, "  alpha_beta() test failed."
      else
            print *, "  alpha_beta() test succeeded."
      end if 
      print *, " "
      !===============================================================
      !///////////////////////////////////////////////////////////////
      !===============================================================
      print *, "ACTION_UTIL TESTING..."  
      !===============================================================
      !==== hv_diagonals() ===========================================
      print *, "Testing hv_diagonals() subroutine... "
      hv_d_ciDim=1000
      !Allocate arrays
      ALLOCATE(hv_d_D(hv_d_ciDim))
      ALLOCATE(hv_d_IV(hv_d_ciDim))
      ALLOCATE(hv_d_OV(hv_d_ciDim))
      ALLOCATE(hv_d_OV_test(hv_d_ciDim))
      do i=1,hv_d_ciDim
            hv_d_D=REAL(i)
            hv_d_OV_test=REAL(i)
      end do
      hv_d_IV=1d0
      hv_d_OV=0d0
      CALL hv_diagonals(hv_d_IV,hv_d_D,hv_d_ciDim,hv_d_OV)
      hv_d_test=0
      do i=1,hv_d_ciDim
            if(hv_d_OV(i).NE.hv_d_OV_test(i))then
                  hv_d_test=hv_d_test+1
            end if
      end do
      if(hv_d_test.NE.0)then
            print *, "  hv_diagonals() test failed..."
      else
            print *, "  hv_diagonals() test succeeded..."
      end if
      DEALLOCATE(hv_d_D,hv_d_IV,hv_d_OV,hv_D_OV_test)
      !===============================================================
      !==== SRepInfo() ===============================================
      print *, "Testing SRepInfo() subroutine... "
      srI_sLen=3
      srI_IndxSwp=2
      srI_NewOrb=5
      srI_Orbitals=6
      srI_test_Eps=-1
      srI_test_Indx=6
      ALLOCATE(srI_String(srI_sLen))
      do i=1,srI_sLen
            srI_String(i)=i
      end do
      srI_String(srI_IndxSwp)=srI_NewOrb
      CALL SRepInfo(srI_String,srI_sLen,srI_NewOrb,srI_IndxSwp,srI_Orbitals,&
            srI_Eps,srI_Indx)
      if(srI_Eps.NE.srI_test_Eps)then
            print *, "  SRepInfo() test of Eps  value failed... "
      else
            print *, "  SRepInfo() test of Eps  value succeeded... "
      end if
      if(srI_Indx.NE.srI_test_Indx)then
            print *, "  SRepInfo() test of Indx value failed... "
      else
            print *, "  SRepInfo() test of Indx value succeeded..."
      end if
      DEALLOCATE(srI_String)
      !===============================================================
      !==== DRepInfo() ===============================================
      print *, "Testing DRepInfo() subroutine... "
      drI_sLen=3
      drI_Orbitals=6
      drI_no1=6
      drI_no2=5
      drI_ind1=1
      drI_ind2=2
      drI_test_eps=-1
      drI_test_indx=19
      ALLOCATE(drI_string(drI_sLen))
      do i=1,drI_sLen
            drI_string(i)=i
      end do
      CALL DRepInfo(drI_string,drI_sLen,drI_no1,drI_ind1,drI_no2,drI_ind2,&
            drI_Orbitals,drI_eps1,drI_indx)
      if(drI_eps1.NE.drI_test_eps)then
            print *, "  DRepInfo() test of Eps1 value failed... "
      else
            print *, "  DRepInfo() test of Eps1 value succeeded..."
      end if
      if(drI_indx.NE.drI_test_indx)then
            print *, "  DRepInfo() test of Indx value failed..."
      else
            print *, "  DRepInfo() test of Indx value succeeded..."
      end if
      DEALLOCATE(drI_string)
      !===============================================================
      !==== int_search() =============================================
      print *, "Testing int_search() subroutine... "
      int_search_length=100
      int_search_test=0 !Should not find element
      int_search_find=900
      ALLOCATE(int_search_list(int_search_length))
      do i=1,100
            int_search_list(i)=i*3
      end do
      CALL int_search2(int_search_find,int_search_list,int_search_length,int_search_out)
      if(int_search_out.NE.int_search_test)then
            print *, "  int_search() test failed... "
            print *, int_search_out, int_search_test
      else
            print *, "  int_search() test succeeded... "
      end if
      DEALLOCATE(int_search_list)
END PROGRAM
