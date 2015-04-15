ciwfmt.f
ciwfmt part=1 of 1.  formatted integral file utility program
cversion=1.0 last modified: 04-may-92
cversion 5.0
c
      subroutine iwfmt( moints1, moints2, type1, bsfcns, moflnm, 
     &  max1e, max2e)
c
c  read an integral file and write a formatted version.
c
c  version 1.0b5 04-may-92
c
c  version history:
c  04-may-92 minor ftnchek-related cleanup. -rls
c  21-oct-90 written by ron shepard, based on program istat.
c  istat history:
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  06-oct-90 fsplit=2 capability added. -rls
c  11-aug-89 siftyp() call added. -rls
c  08-aug-89 SIFS version. -rls
c  written by ron shepard.
c
       implicit none
      integer aoints, buf, i, ibitv, idummy, ierr, ilab, itotal,
     & lenbuf, nbft, nenrgy, ninfo, nlist, nmap, nmax, nsym, ntitle,
     & val, nin, type1,bsfcns, max1e, max2e
      character*80 title(20)
      character*20 moflnm
      integer nbpsy(8)
      integer choice
c
      integer   nbfmxp,     nengmx,    nmapmx
      parameter(nbfmxp=511, nengmx=20, nmapmx=10)
c
      integer ietype(nengmx), info(10)
      integer imtype(nmapmx), map(nbfmxp*nmapmx)
      real*8 energy(nengmx)
      character*4 slabel(8)
      character*8 bfnlab(nbfmxp)
      character*60 fname
c
      integer   lencor
      parameter(lencor=200000)
      real*8 core(lencor)
      real*8 moints1(max1e)
      real*8 moints2(max2e)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
*@ifndef nonamelist
      namelist /input/ fname
*@endif
c
      integer  forbyt, atebyt
      external forbyt, atebyt
c
      nlist = 6
      nin = 5
c     open( unit=nin,   file='iwfmtin',  status='unknown'     )
      aoints= 10
c
*@ifdef future
*C  open nlist here if necessary.
*@else
c  assume preconnected unit.
*@endif
c
      call ibummr(nlist)
c
c     fname = 'aoints'
c      write(0,*) 'filename 1'
      fname = moflnm
c     read( nin, input, end=3 )
c3     continue
      call trnfln( 1, fname )
      open(unit=aoints,file=fname,form='unformatted',status='old')
c
c     # read the header info.
c
      call sifrh1( aoints, ntitle, nsym, nbft,
     & ninfo, nenrgy, nmap, ierr )
      if ( ierr .ne. 0 ) then
         call bummer('iwfmt: ierr=',ierr,faterr)
      elseif ( ntitle .gt. 20 ) then
         call bummer('iwfmt: ntitle=',ntitle,faterr)
      elseif ( nbft .gt. nbfmxp ) then
         call bummer('iwfmt: nbft=',nbft,faterr)
      elseif ( ninfo .gt. 10 ) then
         call bummer('iwfmt: ninfo=',ninfo,faterr)
      elseif ( nenrgy .gt. nengmx ) then
         call bummer('iwfmt: nenrgy=',nenrgy,faterr)
      elseif( nmap .gt. nmapmx ) then
         call bummer('iwfmt: nmap=',nmap,faterr)
      endif
c
c      write(nlist,6100) ntitle, nsym, nbft, ninfo, nenrgy, nmap
6100  format(1x,10i8)
c
      call sifrh2(
     & aoints, ntitle, nsym,   nbft,
     & ninfo,  nenrgy, nmap,   title,
     & nbpsy,  slabel, info,   bfnlab,
     & ietype, energy, imtype, map,
     & ierr )
      if ( ierr.ne.0 ) call bummer('iwfmt: ierr=',ierr,faterr)
c
c      write(nlist,6020) (title(i),i=1,ntitle)
6020  format(1x,a)
c
c      write(nlist,6100) (nbpsy(i),i=1,nsym)
c
c      write(nlist,6030) (slabel(i),i=1,nsym)
6030  format(1x,8a5)
c
c      write(nlist,6100) (info(i),i=1,ninfo)
c
c      write(nlist,6040) (bfnlab(i),i=1,nbft)
6040  format(1x,8a10)
c
c      write(nlist,6100) (ietype(i),i=1,nenrgy)
c
c      write(nlist,6050) (energy(i),i=1,nenrgy)
6050  format(1x,1p4e20.12)
c
c      if ( nmap .gt. 0 ) then
c         write(nlist,6100) ( imtype(i), i = 1, nmap )
c         write(nlist,6060) ( map(i), i = 1, nmap*nbft )
c      endif
6060  format(1x,20i4)
c
      if ( (info(1) .ne. 1) .and. (info(1) .ne. 2) ) then
         call bummer('iwfmt: fsplit=',info(1),faterr)
      endif
c
c     # print the 1-e integrals.
c     # core=>1:buffer,2:ilab1,3:val1,4:ibitv
c
      lenbuf = info(2)
      nmax   = info(3)
c
      buf    = 1
      ilab   = buf   + atebyt( lenbuf )
      val    = ilab  + forbyt( 2*nmax )
      ibitv  = val   + atebyt( nmax )
      itotal = ibitv + forbyt( ((nmax+63)/64)*64 ) -1
      if ( itotal .gt. lencor ) then
         call bummer( 'iwfmt: prt1e itotal=',itotal,faterr)
      endif
c
c      write(0,*) 'Print 1e integral record structure (1)'
c      write(0,*) 'Print 1e integrals in array format (2)'
      choice = 1
c      if (choice.eq.1) then 
      call prt1e(
     & nlist,      aoints,    info,        core(buf),
     & core(ilab), core(val), core(ibitv), type1, max1e, moints1 )
c      elseif (choice.eq.2) then
c      call prtarr1e(
c     . nlist,      aoints,    info,    core(1), lencor,
c     . nbpsy, nsym) 
c      endif 

c     use a modified sifrsh version for printing the arrays
c     and siftype 
c

c
c     # open the 2-e file.
c     fname = 'aoints2'
c      write(0,*) 'filename 2'
      fname = moflnm
      call trnfln( 1, fname )
      call sifo2f( aoints, aoints, fname, info, idummy, ierr )
      if ( ierr .ne. 0 ) then
         call bummer('iwfmt: from sifo2f, ierr=',ierr,faterr)
      endif
c
c     # print the 2-e integrals.
c     # core=>1:buffer,2:ilab1,3:val1,4:ibitv
c
      lenbuf = info(4)
      nmax   = info(5)
c
      buf    = 1
      ilab   = buf   + atebyt( lenbuf )
      val    = ilab  + forbyt( 4*nmax )
      ibitv  = val   + atebyt( nmax )
      itotal = ibitv + forbyt( ((nmax+63)/64)*64 ) -1
      if ( itotal .gt. lencor ) then
         call bummer( 'iwfmt: prt2e itotal=',itotal,faterr)
      endif
c
c      print *, "CALLING PRT2E"
      call prt2e(
     & nlist,      aoints,    info,        core(buf),
     & core(ilab), core(val), core(ibitv), max2e, moints2 )
c
c      print *, "CALLING SIFC2F"
      call sifc2f( aoints, info, ierr )
      if ( ierr .ne. 0 ) then
         call bummer('iwfmt: from sifc2f, ierr=',ierr,faterr)
      endif
c      print *, "CLOSING AOINTS"
      close (unit = aoints)
c
c      print *, "CALLING BUMMER"
      call bummer('normal termination',0,3)
c      stop "SERIOUSLY, STOP!"
c      print *, "END OF PROGRAM"
      return
      end
c deck prt1e
      subroutine prt1e(
     & nlist,  ntape,  info,   buf,
     & ilab,   val,    ibitv, type1, max1e, mstrmat )
c
c  read and dump the 1-e integral records.
c
c  21-oct-90 written by ron shepard.
c
       implicit none
      integer   nipv,   msame,   nmsame,   nomore
      parameter(nipv=2, msame=0, nmsame=1, nomore=2)
c
      integer nlist, ntape, type1
      real*8 buf(*), val(*)
      integer info(*), ilab(nipv,*), ibitv(*)
c
      integer ibvtyp, ierr, ifmt, itypea, itypeb, last, num,
     & num1, i, j, max1e
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer   iretbv
      parameter(iretbv=-1)
c
      real*8 fcore
      real*8 moints(max1e)
      real*8 mstrmat(max1e)
c
      last = msame
c CLM: num1 is previous num value. 
c      for creation of moints1
      num1 = 0
      do 300 i=1,max1e
        mstrmat(i) = 0d0
300   continue
c
200   if ( last .ne. nomore ) then
c
         call sifrd1( ntape, info, nipv, iretbv,
     &    buf, num, last, itypea,
     &    itypeb, ibvtyp, val,
     &    ilab, fcore, ibitv, ierr )
c
         if ( ierr.ne.0 ) call bummer(
     &    'iwfmt: sifrd1 ierr=',ierr,faterr)
         if ( type1 .eq. 1 ) then
             if ( ibvtyp .eq. 0 .and. itypea .eq. 0 
     &             .and. itypeb .eq. 2
     &             .and. last .eq. 1 .and. nipv .eq. 2 ) then 
               call sifuw1(
     &         info,   nipv,   num,    last,
     &         itypea, itypeb,   ibvtyp,
     &         val,    ilab,   fcore,  ibitv,
     &         nlist,  ierr, max1e, moints  )
               do 202 i=1,max1e
                   mstrmat(i) = mstrmat(i) + moints(i)
202            continue
               num1 = num
             else if ( ibvtyp .eq. 0 .and. itypea .eq. 0 
     &             .and. itypeb .eq. 1
     &             .and. last .eq. 1 .and. nipv .eq. 2 ) then
               call sifuw1(
     &         info,   nipv,   num,    last,
     &         itypea, itypeb,   ibvtyp,
     &         val,    ilab,   fcore,  ibitv,
     &         nlist,  ierr, max1e, moints  )
               do 203 i=1,max1e
                   mstrmat(i) = mstrmat(i) + moints(i)
203            continue
             else
             end if
             if ( ierr .ne. 0 ) call bummer(
     &       'iwfmt: siffw1 ierr=',ierr,faterr)
c
         end if
         goto 200
       end if
c
      return
      end
c deck prt2e
      subroutine prt2e(
     & nlist,  ntape,  info,   buf,
     & ilab,   val,    ibitv, max2e, mstrmat )
c
c  read and dump the 2-e integral records.
c
c  21-oct-90 written by ron shepard.
c
       implicit none
      integer   nipv,   msame,   nmsame,   nomore
      parameter(nipv=4, msame=0, nmsame=1, nomore=2)
c
      integer nlist, ntape, max2e
      real*8 buf(*), val(*)
      integer info(*), ilab(nipv,*), ibitv(*)
c
      integer ibvtyp, ierr, ifmt, itypea, itypeb, last, num,
     & num1, i, j, k, l
      integer   iretbv
      parameter(iretbv=-1)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      real*8  moints(max2e)
      real*8  mstrmat(max2e)
      last = msame
      do 210 i=1,max2e
        mstrmat(i) = 0d0
210   continue
200   if ( last .ne. nomore ) then
c
         call sifrd2(
     &    ntape,  info,   nipv,   iretbv,
     &    buf,    num,    last,   itypea,
     &    itypeb,    ibvtyp, val,
     &    ilab,   ibitv,  ierr )
c
         if ( ierr.ne.0 ) call bummer(
     &    'iwfmt: sifrd2 ierr=',ierr,faterr)
c
c        print *,"CALLING SIFUW2"
         call sifuw2(
     &    info, nipv, num, last,
     &    itypea, itypeb,  ibvtyp,
     &    val, ilab, ibitv, nlist,
     &    ierr, max2e, moints )
c
         do 201 i=1, max2e
                 mstrmat(i) = mstrmat(i) + moints(i)
201      continue
         if ( ierr.ne.0 ) call bummer(
     &    'iwfmt: siffw2 ierr=',ierr,faterr)
c
         goto 200
c         print *, "IF LOOP"
      end if
c
c      print *,"END"
      return
      end

      subroutine  prtarr1e( nlist, aoints, info,core, lencor,
     .                      nbpsy,nsym)
      implicit none
      integer nlist,aoints,info(6),lencor,nsym,nbpsy(nsym)
      real*8  core(*),fcore
      integer  icd(20),itypea,maxtypeb(0:2), nntot,isym,itypeb
      integer  nttot,symoff(8,8),mapin(511),symb(511),kntin(36)
      integer  lasta,lastb,nrec,ierr,ntot,i,btypes(0:41),last
      integer  choice
      character*8  inttype
      nttot=0
      nntot=0
      ntot=0 
 100  write(0,*) ' printing in fixed point format (1)'
      write(0,*) ' printing in floating point format (2)'
      read(5,*) choice
      if (choice.ne.1 .and. choice.ne.2) goto 100
      do isym=1,nsym
          nttot=nttot+nbpsy(isym)*(nbpsy(isym)+1)/2
          nntot=nntot+nbpsy(isym)*nbpsy(isym)
          ntot=ntot+nbpsy(isym)
      enddo
      do i=1,ntot
        mapin(i)=i
      enddo 
      maxtypeb(0)=9 
      maxtypeb(1)=41
      maxtypeb(2)=9 
c
c     1:  buffer (info(2))
c     2:  values (info(3))
c     3:  labelsi(2,info(3))
c     4:  array  (nntot)
c
      icd(1) = 1   
      icd(2) = icd(1)+info(2)
      icd(3) = icd(2)+info(3)
      icd(4) = icd(3)+2*info(3)
      icd(5) = icd(4)+ nttot
      icd(6) = icd(5)+ nntot
      if (icd(6).gt.lencor) 
     .   call bummer('insufficient memory in pretarr1e',icd(5),2)

      do itypea=0,0
        do itypeb=0,maxtypeb(itypea)
         rewind(aoints)
         call sifskh(aoints,ierr)
         if (ierr.ne.0) call bummer('sifskh failed',ierr,2)
         call wzero(nntot,core(icd(5)),1) 
         fcore=0
         btypes(0:41)=0
         btypes(itypeb)=1
         call sifr1n(aoints,info,itypea,itypeb,btypes,
     .               core(icd(1)),core(icd(2)),core(icd(3)),
     .               nsym,nbpsy,symoff,mapin,nttot,core(icd(5)),
     .               fcore,symb,kntin,lasta,lastb,last,nrec,ierr)

        if (ierr.eq.-4) cycle
        if (ierr.ne.0) call bummer('integral read failure',0,2)
        call siftyp(itypea,itypeb,inttype)
        call plblks(' printing 1e integrals: '//inttype,core(icd(5)),
     .                nsym,nbpsy,'AO',2,nlist)
        enddo
      enddo  
      return
      end 
c deck sifuw1
      subroutine sifuw1(
     & info,    nipv,    num,     last,
     & itypea,  itypeb,  ibvtyp,
     & values,  labels,  fcore,   ibitv,
     & nlist,   ierr, max1e, moints )
c
c  UNformatted write of a 1-e buffer.
c
c  input:
c  info(*) = info array for this file.
c  nipv = number of integers per value.
c       = 1 two orbital labels are packed in each labels(*) entry.
c       = 2 one orbital label is packed in each labels(*) entry.
c  num = number of values to place in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format to use for the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = integral values.
c  labels(1:nipv,1:num) = integral labels
c           note: if ifmt=0, then as many as ((nipv*n2max+3)/4)*4
c                 elements of labels(*) are referenced.
c  fcore = frozen core contribution.
c  ibitv(*) = unpacked bit vector (referenced only if ibvtyp.ne.0).
c             note: as many as ((n2max+63)/64)*64 elements of this
c                   array are referenced.
c  ierr = error return code. 0 for normal return.
c
c  24-apr-92 nipv added to the output record. -rls
c  16-aug-89 num=0 short return added.
c  26-jun-89 written by ron shepard.
c
       implicit logical(a-z)
      integer  nlist,  nipv,   num,    itypea, itypeb,
     & last,   ifmt,   ibvtyp, ierr, max1e
      real*8   values(*),      fcore
      integer  info(*),        labels(nipv,*), ibitv(*)
c
      real*8   moints(max1e)
      integer  l1rec,  n1max,  lenpl,  lenbv,  l1recx,
     & lab1,   ifmtv,  i, j, h , l,Ind2Val
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
      character*20 fstring(2)
      data fstring /'(1x,1pe20.12, i7  )',
     .              '(1x,1pe20.12, 2i4 )'/
c
6021  format(1x,1pe20.12, i7  )
6022  format(1x,1pe20.12, 2i4 )
      ierr  = 0
      l1rec = info(2)
      n1max = info(3)
      ifmt  = info(6)
c
c     # check parameters for consistency...
c
      if(num.gt.n1max)then
         call bummer('siffw1: num=',num,wrnerr)
         ierr = -1
         return
      elseif(itypea.lt.0 .or. itypea.gt.7)then
         call bummer('siffw1: itypea=',itypea,wrnerr)
         ierr = -1
         return
      elseif(itypeb.lt.0 .or. itypeb.gt.1023)then
         call bummer('siffw1: itypeb=',itypeb,wrnerr)
         ierr = -1
         return
      elseif(last.lt.0 .or. last.gt.3)then
         call bummer('siffw1: last=',last,wrnerr)
         ierr = -1
         return
      elseif(ibvtyp.lt.0 .or. ibvtyp.gt.7)then
         call bummer('siffw1: ibvtyp=',ibvtyp,wrnerr)
         ierr = -1
         return
      endif
c
      if(ifmt.eq.0)then
         lenpl=(num+3)/4
      elseif(ifmt.eq.1)then
         lenpl=(num+1)/2
      else
         call bummer('siffw1: ifmt=',ifmt,wrnerr)
         ierr = -1
         return
      endif
      if(ibvtyp.ne.0)then
         lenbv=(n1max+63)/64
      else
         lenbv=0
      endif
      l1recx=(2+num+lenpl+lenbv)
      if( l1recx .gt. l1rec )then
         call bummer('siffw1: (l1recx-l1rec)=',(l1recx-l1rec),0)
c        ierr = -1
c        return
      endif
c
      lab1 = num + 4
c
c     # write out the dword information.
c
c     write(nlist,6010)
c     & num,    lab1,   ibvtyp, itypea,
c     & itypeb, ifmt,   last,   nipv
6010  format(1x,8i7)
c
      if ( nipv .eq. 1 ) then
         ifmtv=nipv
      elseif( nipv .eq. 2 ) then
         ifmtv=nipv
      else
         call bummer('siffw1: nipv=',nipv,wrnerr)
         ierr = -1
         return
      endif
c
c      write(nlist,6021) fcore
c
c      do 9 h = 1, num
c        do 10 i = 1, bsfcns
c          do 11 j=1, i
c            if (labels(1,h) == i .and. labels(2, h) == j) then 
c              moints(i) = values(h)
c              moints(j,i) = values(h)
c            end if
c11        continue
c10      continue 
c9     continue
			 do 8 l = 1, max1e
			   moints(l) = 0d0
8      continue
       do 9 h = 1, num
         do 10 i = 1, max1e
           if ( Ind2Val(labels(1,h),labels(2,h)) == i ) then
             moints(i) = values(h)
           end if
10       continue
9      continue
c
c      if ( ibvtyp .ne. 0 ) then
c         write(nlist,6030) ( ibitv(i), i = 1, num )
c      end if
6030  format(1x,20i2)
c
      return
      end
c deck sifuw2
      subroutine sifuw2(
     & info,    nipv,    num,     last,
     & itypea,  itypeb,  ibvtyp,
     & values,  labels,  ibitv,   nlist,
     & ierr, max2e, moints )
c
c  UNformatted write of a 2-e buffer.
c
c  input:
c  info(*) = info array for this file.
c  nipv = number of integers per value.
c       = 1 four orbital labels are packed in each labels(*) entry.
c       = 2 two orbital labels are packed in each labels(*) entry.
c       = 4 one orbital label is stored in each labels(*) entry.
c  num = number of values to place in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format to use for the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = integral values.
c  labels(1:nipv,1:num) = integral labels
c           note: if ifmt=0, then as many as ((nipv*n2max+7)/8)*8
c                 elements of labels(*) are referenced.
c  ibitv(*) = unpacked bit vector (referenced only if ibvtyp.ne.0).
c             note: as many as ((n2max+63)/64)*64 elements of this
c                   array are referenced.
c  ierr = error return code. 0 for normal return.
c
c  24-apr-92 nipv added to the output record. -rls
c  16-aug-89 num=0 short return added.
c  26-jun-89 written by ron shepard.
c
       implicit logical(a-z)
      integer  nipv,   num,    itypea, itypeb, last,
     & ifmt,   ibvtyp, nlist,  ierr, max2e
      real*8   values(*)
      real*8   moints(max2e)
      integer  info(*), labels(nipv,*), ibitv(*)
c
      integer  l2rec,  n2max,  lenpl,  lenbv,  l2recx,
     & lab1,   ifmtv,  h, i, j, k, l, m,Index2E
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
      character*20 fstring(3) 
      data fstring /'(1x,1pe20.12, i11 )',
     .              '(1x,1pe20.12, 2i7 )',
     .              '(1x,1pe20.12, 4i4 )'/
6021  format(1x,1pe20.12, i11 )
6022  format(1x,1pe20.12, 2i7 )
6024  format(1x,1pe20.12, 4i4 )
c
      ierr  = 0
      l2rec = info(4)
      n2max = info(5)
      ifmt  = info(6)
c
c     # check parameters for consistency...
c
      if(num.gt.n2max)then
         call bummer('siffw2: num=',num,wrnerr)
         ierr = -1
         return
      elseif(itypea.lt.0 .or. itypea.gt.7)then
         call bummer('siffw2: itypea=',itypea,wrnerr)
         ierr = -1
         return
      elseif(itypeb.lt.0 .or. itypeb.gt.1023)then
         call bummer('siffw2: itypeb=',itypeb,wrnerr)
         ierr = -1
         return
      elseif(last.lt.0 .or. last.gt.3)then
         call bummer('siffw2: last=',last,wrnerr)
         ierr = -1
         return
      elseif(ibvtyp.lt.0 .or. ibvtyp.gt.7)then
         call bummer('siffw2: ibvtyp=',ibvtyp,wrnerr)
         ierr = -1
         return
      endif
c
      if(ifmt.eq.0)then
         lenpl=(num+1)/2
      elseif(ifmt.eq.1)then
         lenpl=num
      else
         call bummer('siffw2: ifmt=',ifmt,wrnerr)
         ierr = -1
         return
      endif
      if(ibvtyp.ne.0)then
         lenbv=(n2max+63)/64
      else
         lenbv=0
      endif
      l2recx=(1+num+lenpl+lenbv)
      if( l2recx .gt. l2rec )then
         call bummer('siffw2: (l2recx-l2rec)=',(l2recx-l2rec),0)
c        ierr = -1
c        return
      endif
c
      lab1 = num + 2
c
c     # write out the dword information.
c
      do 20 i = 1, max2e
        moints(i) = 0d0
20    continue
c      write(nlist,6010)
c     & num,    lab1,   ibvtyp, itypea,
c     & itypeb, ifmt,   last,   nipv
6010  format(1x,8i7)
c
      if ( nipv .eq. 1 ) then
         ifmtv=1
      elseif( nipv .eq. 2 ) then
         ifmtv=2
      elseif( nipv .eq. 4 ) then
         ifmtv=3
      else
         call bummer('siffw2: nipv=',nipv,wrnerr)
         ierr = -1
         return
      endif
c      print *, num
c      do 10 i = 1, max2e
      do 11 h = 1, num
         i = Index2E(labels(1,h),labels(2,h),labels(3,h),labels(4,h))
         moints(i) = values(h)  
c            go to 10
c          end if
c11      continue
11    continue
c
c      do 10 h = 1, num
c        do 11 i = 1, bsfcns
c          do 12 j = 1, bsfcns
c            do 13 k = 1, bsfcns
c              do 14 l = 1, bsfcns
c                if (labels(1,h) == i .and. labels(2,h) == j .and.
c     &              labels(3,h) == k .and. labels(4,h) == l ) then
c  # All 4 permutations of i,j,k,l
c                  moints(i,j,k,l) = values(h)
c                  moints(j,i,k,l) = values(h)
c                  moints(j,i,l,k) = values(h)
c                  moints(i,j,l,k) = values(h)
c                endif
c14            continue
c13          continue
c12        continue
c11       continue
c10    continue
c
c      if(ibvtyp.ne.0)then
c         write(nlist,6030) ( ibitv(i), i = 1, num )
c      endif
6030  format(1x,20i2)
c
      return
c      print *, "End of sifuw2"
      end
c deck tsifrd1
      subroutine tsifrd1(
     & aoints,  info,    nipv,    iretbv,
     & buffer,  num,     last,    itypea,
     & itypeb,  ibvtyp,  values,
     & labels,  fcore,   ibitv,   ierr )
c
c  read and decode a 1-e integral record.
c
c  input:
c  aoints = input file unit number.
c  info(*) = info array for this file.
c  nipv = number of integers per value to be returned.
c       = 0 only unpack dword.  values(*), labels(*), and ibitv(*)
c           are not referenced.
c       = 1 return two orbital labels packed in each labels(*) entry.
c       = 2 return one orbital label in each labels(*) entry.
c  iretbv = bit vector request type.
c     if ( iretbv=0 ) then
c         null request, don't return ibitv(*).
c     elseif ( iretbv=ibvtyp ) then
c         request return of the bit-vector of type iretbv.
c     elseif ( iretbv=-1 .and. ibvtyp<>0 ) then
c         return any type of bit-vector that is on the record.
c     else
c        error. requested bit-vector is not available in buffer(*).
c     endif
c  buffer(1:l1rec) = packed  buffer.
c
c  output:
c  num = actual number of values in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format of the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = values (referenced only if nipv.ne.0).
c  labels(1:nipv,1:num) = integral labels
c           (referenced only if nipv.ne.0).
c           note: if ifmt=0, then as many as ((nipv*n2max+7)/8)*8
c                 elements of labels(*) are referenced.
c  fcore = frozen core contribution.
c  ibitv(*) = unpacked bit vector (referenced only if iretbv.ne.0).
c             note: as many as ((n1max+63)/64)*64 elements of this
c                   array are referenced.
c  ierr = return code. 0 for normal return.
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  26-jun-89 written by ron shepard.
c
       implicit logical(a-z)
      integer  aoints, nipv,   iretbv, num,   last,
     & itypea, itypeb,    ibvtyp, ierr
      integer  info(*),        labels(*),     ibitv(*)
      real*8   buffer(*),      values(*),     fcore
c
      integer  l1rec,  n1max
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr  = 0
      l1rec = info(2)
      n1max = info(3)
c
c     # read the input file...
c
      call seqrbf ( aoints, buffer, l1rec )
c     # seqrbf() does not return error codes.
      ierr = 0
c
c     # unpack the buffer...
c
      call tsifd1(
     & info,  nipv,   iretbv, buffer,
     & num,   last,   itypea, itypeb,
     & ibvtyp, values, labels,
     & fcore, ibitv,  ierr )
c
c      print *, num
      return
      end
c deck tsifd1
      subroutine tsifd1(
     & info,   nipv,    iretbv,  buffer,
     & num,    last,    itypea,  itypeb,
     & ibvtyp,  values,  labels,
     & fcore,  ibitv,   ierr )
c
c  decode a 1-e buffer.
c  buffer has the form:
c    dword // packed_values(*) // fcore // packed_labels(*) //
c          //packed_bit_vector(*)
c
c  input:
c  info(*) = info array for this file.
c  nipv = number of integers per value to be returned.
c       = 0 only unpack dword.  values(*), labels(*), and ibitv(*)
c           are not referenced.
c       = 1 return two orbital labels packed in each labels(*) entry.
c       = 2 return one orbital label in each labels(*) entry.
c  iretbv = bit vector request type.
c     if ( iretbv=0 ) then
c         null request, don't return ibitv(*).
c     elseif ( iretbv=ibvtyp ) then
c         request return of the bit-vector of type iretbv.
c     elseif ( iretbv=-1 .and. ibvtyp<>0 ) then
c         return any type of bit-vector that is on the record.
c     else
c        error. requested bit-vector is not available in buffer(*).
c     endif
c  buffer(1:l1rec) = packed  buffer.
c
c  output:
c  num = actual number of values in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = values (referenced only if nipv.ne.0).
c  labels(1:nipv,1:num) = integral labels
c           (referenced only if nipv.ne.0).
c           note: if ifmt=0, then as many as ((nipv*n1max+7)/8)*8
c                 elements of labels(*) are referenced.
c  fcore = frozen core contribution.
c  ibitv(*) = unpacked bit vector (referenced only if iretbv.ne.0).
c             note: as many as ((n1max+63)/64)*64 elements of this
c                   array are referenced.
c  ierr = error return code. 0 for normal return.
c
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  26-jun-89 written by ron shepard.
c
       implicit logical(a-z)
      integer  nipv,   iretbv, num,    itypea, itypeb,
     & last,   ifmt,   ibvtyp, ierr
      real*8 buffer(*), values(*), fcore
      integer info(*), labels(*), ibitv(*)
c
      integer  l1rec,  n1max,  lab1,   lenpl,  lenbv,
     & l1recx, nuw,    bv1
      integer  unpack(4)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr  = 0
c
      l1rec = info(2)
      n1max = info(3)
c
c     # unpack dword...
c
      call ulab16( buffer(1), unpack, 4 )
c      print *, "From ulab"
c      print *, num
      num    = unpack(1)
      lab1   = unpack(2)
      last   = mod(unpack(4),         4)
      ifmt   = mod(unpack(4)/2**2,    8)
      itypeb = mod(unpack(3),      1024)
      itypea = mod(unpack(3)/2**10,   8)
      ibvtyp = mod(unpack(3)/2**13,   8)
c
c     ifmt is now double checked
c
       if (ifmt.ne.info(6)) then
        call bummer(' sifd1: ifmt.ne.info(6)',0,2)
       endif
c
c     # if nipv=0 then only dword is unpacked...
c
      if ( nipv .eq. 0 ) return
c
      if ( ifmt .eq. 0 ) then
         lenpl=(num+3)/4
      elseif ( ifmt .eq. 1 ) then
c   # now 10bit packing
c        lenpl=(num+1)/2
         lenpl=(num/16)*5
      else
c        # illegal ifmt.
         ierr = -1
         return
      endif
      if ( ibvtyp .ne. 0 ) then
         lenbv=(n1max+63)/64
      else
         lenbv=0
      endif
      l1recx=(2+num+lenpl+lenbv)
      if ( l1recx .gt. l1rec ) then
c        # inconsistent l1rec.
         ierr = -2
         return
      endif
c
c     # unpack/copy the values(*)...
c
      call dcopy_wr( num, buffer(2), 1,  values, 1 )
c
      fcore = buffer( num + 2 )
c
c     # unpack the labels(*)...
c
      if ( ifmt .eq. 0 ) then
c        # 8-bit packing of orbital labels.
         if ( nipv .eq. 1 ) then
c           # 1 integer/value output.
            nuw=num
            call ulab16( buffer(lab1), labels, nuw )
         elseif ( nipv .eq. 2 ) then
c           # 2 integers/value output.
            nuw=2*num
            call ulab8( buffer(lab1), labels, nuw )
         else
c           # illegal nipv.
            ierr = -3
            return
         endif
      elseif ( ifmt .eq. 1 ) then
c        # 10-bit packing of orbital labels.
         if ( nipv .eq. 1 ) then
c           # 1 integer/value output.
            nuw=num
            call bummer('sifd1: ifmt=1 && nipv=1 unsupported',0,2)
         elseif ( nipv .eq. 2 ) then
c           # 2 integers/value output.
            nuw=2*num
            call ulab10( buffer(lab1), labels, nuw )
         else
c           # illegal nipv.
            ierr = -3
            return
         endif
      endif
c
      if ( iretbv .eq. 0 ) then
c        # ignore bit-vector processing.
         continue
      elseif ( (iretbv .eq. ibvtyp) .or.
     &    ( (iretbv .eq. -1) .and. (ibvtyp .ne. 0 ) ) ) then
c        # unpack the bit vector from the end of buffer(*).
         bv1=l1rec+1-lenbv
         nuw=num
         call ulab1( buffer(bv1), ibitv, nuw )
      elseif ( iretbv .eq. -1 ) then
c        # general bitvector request with ibvtyp=0.  not an error.
         continue
      else
c        # inconsistent ibvtyp.
         ierr = -4
         return
      endif
c
      return
      end
c
      integer function Ind2Val(i,j)
         integer :: i, j
c # ind2val = (i-1)i/2 + j
         if ( i .ge. j ) then
	       Ind2Val = (i - 1)*i/2 + j
         else
	       Ind2Val = (j - 1)*j/2 + i
	     end if
      end function Ind2Val
      integer function Index2E(i,j,k,l)
         implicit none
         integer :: i, j, k, l, Ind2Val
c
         Index2E = Ind2val(Ind2Val(i,j),Ind2Val(k,l))
      end function
