
        real*8 FUNCTION DNRM2_WR ( N, X, INCX )
         implicit none
         real*8 dnrm2,X
*@if defined (int64) && defined (noblas64)
*         integer*8 n,incx
*         integer*4 ln,lincx
*         integer*8 maxval
*         parameter (maxval=2**31-1)
*
*         if (n.gt.maxval) then
*           call bummer('dnrm2_wr, n exceeds maxval',n,2)
*         endif
*         ln=n
*         lincx=incx 
*         dnrm2_wr = dnrm2(ln,x,lincx)
*@else
         integer*4 n,incx
         dnrm2_wr= dnrm2(n,x,incx)
*@endif 
        return
        end




      SUBROUTINE DGEMM_WR ( TRANSA, TRANSB, M, N, K, ALPHA, 
     .    A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      implicit none
      CHARACTER*1        TRANSA, TRANSB
      real*8   ALPHA, BETA
c     .. Array Arguments ..

*@if defined (int64) && defined (noblas64) 
*      INTEGER*8            M, N, K, LDA, LDB, LDC
*      INTEGER*4            LM, LN, LK, LLDA, LLDB, LLDC
*      real*8   A( LDA, * ), B( LDB, * ), C( LDC, * )
*      
*       lm=m
*       ln=n
*       lk=k
*       llda=lda
*       lldb=ldb
*       lldc=ldc
*      
*       call dgemm(TRANSA, TRANSB, lM, lN, lK, ALPHA, A, lLDA, B, lLDB,
*     .    BETA, C, lLDC )
* 
*@else
      INTEGER*4            M, N, K, LDA, LDB, LDC
       real*8   A( LDA, * ), B( LDB, * ), C( LDC, * )
       call dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     .    BETA, C, LDC )
*@endif
       
       return
       end

    

      SUBROUTINE DGEMV_WR ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
c     .. Scalar Arguments ..
      implicit none
      real*8   ALPHA, BETA
      CHARACTER*1        TRANS
c     .. Array Arguments ..
 
*@if defined (int64) && defined (noblas64) 
*      INTEGER*8            INCX, INCY, LDA, M, N
*      INTEGER*4            lINCX, lINCY, lLDA, lM, lN
*      real*8   A( LDA, * ), X( * ), Y( * )
*       lincx=incx
*       lincy=incy
*       llda = lda
*       lm=m
*       ln=n
*       call  DGEMV ( TRANS, lM, lN, ALPHA, A, lLDA, X, lINCX,
*     $                   BETA, Y, lINCY )
*@else
       INTEGER*4            INCX, INCY, LDA, M, N
       real*8   A( LDA, * ), X( * ), Y( * )
       call  DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*@endif
       return
       end


      real*8 function ddot_wr(n,dx,incx,dy,incy)
      implicit none
      real*8 ddot, dx(*),dy(*),dtemp
*@if defined (int64) && defined (noblas64) 
*      integer*8 n,incx,incy
*      integer*4 ln,lincx,lincy
*       
*      ln=n
*      lincx=incx
*      lincy=incy 
*      ddot_wr=ddot(ln,dx,lincx,dy,lincy)
*@else
      integer*4 n,incx,incy
      ddot_wr=ddot(n,dx,incx,dy,incy)
*@endif
      return
      end



      subroutine  dcopy_wr(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 dx(*),dy(*)
*@if defined (int64) && defined (noblas64) 
*      integer*8 n,incx,incy
*      integer*4 ln,lincx,lincy
*      ln=n
*      lincx=incx
*      lincy=incy 
*       call  dcopy(ln,dx,lincx,dy,lincy)
*@else
      integer*4 n,incx,incy
       call  dcopy(n,dx,incx,dy,incy)
*@endif
      return
      end



      subroutine  dscal_wr(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 da,dx(*)
*@if defined (int64) && defined (noblas64)
*      integer*8 n,incx
*      integer*4 ln,lincx
*      ln=n
*      lincx=incx
*      call dscal(ln,da,dx,lincx)
*@else
      integer*4 n,incx
      call dscal(n,da,dx,incx)
*@endif

      return
      end






      SUBROUTINE DGER_WR  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
c     .. Scalar Arguments ..
      implicit none
      real*8   ALPHA
c     .. Array Arguments ..
*@if defined (int64) && defined (noblas64)
*      INTEGER*8          INCX, INCY, LDA, M, N
*      INTEGER*4          lINCX, lINCY, lLDA, lM, lN
*      real*8   A( LDA, * ), X( * ), Y( * )
*       lincx=incx
*       lincy=incy 
*       llda = lda
*       lm = m
*       ln = n
*       call dger( lM, lN, ALPHA, X, lINCX, Y, lINCY, A, lLDA )
*     
*@else
      INTEGER*4          INCX, INCY, LDA, M, N
       real*8   A( LDA, * ), X( * ), Y( * )
      call dger( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*@endif
      return
      end

      real*8 function dasum_wr(n,dx,incx)
      implicit none
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      real*8 dasum,dx(*),dtemp
*@if defined (int64) && defined (noblas64)
*      integer*8 n,incx
*      integer*4 ln,lincx
*       ln=n
*       lincx=incx
*      dasum_wr=dasum(ln,dx,lincx)
*@else
      integer*4 n,incx
      dasum_wr=dasum(n,dx,incx)
*@endif
      return
      end



      subroutine daxpy_wr(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 dx(*),dy(*),da
*@if defined (int64) && defined (noblas64)
*      integer*8 n,incx,incy
*      integer*4 ln,lincx,lincy
*       ln=n
*       lincx=incx
*       lincy=incy
*       call daxpy(ln,da,dx,lincx,dy,lincy)
*@else
      integer*4 n,incx,incy
      call daxpy(n,da,dx,incx,dy,incy)
*@endif
      return
      end

c


      subroutine  dswap_wr (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none
      real*8 dx(*),dy(*),dtemp
*@if defined (int64) && defined (noblas64)
*      integer*8 n,incx,incy
*      integer*4 ln,lincx,lincy
*       ln=n
*       lincx=incx
*       lincy=incy
*       call dswap(ln,dx,lincx,dy,lincy)
*@else
      integer*4 n,incx,incy
      call dswap (n,dx,incx,dy,incy)
*@endif
      return
      end


*@if defined (int64) 
*      integer*8 function idamax_wr(n,dx,incx)
*@else
      integer*4 function idamax_wr(n,dx,incx)
*@endif
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      real*8 dx(*),dmax
*@if defined (int64) && defined (noblas64)
*      integer*8 n, incx
*      integer*4 ln, lidamax,lincx
*       ln=n
*       lincx=incx
*       lidamax=idamax(ln,dx,lincx)
*       idamax_wr=lidamax
*@else
      integer*4 n, incx
       idamax_wr=idamax(n,dx,incx)
*@endif
      return
      end


      SUBROUTINE DSPR_WR  ( UPLO, N, ALPHA, X, INCX, AP )
       implicit none
c     .. Scalar Arguments ..
      real*8   ALPHA
c     .. Array Arguments ..
      real*8   AP( * ), X( * )
      CHARACTER*1        UPLO
*@if defined (int64) && defined (noblas64)
*      INTEGER*8          INCX, N
*      INTEGER*4         LINCX,LN
*      LINCX=INCX
*      LN=N
*c
*      call DSPR  ( UPLO, lN, ALPHA, X, lINCX, AP )
*@else
      INTEGER*4          INCX, N
      call DSPR  ( UPLO, N, ALPHA, X, INCX, AP )
*@endif
      return
      end

