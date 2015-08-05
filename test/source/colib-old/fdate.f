      function fdate()
      character*24 fdate
if rs6000 .and. (.not. extname)
      character*24 fdate_
      external     fdate_
c     # call the fortran library routine
      fdate = fdate_()
mdc*else
c     # call the C function with the argument exposed.
      call c_fdate( fdate )
mdc*endif
      return
      end
