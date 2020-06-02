program read_netcdf

!-- a generic netcdf file reader that uses the netcdf fortran 90 interface.
!-- use compncdf to compile, and be sure that the following symbolic link
!-- has been put into place:
!--   ln -s /usr/local/absoft/include/netcdf.mod .

!-- For more info, see:  
!--   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90/

!-- Mark Branson, 2/2007

  use netcdf
  implicit none

  integer :: i, n, k      ! loop counters
  integer :: ncid, ndims, nvars, nglobalatts, unlimdimid
  integer :: dimlen, attlen, att_type
  integer :: vardims, varnatts, vartype
  integer, dimension(nf90_max_var_dims) :: vardimids
  integer, dimension(:), allocatable :: dimsize
  real, dimension(:), allocatable :: float_1din
  real, dimension(:,:), allocatable :: float_2din
  logical :: exans
  character(len=80) :: infile, dname, attname, attvalue, varname
  character(len=20) :: vtype
  character(len=80), dimension(:), allocatable :: dimname

!-- ensure that the chosen input file actually exists
  !write(*,*) 'Enter a netcdf file:'
  !read(*,'(a80)') infile
  infile = 'psib_199707p001.pbp2.nc'
  !infile = 'seminar/simple_xy.nc'
  inquire(file=trim(infile),exist=exans)
  if (exans) then
    write(*,*) '*** Info for ',trim(infile),' ***'
    
!-- open the netcdf input file
    call check( nf90_open(infile, NF90_NOWRITE, ncid) )
  
!-- get some general information:  
!--   1) number of dimensions
!--   2) number of variables
!--   3) number of global attributes
!--   4) id of unlimited dimension

    call check( nf90_inquire(ncid, nDims, nVars, nGlobalAtts, unlimdimid) )

    write(*,*) '# of dimensions = ',ndims
    write(*,*) '# of variables = ',nvars
    write(*,*) '# of global attributes = ',nglobalatts
    if (unlimdimid == -1) then
      write(*,*) 'There is no unlimited dimension defined for this file.'
    else
      write(*,*) 'id of unlimited dimension = ',unlimdimid
    endif
    write(*,*)
    
    allocate(dimname(ndims),dimsize(ndims))
    
!-- print out the global attributes
    do k = 1,nglobalatts
       call check( nf90_inq_attname(ncid, nf90_global, k, attname) )
       write(*,*) '--> global attribute name = ',trim(attname)
       call check( nf90_inquire_attribute(ncid, nf90_global, attname,  &
                    xtype = att_type, len = attlen) )
       write(*,*) '--> len = ',attlen, ' type = ',att_type

       do i = 1,80
         attvalue(i:i) = ' '  ! get goofy output if we don't do this
       enddo
       if (att_type == NF90_CHAR) then
         call check( nf90_get_att(ncid, nf90_global, attname, attvalue) )
         write(*,*) '                     value = ',trim(attvalue)
       endif
    enddo

!-- print out the dimension info
    do n = 1,ndims
       call check( nf90_inquire_dimension(ncid, n, dname, dimlen) )
       write(*,*) '>>> dimension name = ',trim(dname), ' length = ',dimlen
       dimname(n) = dname; dimsize(n) = dimlen
    enddo
    write(*,*)

!-- print out the variable info
    do n = 1,nvars
       call check( nf90_inquire_variable(ncid, n, varname, vartype,  &
                    vardims, vardimids, varnatts) )
       !write(*,*,advance='no') '>>> variable name = ',varname, '('
       write(*,'(3a)',advance='no') '>>> variable name = ',trim(varname), '('
       do i = 1,vardims
         if (i < vardims) then
           write(*,'(2a)',advance='no') trim(dimname(vardimids(i))), ','
         else
           write(*,*) trim(dimname(vardimids(i))), ')'
         endif
       enddo
       if (vartype == NF90_BYTE) then
         vtype = 'BYTE'
       elseif (vartype == NF90_CHAR) then
         vtype = 'CHARACTER'
       elseif (vartype == NF90_SHORT) then
         vtype = 'SHORT'
       elseif (vartype == NF90_INT) then
         vtype = 'INTEGER'
       elseif (vartype == NF90_FLOAT) then
         vtype = 'FLOAT'
       elseif (vartype == NF90_DOUBLE) then
         vtype = 'DOUBLE'
       endif
       write(*,*) '       type = ',trim(vtype)
       do k = 1,varnatts
          call check( nf90_inq_attname(ncid, n, k, attname) )
          write(*,*) '       *attribute = ',attname
       enddo
!-- sample:  grab the data if the variable is a 1d float and save it
!-- *** alternative to using dynamic arrays:  use ncdump on the file
!--     before you read it so you know the size of all of the dimensions
       if (vartype == nf90_float .and. vardims == 1) then
         allocate(float_1din(dimsize(vardimids(1))))
         call check( nf90_get_var(ncid, n, float_1din) )
         write(*,*) 'float_1din = ',float_1din
         deallocate(float_1din)
       endif
    enddo
  else 
    write(*,*) 'The selected netcdf file ',trim(infile),' could not be opened.'
  endif

  deallocate(dimname,dimsize)
  
  contains

  subroutine check(status)

    integer, intent(in) :: status
    
    if (status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    endif
    
  end subroutine check  
  
end program read_netcdf
