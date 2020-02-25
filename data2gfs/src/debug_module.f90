module debug_module
  use type_module

contains

  subroutine write_grads(unit, fname, data, idrt)
    integer(kind=i4b), intent(in) :: unit, idrt
    character(len=*), intent(in)  :: fname
    real(kind=sp), intent(in)     :: data(:,:,:)
    integer(kind=i4b) :: xn,yn,zn
    real(kind=sp) :: slat(size(data,2)), wlat(size(data,2))

    xn = size(data,1)
    yn = size(data,2)
    zn = size(data,3)

    open(unit,file=fname,recl=4*xn*yn*zn,form='unformatted',access='direct')
    write(unit,rec=1) data(:,:,:)
    close(unit)

    open(unit,file=fname//'.ctl')
    write(unit,'("dset ",a)') fname
    write(unit,'("options yrev")')
    write(unit,'("undef -9.99E+33")')
    write(unit,'("title debug_module")')
    write(unit,'("xdef",i6," linear",2f12.6)') xn,0.d0,360.d0/xn
    if(idrt.eq.0) then
       write(unit,'("ydef",i6," linear",2f12.6)')&
            yn,-90.d0,180.d0/(yn-1)
    elseif(idrt.eq.256) then
       write(unit,'("ydef",i6," linear",2f12.6)')&
            yn,-90.d0*(yn-1)/yn,180.d0/yn
    elseif(idrt.eq.4) then
       call splat(idrt,yn,slat,wlat)
       write(unit,'("ydef",i6," levels")') yn
       write(unit,'(5f12.6)') 180.d0/acos(-1.d0)*asin(dble(slat(yn:1:-1)))
    endif
    write(unit,'("zdef",i6," linear 1 1")') zn
    write(unit, '(a)') 'tdef 1 linear 00z01jan1900 1yr'
    write(unit,'("vars",i6)') 1
    write(unit,'("var  ",i3," 99 **")') zn
    write(unit,'("endvars")')
    close(unit)

  end subroutine write_grads

end module debug_module
