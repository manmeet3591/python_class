module param_module
  use type_module
  implicit none
   
  !! file unit
  integer(kind=i4b), parameter, public :: lunm=11, lugb=12, lugbo=13, luo3=14, luvc=15, lugboi=16, lutb=17, luso=21

  !! GFS sigma file parameter
  !! These parameters depend on GFS version and experimental setting
  integer(kind=i4b),parameter,public :: nvcoord=2, & ! number of vertical coordinate parameter
       &                                idsl=2,    & ! vertical coordinate type (2:hyblid sigma-p)
       &                                idvc=2,    & ! how to define model level (2:mean)
       &                                ntrac=3,   & ! number of tracer
       &                                idvt=21,   & ! tracer variable id (21:vapor, ozone, cloud content)
       &                                idvm=0       ! mass variable id (0: default in GFS)

end module param_module
