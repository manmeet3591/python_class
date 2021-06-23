;***************************************************************
; hdf4sds_3.ncl
;
; Concepts illustrated:
;   - Reading multiple files with data from Trmm satellite
;   - Binning and summing the data from each pass [swath]
;   - Averaging binned data
;   - Writing data to a NetCDF file using the easy but inefficient method
;***************************************************************

;*****************Load Libraries ************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
     
;********************************************************************
; Utility functions to get (1) date information (2) facilitate netCDF 
;********************************************************************
function parseFileName ( fNam:string )
; parse the file names and extract data information
;           1         2
; 01234567890123456789012345678
; MYDATML2.A2005364.1820.005.nc
local onNam, n, output, tmp_c
begin
   nNam = dimsizes( fNam )

   if (nNam.eq.1) then
       output    = new( 4,integer,"No_FillValue")
       tmp_c     = stringtochar(fNam)
       
       output(0) = stringtointeger((/tmp_c(10:13)/))   ; YYYY 
       output(1) = stringtointeger((/tmp_c(14:16)/))   ; MM
       output(2) = stringtointeger((/tmp_c(18:19)/))   ; HH
       output(3) = stringtointeger((/tmp_c(20:21)/))   ; MN        
   else
       output = new( (/nNam,4/),integer,"No_FillValue")
       do n=0,nNam-1
          tmp_c       = stringtochar(fNam(n))       
          output(n,0) = stringtointeger((/tmp_c(10:13)/))   ; YYYY 
          output(n,1) = stringtointeger((/tmp_c(14:16)/))   ; MM
          output(n,2) = stringtointeger((/tmp_c(18:19)/))   ; HH
          output(n,3) = stringtointeger((/tmp_c(20:21)/))   ; MN        
       end do
   end if
   
   return (output)   
end
;---------------------------------------------------
function create3d (x[*][*]:numeric, time[*]:numeric)
; convert a 2-dimensional to a 3-dimensional variable
; with time as a coordinate variable.
local dimx, ntim, nlat, mlon, x3d
begin
  dimx    = dimsizes(x)

  ntim    = dimsizes(time)
  nlat    = dimx(0)
  mlon    = dimx(1)
          ; perform expansion
  x3d     = conform_dims( (/ntim,nlat,mlon/), x, (/1,2/))
  copy_VarAtts(x, x3d)
  x3d!0   = "time"
  x3d&time=  time
  copy_VarCoords(x, x3d(0,:,:) )
  return(x3d)
end

;**************************************************************
;                     MAIN
;**************************************************************
begin
  vNam  = "Cloud_Water_Path_1621"    
  vNam  = "nearSurfZ"    
  vNam  = "correctZFactor"    

  PLOT  = True
  netCDF= True

;*****************************************************************
; binned grid definition     
; The grid must be constant spacing in lat and in lon.
; The spacing in the lat and lon directions may be different.
; so: 2x2, 3x5, etc
;*****************************************************************

  nlat  = 480                       ; grid coordinates
  mlon  = 1440   
	nlev  = 80  
  lat   = latGlobeFo(nlat,"lat","latitude","degrees_north")
  lat   = lat(::-1)         ; North to South
  lon   = lonGlobeFo(mlon,"lon","longitude","degrees_east")
  lon   = (/ lon - 180. /)  ; subtract 180 from all values 
  lon&lon = lon             ; update coordinates
	lev=new(80,float)
;  printVarSummary(lev)

      tunits  = "hours since 1990-1-1 00:00:0.0"
      time    = new((/1/),integer)
			time    = 2005
      time!0  = "time"
      time@units = tunits 

	do k=0,nlev-1
	
	lev(k) = 0.25 * (k+1)

	end do

;	print(lev)
;exit
;end
; to add lev

;*****************************************************************
; Variables to hold binned quantities
;*****************************************************************
  GBIN  = new ( (/nlat,mlon/), float ) 
  GKNT  = new ( (/nlat,mlon/), integer ) 
	Z     = new ( (/nlat,mlon,nlev,1/),float )
	Z = 0.0
;*****************************************************************
; File info                          
;*****************************************************************

  diri = "./"
  fili = systemfunc("cd "+diri+" ; ls 2A25*.HDF")
  nfil = dimsizes( fili )
  ;print("nfil="+nfil)  
 

      diro    = "./"                   ; output directory
      filo    = "TrmmBin_bin_sum"    ; output file name
      fout    = diro+filo+".nc"

      system("/bin/rm -f "+fout)       ; remove any pre-existing file
      ncdf    = addfile(fout ,"c")     ; open output netCDF file

      filedimdef(ncdf,"time",-1,True)  ; make time UNLIMITED dimension
      ncdf->time     = time

;*****************************************************************
; Loop over the files: Sum the quantities
;*****************************************************************


		do k=0,nlev-1
  vNam  = "correctZFactor"    
	  GBIN  = 0.0                      ; initialization
	  GKNT  = 0
			
		do nf=0,nfil-1
     print(nf+"   "+fili(nf))
     f       = addfile(diri+fili(nf), "r")
                                     ; read data
     lat2d   =  f->Latitude 
     lon2d   =  f->Longitude
     x       =  f->$vNam$   
     nx      = product(dimsizes(x))
     bin_sum(GBIN,GKNT,lon,lat   \
            ,ndtooned(lon2d), ndtooned(lat2d),ndtooned(x(:,:,k)) )

     delete(  x  )                   ; size may change
     delete(lat2d)  
     delete(lon2d)
		end do  
;*****************************************************************
; Perform averaging
;*****************************************************************
                                    ; avoid division by 0.0
  GKNT     = where(GKNT.eq.0 , GKNT@_FillValue, GKNT)

  GBIN     = GBIN/GKNT
 
  GBIN!0   = "lat"
  GBIN!1   = "lon"
  GBIN&lat =  lat
  GBIN&lon =  lon

  Z!0   = "lat"
  Z!1   = "lon"
  Z!2   = "lev"
	Z!3   = "time"
  Z&lat =  lat
  Z&lon =  lon
  Z&lev =  lev
	Z&time = time
 copy_VarCoords(GBIN, GKNT)        ; copy coords

  if (isfilevaratt(f, vNam, "long_name")) then
      GBIN@long_name = "BINNED: "+vNam
      GKNT@long_name = "BINNED COUNT: "+vNam
  end if

  if (isfilevaratt(f, vNam, "units")) then
      GBIN@units     = f->$vNam$@units
  end if


;Z(:,:,k,:) = create3d(GBIN,time)

			vNam = vNam + k
      ncdf->$vNam$   = create3d(GBIN,time) 

            
end do

end
