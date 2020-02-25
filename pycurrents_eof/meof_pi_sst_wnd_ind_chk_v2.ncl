;*********************************************
load "/usr/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "/usr/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "/usr/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/usr/lib/ncarg/nclscripts/csm/shea_util.ncl"

f1=addfile("sst_yrmean_1951_2124_remap_ind.nc","r")
f2=addfile("esm_pi_taux_1951_2124_ind.nc","r")
f3=addfile("esm_pi_tauy_1951_2124_ind.nc","r")



sst=f1->SST
taux=f2->TAUX
tauy=f3->TAUY

lat=f1->lat
lon=f1->lon
TIME=f1->time

mlat=dimsizes(lat)
mlon=dimsizes(lon)
ntim=dimsizes(TIME)


sst1=doubletofloat(f1->SST(:,:,:))
taux1=doubletofloat(f2->TAUX(:,:,:))
tauy1=doubletofloat(f3->TAUY(:,:,:))


copy_VarCoords(sst,sst1)
copy_VarCoords(taux,taux1)
copy_VarCoords(tauy,tauy1)


mean_sst  = dim_avg_n(sst1,0)          ; (time)
mean_taux = dim_avg_n(taux1,0)          ; (time)
mean_tauy = dim_avg_n(tauy1,0)          ; (time)

var_sst = dim_variance_n_Wrap(sst1,0)          ; (time)
var_taux = dim_variance_n_Wrap(taux1,0)          ; (time)
var_tauy = dim_variance_n_Wrap(tauy1,0)          ; (time)



sst2=sst1
taux2=taux1
tauy2=tauy1

   do j=0,39
	do i=0,70
	  sst2(:,j,i)  = (sst1(:,j,i)  - mean_sst(j,i))/sqrt(var_sst(j,i))          ; (lon,time)
	  taux2(:,j,i) = (taux1(:,j,i) - mean_taux(j,i))/sqrt(var_taux(j,i))          ; (lon,time)
	  tauy2(:,j,i) = (tauy1(:,j,i) - mean_tauy(j,i))/sqrt(var_tauy(j,i))          ; (lon,time)
	end do
   end do

sst2_n=sst2(lat|:,lon|:,time|:)
taux2_n=taux2(lat|:,lon|:,time|:)
tauy2_n=tauy2(lat|:,lon|:,time|:)

;printVarSummary(sst2_n)
;printVarSummary(taux2_n)
;printVarSummary(tauy2_n)

 cdata = new ((/mlat*mlon*3,ntim/),typeof(sst2_n), getFillValue(sst2_n))
; printVarSummary(cdata)


sst2_nr = reshape(sst2_n,(/mlat*mlon,ntim/))
taux2_nr = reshape(taux2_n,(/mlat*mlon,ntim/))
tauy2_nr = reshape(tauy2_n,(/mlat*mlon,ntim/))

latmlon = mlat*mlon

do mt=0,ntim-1
   cdata(0:latmlon-1,:) = sst2_nr
   cdata(latmlon:2*latmlon-1,:) = taux2_nr
   cdata(2*latmlon:3*latmlon-1,:) = tauy2_nr
end do


;printVarSummary(cdata)
;
neof=3
eof_cdata    = eofunc(cdata,neof, False)     ; (neof,3*mlon)
printVarSummary(eof_cdata)

eof_sst = eof_cdata(:,0:latmlon-1) 
eof_taux = eof_cdata(:,latmlon:2*latmlon-1) 
eof_tauy = eof_cdata(:,2*latmlon:3*latmlon-1) 


eof_cdata_sstr  = reshape(eof_sst,(/neof,mlat,mlon/))
eof_cdata_tauxr = reshape(eof_taux,(/neof,mlat,mlon/))
eof_cdata_tauyr = reshape(eof_tauy,(/neof,mlat,mlon/))
;;printVarSummary(eof_cdatar)
;
;;
;eof_ts = eofunc_ts(cdata,eof_cdata,False)  ; (neof,time)
eof_ts_sst  = eofunc_ts(sst2_nr,reshape(eof_cdata_sstr,(/neof,latmlon/)),False)  ; (neof,time)
eof_ts_taux = eofunc_ts(taux2_nr,reshape(eof_cdata_tauxr,(/neof,latmlon/)),False)  ; (neof,time)
eof_ts_tauy = eofunc_ts(tauy2_nr,reshape(eof_cdata_tauyr,(/neof,latmlon/)),False)  ; (neof,time)
printVarSummary(eof_ts_sst)

eof_cdata_sstr!1="lat"
eof_cdata_sstr&lat=lat
eof_cdata_sstr!2="lon"
eof_cdata_sstr&lon=lon
;printVarSummary(eof_cdata_sstr)

eof_cdata_tauxr!1="lat"
eof_cdata_tauxr&lat=lat
eof_cdata_tauxr!2="lon"
eof_cdata_tauxr&lon=lon

eof_cdata_tauyr!1="lat"
eof_cdata_tauyr&lat=lat
eof_cdata_tauyr!2="lon"
eof_cdata_tauyr&lon=lon



eof_ts_sst!0="neof"
;eof_ts_sst&neof=neofa
eof_ts_sst!1="time"
eof_ts_sst&time=TIME
;printVarSummary(eof_ts_sst)

eof_ts_taux!0="neof"
;eof_ts_taux&neof=neofa
eof_ts_taux!1="time"
eof_ts_taux&time=TIME
;printVarSummary(eof_ts_sst)

eof_ts_tauy!0="neof"
;eof_ts_tauy&neof=neofa
eof_ts_tauy!1="time"
eof_ts_tauy&time=TIME
;printVarSummary(eof_ts_sst)


diro="./"
filo="esm_eof.nc"
system("/bin/rm -f "+diro+filo)      ; rm any pre-existing file
setfileoption("nc","Format","LargeFile")
f = addfile(diro+filo, "c")
f->eof_ts_sst=eof_ts_sst
f->eof_ts_taux=eof_ts_taux
f->eof_ts_tauy=eof_ts_tauy
f->eof_cdata_sst =eof_cdata_sstr
f->eof_cdata_taux =eof_cdata_tauxr
f->eof_cdata_tauy =eof_cdata_tauyr


