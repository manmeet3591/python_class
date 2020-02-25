#!/bin/bash
#----------------------------------------------
# SCRIPT ID:2
# for DPAC GANAL binary file
#
# 注意:250hPa以上のTTDは定数。雲水70hPa以上は未定義。
#      ganal_to_bin.pyを使って雲水70hPa以上はゼロにする。オゾンはゼロ。
#      オゾンはo3clim=.true.にすると気候値のプロファイルを東西一様に与える。
#----------------------------------------------

export IDATE=2009093012
export JCAP=190

. ../configure

LEVTYP=pl

year=`echo $IDATE | cut -c 1-4`
month=`echo $IDATE | cut -c 5-6`
day=`echo $IDATE | cut -c 7-8`
hour=`echo $IDATE | cut -c 9-10`

[ ! -d $DATA ] && mkdir -p $DATA
PWD0=$(pwd)
cd $DATA

[ ! -d ${COMOUT}/ganal/${IDATE} ]&&mkdir -p ${COMOUT}/ganal/${IDATE}

ln -sf /mnt/climf/jma/ganal/${year}/${year}${month}/${IDATE}00 input.bin
ln -sf ${OROGRAPHY} input.orography.grib
ln -sf ${O3CLIM} input.o3clim.txt
ln -sf ${COMOUT}/ganal/${IDATE}/siganl.${IDATE}.t${JCAP} output.siganl

cat <<EOF > data2gfs_namelist
 &gfsparam
  jcap=${JCAP},lonr=${LONB},latr=${LATB},levs=${LEVS},o3clim=.true.,idusr=2
 /
 &dataparam
  in=288,jn=145,kn=23,yr=${year},mn=${month},dy=${day},hr=${hour},
  levtyp=${LEVTYP},fhour=0.,form='grads',qvar='ttd',
  yrev=.true.,big_endian=.true.,undef=9e+20
  recs=2,74,26,50,98,-1,146
 /
 &ipparam
  iptype=4,maxwv=-1
 /
 &level
  levels=1000,925,850,700,600,500,400,300,250,200,150,100, 70, 50, 30,
           20, 10,  7,  5,  3,  2,  1,0.4
 /
 &end
EOF

ulimit -s unlimited
$PGM
rm -f -r $DATA


