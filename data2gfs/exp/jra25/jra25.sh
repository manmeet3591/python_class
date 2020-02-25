#!/bin/bash
#----------------------------------------------
# SCRIPT ID:3
# for JRA-25 pressure level grib data
#
# 注意:JRA25はオゾン混合比、70hPa以上の比湿、雲水が含まれていない。
#      これらは読み込み時にゼロに設定されてから鉛直内挿される。
#
#      オゾンはo3clim=.true.にすると気候値のプロファイルを東西一様に与える。
#----------------------------------------------
export IDATE=2009093012
export JCAP=190
export COMOUT=$(pwd)
. ../configure

LEVTYP=pl

year=`echo $IDATE | cut -c 1-4`
month=`echo $IDATE | cut -c 5-6`
day=`echo $IDATE | cut -c 7-8`
hour=`echo $IDATE | cut -c 9-10`

[ ! -d $DATA ]&&mkdir -p $DATA
[ ! -d $COMOUT/jra25/${IDATE} ] && mkdir -p $COMOUT/jra25/${IDATE}

PWD0=$(pwd)
cd $DATA

ln -sf /mnt/clim6/jra25_grib/1.25x1.25/anal_p/${year}${month}/anl_p.${IDATE} input.grib
ln -sf ${OROGRAPHY} input.orography.grib
ln -sf ${O3CLIM} input.o3clim.txt
ln -sf $COMOUT/jra25/${IDATE}/siganl.${IDATE}.t${JCAP} output.siganl
ln -sf $PWD0/jma.grib1.table input.grib1.table

cat <<EOF > data2gfs_namelist
 &gfsparam
  jcap=${JCAP},lonr=${LONB},latr=${LATB},levs=${LEVS},idusr=3,o3clim=.true.
 /
 &dataparam
  in=288,jn=145,kn=23,yr=${year},mn=${month},dy=${day},hr=${hour},
  qvar='qq ',levtyp=${LEVTYP},fhour=0.,form='grib1'
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

cd $pwd

