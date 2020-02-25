#!/bin/bash
#----------------------------------------------------
# SCRIPT ID: 1
# YOTC lon/lat hybrid level grib data
#----------------------------------------------------

export IDATE=2009093012
export JCAP=190

. ../configure

LEVTYP=hl

year=`echo $IDATE | cut -c 1-4`
month=`echo $IDATE | cut -c 5-6`
day=`echo $IDATE | cut -c 7-8`
hour=`echo $IDATE | cut -c 9-10`

[ ! -d $DATA ]&&mkdir -p $DATA
PWD0=$(pwd)
cd $DATA

if [ ! -d $COMOUT/ecmwf/${IDATE} ] ; then
    mkdir -p $COMOUT/ecmwf/${IDATE}
fi

ln -sf /mnt/drobo1/Public/YOTC/0.5x0.5/anal_hyblev/${year}${month}/hyblev.${IDATE}.grib input.grib
ln -sf /mnt/drobo1/Public/YOTC/0.5x0.5/fix/geopotential.grib input.geo.grib
ln -sf ${OROGRAPHY} input.orography.grib
ln -sf $COMOUT/ecmwf/${IDATE}/siganl.${IDATE}.t${JCAP} output.siganl
ln -sf $PWD0/ecmwf_hyblev.l91.txt input.vcoord.txt
ln -sf $PWD0/ecmwf.grib1.table input.grib1.table

cat <<EOF > data2gfs_namelist
 &gfsparam
  jcap=${JCAP},lonr=${LONB},latr=${LATB},levs=${LEVS},idusr=1
 /
 &dataparam
  in=720,jn=361,kn=91,yr=${year},mn=${month},dy=${day},hr=${hour},
  qvar='qq ',levtyp=${LEVTYP},fhour=0.,form='grib1'
 /
 &ipparam
  iptype=4,maxwv=-1
 /
 &level
  levels=91,90,89,88,87,86,85,84,83,82,81,80,79,78,77,76,75,74,73,72,
         71,70,69,68,67,66,65,64,63,62,61,60,59,58,57,56,55,54,53,52,
         51,50,49,48,47,46,45,44,43,42,41,40,39,38,37,36,35,34,33,32,
         31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,
         11,10, 9, 8, 7, 6, 5, 4, 3, 2, 1
 /
 &end
EOF

ulimit -s unlimited
$PGM
rm -f -r $DATA


