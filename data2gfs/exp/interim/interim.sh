#!/bin/bash
#----------------------------------------------------
# SCRIPT ID: 5
# YOTC lon/lat pressure level grib data
#----------------------------------------------------

export IDATE=2013060812
export JCAP=382

. ../configure

LEVTYP=pl

year=`echo $IDATE | cut -c 1-4`
month=`echo $IDATE | cut -c 5-6`
day=`echo $IDATE | cut -c 7-8`
hour=`echo $IDATE | cut -c 9-10`
echo $(pwd)
[ ! -d $DATA ]&&mkdir -p $DATA
PWD0=$(pwd)
cd $DATA

${PWD0}/download_interim_plev.py ${IDATE} 0.5 $DATA/input.grib

if [ ! -d $COMOUT/erai/${IDATE} ] ; then
    mkdir -p $COMOUT/erai/${IDATE}
fi

#ln -sf ERAINTERIM/0.5x0.5/anal_plev/${year}${month}/plev.${IDATE}.grib input.grib
ln -sf ${OROGRAPHY} input.orography.grib
ln -sf $COMOUT/erai/${IDATE}/siganl.${IDATE}.t${JCAP} output.siganl
ln -sf $PWD0/ecmwf.grib1.table input.grib1.table

cat <<EOF > data2gfs_namelist
 &gfsparam
  jcap=${JCAP},lonr=${LONB},latr=${LATB},levs=${LEVS},idusr=4
 /
 &dataparam
  in=720,jn=181,kn=37,yr=${year},mn=${month},dy=${day},hr=${hour},
  qvar='qq ',levtyp=${LEVTYP},fhour=0.,form='grib1'
 /
 &ipparam
  iptype=4,maxwv=-1
 /
 &level
  levels=1000,975,950,925,900,875,850,825,800,775,
         750, 700,650,600,550,500,450,400,350,300,
         250, 225,200,175,150,125,100,70, 50, 30,
         20, 10, 7,  5,  3,  2,  1
 /
 &end
EOF

ulimit -s unlimited
$PGM
rm -f -r $DATA


