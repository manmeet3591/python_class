#!/bin/bash
#----------------------------------------------
# SCRIPT ID:6
# for ALERA2 sigma面バイナリ
#
# 注意:オゾンはo3clim=.true.にすると気候値のプロファイルを東西一様に与える。
#----------------------------------------------

export IDATE=2009093012
export JCAP=190
export COMOUT=$(pwd)
. ../configure

LEVTYP=sl
MEMBER=${MEMBER:-mean}

year=`echo $IDATE | cut -c 1-4`
month=`echo $IDATE | cut -c 5-6`
day=`echo $IDATE | cut -c 7-8`
hour=`echo $IDATE | cut -c 9-10`

if [ ! -d $DATA ] ; then 
    mkdir -p $DATA
fi
PWD0=$(pwd)
cd $DATA

if [ ! -d ${COMOUT}/alera2/${IDATE} ]; then
    mkdir -p ${COMOUT}/alera2/${IDATE}
fi

ln -sf /mnt/drobo1/Public/alera2/stream2008/slev/anal/${MEMBER}/${IDATE}.grd input.bin
ln -sf /mnt/drobo1/Public/alera2/geo/Grads.grd input.geo.bin
ln -sf ${OROGRAPHY} input.orography.grib
ln -sf ${O3CLIM} input.o3clim.txt
ln -sf $PWD0/alera2_sigmalev.l48.txt input.vcoord.txt
if [ $MEMBER = 'mean' ]; then
    MN='c00'
else
    MN=p${MEMBER:1:2}
fi
ln -sf ${COMOUT}/alera2/${IDATE}/siganl${MN}.${IDATE}.t${JCAP} output.siganl

cat <<EOF > data2gfs_namelist
 &gfsparam
  jcap=${JCAP},lonr=${LONB},latr=${LATB},levs=${LEVS},o3clim=.true.,idusr=6
 /
 &dataparam
  in=360,jn=180,kn=48,yr=${year},mn=${month},dy=${day},hr=${hour},
  levtyp=${LEVTYP},fhour=0.,form='grads',qvar='qq ',
  yrev=.false.,big_endian=.true.,undef=9e+20,zvar='hs',grid='gg',
  recs=-1,97,1,49,145,-1,193,psrec=241
 /
 &ipparam
  iptype=4,maxwv=-1
 /
 &level
 /
 &end
EOF

ulimit -s unlimited
$PGM

rm -f -r $DATA

cd $PWD0

