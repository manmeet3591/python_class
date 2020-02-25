#!/bin/sh
###############################################################
#
#   AUTHOR:    Gilbert - W/NP11
#
#   DATE:      01/11/1999
#
#   PURPOSE:   This script uses the make utility to update the bacio 
#              archive libraries.
#
#
###############################################################
#
#
#
#     Remove make file, if it exists.  May need a new make file
#
if [ -f make.bacio ] ;  then
  rm -f make.bacio
fi
#
#     Generate a make file ( make.bacio) from this HERE file.
#
cat > make.bacio << EOF
SHELL=/bin/sh

\$(LIB):	\$(LIB)( bacio.o baciof.o bafrio.o byteswap.o chk_endianc.o)

\$(LIB)(bacio.o):       bacio.c clib.h
	\${CCMP} -c \$(CFLAGS) bacio.c
	ar -rv \$(AFLAGS) \$(LIB) bacio.o

\$(LIB)(baciof.o):   baciof.f
	\${FCMP} -c \$(FFLAGS) baciof.f
	ar -rv \$(AFLAGS) \$(LIB) baciof.o 

\$(LIB)(bafrio.o):   bafrio.f
	\${FCMP} -c \$(FFLAGS) bafrio.f
	ar -rv \$(AFLAGS) \$(LIB) bafrio.o 

\$(LIB)(byteswap.o):       byteswap.c 
	\${CCMP} -c \$(CFLAGS) byteswap.c
	ar -rv \$(AFLAGS) \$(LIB) byteswap.o

\$(LIB)(chk_endianc.o):       chk_endianc.f 
	\${FCMP} -c \$(FFLAGS) chk_endianc.f
	ar -rv \$(AFLAGS) \$(LIB) chk_endianc.o
	rm -f baciof.o bafrio.o bacio.o *.mod byteswap.o chk_endianc.o

EOF

#
export FCMP=${1:-gfortran}
export CCMP=${2:-gcc}
#
#     Update 4-byte version of libbacio_4.a
#
export LIB="../libbacio_4.a"

export FFLAGS=" -O2 -fPIC"
export AFLAGS=" "
# -I/usr/include/malloc needed for macos x
export CFLAGS=" -O2 -DUNDERSCORE -DLINUX -fPIC -I/usr/include/malloc"
make -f make.bacio
rm -f make.bacio
