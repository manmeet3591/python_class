#!/bin/bash
#combine_lev.ncl  jjas_2007.sh  jjas_2009.sh  jjas_2011.sh  jjas.sh     swath2grid.ncl
for year in {2006..2011}
do
		for dd in {153..275}
			do
				cp swath2grid.ncl combine_lev.ncl ${year}/${dd}/.
			done
done
