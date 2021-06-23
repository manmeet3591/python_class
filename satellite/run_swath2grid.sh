#!/bin/bash
for year in {2006..2011}
do
		for dd in {153..275}
			do
				cd ${year}/${dd}/
					tmpfile=$(mktemp --suffix=".sh")
					echo "ncl swath2grid.ncl" > $tmpfile
					echo "ncl combine_lev.ncl" >> $tmpfile
					bsub -q "cccr-res" -W 10:00 -J "swath_${year}_${dd}" -o "swath.out" < $tmpfile
				cd -
			done
done
