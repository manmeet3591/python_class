#!/bin/bash
source activate py37_
for file in /home/cccr/msingh/dev_lab/ndvi3g_geo_v1_????_????.nc4 
	do
		echo $file
		python preprocess_ndvi.py $file 
	done

