#!/bin/bash
for file in temp_ocean_????.nc 

do
tmpfile=$(mktemp --suffix=".sh")
echo $file
ferretscript="extract_d20.jnl"

rm -f $ferretscript
cat >> $ferretscript << EOF
set memory/size=60000
use "$file"
sh da
let d20 = temp[z=@loc:20]
save/file=d20_${file} d20

EOF
echo "/iitm2/cccr-res/msingh/anaconda3/envs/FERRET/bin/pyferret -nojnl -script $ferretscript" >> $tmpfile
#pyferret -nojnl -script $ferretscript
bsub -q "cccr-res" -W 10:00 -J "$file" -o "d20.out" < $tmpfile

done
