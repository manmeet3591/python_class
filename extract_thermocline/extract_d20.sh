#!/bin/bash
for file in temp_ocean_????.nc 
do
echo $file
ferretscript="extract_d20.jnl"

rm -f $ferretscript
cat >> $ferretscript << EOF

use "$file"
sh da
let d20 = temp[z=@loc:20]
save/file=d20_${file} d20

EOF

#pyferret -nojnl -script $ferretscript


done
