# Running parallel jobs on Pratyush HPC

#!/bin/bash
  
for file in  2m_temperature_2018_1.40625deg_12_31.nc
do
    tmpfile=$(mktemp --suffix=".sh")
    echo "#!/bin/sh" > $tmpfile
    echo "#PBS -N test_job" >> $tmpfile
    echo "#PBS -q research" >> $tmpfile
    echo "#PBS -l select=1:ncpus=1:vntype=cray_compute" >> $tmpfile
    echo "#PBS -l walltime=500:00:00" >> $tmpfile
    echo "cd /lus/dal/cccr_rnd/manmeet/AI_IITM/WeatherBench/data/dataserv.ub.tum.de/1.40625deg/2m_temperature" >> $tmpfile
  #  echo "cd $PBS_O_WORKDIR" >> $tmpfile
#    echo "module load cudatoolkit"
    echo "aprun  -n 1 -N 1 /lus/dal/cccr_rnd/manmeet/anaconda3/envs/py36/bin/python /lus/dal/cccr_rnd/manmeet/AI_IITM/WeatherBench/data/dataserv.ub.tum.de/1.40625deg/2m_temperature/preprocess_and_to_cs512.py $file  > /lus/dal/cccr_rnd/manmeet/AI_IITM/WeatherBench/data/dataserv.ub.tum.de/1.40625deg/2m_temperature/output.log" >> $tmpfile
    qsub  $tmpfile
done
