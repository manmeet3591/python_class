#!/bin/tcsh
#PBS -A NTDD0004
#PBS -N ensSumPop
#PBS -q regular
#PBS -l select=4:ncpus=3:mpiprocs=3
#PBS -l walltime=0:45:00
#PBS -j oe
#PBS -M abaker@ucar.edu


mpiexec_mpt python pyEnsSumPop.py --verbose --indir  /glade/p/cisl/iowa/pop_verification/cesm2_0_beta10/ensembles --sumfile pop.ens.sum.cesm2.0.nc --tslice 0 --nyear 1 --nmonth 12 --esize 40 --jsonfile pop_ensemble.json  --mpi_enable --mach cheyenne --compset G --tag cesm2.0 --res T62_g17
