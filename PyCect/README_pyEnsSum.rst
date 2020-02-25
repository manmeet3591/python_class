===============
README.pyEnsSum
===============

CESM-ECT (CESM Ensemble Consistency Test) is a suite of tests to 
determine whether a new simulation run(s) (e.g., from a new machine, 
compiler, etc.) is statistically distinguishable from an accepted 
ensemble.  The verification tools in the CESM-ECT suite are:

CAM-ECT - detects issues in CAM and CLM (12 month runs)
UF-CAM-ECT - detects issues in CAM and CLM (9 time step runs)
POP-ECT - detects issues in POP and CICE (12 month runs)

All of the above tools require an ensemble summary file (which contains
statistics describing the ensemble distribution). 

Ensemble summary files for existing CESM tags for CAM-ECT, UF-CAM-ECT, 
and POP-ECT that were created by CSEG are located (respectively) in the 
CESM input data directories:

$CESMDATAROOT/inputdata/validation/ensembles
$CESMDATAROOT/inputdata/validation/uf_ensembles
$CESMDATAROOT/inputdata/validation/pop_ensembles

Alternatively, *this package* can be used to create a summary file for CAM-ECT or
UF-CAM-ECT, given the location of appropriate ensemble history files (which should 
be generated in CIME via $CIME/tools/statistical_ensemble_test/ensemble.py).

(Note: to generate a summary file for POP-ECT, you must use pyEnsSumPop.py,
which has its own corresponding README file.)

:AUTHORS: Haiying Xu, Allison Baker
:COPYRIGHT: See the document entitled LICENSE.txt

Send questions and comments to Haiying Xu (haiyingx@ucar.edu) or Allison 
Baker (abaker@ucar.edu).

Github location:  https://github.com/NCAR/PyCECT/releases


This package includes:  
----------------------
     	pyEnsSum.py             
                            A script that generates an ensemble summary file 
     		            from a collection of CESM output files.

        pyEnsLib.py     
                            Library python script used by pyEnsSum.py

        pyEnsSum_test.sh        
                            Example PBS script to submit pyEnsSum.py to Cheyenne

        excluded_varlist.json
                            An example of a variable list that will be excluded from
                            reading and processing

        included_varlist.json
                            An example of a variable list that will be included for
                            reading and processing

	empty_excluded.json
	                   An empty exclude variable list, useful as a template
			   for listing variables to exclude


To use pyEnsSum: 
----------------------------------------------------------------------------
       Note: compatible with python 2.7 or greater (same as CESM) - NOT python 3

       1) For example, on NCAR's Cheyenne machine:

	  module load python/2.7.14
	  ncar_pylib
	  qsub test_pyEnsSum.sh


      2) Otherwise you need these packages:
       - module load python 
       - module load numpy
       - module load scipy
       - module load netCDF4
       - module load asaptools (https://github.com/NCAR/ASAPPyTools)
       - module load mpi4py
      
 
To see all options (and defaults):
----------------------------------
       python pyEnsSum.py -h

       PyCECT> python pyEnsSum.py -h

       Creates the summary file for an ensemble of CAM data. 

       ------------------------
       Args for pyEnsSum : 
       ------------------------
       pyEnsSum.py
       -h                   : prints out this usage message
       --verbose            : prints out in verbose mode (off by default)
       --sumfile <ofile>    : the output summary data file (default = ens.summary.nc)
       --indir <path>       : directory containing all of the ensemble runs (default = ./)
       --esize  <num>       : Number of ensemble members (default = 350)
       --tag <name>         : Tag name used in metadata (default = cesm2_0)
       --compset <name>     : Compset used in metadata (default = F2000climo)
       --res <name>         : Resolution used in metadata, (default = f19_f19)
       --tslice <num>       : the index into the time dimension, (default = 1)
       --mach <num>         : Machine name used in the metadata, (default = cheyenne)
       --jsonfile <fname>   : Jsonfile to provide that a list of variables that will 
                              be excluded or included  (default = exclude_empty.json)
       --mpi_disable        : Disable mpi mode to run in serial (off by default)
       --fIndex <num>       : Use this to start at ensemble member <num> instead of 000 (so 
                              ensembles with numbers less than <num> are excluded from summary file) 
   




Notes:
------

       1) CAM-ECT uses yearly average files, which by default (in the ensemble.py
	  generation script in CIME) also contains the initial conditions.  Therefore, 
	  one typically needs to set tslice=1 to use the yearly average (because 
	  tslice ==0 is the initial conditions.)

       2) UF-CAM-ECT uses timestep nine.  By default (in the ensemble.py
          generation script in CIME) the ouput file also contains the initial conditions.
	  Therefore, one typically needs to set tslice=1 to use time step nine (because
          tslice ==0 is the initial conditions.)

       3) There is no need to indicate UF-CAM-ECT vs. CAM-ECT to this routine.  It 
	  simply creates statistics for the supplied history files at the specified
	  time slice. For example, if you want to look at monthly files, simply 
	  supply their location.  Monthly files typically do not contain an initial 
	  condition and would requiree tslice=0.

       4) Esize (the ensemble size) can be less than or equal to the number of files 
	  in "--indir" (but ensembles 000-(esize-1) will be included unless fIndex
	  is specified.  UF-CAM-ECT typically uses at least 350 members (the default),
	  whereas CAM-ECT does not require as many.

       5) Note that --res, --tag, --compset, and --mach only affect the metadata 
	  in the summary file.

       6) When running in parallel, the recommended number of cores to use is one 
	  for each 3D variable. The default is to run in paralalel (recommended).

       7) You must specify a json file that indicates variables in the ensemble 
	  output files that you want to include or exclude from the summary file
	  statistics (see example json files).  We recommend excluding variables, as
	  this is typically less work and pyEnsSum will let you know if you have not
	  listed variables that need to be excluded (see next note).  Keep in mind that
	  you must have *fewer* variables included than ensemble members.

       8) IMPORTANT: If there are variables that need to be excluded (that are not in 
	  the .json file  already), pyEnsSum will exit early and provide a list of the
	  variables to exclude in the output.  These should be added to your exclude
	  variable list  (or removed from an include list), and then pyEnsSum can
	  be re-run.  Note that additional problematic variables may be found by 
	  pyEnsSum as variables are detected in three stages. (First any variables that 
	  are constant across the ensemble are identified.  Once these are removed, 
	  linearly dependant variables are indentified for removal. Finally, variables
	  that are not constant but have very few unique values are identified.)


Example for generating summary files:
--------------------------------------
	To generate a summary file for 350 UF-CAM-ECT simulations runs (time step nine), 
       	 
           we specify the size (this is optional since 350 is the default) and data location:
	    --esize 350
	    --indir /glade/p/cisl/iowa/verification/cesm2.1.0-rc.01/uf/

           We also specify the name of file to create for the summary:
 	    --sumfile ens.sum.cesm2.1.nc 

	   Since the ensemble files contain the intial conditions  as well
	   as the values at time step 9 (this is optional as 1 is the 
	   default), we set
	    --tslice 1 
	   
	   We also specify the CESM tag, compset and resolution of our ensemble data so 
	   that it can be written to the metadata of the summary file (ens.sum.cesm2.1.nc):
	    --tag cesm2.1.0 --compset F2000climo --res f19_f19 

           We can exclude or include some variables from the analysis by specifying them 
	   in a json file:
            --jsonfile excluded_varlist.json

	   This yields the following command for your job submission script:


	   python pyEnsSum.py  --esize 350 --indir /glade/p/cisl/iowa/verification/cesm2.1.0-rc.01/uf/ --tslice 1 --sumfile ens.sum.cesm2.1.0.nc --tag cesm2.1.0 --compset F2000climo --res f19_f19 --jsonfile excluded_varlist.json


	   
