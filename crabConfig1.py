[CRAB]

jobtype = cmssw
scheduler = lsf
#scheduler = condor
use_server = 0

[CMSSW]

### The data you want to access (to be found on DBS)
datasetpath=none

pset=none

number_of_jobs=2
events_per_job=1

### used to define random numbers
first_lumi=1

### The output files (comma separated list)
output_file=fitresults.root, plotsFitResults.root

[USER]
script_exe=scripts/RunGRID.sh

additional_input_files=topmassforgrid.tar.gz

ui_working_dir = crab_20150121_220test

return_data = 1
copy_data = 0

# remove if working at lpc
[LSF]
queue = 1nd
resource = "type=SLC6_64"
