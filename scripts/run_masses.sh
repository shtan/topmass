#! /bin/bash

cd /uscms/home/nmirman/CMSSW_5_3_5/src
eval `scramv1 runtime -sh`

export WORKING_DIR=$_CONDOR_SCRATCH_DIR

# mc top masses
mass=(161.5 163.5 166.5 169.5 172.5 175.5 178.5 181.5)

iter=$(($1/100))
#echo ${mass[$iter]} $lbound $rbound

cd /uscms/home/nmirman/topmass/
./DoFit --run_number $1 --bootstrap --fit --masspnt ${mass[$iter]}
#./DoFit --run_number $1 --bootstrap --fit --masspnt ${mass[$iter]} --gnorm 30000
#./DoFit --run_number $1 --bootstrap --fit --masspnt ${mass[$iter]} --lbnd 200 --rbnd 300
#./DoFit --fit --masspnt ${mass[$iter]} --lbnd 0 --rbnd 30
#./DoFit --fit --masspnt ${mass[$1]} --lmbl 13 --lmt 32
