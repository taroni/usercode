#!/bin/tcsh

cd /afs/cern.ch/work/t/taroni/private/ECAL/CMSSW_7_4_2_patch1/src/EcalPFG/Scripts/TriggerAnalysis/

eval `scramv1 runtime -csh`

./makeTrigPrimPlotsOld.sh -r 246926

exit
