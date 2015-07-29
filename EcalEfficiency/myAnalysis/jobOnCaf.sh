#!/bin/bash 

cd /afs/cern.ch/work/t/taroni/private/ECAL/EcalEleStudy/CMSSW_5_3_28/src/L1Studies/EGamma/Macros/myAnalysis/
eval `scram runtime -sh`
make runAll

exit
