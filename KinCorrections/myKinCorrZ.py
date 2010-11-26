#!/usr/bin/env python
import setupSUSY
from libFrameworkSUSY import *

from libWCharge       import *
from libSSDL      import *

from icf.core import PSet, Analysis
from icf.config import defaultConfig

import icf.utils as utils

import montecarlo.TTBarTauola as TT
import montecarlo.WJets_Madgraph as W
import montecarlo.Zjets_madgraph as Z


conf = defaultConfig.copy()
conf.Ntuple.Jets.Prefix="ak5JetPF"
conf.Common.Jets.EtCut = 20.0
conf.Common.Jets.EtaCut = 5.0
conf.Common.Muons.EtaCut = 2.4
conf.Common.Muons.PtCut = 15.0
conf.Common.Muons.TrkIsoCut = -1.
conf.Common.Muons.CombIsoCut = 0.1
conf.Common.Electrons.PtCut = 15.0
conf.Common.Electrons.EtaCut = 2.4
conf.Common.Photons.EtCut=30.
conf.Common.Electrons.ApplyID=False

# Create the analysis
a = Analysis("KinCorrection")

tree = Tree("Main")
Trigger=PSet(
       McAnal =False,
       MSugraScan = 0.,
       TriggerBits=PSet(
         bit1=""
        )
    )
ssdlTrigg=SSDLTrigger("Trigger",Trigger.ps())



ZeroMuons = OP_NumComMuons("==", 0)

KinCorrParPlus = PSet(
	zMassMin = 71.,
	zMassMax = 111.,
        BarrelEtaMax = 1.4442,
	EndCapEtaMin = 1.56,
	EndCapEtaMax = 2.5,
        MinElecRescaledPt=0.,
	c_ErNu_pt = 0.,
        wID=+24., ##24 W+, -24 W-
        mW = 80.398,
        mZ = 90.188,       
        ElecET     = 20.0,
        ChCheck=True,
        ConvCheck=True,
        EleVeto=True,
        ElePtVeto=15.0,
        CorrEEMisalig = True,
        WorkingPoint = 80,
        DataSet = False, #True = W, False =Z
        DataType = False, #True = LHC Data, False=MonteCarlo
	ReWeightHistos	= False,
	ReWeightCorrectionFile ="KinematicReWeight.root"
)
myKinCorrPlus=KinCorrection("KinCorrectionPlus_WP80",KinCorrParPlus.ps())

KinCorrParMinus = PSet(
	zMassMin = 71.,
	zMassMax = 111.,
        BarrelEtaMax = 1.4442,
	EndCapEtaMin = 1.56,
	EndCapEtaMax = 2.5,
        MinElecRescaledPt=0.,
	c_ErNu_pt = 0.,
        wID=-24.,
        mW = 80.398,
        mZ = 90.188,       
        ElecET     = 20.0,
        ChCheck=True,
        ConvCheck=True,
        EleVeto=True,
        ElePtVeto=15.0,
        CorrEEMisalig = True,
        WorkingPoint = 80,
        DataSet = False, #True = W, False =Z
        DataType = False, #True = LHC Data, False=MonteCarlo
	ReWeightHistos	= False,
	ReWeightCorrectionFile ="KinematicReWeight.root"
)

myKinCorrMinus=KinCorrection("KinCorrectionMinus_WP80",KinCorrParMinus.ps())
###------------------
KinCorrParPlus = PSet(
	zMassMin = 71.,
	zMassMax = 111.,
        BarrelEtaMax = 1.4442,
	EndCapEtaMin = 1.56,
	EndCapEtaMax = 2.5,
        MinElecRescaledPt=0.,
	c_ErNu_pt = 0.,
        wID=+24., ##24 W+, -24 W-
        mW = 80.398,
        mZ = 90.188,       
        ElecET     = 20.0,
        ChCheck=True,
        ConvCheck=True,
        EleVeto=True,
        ElePtVeto=15.0,
        CorrEEMisalig = True,
        WorkingPoint = 80,
        DataSet = False, #True = W, False =Z
        DataType = False, #True = LHC Data, False=MonteCarlo
	ReWeightHistos	= True,
	ReWeightCorrectionFile ="KinematicReWeight.root"
)
myKinCorrPlusReweight=KinCorrection("KinCorrectionPlus_Reweight_WP80",KinCorrParPlus.ps())

KinCorrParMinus = PSet(
	zMassMin = 71.,
	zMassMax = 111.,
        BarrelEtaMax = 1.4442,
	EndCapEtaMin = 1.56,
	EndCapEtaMax = 2.5,
        MinElecRescaledPt=0.,
	c_ErNu_pt = 0.,
        wID=-24.,
        mW = 80.398,
        mZ = 90.188,       
        ElecET     = 20.0,
        ChCheck=True,
        ConvCheck=True,
        EleVeto=True,
        ElePtVeto=15.0,
        CorrEEMisalig = True,
        WorkingPoint = 80,
        DataSet = False, #True = W, False =Z
        DataType = False, #True = LHC Data, False=MonteCarlo
	ReWeightHistos	= True,
	ReWeightCorrectionFile ="KinematicReWeight.root"
)
myKinCorrMinusReweight=KinCorrection("KinCorrectionMinus_Reweight_WP80",KinCorrParMinus.ps())

# Add the tree to our analysis
a += tree
##tree.Attach(ssdlTrigg)
tree.Attach(ZeroMuons)
tree.TAttach(ZeroMuons, myKinCorrPlus)
tree.TAttach(ZeroMuons, myKinCorrMinus)
#tree.TAttach(ZeroMuons, myKinCorrPlusReweight)
#tree.TAttach(ZeroMuons, myKinCorrMinusReweight)

## Z.zjets_madgraph_vols.LastEntry=1084920   ###METTERE FALSE NEL DATASET TYPE
#Z.zjets_madgraph_vols.LastEntry=10000  ###METTERE FALSE NEL DATASET TYPE
samples=[
           Z.zjets_madgraph_vols,
           ]




a.Run(".", conf, samples)
