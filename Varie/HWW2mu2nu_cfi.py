import FWCore.ParameterSet.Config as cms

source = cms.Source("EmptySource")

generator = cms.EDFilter("Pythia6GeneratorFilter",
    comEnergy = cms.double(14000.),
    PythiaParameters = cms.PSet(
        #
        # Default cards for minimum bias events (unfiltered)
        # Name of the set is "pythiaMinBias"
        #include "IOMC/GeneratorInterface/test/pythiaMinBias.cfg"
        #
        # User cards - name is "myParameters"
        # Pythia's random generator initialization 
        myParameters = cms.vstring('MSEL=0 ! Users defined processes', 
            'MSUB(102)=1 ! Define the process : gg -> H', 
            'PMAS(23,1)=91.188 ! Z mass', 
            'PMAS(24,1)=80.450 ! W mass', 
            'PMAS(25,1)=140.00 ! H mass', 
            'MDME(210,1)=0 ! Switch off Higgs decay channels', 
            'MDME(211,1)=0', 
            'MDME(212,1)=0', 
            'MDME(213,1)=0', 
            'MDME(214,1)=0', 
            'MDME(215,1)=0', 
            'MDME(216,1)=0', 
            'MDME(217,1)=0', 
            'MDME(218,1)=0', 
            'MDME(219,1)=0', 
            'MDME(220,1)=0', 
            'MDME(221,1)=0', 
            'MDME(222,1)=0', 
            'MDME(223,1)=0', 
            'MDME(224,1)=0', 
            'MDME(225,1)=0 ! H -> ZZ switched off', 
            'MDME(226,1)=1 ! H -> WW switched on', 
            'MDME(190,1)=0           ',
            'MDME(191,1)=0           ',
            'MDME(192,1)=0           ',
            'MDME(193,1)=0           ',
            'MDME(194,1)=0           ',
            'MDME(195,1)=0           ',
            'MDME(196,1)=0           ',
            'MDME(197,1)=0           ',
            'MDME(198,1)=0           ',
            'MDME(199,1)=0           ',
            'MDME(200,1)=0           ',
            'MDME(201,1)=0           ',
            'MDME(202,1)=0           ',
            'MDME(203,1)=0           ',
            'MDME(204,1)=0           ',
            'MDME(205,1)=0           ',
            'MDME(206,1)=0           ',
            'MDME(207,1)=1           !W decay to mu nu',
            'MDME(208,1)=0           ',
            'MDME(209,1)=0           ',                                       
            'MSTJ(22)=2   ! Do not decay unstable particles', 
            'PARJ(71)=10. ! with c*tau > cTauMin (in mm) in PYTHIA'),
        # This is a vector of ParameterSet names to be read, in this order
        # The first two are in the include files below
        # The last one are simply my additional parameters
        parameterSets = cms.vstring('myParameters')
    )
)

ProductionFilterSequence = cms.Sequence(generator)
