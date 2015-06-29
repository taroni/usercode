// -*- C++ -*-
//
// Class:      EcalTPGAnalyzer
//
//
// Original Author:  Pascal Paganini
//
//


// system include files
#include <memory>
#include <utility>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"

#include "DataFormats/DetId/interface/DetId.h"

#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMask.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMaskAlgoTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMaskTechTrigRcd.h"

#include "CondFormats/DataRecord/interface/EcalTPGTowerStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalTPGTowerStatus.h"




//L1 Trigger
#include "L1Trigger/L1ExtraFromDigis/interface/L1ExtraParticleMapProd.h"

//RCT
#include "CondFormats/L1TObjects/interface/L1RCTChannelMask.h"
#include "CondFormats/DataRecord/interface/L1RCTChannelMaskRcd.h"

#include "EcalTPGAnalyzer.h"


//#include "DQM/EcalCommon/interface/Numbers.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"


// Trip Masking info
#include "CondFormats/EcalObjects/interface/EcalTPGStripStatus.h"
#include "CondFormats/DataRecord/interface/EcalTPGStripStatusRcd.h"


#include <TMath.h>
#include <sstream>

using namespace edm;
class CaloSubdetectorGeometry;

EcalTPGAnalyzer::EcalTPGAnalyzer(const edm::ParameterSet&  iConfig)
{
    noL1Info = false;

    tpCollection_ = iConfig.getParameter<edm::InputTag>("TPCollection") ;
    tpEmulatorCollection_ = iConfig.getParameter<edm::InputTag>("TPEmulatorCollection") ;
    digiCollectionEB_ = iConfig.getParameter<edm::InputTag>("DigiCollectionEB") ;
    digiCollectionEE_ = iConfig.getParameter<edm::InputTag>("DigiCollectionEE") ;
    gtRecordCollectionTag_ = iConfig.getParameter<std::string>("GTRecordCollection") ;
    EcalRecHitCollectionEB_ = iConfig.getParameter<edm::InputTag>("EcalRecHitCollectionEB") ;
    EcalRecHitCollectionEE_ = iConfig.getParameter<edm::InputTag>("EcalRecHitCollectionEE") ;
  
    allowTP_ = iConfig.getParameter<bool>("ReadTriggerPrimitives");
    useEE_ = iConfig.getParameter<bool>("UseEndCap");
    skipWritingEndcapDigis_ = iConfig.getParameter<bool>("skipWritingEndcapDigi");
    keepOnlyTowersAboveZero_ = iConfig.getParameter<bool>("keepOnlyTowersAboveZero") ;
    print_ = iConfig.getParameter<bool>("Print");
    L1print_ = iConfig.getParameter<bool>("L1Print");
   
    // file
    file_ = new TFile("ECALTPGtree.root","RECREATE");
    file_->cd() ;
  
    // trees
    tree_ = new TTree( "EcalTPGAnalysis","EcalTPGAnalysis" );
    tree_->Branch("runNb",&treeVariables_.runNb,"runNb/i"); //
    tree_->Branch("evtNb",&treeVariables_.evtNb,"evtNb/l"); //
    tree_->Branch("bxNb",&treeVariables_.bxNb,"bxNb/i"); //
    tree_->Branch("bxGT",&treeVariables_.bxGT,"bxGT/i"); //
    tree_->Branch("orbitNb",&treeVariables_.orbitNb,"orbitNb/l"); //
    tree_->Branch("lumiBlock",&treeVariables_.lumiBlock,"lumiBlock/i"); //
    tree_->Branch("timeStamp",&treeVariables_.timeStamp,"timeStamp/d"); //
    tree_->Branch("nbOfActiveTriggers",&treeVariables_.nbOfActiveTriggers,"nbOfActiveTriggers/i"); //
    tree_->Branch("activeTriggers",treeVariables_.activeTriggers,"activeTriggers[nbOfActiveTriggers]/I"); //
    tree_->Branch("nbOfActiveTechTriggers",&treeVariables_.nbOfActiveTechTriggers,"nbOfActiveTechTriggers/i"); //
    tree_->Branch("activeTechTriggers",treeVariables_.activeTechTriggers,"activeTechTriggers[nbOfActiveTechTriggers]/I");
    tree_->Branch("activeTriggersBX",treeVariables_.activeTriggersBX,"activeTriggersBX[128]/I"); //  

    tree_->Branch("nbOfTowers",&treeVariables_.nbOfTowers,"nbOfTowers/i"); //

    tree_->Branch("ieta", treeVariables_.ieta,"ieta[nbOfTowers]/I");//
    tree_->Branch("iphi", treeVariables_.iphi,"iphi[nbOfTowers]/I");//
    tree_->Branch("nbOfXtals", treeVariables_.nbOfXtals,"nbOfXtals[nbOfTowers]/I");//
    tree_->Branch("rawTPData", treeVariables_.rawTPData,"rawTPData[nbOfTowers]/I");//
    tree_->Branch("rawTPEmul1", treeVariables_.rawTPEmul1,"rawTPEmul1[nbOfTowers]/I");//
    tree_->Branch("rawTPEmul2", treeVariables_.rawTPEmul2,"rawTPEmul2[nbOfTowers]/I");//
    tree_->Branch("rawTPEmul3", treeVariables_.rawTPEmul3,"rawTPEmul3[nbOfTowers]/I");//
    tree_->Branch("rawTPEmul4", treeVariables_.rawTPEmul4,"rawTPEmul4[nbOfTowers]/I");//
    tree_->Branch("rawTPEmul5", treeVariables_.rawTPEmul5,"rawTPEmul5[nbOfTowers]/I");//
    tree_->Branch("crystNb", treeVariables_.crystNb,"crystNb[nbOfTowers]/I");//
    tree_->Branch("eRec", treeVariables_.eRec,"eRec[nbOfTowers]/F");//
    tree_->Branch("ttFlag", treeVariables_.ttFlag,"ttFlag[nbOfTowers]/I");//
    tree_->Branch("sevlv", treeVariables_.sevlv,"sevlv[nbOfTowers]/I");//
    tree_->Branch("spike", treeVariables_.spike,"spike[nbOfTowers]/I");//
    tree_->Branch("rawTPEmulsFGVB1", treeVariables_.rawTPEmulsFGVB1,"rawTPEmulsFGVB1[nbOfTowers]/I");//
    tree_->Branch("rawTPEmulsFGVB2", treeVariables_.rawTPEmulsFGVB2,"rawTPEmulsFGVB2[nbOfTowers]/I");//
    tree_->Branch("rawTPEmulsFGVB3", treeVariables_.rawTPEmulsFGVB3,"rawTPEmulsFGVB3[nbOfTowers]/I");//
    tree_->Branch("rawTPEmulsFGVB4", treeVariables_.rawTPEmulsFGVB4,"rawTPEmulsFGVB4[nbOfTowers]/I");//
    tree_->Branch("rawTPEmulsFGVB5", treeVariables_.rawTPEmulsFGVB5,"rawTPEmulsFGVB5[nbOfTowers]/I");//
    tree_->Branch("rawTPEmulttFlag1", treeVariables_.rawTPEmulttFlag1,"rawTPEmulttFlag1[nbOfTowers]/I");//
    tree_->Branch("rawTPEmulttFlag2", treeVariables_.rawTPEmulttFlag2,"rawTPEmulttFlag2[nbOfTowers]/I");//
    tree_->Branch("rawTPEmulttFlag3", treeVariables_.rawTPEmulttFlag3,"rawTPEmulttFlag3[nbOfTowers]/I");//
    tree_->Branch("rawTPEmulttFlag4", treeVariables_.rawTPEmulttFlag4,"rawTPEmulttFlag4[nbOfTowers]/I");//
    tree_->Branch("rawTPEmulttFlag5", treeVariables_.rawTPEmulttFlag5,"rawTPEmulttFlag5[nbOfTowers]/I");//
  
    treeShape_ = new TTree( "EcalShape","EcalShape" );
    treeShape_->Branch("evtNb",&treeVariablesShape_.evtNb,"evtNb/l"); //
    treeShape_->Branch("ieta", &treeVariablesShape_.ieta,"ieta/i");//
    treeShape_->Branch("iphi", &treeVariablesShape_.iphi,"iphi/i");//
    treeShape_->Branch("ixXtal", &treeVariablesShape_.ixXtal,"ixXtal/i");//
    treeShape_->Branch("iyXtal", &treeVariablesShape_.iyXtal,"iyXtal/i");//
    treeShape_->Branch("TCCid", &treeVariablesShape_.TCCid,"TCCid/i");//
    treeShape_->Branch("TowerInTCC", &treeVariablesShape_.TowerInTCC,"TowerInTCC/i");//
    treeShape_->Branch("strip", &treeVariablesShape_.strip,"strip/i");//
    treeShape_->Branch("nbOfSamples", &treeVariablesShape_.nbOfSamples,"nbOfSamples/i");//
    treeShape_->Branch("samp", treeVariablesShape_.samp,"samp[nbOfSamples]/I");//
  
  
    treeAux = new TTree("treeAux","treeAux");
    treeAux->Branch("runNb",&ecalAux_.runNb,"runNb/i"); //
    treeAux->Branch("nMaskedCh",&ecalAux_.nMaskedCh,"nMaskedCh/i"); //
    treeAux->Branch("iMaskedTTeta",ecalAux_.iMaskedTTeta,"iMaskedTTeta[nMaskedCh]/I"); //
    treeAux->Branch("iMaskedTTphi",ecalAux_.iMaskedTTphi,"iMaskedTTeta[nMaskedCh]/I"); //

    treeAux->Branch("nMaskedChannelsFromStrips",&ecalAux_.nMaskedChannelsFromStrips,"nMaskedChannelsFromStrips/i"); //
    treeAux->Branch("iMaskedChannelsFromStripsX",ecalAux_.iMaskedChannelsFromStripsX,"iMaskedChannelsFromStripsX[nMaskedChannelsFromStrips]/I"); //
    treeAux->Branch("iMaskedChannelsFromStripsY",ecalAux_.iMaskedChannelsFromStripsY,"iMaskedChannelsFromStripsY[nMaskedChannelsFromStrips]/I"); //
    treeAux->Branch("iMaskedChannelsFromStripsZ",ecalAux_.iMaskedChannelsFromStripsZ,"iMaskedChannelsFromStripsZ[nMaskedChannelsFromStrips]/I"); //
        
    //L1 tree branches
    tree_->Branch("nbOfL1IsoCands",&treeVariables_.nbOfL1IsoCands,"nbOfL1IsoCands/i"); //
    tree_->Branch("L1IsoIeta", treeVariables_.L1IsoIeta,"L1IsoIeta[nbOfL1IsoCands]/I");//
    tree_->Branch("L1IsoIphi", treeVariables_.L1IsoIphi,"L1IsoIphi[nbOfL1IsoCands]/I");//
    tree_->Branch("L1IsoRank", treeVariables_.L1IsoRank,"L1IsoRank[nbOfL1IsoCands]/I");//
  
    tree_->Branch("nbOfL1NonisoCands",&treeVariables_.nbOfL1NonisoCands,"nbOfL1NonisoCands/i"); //
    tree_->Branch("L1NonisoIeta", treeVariables_.L1NonisoIeta,"L1NonisoIeta[nbOfL1NonisoCands]/I");//
    tree_->Branch("L1NonisoIphi", treeVariables_.L1NonisoIphi,"L1NonisoIphi[nbOfL1NonisoCands]/I");//
    tree_->Branch("L1NonisoRank", treeVariables_.L1NonisoRank,"L1NonisoRank[nbOfL1NonisoCands]/I");//

    tree_->Branch("nMaskedRCT",&treeVariables_.nMaskedRCT,"nMaskedRCT/i"); //
    tree_->Branch("iMaskedRCTeta",treeVariables_.iMaskedRCTeta,"iMaskedRCTeta[nMaskedRCT]/I"); //
    tree_->Branch("iMaskedRCTcrate",treeVariables_.iMaskedRCTcrate,"iMaskedRCTcrate[nMaskedRCT]/I"); //
    tree_->Branch("iMaskedRCTphi",treeVariables_.iMaskedRCTphi,"iMaskedRCTphi[nMaskedRCT]/I"); //

    // L1 pre-firing tree branches

    tree_->Branch("nbOfL1preIsoCands",&treeVariables_.nbOfL1preIsoCands,"nbOfL1preIsoCands/i"); //
    tree_->Branch("L1preIsoIeta", treeVariables_.L1preIsoIeta,"L1preIsoIeta[nbOfL1preIsoCands]/I");//
    tree_->Branch("L1preIsoIphi", treeVariables_.L1preIsoIphi,"L1preIsoIphi[nbOfL1preIsoCands]/I");//
    tree_->Branch("L1preIsoRank", treeVariables_.L1preIsoRank,"L1preIsoRank[nbOfL1preIsoCands]/I");//
  
    tree_->Branch("nbOfL1preNonisoCands",&treeVariables_.nbOfL1preNonisoCands,"nbOfL1preNonisoCands/i"); //
    tree_->Branch("L1preNonisoIeta", treeVariables_.L1preNonisoIeta,"L1preNonisoIeta[nbOfL1preNonisoCands]/I");//
    tree_->Branch("L1preNonisoIphi", treeVariables_.L1preNonisoIphi,"L1preNonisoIphi[nbOfL1preNonisoCands]/I");//
    tree_->Branch("L1preNonisoRank", treeVariables_.L1preNonisoRank,"L1preNonisoRank[nbOfL1preNonisoCands]/I");//


 // ontime firing tree branches
 /*
     tree_->Branch("nbOfL1ontimeIsoCands",&treeVariables_.nbOfL1ontimeIsoCands,"nbOfL1ontimeIsoCands/i"); //
     tree_->Branch("L1ontimeIsoIeta", treeVariables_.L1ontimeIsoIeta,"L1ontimeIsoIeta[nbOfL1ontimeIsoCands]/I");//
     tree_->Branch("L1ontimeIsoIphi", treeVariables_.L1ontimeIsoIphi,"L1ontimeIsoIphi[nbOfL1ontimeIsoCands]/I");//
     tree_->Branch("L1ontimeIsoRank", treeVariables_.L1ontimeIsoRank,"L1ontimeIsoRank[nbOfL1ontimeIsoCands]/I");//
   
     tree_->Branch("nbOfL1ontimeNonisoCands",&treeVariables_.nbOfL1ontimeNonisoCands,"nbOfL1ontimeNonisoCands/i"); //
     tree_->Branch("L1ontimeNonisoIeta", treeVariables_.L1ontimeNonisoIeta,"L1ontimeNonisoIeta[nbOfL1ontimeNonisoCands]/I");//
     tree_->Branch("L1ontimeNonisoIphi", treeVariables_.L1ontimeNonisoIphi,"L1ontimeNonisoIphi[nbOfL1ontimeNonisoCands]/I");//
     tree_->Branch("L1ontimeNonisoRank", treeVariables_.L1ontimeNonisoRank,"L1ontimeNonisoRank[nbOfL1ontimeNonisoCands]/I");//
 */

    // L1 post-firing tree branches
    tree_->Branch("nbOfL1postIsoCands",&treeVariables_.nbOfL1postIsoCands,"nbOfL1postIsoCands/i"); //
    tree_->Branch("L1postIsoIeta", treeVariables_.L1postIsoIeta,"L1postIsoIeta[nbOfL1postIsoCands]/I");//
    tree_->Branch("L1postIsoIphi", treeVariables_.L1postIsoIphi,"L1postIsoIphi[nbOfL1postIsoCands]/I");//
    tree_->Branch("L1postIsoRank", treeVariables_.L1postIsoRank,"L1postIsoRank[nbOfL1postIsoCands]/I");//
  
    tree_->Branch("nbOfL1postNonisoCands",&treeVariables_.nbOfL1postNonisoCands,"nbOfL1postNonisoCands/i"); //
    tree_->Branch("L1postNonisoIeta", treeVariables_.L1postNonisoIeta,"L1postNonisoIeta[nbOfL1postNonisoCands]/I");//
    tree_->Branch("L1postNonisoIphi", treeVariables_.L1postNonisoIphi,"L1postNonisoIphi[nbOfL1postNonisoCands]/I");//
    tree_->Branch("L1postNonisoRank", treeVariables_.L1postNonisoRank,"L1postNonisoRank[nbOfL1postNonisoCands]/I");//
    //


    //   ADD the branches from here http://azabi.web.cern.ch/azabi/TRIGGER/CODE/SimpleNtpleCustom.cc
    tree_->Branch("trig_tower_adc",  &treeVariables_.twrADC,  "twrADC[nbOfTowers]/I");
    tree_->Branch("trig_tower_sFGVB",  &treeVariables_.sFGVB,  "sFGVB[nbOfTowers]/I");

  

}


EcalTPGAnalyzer::~EcalTPGAnalyzer()
{
    file_->cd();
    tree_->Write();
    treeShape_->Write();
    treeAux->Write();
    file_->Close();
}

void EcalTPGAnalyzer::beginJob()
{
    myevt = 0;
    std::cout << " beginJob myevt = " << myevt << std::endl;
}

void EcalTPGAnalyzer::beginRun()
{
    myevt = 0;
    std::cout << " beginRun myevt = " << myevt << std::endl;

}


void EcalTPGAnalyzer::analyze(const edm::Event& iEvent, const  edm::EventSetup & iSetup)
{


    using namespace edm;
    using namespace std;

    myevt++;
    // cout << " analyze myevt = " << myevt << endl;

    if ( myevt == 1)
    {

        // electronics mapping
        ESHandle< EcalElectronicsMapping > ecalmapping;
        iSetup.get< EcalMappingRcd >().get(ecalmapping);
        theMapping_ = ecalmapping.product();
  
        // Reading strip statuses
        try
        {
           edm::ESHandle<EcalTPGStripStatus> theEcalTPGStripStatus_handle;
           iSetup.get<EcalTPGStripStatusRcd>().get(theEcalTPGStripStatus_handle);
           const EcalTPGStripStatus * ecaltpgStripStatus = theEcalTPGStripStatus_handle.product();
           const EcalTPGStripStatusMap &stripMap = ecaltpgStripStatus->getMap();
           EcalTPGStripStatusMapIterator itSt;
           ecalAux_.nMaskedChannelsFromStrips = 0;
           for (itSt = stripMap.begin();itSt != stripMap.end();++itSt) {
              if(itSt->second > 0) {
                 // lets decode the ID
                 int strip = itSt->first/8;
                 int pseudostrip = strip & 0x7;
                 strip /= 8;
                 int tt = strip & 0x7F;
                 strip /= 128;
                 int tccid = strip & 0x7F;
                 std::cout << "Bad strip ID = " << itSt->first
                 << " TCC " << tccid << " TT " << tt << " ST " << pseudostrip
                 << ", status = " << itSt->second << std::endl;
                 std::vector<DetId> xtalId =
                 theMapping_->pseudoStripConstituents(tccid, tt, pseudostrip);
                 std::vector<DetId>::iterator it;
                 for(it = xtalId.begin(); it != xtalId.end(); it++) {
                    EEDetId eeid = *it;
                    std::cout << eeid << " ";
                    // int iz = eeid.zside();
                    // int ix = eeid.ix();
                    // int iy = eeid.iy();
                    ecalAux_.iMaskedChannelsFromStripsX[ecalAux_.nMaskedChannelsFromStrips]=eeid.ix();
                    ecalAux_.iMaskedChannelsFromStripsY[ecalAux_.nMaskedChannelsFromStrips]=eeid.iy();
                    ecalAux_.iMaskedChannelsFromStripsZ[ecalAux_.nMaskedChannelsFromStrips]=eeid.zside();
                    ecalAux_.nMaskedChannelsFromStrips++;
                 }
                 std::cout << std::endl;
              }
           }
        }
        catch (cms::Exception& exception)
        {
           std::cout <<  "   ===> EcalTPGStripStatusRcd not found => skip strip masking info" << std::endl;
        }
        
        // Reading TPG tower statuses for first event (masked/not masked information)        
        edm::ESHandle<EcalTPGTowerStatus> theEcalTPGTowerStatus_handle;
        iSetup.get<EcalTPGTowerStatusRcd>().get(theEcalTPGTowerStatus_handle);
        const EcalTPGTowerStatus * ecaltpgTowerStatus=theEcalTPGTowerStatus_handle.product();

        const EcalTPGTowerStatusMap &towerMap=ecaltpgTowerStatus->getMap();
        EcalTPGTowerStatusMapIterator  it;
        

        uint nMaskedChannels = 0;
        for (it=towerMap.begin();it!=towerMap.end();++it) {
            if ((*it).second > 0) // if status not equals 0 then channle was massked 
            {
                EcalTrigTowerDetId  ttId((*it).first);
                ecalAux_.iMaskedTTeta[nMaskedChannels] = ttId.ieta();
                ecalAux_.iMaskedTTphi[nMaskedChannels] = ttId.iphi();
                nMaskedChannels++;
            }
        }

        ecalAux_.nMaskedCh = nMaskedChannels;

        // Storing eta and phi of masked TT to 'treeAux' root tree
        treeAux->Fill() ;



        // geometry
        ESHandle<CaloGeometry> theGeometry;
        ESHandle<CaloSubdetectorGeometry> theEndcapGeometry_handle, theBarrelGeometry_handle;
  
        iSetup.get<CaloGeometryRecord>().get( theGeometry );
        iSetup.get<EcalEndcapGeometryRecord>().get("EcalEndcap",theEndcapGeometry_handle);
        iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel",theBarrelGeometry_handle);
  
        iSetup.get<IdealGeometryRecord>().get(eTTmap_);
        theEndcapGeometry_ = &(*theEndcapGeometry_handle);
        theBarrelGeometry_ = &(*theBarrelGeometry_handle);
  
//        Numbers::initGeometry(iSetup, false);
  


        //Adding RCT mask
        // list of RCT channels to mask
        vector<int> RCTMaskedEta;
        vector<int> RCTMaskedPhi;
        edm::ESHandle<L1RCTChannelMask> channelMask;
        iSetup.get<L1RCTChannelMaskRcd>().get(channelMask);
        const L1RCTChannelMask* cEs = channelMask.product();
        uint nMaskedRCT = 0;
        for(int i = 0; i< 18; i++)
            for(int j =0; j< 2; j++)
                for(int k =0; k<28; k++)
                    if(cEs->ecalMask[i][j][k]){
                        //cout << "ECAL masked channel: RCT crate " << i << " iphi " << j <<" ieta " <<k <<endl; 
                        treeVariables_.iMaskedRCTeta[nMaskedRCT]=k;
                        treeVariables_.iMaskedRCTphi[nMaskedRCT]=j;
                        treeVariables_.iMaskedRCTcrate[nMaskedRCT]=i;
                        nMaskedRCT++;
                    }
        treeVariables_.nMaskedRCT = nMaskedRCT;

    }


    map<EcalTrigTowerDetId, towerEner> mapTower ;
    map<EcalTrigTowerDetId, towerEner>::iterator itTT ;
    map<EcalTrigTowerDetId, EBDataFrame> mapShapeEB ;
    map<EcalTrigTowerDetId, EBDataFrame>::iterator itEB ;
    map<EcalTrigTowerDetId, EEDataFrame> mapShapeEE ;
    map<EcalTrigTowerDetId, EEDataFrame> ::iterator itEE ; 


    ///////////////////////////
    // get Evts info
    ///////////////////////////

    treeVariables_.runNb = iEvent.id().run() ;
    treeVariables_.evtNb = iEvent.id().event() ;
    treeVariables_.bxNb = iEvent.bunchCrossing() ;
    treeVariables_.orbitNb = iEvent.orbitNumber() ;
    treeVariables_.lumiBlock = iEvent.luminosityBlock();

    treeVariablesShape_.evtNb = iEvent.id().event() ;


  
    //===================TIMESTAMP INFORMTION=================================
    // Code added to get the time in seconds
    double startTime_= 0.0;
    startTime_ += 1215107133.;
    startTime_ += 33470900.; //start time of run (empirically)
    unsigned int  timeStampLow = ( 0xFFFFFFFF & iEvent.time().value() );
    unsigned int  timeStampHigh = ( iEvent.time().value() >> 32 );
    double rawEventTime = ( double)(timeStampHigh)+((double )(timeStampLow)/1000000.);
    double eventTime = rawEventTime - startTime_; //Notice the subtraction of the "starttime"
    treeVariables_.timeStamp = eventTime;

 
    // print_ = true;
    if (print_)
      std::cout<<"==========="<<iEvent.id()<<"   bxNb:"<<iEvent.bunchCrossing()<<" orbit:"<<iEvent.orbitNumber()<<std::endl ;


    if (!noL1Info) {
    ///////////////////////////
    // get L1 candidate
    ///////////////////////////
      edm::Handle< l1extra::L1EmParticleCollection > emNonisolColl ;
      iEvent.getByLabel("l1extraParticles","NonIsolated", emNonisolColl ) ;
      edm::Handle< l1extra::L1EmParticleCollection > emIsolColl ;
      iEvent.getByLabel("l1extraParticles","Isolated", emIsolColl ) ;
      

      if (L1print_)  std::cout << "iso size " << emIsolColl->size() << std::endl;
      treeVariables_.nbOfL1IsoCands = emIsolColl->size();
      int isocounter = 0;
      for( l1extra::L1EmParticleCollection::const_iterator emItr = emIsolColl->begin(); emItr != emIsolColl->end() ;++emItr){
	
        treeVariables_.L1IsoIphi[isocounter] = emItr->gctEmCand()->regionId().iphi();
        treeVariables_.L1IsoIeta[isocounter] = emItr->gctEmCand()->regionId().ieta();
        treeVariables_.L1IsoRank[isocounter] = emItr->gctEmCand()->rank();
        isocounter++;
        if (L1print_)  cout<< " L1 Isolated EM object number, type, id, e, et, eta, phi"<<
			 emIsolColl->size()<<" "<<emItr->type()<<" "<<emItr->pdgId()<<" "<<emItr->energy()<<" "<<emItr->et()<<" "<<emItr->eta()<<" "<<emItr->phi()<<endl;
	
	
      }

      if (L1print_) std::cout << "noniso size " << emNonisolColl->size() << std::endl;
      treeVariables_.nbOfL1NonisoCands = emNonisolColl->size();
      int nonisocounter = 0;
      for( l1extra::L1EmParticleCollection::const_iterator emItr = emNonisolColl->begin(); emItr != emNonisolColl->end() ;++emItr){
	
        treeVariables_.L1NonisoIphi[nonisocounter] = emItr->gctEmCand()->regionId().iphi();
        treeVariables_.L1NonisoIeta[nonisocounter] = emItr->gctEmCand()->regionId().ieta();
        treeVariables_.L1NonisoRank[nonisocounter] = emItr->gctEmCand()->rank();
        nonisocounter++;
        if (L1print_) cout<< " L1 Nonisolated EM object number, type, id, e, et, eta, phi"<<
			emIsolColl->size()<<" "<<emItr->type()<<" "<<emItr->pdgId()<<" "<<emItr->energy()<<" "<<emItr->et()<<" "<<emItr->eta()<<" "<<emItr->phi()<<endl;
      }
      
      
    }
    ///////////////////////////
    // get L1 Trigger info
    ///////////////////////////
  
    edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
    iEvent.getByLabel( edm::InputTag(gtRecordCollectionTag_), gtRecord);
    DecisionWord dWord = gtRecord->decisionWord();   // this will get the decision word *before* masking disabled bits

    edm::ESHandle< L1GtTriggerMask > l1GtTmAlgo;
    iSetup.get< L1GtTriggerMaskAlgoTrigRcd >().get( l1GtTmAlgo );        
    std::vector<unsigned int> triggerMaskAlgoTrig = l1GtTmAlgo.product()->gtTriggerMask();

    L1GtfeWord gtfeWord = gtRecord->gtfeWord();
    treeVariables_.bxGT = gtfeWord.bxNr();


    int hashedTr[128];
    for ( int i = 0; i < 128; ++i )
    {
        hashedTr[i]=0;
        treeVariables_.activeTriggersBX[i]= 0; 
        treeVariables_.activeTriggers[i]= 0;
    }
    for(int iebx=0; iebx<=4; iebx++) {
        DecisionWord gtDecisionWord = gtRecord->decisionWord(iebx-2);
      
        for ( int i = 0; i < 128; ++i ) {
            if ( gtDecisionWord[i] ) {
                int coef = 0;
                if (iebx == 0 ) coef = 1;
                if (iebx == 1 ) coef = 10;
                if (iebx == 2 ) coef = 100;
                if (iebx == 3 ) coef = 1000;
                if (iebx == 4 ) coef = 10000;;
                hashedTr[i]+= coef;
                treeVariables_.activeTriggersBX[i]= hashedTr[i];
            }
        }
    }



    // apply masks on algo    
    int iDaq = 0;
    int iBit = -1;
    treeVariables_.nbOfActiveTriggers = 0 ;
    for (std::vector<bool>::iterator itBit = dWord.begin(); itBit != dWord.end(); ++itBit) {        
        iBit++;
        int maskBit = triggerMaskAlgoTrig[iBit] & (1 << iDaq);
        if (maskBit) *itBit = false;
        if (*itBit) {
            treeVariables_.activeTriggers[treeVariables_.nbOfActiveTriggers] = iBit ;
            treeVariables_.nbOfActiveTriggers++ ;      
        }
    }
    
    treeVariables_.nbOfActiveTechTriggers = 0;
    TechnicalTriggerWord tw = gtRecord->technicalTriggerWord();
    if ( ! tw.empty() ) {
        // loop over dec. bit to get total rate (no overlap)
        for ( int itechbit = 0; itechbit < 64; ++itechbit ) {
   
            treeVariables_.activeTechTriggers[treeVariables_.nbOfActiveTechTriggers] = 0; // ADD THIS 
   
            if ( tw[itechbit] ){
                treeVariables_.activeTechTriggers[treeVariables_.nbOfActiveTechTriggers] = itechbit;
                treeVariables_.nbOfActiveTechTriggers++ ;
            }
   
        }
    }


    //PRE-FIRING INFO
    const L1GtPsbWord psb = gtRecord->gtPsbWord(0xbb0d, -1);
    //psb.print(cout); 
    vector<int> psbel;
    psbel.push_back(psb.aData(4));
    psbel.push_back(psb.aData(5));
    psbel.push_back(psb.bData(4));
    psbel.push_back(psb.bData(5));
    int isocounter1 = 0;
    treeVariables_.nbOfL1preNonisoCands = psbel.size();
    std::vector<int>::const_iterator ipsbel;
    for(ipsbel=psbel.begin(); ipsbel!=psbel.end(); ipsbel++) {
        //int rank = (*ipsbel)&0x3f;
        //double pt = emScale->et(rank); 
        //if(rank>0 && pt>0.6) {// firing EG1
        int iEta = int(((*ipsbel)>>6)&7);
        int sign = ( ((*ipsbel>>9)&1) ? -1. : 1. ); 
        int regionEtaRec;
        if(sign > 0) regionEtaRec = iEta + 11;
        if(sign < 0) regionEtaRec = 10 - iEta;
        treeVariables_.L1preNonisoIeta[isocounter1] = regionEtaRec;
        treeVariables_.L1preNonisoIphi[isocounter1] = int(((*ipsbel)>>10)&0x1f);
        treeVariables_.L1preNonisoRank[isocounter1] = (*ipsbel)&0x3f;
        isocounter1++;
        //}
    }//loop Noniso
    //treeVariables_.nbOfL1preNonisoCands = isocounter1;

    psbel.clear();
    psbel.push_back(psb.aData(6));
    psbel.push_back(psb.aData(7));
    psbel.push_back(psb.bData(6));
    psbel.push_back(psb.bData(7));
    isocounter1 = 0;
    treeVariables_.nbOfL1preIsoCands = psbel.size();
    for(ipsbel=psbel.begin(); ipsbel!=psbel.end(); ipsbel++) {
        //int rank = (*ipsbel)&0x3f;
        //double pt = emScale->et(rank); 
        //if(rank>0 && pt>0.6) {// firing EG1
        int iEta = int(((*ipsbel)>>6)&7);
        int sign = ( ((*ipsbel>>9)&1) ? -1. : 1. ); 
        int regionEtaRec;
        if(sign > 0) regionEtaRec = iEta + 11;
        if(sign < 0) regionEtaRec = 10 - iEta;
        treeVariables_.L1preIsoIeta[isocounter1] = regionEtaRec;
        treeVariables_.L1preIsoIphi[isocounter1] = int(((*ipsbel)>>10)&0x1f);
        treeVariables_.L1preIsoRank[isocounter1] = (*ipsbel)&0x3f;
        isocounter1++;
        //}
    }//loop Iso
    //treeVariables_.nbOfL1preIsoCands = isocounter1;

 
     //ontime fire INFO
   
     const L1GtPsbWord psb1 = gtRecord->gtPsbWord(0xbb0d,0);
     //psb.print(cout); 
     vector<int> psbel1;
     psbel1.push_back(psb1.aData(4));
     psbel1.push_back(psb1.aData(5));
     psbel1.push_back(psb1.bData(4));
     psbel1.push_back(psb1.bData(5));
     isocounter1 = 0;
     treeVariables_.nbOfL1NonisoCands = psbel1.size();
     std::vector<int>::const_iterator ipsbel1;
     for(ipsbel1=psbel1.begin(); ipsbel1!=psbel1.end(); ipsbel1++) {
         //int rank = (*ipsbel)&0x3f;
         //double pt = emScale->et(rank); 
         //if(rank>0 && pt>0.6) {// firing EG1
         int iEta = int(((*ipsbel1)>>6)&7);
         int sign = ( ((*ipsbel1>>9)&1) ? -1. : 1. ); 
         int regionEtaRec;
         if(sign > 0) regionEtaRec = iEta + 11;
         if(sign < 0) regionEtaRec = 10 - iEta;
         treeVariables_.L1NonisoIeta[isocounter1] = regionEtaRec;
         treeVariables_.L1NonisoIphi[isocounter1] = int(((*ipsbel1)>>10)&0x1f);
         treeVariables_.L1NonisoRank[isocounter1] = (*ipsbel1)&0x3f;
         isocounter1++;
         //}
     }//loop Noniso
     //treeVariables_.nbOfL1preNonisoCands = isocounter1;
 
     psbel1.clear();
     psbel1.push_back(psb1.aData(6));
     psbel1.push_back(psb1.aData(7));
     psbel1.push_back(psb1.bData(6));
     psbel1.push_back(psb1.bData(7));
     isocounter1 = 0;
     treeVariables_.nbOfL1IsoCands = psbel1.size();
     for(ipsbel1=psbel1.begin(); ipsbel1!=psbel1.end(); ipsbel1++) {
         //int rank = (*ipsbel)&0x3f;
         //double pt = emScale->et(rank); 
         //if(rank>0 && pt>0.6) {// firing EG1
         int iEta = int(((*ipsbel1)>>6)&7);
         int sign = ( ((*ipsbel1>>9)&1) ? -1. : 1. ); 
         int regionEtaRec;
         if(sign > 0) regionEtaRec = iEta + 11;
         if(sign < 0) regionEtaRec = 10 - iEta;
         treeVariables_.L1IsoIeta[isocounter1] = regionEtaRec;
         treeVariables_.L1IsoIphi[isocounter1] = int(((*ipsbel1)>>10)&0x1f);
         treeVariables_.L1IsoRank[isocounter1] = (*ipsbel1)&0x3f;
         isocounter1++;
         //}
     }//loop Iso
     //treeVariables_.nbOfL1preIsoCands = isocounter1;

    //POST-FIRING INFO
    const L1GtPsbWord psb2 = gtRecord->gtPsbWord(0xbb0d, 1);
    //psb.print(cout); 
    vector<int> psbel2;
    psbel2.push_back(psb2.aData(4));
    psbel2.push_back(psb2.aData(5));
    psbel2.push_back(psb2.bData(4));
    psbel2.push_back(psb2.bData(5));
    treeVariables_.nbOfL1postNonisoCands = psbel2.size();
    isocounter1 = 0;
    std::vector<int>::const_iterator ipsbel2;
    for(ipsbel2=psbel2.begin(); ipsbel2!=psbel2.end(); ipsbel2++) {
        //int rank = (*ipsbel2)&0x3f;
        //double pt = emScale->et(rank); 
        //if(rank>0 && pt>0.6) {// firing EG1
        int iEta = int(((*ipsbel2)>>6)&7);
        int sign = ( ((*ipsbel2>>9)&1) ? -1. : 1. ); 
        int regionEtaRec;
        if(sign > 0) regionEtaRec = iEta + 11;
        if(sign < 0) regionEtaRec = 10 - iEta;
        treeVariables_.L1postNonisoIeta[isocounter1] = regionEtaRec;
        treeVariables_.L1postNonisoIphi[isocounter1] = int(((*ipsbel2)>>10)&0x1f);
        treeVariables_.L1postNonisoRank[isocounter1] = (*ipsbel2)&0x3f;
        isocounter1++;
        //}
    }//loop Noniso
    //treeVariables_.nbOfL1postNonisoCands = isocounter1;

    psbel2.clear();
    psbel2.push_back(psb2.aData(6));
    psbel2.push_back(psb2.aData(7));
    psbel2.push_back(psb2.bData(6));
    psbel2.push_back(psb2.bData(7));
    treeVariables_.nbOfL1postIsoCands = psbel2.size();
    isocounter1 = 0;
    for(ipsbel2=psbel2.begin(); ipsbel2!=psbel2.end(); ipsbel2++) {
        //int rank = (*ipsbel2)&0x3f;
        //double pt = emScale->et(rank); 
        //if(rank>0 && pt>0.6) {// firing EG1
        int iEta = int(((*ipsbel2)>>6)&7);
        int sign = ( ((*ipsbel2>>9)&1) ? -1. : 1. ); 
        int regionEtaRec;
        if(sign > 0) regionEtaRec = iEta + 11;
        if(sign < 0) regionEtaRec = 10 - iEta;
        treeVariables_.L1postIsoIeta[isocounter1] = regionEtaRec;
        treeVariables_.L1postIsoIphi[isocounter1] = int(((*ipsbel2)>>10)&0x1f);
        treeVariables_.L1postIsoRank[isocounter1] = (*ipsbel2)&0x3f;
        isocounter1++;
        //}
    }//loop Iso
    //treeVariables_.nbOfL1postIsoCands = isocounter1;
 





 

    /////////////////////////// 
    // Get TP data  
    ///////////////////////////

    edm::Handle<EcalTrigPrimDigiCollection> tp;
    iEvent.getByLabel(tpCollection_,tp);
    if (print_) std::cout<<"TP collection size="<<tp.product()->size()<<std::endl ;
    
    for (unsigned int i=0;i<tp.product()->size();i++) {
        EcalTriggerPrimitiveDigi d = (*(tp.product()))[i];
        const EcalTrigTowerDetId TPtowid= d.id();
        towerEner tE ;
        tE.iphi_ = TPtowid.iphi() ;
        tE.ieta_ = TPtowid.ieta() ;
	tE.ttFlag_ = d[0].ttFlag();
        tE.tpgADC_ = (d[0].raw()&0xfff) ;
	tE.twrADC = (d[0].raw()&0xff) ;
	// if ((d[0].raw()&0xfff)!=0 ||  (d[0].raw()&0xff)!=0){
	//   if (fabs(TPtowid.ieta()) < 18 && d[0].sFGVB()==0) {
	//       std::cout << i  << " evt:" << myevt<<  " adc size: " << d[0].raw() <<  " tpgADC: " << (d[0].raw()&0xfff) << " twrADC: " << (d[0].raw()&0xff)  <<  " sFGVB: " << d[0].sFGVB() <<   " eta: " << TPtowid.ieta() << endl; 
	//     }
	// }
	tE.sFGVB = (d[0].sFGVB());


        mapTower[TPtowid] = tE ;
    }

    ///////////////////////////
    // Get Emulators TP
    ///////////////////////////

    edm::Handle<EcalTrigPrimDigiCollection> tpEmul ;
    iEvent.getByLabel(tpEmulatorCollection_, tpEmul);
    if (print_) std::cout<<"TPEmulator collection size="<<tpEmul.product()->size()<<std::endl ;
    
    for (unsigned int i=0;i<tpEmul.product()->size();i++) {
        EcalTriggerPrimitiveDigi d = (*(tpEmul.product()))[i];
        const EcalTrigTowerDetId TPtowid= d.id();
        itTT = mapTower.find(TPtowid) ;
        if (itTT != mapTower.end())
	  for (int j=0 ; j<5 ; j++) {
	    (itTT->second).tpgEmul_[j] = (d[j].raw()&0xfff) ;
	    (itTT->second).tpgEmulFlag_[j] = d[j].ttFlag();
	    (itTT->second).tpgEmulsFGVB_[j] = d[j].sFGVB();
	  }
    }

    ///////////////////////////
    // Get nb of crystals read out
    ///////////////////////////

    // Get EB xtal digi inputs
    edm::Handle<EBDigiCollection> digiEB;
    iEvent.getByLabel(digiCollectionEB_, digiEB);
    EBDataFrame dfMaxEB ;
    if(print_) cout << " EB digi collection size " << digiEB.product()->size()  << endl ;
    for (unsigned int i=0;i<digiEB.product()->size();i++) {
        const EBDataFrame & df = (*(digiEB.product()))[i];    
        const EBDetId & id = df.id();
        const EcalTrigTowerDetId towid = id.tower();
        itTT = mapTower.find(towid) ;
        if (itTT != mapTower.end()) {
            (itTT->second).nbXtal_++ ;
            bool fill(false) ;
            if (((itTT->second).tpgADC_ & 0xff)>0) fill = true ;   
	    for (int j=0 ; j<5 ; j++) if (((itTT->second).tpgEmul_[j] & 0xff)>8) fill = true ;			
            if (fill) {
	      if(print_) cout<<"TP="<<((itTT->second).tpgADC_ & 0xff)<<" eta="<<towid.ieta()<<" phi="<<towid.iphi()<<endl ;
	      if (print_) for (int j=0 ; j<5 ; j++) if (((itTT->second).tpgEmul_[j] & 0xff)>8) cout << "tp emul "<<  j << " " << ((itTT->second).tpgEmul_[j] & 0xff)<< endl;
                treeVariablesShape_.ieta = towid.ieta() ;
                treeVariablesShape_.iphi = towid.iphi() ;
                treeVariablesShape_.ixXtal = id.iphi() ;
                treeVariablesShape_.iyXtal = id.ieta() ;
                treeVariablesShape_.TCCid = theMapping_->TCCid(towid);
                treeVariablesShape_.TowerInTCC = theMapping_->iTT(towid);
                const EcalTriggerElectronicsId elId = theMapping_->getTriggerElectronicsId(id) ;
                treeVariablesShape_.strip = elId.pseudoStripId() ;
                treeVariablesShape_.nbOfSamples = df.size() ;
                for (int s=0 ; s<df.size() ; s++) treeVariablesShape_.samp[s] = df[s].adc() ; 
                treeShape_->Fill() ;
            }
        }
    }


    if (useEE_) {
        // Get EE xtal digi inputs
        edm::Handle<EEDigiCollection> digiEE;
        iEvent.getByLabel(digiCollectionEE_, digiEE);
        EEDataFrame dfMaxEE ;
	if(print_) cout << " EE digi collection size " << digiEE.product()->size()  << endl ;
        for (unsigned int i=0;i<digiEE.product()->size();i++) {
            const EEDataFrame & df = (*(digiEE.product()))[i];
            const EEDetId & id = df.id();
            const EcalTrigTowerDetId towid = (*eTTmap_).towerOf(id);
            itTT = mapTower.find(towid) ;
            if (itTT != mapTower.end()) {
                (itTT->second).nbXtal_++ ;
                if (!skipWritingEndcapDigis_)
                {
                    bool fill(false) ;
                    if (((itTT->second).tpgADC_ & 0xff)>0) fill = true ;
                    for (int j=0 ; j<5 ; j++) if (((itTT->second).tpgEmul_[j] & 0xff)>0) fill = true ;
                    for (int s=1 ; s<df.size() ; s++) if (df[s].adc()-df[0].adc()>15) fill = true ;
                    if (fill) {
		      if(print_) cout<<"TP="<<((itTT->second).tpgADC_ & 0xff)<<" eta="<<towid.ieta()<<" phi="<<towid.iphi()<<endl ;
                        treeVariablesShape_.ieta = towid.ieta() ;
                        treeVariablesShape_.iphi = towid.iphi() ;
                        treeVariablesShape_.ixXtal = id.ix() ;
                        treeVariablesShape_.iyXtal = id.iy() ;
			const EcalTriggerElectronicsId elId = theMapping_->getTriggerElectronicsId(id) ;
                        treeVariablesShape_.TCCid = theMapping_->TCCid(towid);
                        treeVariablesShape_.TowerInTCC = theMapping_->iTT(towid);
                        treeVariablesShape_.strip = elId.pseudoStripId() ;
                        treeVariablesShape_.nbOfSamples = df.size() ;
                        for (int s=0 ; s<df.size() ; s++) treeVariablesShape_.samp[s] = df[s].adc() ; 
                        treeShape_->Fill() ;
                    }
                }
            }
        }  // loop over EE digi
    } // UseEE_

    ///////////////////////////
    // Get rechits and spikes
    ///////////////////////////

    edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
    iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);


    
    // channel status
    edm::ESHandle<EcalChannelStatus> pChannelStatus;
    iSetup.get<EcalChannelStatusRcd>().get(pChannelStatus);
  

    //  const EcalChannelStatus *chStatus = pChannelStatus.product();
    //const EcalRecHit * rh; 
   
    // Get EB rechits
    edm::Handle<EcalRecHitCollection> rechitsEB; 
    iEvent.getByLabel(EcalRecHitCollectionEB_, rechitsEB) ;
    //std::cout << " rechitsEB size " << rechitsEB.product()->size() << std::endl;
    float maxRecHitEnergy = 0. ;
    if (rechitsEB.product()->size()!=0) {
      for ( EcalRecHitCollection::const_iterator rechitItr = rechitsEB->begin(); rechitItr != rechitsEB->end(); ++rechitItr ) {   
            EBDetId id = rechitItr->id(); 
            const EcalTrigTowerDetId towid = id.tower();
            itTT = mapTower.find(towid) ;

            if (itTT != mapTower.end()) {
	     
	      double theta = theBarrelGeometry_->getGeometry(id)->getPosition().theta() ;
                (itTT->second).eRec_ += rechitItr->energy()*sin(theta) ;
		if (maxRecHitEnergy < rechitItr->energy()*sin(theta) && rechitItr->energy()*sin(theta) > 1. ){
		  (itTT->second).sevlv_ = sevlv->severityLevel(id, *rechitsEB); 

		}
		//(itTT->second).sevlv = sevlv->severityLevel(*rh); 
		//cout << "severity level barrel " << sevlv->severityLevel(id, *rechitsEB) << endl; 
		(itTT->second).crystNb_++;
            }

          //  uint32_t sev = sevlv->severityLevel( id, *rechitsEB);
            //sevlv->severityLevel( id, *rechitsEB);
	    
            /*
              if (sev == EcalSeverityLevelAlgo::kWeird) {
              itTT = mapTower.find(towid) ;
              if (itTT != mapTower.end()) {
              (itTT->second).spike_++ ;  
              }
              }
            */
        }
    }

    // Get EE rechits
    edm::Handle<EcalRecHitCollection> rechitsEE; 
    if (iEvent.getByLabel(EcalRecHitCollectionEE_, rechitsEE) ) {
      
        for ( EcalRecHitCollection::const_iterator rechitItr = rechitsEE->begin(); rechitItr != rechitsEE->end(); ++rechitItr ) {   
            EEDetId id = rechitItr->id();
            const EcalTrigTowerDetId towid = (*eTTmap_).towerOf(id);
            itTT = mapTower.find(towid) ;
            if (itTT != mapTower.end()) {
                double theta = theEndcapGeometry_->getGeometry(id)->getPosition().theta() ;
                (itTT->second).eRec_ += rechitItr->energy()*sin(theta) ;
		//rh = &*rechitItr;
		(itTT->second).sevlv_ = sevlv->severityLevel(id, *rechitsEE); 
		//cout << "severity level endcap " << sevlv->severityLevel(id, *rechitsEE) << endl;
            }
        }
    }




    ///////////////////////////  
    // fill tree
    ///////////////////////////  
    
    int towerNb = 0 ;
    for (itTT = mapTower.begin() ; itTT != mapTower.end() ; ++itTT) {
       
        // select only non zero towers
        bool fill(true) ;
        bool nonZeroEmul(false) ;
        for (int i=0 ; i<5 ; i++) if (((itTT->second).tpgEmul_[i]&0xff) > 0) nonZeroEmul = true ;
        if ( keepOnlyTowersAboveZero_ && ((itTT->second).tpgADC_&0xff) <= 0 && (!nonZeroEmul) ) fill = false ;
        if (print_ && fill) {
            std::cout<<"ieta="<<(itTT->second).ieta_<<" "<<(itTT->second).iphi_<<" tp="<<((itTT->second).tpgADC_&0xff)<<" tpEmul=" ;
            for (int i=0 ; i<5 ; i++) std::cout<<((itTT->second).tpgEmul_[i]&0xff)<<" " ;
            std::cout<<" nbXtal="<<(itTT->second).nbXtal_ ;
            std::cout<<std::endl ;      
        }
        if (fill) {
            treeVariables_.ieta[towerNb] = (itTT->second).ieta_ ;
            treeVariables_.iphi[towerNb] = (itTT->second).iphi_ ;
            treeVariables_.nbOfXtals[towerNb] = (itTT->second).nbXtal_ ;
            treeVariables_.rawTPData[towerNb] = (itTT->second).tpgADC_ ;
            treeVariables_.rawTPEmul1[towerNb] = (itTT->second).tpgEmul_[0] ;
            treeVariables_.rawTPEmul2[towerNb] = (itTT->second).tpgEmul_[1] ;
            treeVariables_.rawTPEmul3[towerNb] = (itTT->second).tpgEmul_[2] ;
            treeVariables_.rawTPEmul4[towerNb] = (itTT->second).tpgEmul_[3] ;
            treeVariables_.rawTPEmul5[towerNb] = (itTT->second).tpgEmul_[4] ;
            treeVariables_.rawTPEmulttFlag1[towerNb] = (itTT->second).tpgEmulFlag_[0] ;
            treeVariables_.rawTPEmulttFlag2[towerNb] = (itTT->second).tpgEmulFlag_[1] ;
            treeVariables_.rawTPEmulttFlag3[towerNb] = (itTT->second).tpgEmulFlag_[2] ;
            treeVariables_.rawTPEmulttFlag4[towerNb] = (itTT->second).tpgEmulFlag_[3] ;
            treeVariables_.rawTPEmulttFlag5[towerNb] = (itTT->second).tpgEmulFlag_[4] ;
            treeVariables_.rawTPEmulsFGVB1[towerNb] = (itTT->second).tpgEmulsFGVB_[0] ;
            treeVariables_.rawTPEmulsFGVB2[towerNb] = (itTT->second).tpgEmulsFGVB_[1] ;
            treeVariables_.rawTPEmulsFGVB3[towerNb] = (itTT->second).tpgEmulsFGVB_[2] ;
            treeVariables_.rawTPEmulsFGVB4[towerNb] = (itTT->second).tpgEmulsFGVB_[3] ;
            treeVariables_.rawTPEmulsFGVB5[towerNb] = (itTT->second).tpgEmulsFGVB_[4] ;
            treeVariables_.crystNb[towerNb] = (itTT->second).crystNb_ ;
            treeVariables_.eRec[towerNb] = (itTT->second).eRec_ ;
            treeVariables_.sevlv[towerNb] = (itTT->second).sevlv_ ;
            treeVariables_.ttFlag[towerNb] = (itTT->second).ttFlag_ ;
	    
            treeVariables_.spike[towerNb] = (itTT->second).spike_ ;
	    treeVariables_.twrADC[towerNb] =  (itTT->second).twrADC;
	    treeVariables_.sFGVB[towerNb] =  (itTT->second).sFGVB;
	    
            if (abs(treeVariables_.ieta[towerNb])>17) {
                unsigned int maxEmul = 0 ;
                for (int i=0 ; i<5 ; i++) if (((itTT->second).tpgEmul_[i]&0xff) > maxEmul) maxEmul = ((itTT->second).tpgEmul_[i]&0xff) ;
            }
            towerNb++ ;
        }

    }
    
    treeVariables_.nbOfTowers = towerNb ;
    tree_->Fill() ;

}

DEFINE_FWK_MODULE(EcalTPGAnalyzer);
