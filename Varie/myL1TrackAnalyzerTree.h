#ifndef myL1TrackAnalyzerTree_h
#define myL1TrackAnalyzerTree_h

/////////////////////////
//       HEADERS       //
/////////////////////////

////////////////
// CLASS HEADER
// No more necessary in the current "no *.h file" implementation

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h" 
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h" 
#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSMatchedRecHit2DCollection.h"
//
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/Ref.h"
//
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
//
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
#include "SimDataFormats/SLHC/interface/LocalStub.h"
//
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
//
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
//
#include "SimDataFormats/SLHC/interface/L1CaloCluster.h"
//
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
//
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "TrackingTools/GeomPropagators/interface/HelixArbitraryPlaneCrossing.h"
////////////////////////
// FAST SIMULATION STUFF
#include "FastSimulation/Particle/interface/RawParticle.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "SimGeneral/HepPDTRecord/interface/PDTRecord.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimGeneral/HepPDTRecord/interface/PdtEntry.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
//
#include "SLHCUpgradeSimulations/Utilities/interface/StackedTrackerGeometryRecord.h"
#include "SLHCUpgradeSimulations/Utilities/interface/StackedTrackerGeometry.h"
#include "SLHCUpgradeSimulations/Utilities/interface/StackedTrackerDetUnit.h"
#include "SLHCUpgradeSimulations/Utilities/interface/StackedTrackerDetId.h"

////////////////
// PHYSICS TOOLS
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "RecoTracker/TkSeedGenerator/interface/FastHelix.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "RecoTauTag/TauTagTools/interface/GeneratorTau.h"
//
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementVector.h"
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
#include "SLHCUpgradeSimulations/Utilities/interface/constants.h"

#include "SLHCUpgradeSimulations/L1Trigger/interface/CircleFit.h"
#include "SLHCUpgradeSimulations/L1Trigger/interface/LineFit.h"
#include "SimDataFormats/SLHC/interface/L1Track.h"
///mine
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

///////////////
// ROOT HEADERS
#include <TROOT.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH1.h>

//////////////
// STD HEADERS
#include <memory>
#include <string>
#include <iostream>

//////////////
// NAMESPACES
// I hate them, but there is something
// magic behind and beyond ...
using namespace std;
using namespace edm;
using namespace reco;
using namespace cmsUpgrades;
using namespace l1slhc;
using namespace l1extra;

//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

class TTree;
class TFile;
class TH1D;
class TH2D;
class TGraph;
class RectangularPixelTopology;
class TransientInitialStateEstimator;
class MagneticField;
class TrackerGeometry;
class TrajectoryStateOnSurface;
class PTrajectoryStateOnDet;
//
class myL1TrackAnalyzerTree : public edm::EDAnalyzer
{
  /// Public methods
  public:
    /// Constructor/destructor
    explicit myL1TrackAnalyzerTree(const edm::ParameterSet& conf);
    virtual ~myL1TrackAnalyzerTree();
    // Typical methods used on Loops over events
    virtual void beginJob(const edm::EventSetup& es);
    virtual void endJob();
    virtual void analyze(const edm::Event& e, const edm::EventSetup& es);

    /// Some Type definitions
    typedef edm::Ptr< GlobalStub<Ref_PixelDigi_> > GlobalStubPtrType;
    typedef std::set< std::pair<unsigned int , GlobalStubPtrType> > TrackletMap;

    typedef edm::Ptr< Tracklet<Ref_PixelDigi_> > TrackletPtrType;
    typedef edm::Ptr< L1Track<Ref_PixelDigi_> > L1TrackPtrType;

    /// Define the type
    typedef std::pair< int, std::vector< GlobalStub<Ref_PixelDigi_> > > L1Trk_base;


    
  /// Protected methods only internally used
  protected:
    /// Get Pt from a Stub with the 2 points + tangent constraints
    double getPtFromStub( double xx, double yy, float ddxx, float ddyy );
    /// EtaPhi space
    double deltaPhi( double phi1, double phi2 );
    double deltaPhiNP( double phi1, double phi2 );
    double deltaEta( double eta1, double eta2 );
    double deltaR( double phi1, double eta1, double phi2, double eta2 );

    std::pair<double,PSimHit> makeHit( const GeomDetUnit* dU, BaseParticlePropagator* tP,
                                       FreeTrajectoryState hV, int pID, int tkID );

                                       /// Make L1Track with 1 Tracklet and N Stubs
    //std::vector< std::pair< int, std::vector< GlobalStub<Ref_PixelDigi_> > > >
    std::vector< L1Track<Ref_PixelDigi_> >
        makeL1Track( const edm::EventSetup& es,
                     const MagneticField *magnetF,
                     Tracklet_PixelDigi_Collection::const_iterator seed,
                     std::vector< GlobalStub_PixelDigi_Collection::const_iterator > bricks,
                     bool erase_nulls );

    struct predicate_unique {
      bool operator() ( std::pair< int, std::vector< GlobalStub<Ref_PixelDigi_> > > a, std::pair< int, std::vector< GlobalStub<Ref_PixelDigi_> > > b ) {
        if (a.second.size() != b.second.size()) return false;
        else {
          for (unsigned int q=0; q<a.second.size(); q++) {
            if (a.second.at(q).localStub() != b.second.at(q).localStub()) return false;
          }
        }
        return true;
      }
    };

    struct predicate_sort {
      bool operator() ( std::pair< int, std::vector< GlobalStub<Ref_PixelDigi_> > > a, std::pair< int, std::vector< GlobalStub<Ref_PixelDigi_> > > b ) {
        if (a.second.size() != b.second.size())
          return (a.second.size()<b.second.size());
        else {
          for (unsigned int q=0; q<a.second.size(); q++) {
            if (a.second.at(q).localStub() != b.second.at(q).localStub())
              return (a.second.at(q).localStub() < b.second.at(q).localStub());
          }
        }
        return false;
      }
    };


                     
  /// Private methods and variables
  private:

    int PDGCode;
    double lumiZ;
    bool testedGeometry;
    bool beamspot00;
    int seedsuperlayer;
    int numberstubs;
    int howmanytight;

    double probwindow;

    /// TO CHECK ALL EVENTS ARE PROCESSED
    TH1D* h_EvtCnt;

    TH1D *hSimTrkBeforeCutEta;

    TH1D* hSimTrackPt;
    TH1D* hSimTrackPx ;
    TH1D* hSimTrackPy ;
    TH1D* hSimTrackPz ;
    TH1D* hSimTrackEta;
    TH1D* hSimTrackPhi;
    TH1D* hSimTrackvtxx;
    TH1D* hSimTrackvtxy;
    TH1D* hSimTrackvtxz;
    TH1D* hSimTrackq  ;
	      
    TH1D* hL1TrackPt ;
    TH1D* hL1TrackPx ;
    TH1D* hL1TrackPy ;
    TH1D* hL1TrackPz ;
    TH1D* hL1TrackEta;
    TH1D* hL1TrackPhi;
    TH1D* hL1Trackvtxx;
    TH1D* hL1Trackvtxy;
    TH1D* hL1Trackvtxz;
    TH1D* hL1Trackq  ;

    TH1D* hDeltaTrackPt ;
    TH1D* hDeltaTrackPx ;
    TH1D* hDeltaTrackPy ;
    TH1D* hDeltaTrackPz ;
    TH1D* hDeltaTrackEta;
    TH1D* hDeltaTrackPhi;
    TH1D* hDeltaTrackvtxx;
    TH1D* hDeltaTrackvtxy;
    TH1D* hDeltaTrackvtxz;
    TH1D* hDeltaTrackq  ;
  
   TH1D*  hRecoTrackSize ;
  TH1D*  hRecoTrackPt  ;
  TH1D*  hRecoTrackPx ;
  TH1D*  hRecoTrackPy ;
  TH1D*  hRecoTrackPz ; 
  TH1D*  hRecoTrackEta ;
  TH1D*  hRecoTrackPhi ;
  TH1D*  hRecoTrackvtxx;
  TH1D*  hRecoTrackvtxy;
  TH1D*  hRecoTrackvtxz;
  TH1D*  hRecoTrackq   ;
  
  TH1D *hgenMuonPt;

  TH1D* hMuPt    ;
  TH1D* hMuPx    ;
  TH1D* hMuPy    ;
  TH1D* hMuPz    ; 
  TH1D* hMuEta   ; 
  TH1D* hMuPhi   ;
  TH1D* hMuvtxx  ;
  TH1D* hMuvtxy  ;
  TH1D* hMuvtxz  ;
  TH1D* hMuq     ;
  TH1D* hMusize  ;
  
  TH1D* hMuFromHiggsPt    ;
  TH1D* hMuFromHiggsPx    ;
  TH1D* hMuFromHiggsPy    ;
  TH1D* hMuFromHiggsPz    ;
  TH1D* hMuFromHiggsEta   ;
  TH1D* hMuFromHiggsPhi   ;
  TH1D* hMuFromHiggsvtxx  ;
  TH1D* hMuFromHiggsvtxy  ;
  TH1D* hMuFromHiggsvtxz  ;
  TH1D* hMuFromHiggsq     ;
  
    TH2D* hDET_LAYER_r;
    TH2D* hDET_IPHI_p;
    TH2D* hDET_IZ_z;

    TH2D* hL1Trk_e_ST_e;
    TH2D* hL1Trk_p_ST_p;
    TH2D* hL1Trk_Pt_ST_Pt;
    TH2D* hL1Trk_xvtx_ST_xvtx;
    TH2D* hL1Trk_yvtx_ST_yvtx;
    TH2D* hL1Trk_zvtx_ST_zvtx;
    TH2D* hL1Trk_q_ST_q;
    TH2D* hL1Trk_chi2rphi_n;
    TH2D* hL1Trk_chi2rz_n;
    
    TH1D* hL1Trk_FakeRate_n;
    TH1D* hL1Trk_Number_n;
    TH1D* hST_Pt_matched;
    TH1D* hST_Pt;
    TH1D* hST_Pt_matched_any;    

    TH1D* hST_e_Pt10_matched;
    TH1D* hST_e_Pt10;
    TH1D* hST_e_Pt10_matched_any;    

  TH1D* htightNumber ;
  TH2D* hValidL1TracksvsL1Number;
  TH2D* hMatchedL1TracksvsL1Number;

  TH1D * hrecoTrackNumber; 
  TH1D * hIsoPtRel;
  TH1D * hIsoPt;
  TH1D * hL1IsoPt;
  TH1D * hL1IsoPtRel;
  TH1D * hL1size;

    TH1D* hL1TrackNoMatchPt ;
    TH1D* hL1TrackNoMatchPx ;
    TH1D* hL1TrackNoMatchPy ;
    TH1D* hL1TrackNoMatchPz ;
    TH1D* hL1TrackNoMatchEta;
    TH1D* hL1TrackNoMatchPhi;
    TH1D* hL1TrackNoMatchvtxx;
    TH1D* hL1TrackNoMatchvtxy;
    TH1D* hL1TrackNoMatchvtxz;
    TH1D* hL1TrackNoMatchq  ;

    /// Containers of parameters passed by python
    /// configuration file
    edm::ParameterSet config;

    /// Geometry handles etc
    edm::ESHandle<TrackerGeometry> geometryESH;
    const TrackerGeometry*  theGeometry;
    edm::ESHandle<cmsUpgrades::StackedTrackerGeometry> stackedgeometryESH;
    const cmsUpgrades::StackedTrackerGeometry*  theStackedGeometry;


  typedef  struct BRANCH { 
    std::vector<bool> muRecoIsGlobal;
    std::vector<int> simTrkQ, L1TrkQ, recoTrkQ, simTrkPdgId, L1TrkQNoMatch;
    std::vector< pair <int, int > >genWQ, genMuQ, genWPdgId, genMuPdgId;
    std::vector<size_t> simTrkSize, L1TrkSize, recoTrkSize;
    std::vector<bool> recoTrkIsMuon, recoTrkIsHiggsMuon;
    std::vector<double> genHiggsPt, genHiggsPx,genHiggsPy, genHiggsPz, genHiggsEta, genHiggsPhi, genHiggsVtxX,genHiggsVtxY,genHiggsVtxZ,genHiggsMass,genHiggsE,
      simTrkPt, simTrkPx,simTrkPy, simTrkPz, simTrkEta, simTrkPhi, simTrkVtxX,simTrkVtxY,simTrkVtxZ, simTrkId, 
      L1TrkPt, L1TrkPx,L1TrkPy, L1TrkPz, L1TrkEta, L1TrkPhi, L1TrkVtxX,L1TrkVtxY,L1TrkVtxZ,L1TrkId, 
      L1TrkPtNoMatch, L1TrkPxNoMatch,L1TrkPyNoMatch, L1TrkPzNoMatch, L1TrkEtaNoMatch, L1TrkPhiNoMatch, L1TrkVtxXNoMatch,L1TrkVtxYNoMatch,L1TrkVtxZNoMatch,L1TrkIdNoMatch, 
      recoTrkPt, recoTrkPx,recoTrkPy, recoTrkPz, recoTrkEta, recoTrkPhi, recoTrkVtxX,recoTrkVtxY,recoTrkVtxZ,recoTrkId,
      muRecoTrkPt, muRecoTrkPx, muRecoTrkPy, muRecoTrkPz, muRecoTrkEta, muRecoTrkPhi, muRecoTrkVtxX, muRecoTrkVtxY, muRecoTrkVtxZ, muRecoTrkId;
    std::vector< pair <double, double > > genWPt, genWPx,genWPy, genWPz, genWEta, genWPhi, genWVtxX,genWVtxY,genWVtxZ, genWMass, genWE,
      genMuPt, genMuPx,genMuPy, genMuPz, genMuEta, genMuPhi, genMuVtxX,genMuVtxY,genMuVtxZ, genMuE ,
      muHiggsRecoTrkPt, muHiggsRecoTrkPx, muHiggsRecoTrkPy, muHiggsRecoTrkPz, muHiggsRecoTrkEta, muHiggsRecoTrkPhi, muHiggsRecoTrkVtxX, muHiggsRecoTrkVtxY,
      muHiggsRecoTrkVtxZ, muHiggsRecoTrkId;
    std::vector< pair <int, int > > muHiggsRecoIsGlobal;
  };

  struct store_trks {
      int trkid;
      int doublestack;
      int phisector;    
    };
    typedef store_trks StoreTrks;

    struct store_hits {
      int origindoublestack;
      int targetlayer;
      int targetiphisector;
      int targetizsector;
      PSimHit hit;
    };
    typedef store_hits StoreHits;

    struct store_PtEta {
      double pt;
      double eta;
    };
  typedef store_PtEta PtEta;

  BRANCH branch;
  TTree *ntuple;

};
#endif

