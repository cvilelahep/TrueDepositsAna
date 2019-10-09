////////////////////////////////////////////////////////////////////////
//
// \file CAFMaker_module.cc
//
// Chris Marshall's version
// Largely based on historical FDSensOpt/CAFMaker_module.cc
//
///////////////////////////////////////////////////////////////////////

#ifndef TRUEDEPOSITSANA_H
#define TRUEDEPOSITSANA_H

// Generic C++ includes
#include <iostream>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Services/Optional/TFileService.h" 

//#include "Utils/AppInit.h"
//#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
//#include "nusimdata/SimulationBase/MCFlux.h"
//#include "larcoreobj/SummaryData/POTSummary.h"
//#include "dune/FDSensOpt/FDSensOptData/MVASelectPID.h"
//#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
//#include "dune/CVN/func/InteractionType.h"
//#include "dune/CVN/func/Result.h"
//#include "dune/RegCVN/func/RegCVNResult.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
//#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "lardataobj/Simulation/SimChannel.h"

// nutools
//#include "nutools/EventGeneratorBase/GENIE/GENIE2ART.h"

// dunerw stuff
//#include "systematicstools/interface/ISystProviderTool.hh"
//#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
//#include "systematicstools/utility/exceptions.hh"

// root
#include "TTree.h"
#include "TH3F.h"
#include "TRandom3.h"
//#include "TH1D.h"
//#include "TH2D.h"

// pdg
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

// genie
//#include "EVGCore/EventRecord.h"
//#include "GHEP/GHepParticle.h"

//#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
//#include "nutools/RandomUtils/NuRandomService.h"
//#include "CLHEP/Random/RandFlat.h"

#include <Eigen/Dense>

namespace dunemva {

float offset[] = { 0., 305., 5. };
float collarLo[] = {-320., -120., 30. };
float collarHi[] = { 320.,  120., 470.};
float detLo[] = {-350., -150, 0.};
float detHi[] = {+350., +150, 500.};

float containmentCut = 30;

unsigned int evNumToDisplay = 5;
unsigned int nSamplesToDisplay = 30;

unsigned int N = 1000; // 10000 random rotations per event
  
TH3F* display3D(Eigen::Matrix3Xf hitSegments, std::vector<float> energyDeposits,  std::string histName) {
  art::ServiceHandle<art::TFileService> fs;
  std::cout << "got fs"<< std::endl;

  //  TH3F * hXYZ = fs->makeAndRegister<TH3F>((histName+"_XYZ").c_str(),(histName+"_XYZ").c_str(),(histName+"_XYZ").c_str(), "; x [cm]; y [cm]; z [cm]", 200, -500., 500., 200, -500., 500., 200, -250., 750.);
  TH3F * hXYZ = fs->makeAndRegister<TH3F>((histName+"_XYZ").c_str(),(histName+"_XYZ").c_str(),(histName+"_XYZ").c_str(), "; x [cm]; y [cm]; z [cm]", 100, -500., 500., 100, -500., 500., 100, -250., 750.);
  std::cout << "made and registered th3F"<< std::endl;
    //  TH3F * hXYZ = new TH3F((histName+"_XYZ").c_str(), "; x [cm]; y [cm]; z [cm]", 200, -500., 500., 200, -500., 500., 200, -250, 750);

  for (unsigned int i = 0; i < energyDeposits.size(); i++){
    std::cout << "Filling histogram "<< hitSegments(0, i)/10 -offset[0] << " " <<  hitSegments(1, i)/10.-offset[1] << " " <<  hitSegments(2, i)/10.-offset[2] << " " << energyDeposits[i] << std::endl;
    hXYZ->Fill( hitSegments(0, i)/10.-offset[0], hitSegments(1, i)/10.-offset[1], hitSegments(2, i)/10.-offset[2],energyDeposits[i] );
  }

  return hXYZ;
}

bool isFV ( float * vtx){
  bool inDeadRegion = false;
  for( int i = -3; i <= 3; ++i ) {
    // 0.5cm cathode in the middle of each module, plus 0.5cm buffer
    double cathode_center = i*102.1;
    if( vtx[0]/10. - offset[0] > cathode_center - 0.75 && vtx[0]/10. - offset[0] < cathode_center + 0.75 ) inDeadRegion = true;
    
    // 1.6cm dead region between modules (0.5cm module wall and 0.3cm pixel plane, x2)
    // don't worry about outer boundary because events are only generated in active Ar + insides
    double module_boundary = i*102.1 + 51.05;
    if( i <= 2 && vtx[0]/10. - offset[0] > module_boundary - 1.3 && vtx[0]/10. - offset[0] < module_boundary + 1.3 ) inDeadRegion = true;
  }
  for( int i = 1; i <= 4; ++i ) {
    // module boundaries in z are 1.8cm (0.4cm ArCLight plane + 0.5cm module wall, x2)
    // module is 102.1cm wide, but only 101.8cm long due to cathode (0.5cm) being absent in length
    // but ArCLight is 0.1cm thicker than pixel plane so it's 0.3cm difference
    // positions are off-set by 0.6 because I defined 0 to be the upstream edge based on the active volume
    // by inspecting a plot, and aparently missed by 3 mm, but whatever
    // add 8mm = 2 pad buffer due to worse position resolution in spatial dimension z compared to timing direction x
    // so total FV gap will be 1.8 + 2*0.8 = 3.4cm
    double module_boundary = i*101.8 - 0.6;
    if( vtx[2]/10. - offset[2] > module_boundary - 1.7 && vtx[2]/10. - offset[2] < module_boundary + 1.7 ) inDeadRegion = true;
  }
  
  return (
	  abs(vtx[0]/10. - offset[0]) < 300 &&
	  abs(vtx[1]/10. - offset[1]) < 100 &&
	  vtx[2]/10. - offset[2] > 50 &&
	  vtx[2]/10. - offset[2] < 350 &&
	  !inDeadRegion
	  );
}


bool isContained( Eigen::Matrix3Xf hitSegments, std::vector<float> energyDeposits ){
  
  float collarEnergy = 0.;

  for (unsigned int i = 0; i < energyDeposits.size(); i++){
    for (int dim = 0; dim < 3; dim++){
      if ( hitSegments(dim, i)/10.-offset[dim] < collarLo[dim] and hitSegments(dim, i)/10.-offset[dim] > detLo[dim])
	{ 
	  collarEnergy += energyDeposits[i];
	  break;
	}
      else if ( hitSegments(dim, i)/10.-offset[dim] > collarHi[dim] and hitSegments(dim, i)/10.-offset[dim] < detHi[dim])
	{
	  collarEnergy += energyDeposits[i];
	  break;
	}
    }
  }

  return collarEnergy < containmentCut;
	   
}

double totEinActiveVol( Eigen::Matrix3Xf hitSegments, std::vector<float> energyDeposits ){
  
  float totEnergy = 0.;

  for (unsigned int i = 0; i < energyDeposits.size(); i++){
    for (int dim = 0; dim < 3; dim++){
      if ( hitSegments(dim, i)/10.-offset[dim] < detLo[dim]) 
	{
	  break;
	}
      else if ( hitSegments(dim, i)/10.-offset[dim] > detHi[dim]) 
	{
	  break;
	}
      totEnergy += energyDeposits[i];

    }
  }

  return totEnergy;
}


  class TrueDepositsAna : public art::EDAnalyzer {

  public:

    explicit TrueDepositsAna(fhicl::ParameterSet const& pset);
    virtual ~TrueDepositsAna();
    void beginJob() override;
    void endJob() override;
    void beginSubRun(const art::SubRun& sr) override;
    void endSubRun(const art::SubRun& sr) override;
    void reconfigure(fhicl::ParameterSet const& pset) /*override*/;
    void analyze(art::Event const & evt) override;


  private:

    double avgCosDip;
    double avgBeamDir[3];


    std::unique_ptr< TRandom3 > fRandom;

    std::string fHitsModuleLabel;

    std::vector<TH3F*> displays;

    TTree* fTree;

    float ndGeoEff, nd_Ehad, Ev, nd_vtx_x, LepE;

    //    std::vector<float> hadPosX;
    //    std::vector<float> hadPosY;
    //    std::vector<float> hadPosZ;
    //    std::vector<float> hadEdep;


  }; // class TrueDepositsAna


  //------------------------------------------------------------------------------
  TrueDepositsAna::TrueDepositsAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  dunemva::TrueDepositsAna::~TrueDepositsAna(){}

  //------------------------------------------------------------------------------
  void TrueDepositsAna::reconfigure(fhicl::ParameterSet const& pset) 
  {
    fHitsModuleLabel = pset.get<std::string>("HitsModuleLabel");
    
    //    unsigned int seed = 31415;
    fRandom = std::make_unique<TRandom3>();
    
    //    createEngine(seed);
    
    //    art::ServiceHandle< art::NuRandomService > rng;
    //    CLHEP::HepRandomEngine &engine = rng->getEngine();
    //    fRandFlat        = std::make_unique< CLHEP::RandFlat        >(engine);
    
    avgCosDip = -0.09737;
    avgBeamDir[0]  = 0.;
    avgBeamDir[1] = avgCosDip;
    avgBeamDir[2] =   sqrt(1. - pow(avgCosDip,2) );
    

  }


  //------------------------------------------------------------------------------
  void TrueDepositsAna::beginJob()
  {

    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("trueDeposits", "trueDeposits");

    fTree->Branch("ndGeoEff", &ndGeoEff);
    fTree->Branch("Ev", &Ev);
    fTree->Branch("nd_vtx_x", &nd_vtx_x);
    fTree->Branch("nd_Ehad", &nd_Ehad);
    fTree->Branch("LepE", &LepE);
//    fTree->Branch("hadPosX", &hadPosX);
//    fTree->Branch("hadPosY", &hadPosY);
//    fTree->Branch("hadPosZ", &hadPosZ);
//    fTree->Branch("hadEdep", &hadEdep);

  }

  //------------------------------------------------------------------------------
  void TrueDepositsAna::beginSubRun(const art::SubRun& sr)
  {
  }

  //------------------------------------------------------------------------------
  void TrueDepositsAna::analyze(art::Event const & evt)
  {

    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::Handle< std::vector<simb::MCTruth> > mct;
    std::vector< art::Ptr<simb::MCTruth> > truth;
    if( evt.getByLabel("generator", mct) )
      art::fill_ptr_vector(truth, mct);
    else
      mf::LogWarning("TrueDepositsAna") << "No MCTruth.";


    std::map<int,int> primary_pdg; // track ID to PDG code of primary particle
    std::map<int,int> tid_to_mother; // track ID to mother track ID

    double vtx_x, vtx_y, vtx_z;

    for(size_t i=0; i<truth.size(); i++){
      if(i>0){ // multiple interactions per art::Event not currently supported
        mf::LogWarning("TrueDepositsAna") << "Skipping MC truth index " << i;
        continue;
      }
      
      vtx_x     = truth[i]->GetNeutrino().Nu().Vx()*10;// mm for consistency
      vtx_y     = truth[i]->GetNeutrino().Nu().Vy()*10;// mm for consistency
      vtx_z     = truth[i]->GetNeutrino().Nu().Vz()*10;// mm for consistency

      Ev     = truth[i]->GetNeutrino().Nu().E();
      LepE = truth[i]->GetNeutrino().Lepton().Momentum().T();

      const std::vector< const simb::MCParticle* > parts = pi_serv->MCTruthToParticles_Ps(truth[i]);
      for( size_t pp = 0; pp < parts.size(); ++pp ) {
	int tid = parts[pp]->TrackId();
	int mother = parts[pp]->Mother();
	tid_to_mother.emplace(tid, mother);
	if( mother == 0 ) primary_pdg.emplace(tid, parts[pp]->PdgCode());
	int primaryTid = tid;
	while( mother != 0 ) {
	  primaryTid = mother;
	  mother = tid_to_mother[primaryTid]; // find primary
	}
	if( primary_pdg.find(primaryTid) == primary_pdg.end() ) std::cout << "Something is wrong" << std::endl;
	else {
	  primary_pdg.emplace(tid, primary_pdg[primaryTid]);
	}
      }
    } // loop through MC truth i

    // truth-matching hadronic energy to particles
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
    
    std::vector<float> hadPos;
    std::vector<float> hadEdep;

    hadPos.reserve(hitlist.size()*3);
    hadEdep.reserve(hitlist.size());

    double thisEhad = 0.;

    // There might be a better way to do this. Can we access IDE's directly somehow?
    for( size_t i = 0; i < hitlist.size(); ++i ) {
      art::Ptr<recob::Hit> hit = hitlist[i];
    
      if (hit->WireID().Plane != 2) continue; // Do not double count energy deposits

      

      std::vector<sim::IDE> IDEs = bt_serv->HitToAvgSimIDEs(hit);
      for ( sim::IDE IDE : IDEs ){
	int primpdg = primary_pdg[IDE.trackID];
      
	if( abs(primpdg) != 11 && abs(primpdg) != 13 && abs(primpdg) != 15 ) {
	  hadEdep.emplace_back(IDE.energy);
	  hadPos.emplace_back(IDE.x*10.);
	  hadPos.emplace_back(IDE.y*10.);
	  hadPos.emplace_back(IDE.z*10.);
	  //	  std::cout << IDE.energy << " " << IDE.x << " " << IDE.y <<" "<< IDE.z << std::endl;
	  thisEhad+=IDE.energy;
	}
      }
    }

    std::cout << "Ev " << Ev << " EhadDep " << thisEhad << std::endl;

    const int nNonLepHitSegs = hadEdep.size();
    
    Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > hitSegPosOrig(hadPos.data(),3,nNonLepHitSegs,Eigen::OuterStride<>(3)); 

    ndGeoEff = 0.;

    
    // FIRST "MOVE" EVENT TO Random ND POS
    float ndVtx[] = {-1e10, -1e10, -1e10};
    Eigen::Affine3f ndPos (Eigen::Translation3f(Eigen::Vector3f(ndVtx[0]-vtx_x,ndVtx[1]-vtx_y,ndVtx[2]-vtx_z)));
    while( not isFV(ndVtx) ){
      ndVtx[0] = 10*((fRandom->Rndm()-0.5)*2*300 + offset[0]);
      ndVtx[1] = 10*((fRandom->Rndm()-0.5)*2*100 + offset[1]);
      ndVtx[2] = 10*((fRandom->Rndm()-0.5)*2*(350-50)/2. + 200. + offset[2]);
      ndPos = Eigen::Translation3f(Eigen::Vector3f(ndVtx[0]-vtx_x,ndVtx[1]-vtx_y,ndVtx[2]-vtx_z));
    }

    nd_vtx_x = ndVtx[0];

    nd_Ehad = -1.;
    
    if (evt.event() == evNumToDisplay) {
      std::cout << "CALLING DISPLAY3D" << std::endl;
      displays.emplace_back(display3D(ndPos * hitSegPosOrig, hadEdep, std::string("original")));
      std::cout << "emplaced back original th3f" << std::endl;
    }


    ndGeoEff = 0.;
    for (unsigned int rep = 0; rep < N; rep++){
      
      // Randomize direction
      double phi = (fRandom->Rndm()-0.5)*2*TMath::Pi();
      
      // Randomize position Y, Z (keep X)
      float randomVtx[] = {-1e10, -1e10, -1e10};
      while(not isFV(randomVtx) ){
	randomVtx[0] = ndVtx[0];
	randomVtx[1] = 10*((fRandom->Rndm()-0.5)*2*100 + offset[1]);
	randomVtx[2] = 10*((fRandom->Rndm()-0.5)*2*(350-50)/2. + 200. + offset[2]);
      }

      Eigen::Affine3f randomDisplacement(Eigen::Translation3f(Eigen::Vector3f(randomVtx[0]-ndVtx[0],randomVtx[1]-ndVtx[1],randomVtx[2]-ndVtx[2])));
      
      Eigen::Affine3f tThere(Eigen::Translation3f(Eigen::Vector3f(-ndVtx[0],-ndVtx[1],-ndVtx[2])));
      Eigen::Affine3f r = Eigen::Affine3f(Eigen::AngleAxisf(phi, Eigen::Vector3f(avgBeamDir[0], avgBeamDir[1], avgBeamDir[2])));
      Eigen::Affine3f tBack(Eigen::Translation3f(Eigen::Vector3f(ndVtx[0],ndVtx[1],ndVtx[2])));
	
      Eigen::Transform<float,3,Eigen::Affine> m = randomDisplacement * tBack * r * tThere * ndPos ;

      if (isContained(m * hitSegPosOrig, hadEdep)) ndGeoEff += 1.;
      
      if (nd_Ehad < 0.) nd_Ehad =  totEinActiveVol( ndPos * hitSegPosOrig, hadEdep );

      if (evt.event() == evNumToDisplay and rep < nSamplesToDisplay) {
	displays.emplace_back(display3D(m * hitSegPosOrig, hadEdep, ("sample_"+std::to_string(rep)) ));
	std::cout << "emplaced back sample " << rep << std::endl;
      }

	
    }
      
    ndGeoEff /= N;

    
    fTree->Fill();
    std::cout << "filled ttree" << std::endl;
    return;
  }

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
  void TrueDepositsAna::endSubRun(const art::SubRun& sr){
  }

  void TrueDepositsAna::endJob()
  {
    //  fMetaTree->Fill();
  }

  DEFINE_ART_MODULE(TrueDepositsAna)

} // namespace dunemva

#endif // TrueDepositsAna_H
