// *****************************************************
// e+e- ------> Z H ------> q q inv
// Processor for checking Center Mass Energy
//                        --------  Y.Kato
// *****************************************************
#include "EcmCheckProcessor.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/Cluster.h>
#include <UTIL/LCTypedVector.h>
#include <EVENT/Track.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/ParticleID.h>
#include <marlin/Exceptions.h>
#include "UTIL/PIDHandler.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNtupleD.h"
#include "TVector3.h"
#include "TMath.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "Utilities.h"
//#include "HelicityAngle.cc"

using namespace lcio ;
using namespace marlin ;
using namespace std;

using namespace mylib;


EcmCheckProcessor aEcmCheckProcessor ;

EcmCheckProcessor::EcmCheckProcessor() : Processor("EcmCheckProcessor") {
  
  // modify processor description
  _description = "EcmCheckProcessor does whatever it does ..." ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::MCPARTICLE,
			   "InputMCParticlesCollection" , 
			   "Name of the MCParticle collection"  ,
			   _colMCP ,
			   std::string("MCParticle") ) ;
  
  registerProcessorParameter( "CenterOfMassEnergy",
			      "Center of mass energy"  ,
			      _ecm ,
			      double(500) ) ;
  
  registerProcessorParameter( "NHelloAnalysis",
			      "Interval of 'Hello, Analysis!'"  ,
			      _nHello ,
			      int(1000) ) ;

  registerOptionalParameter( "OutputRootFile",
			     "Name of output root file",
			     _outRootFile,
			     std::string("output.root") );
}

void EcmCheckProcessor::init() { 
  
  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  
  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
  stringstream out;
  out << _outRootFile << ends;
  output = new TFile(out.str().data(),"RECREATE");
  hStatAnl = 0;
}

void EcmCheckProcessor::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void EcmCheckProcessor::processEvent( LCEvent * evt ) { 
    
  // this gets called for every event 
  // usually the working horse ...
  _nEvt++;

  // some constants
  static const Double_t fEcm = _ecm;    // center-of-mass energy

  // initialize histo, ntuple...
  if (!hStatAnl) hStatAnl = new TH1D("hStatAnl", "Cut Table", 20, 0, 20);
  Double_t selid = -0.5;
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "No Cuts" << ends;

  TDirectory *last = gDirectory;
  gFile->cd("/");

  if ( (_nEvt%_nHello) == 0 )
    cerr << endl << "Hello, Analysis!" << " No: " << _nEvt << endl;

  // ------------------------------------------------------
  // -- New method to add TLorentzVector Branch in TTree --
  // ------------------------------------------------------

  Int_t nMCP;
  Int_t nOriginP = 0;
  Int_t qpdg[2];
  TLorentzVector qmc[2];
  TLorentzVector lrzZmc;
  TLorentzVector lrzHmc;
  TLorentzVector lrzISRmc[2];
  TLorentzVector lrzEcm = getLorentzEcm(fEcm);
  TLorentzVector lrzEcmmc;
  
  static TTree *hAnl = 0;
  if (!hAnl){
    hAnl = new TTree("hAnl","");
    
    hAnl->Branch("nmcp" ,        &nMCP );    
    hAnl->Branch("norigin" ,     &nOriginP );    
    hAnl->Branch("flvq1mc",      &qpdg[0] );
    hAnl->Branch("flvq2mc",      &qpdg[1] );
    hAnl->Branch("lrzq1mc",      &qmc[0] );
    hAnl->Branch("lrzq2mc",      &qmc[1] );
    hAnl->Branch("lrzZmc" ,      &lrzZmc );
    hAnl->Branch("lrzHmc" ,      &lrzHmc );
    hAnl->Branch("lrzISR1mc",      &lrzISRmc[0] );
    hAnl->Branch("lrzISR2mc",      &lrzISRmc[1] );
    hAnl->Branch("lrzEcm" ,      &lrzEcm );
    hAnl->Branch("lrzqqHisr12" ,      &lrzEcmmc );
  }
  // ------------------------------------------------------

  // ------------------------------------------------
  // -- read out the MCParticles information
  // ------------------------------------------------
  LCCollection *colMC = evt->getCollection(_colMCP);
  // get the truth information
  nMCP = colMC->getNumberOfElements();
  lrzZmc.SetPxPyPzE(0.,0.,0.,0.);
  for (Int_t i=0;i<nMCP;i++) {
    MCParticle *mcPart = dynamic_cast<MCParticle*>(colMC->getElementAt(i));
    Int_t pdg = mcPart->getPDG();
    Int_t ioverlay = mcPart->isOverlay()? 1 : 0;
    Int_t nparents = mcPart->getParents().size();
    Int_t motherpdg = 0;
    if (nparents > 0) {
      MCParticle *mother = mcPart->getParents()[0];
      motherpdg = mother->getPDG();
    }
    if (nparents == 0 && ioverlay == 0) nOriginP++;
    Int_t ndaughters = mcPart->getDaughters().size();
    Int_t daughterpdg = 0;
    if (ndaughters > 0) {
      MCParticle *daughter = mcPart->getDaughters()[0];
      daughterpdg = daughter->getPDG();
    }
    TLorentzVector lortz = getLorentzVector(mcPart);
    //Int_t ileft = mcPart->hasLeftDetector()? 1 : 0;
    if ((pdg > 0 && pdg < 10) && motherpdg == 0 && ioverlay == 0) {
      qpdg[0] = pdg;
      qmc[0]  = lortz;
      lrzZmc += lortz;
    }
    if ((pdg > -10 && pdg < 0) && motherpdg == 0 && ioverlay == 0) {
      qpdg[1] = pdg;
      qmc[1]  = lortz;
      lrzZmc += lortz;
    }
    if (pdg == 25 && motherpdg == 0 && ioverlay == 0) {
      lrzHmc = lortz;
    }
    if (i == 0 && motherpdg == 0) {
      lrzISRmc[0] = lortz;
    }
    if (i == 1 && motherpdg == 0) {
      lrzISRmc[1] = lortz;
    }
  }
  lrzEcmmc = qmc[0] + qmc[1] + lrzHmc + lrzISRmc[0] + lrzISRmc[1];
  
  // ------------------------------------------------------
  hAnl->Fill();

  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
  
  //  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
  //		       << "   in run:  " << evt->getRunNumber() 
  //		       << std::endl ;

  //  _nEvt ++ ;

  last->cd();
}



void EcmCheckProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EcmCheckProcessor::end(){ 

  cerr << "EcmCheckProcessor::end()  " << name() 
       << " processed " << _nEvt << " events in " << _nRun << " runs "
       << endl ;
  //  cerr << endl;
  cerr << "  =============" << endl;
  cerr << "   Cut Summary " << endl;
  cerr << "  =============" << endl;
  cerr << "   ll+4 Jet    " << endl;
  cerr << "  =============" << endl;
  cerr << endl
       << "  -----------------------------------------------------------" << endl
       << "   ID   No.Events    Cut Description                         " << endl
       << "  -----------------------------------------------------------" << endl;
  for (int id=0; id<20 && gCutName[id].str().data()[0]; id++) {
    cerr << "  " << setw( 3) << id
         << "  " << setw(10) << static_cast<int>(hStatAnl->GetBinContent(id+1))
         << "  : " << gCutName[id].str().data() << endl;
  }
  cerr << "  -----------------------------------------------------------" << endl;

  output->Write();
  output->Close();
}
