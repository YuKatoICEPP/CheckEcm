#ifndef EcmCheckProcessor_h
#define EcmCheckProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TLorentzVector.h"
//#include <iostream>
#include <sstream>

class EcmCheckProcessor : public marlin::Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new EcmCheckProcessor ; }  
  
  EcmCheckProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

 protected:

  /** Input collection name.
   */
  std::string _colMCP ;
  std::string _outRootFile;
  
  int _nRun ;
  int _nEvt ;

  double _ecm ;
  int _nHello ;
  
  TFile *output;
  TH1D *hStatAnl;
  std::stringstream gCutName[20];
} ;

#endif



