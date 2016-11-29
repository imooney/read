// STL
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <string>
#include <limits.h>
#include <unistd.h>

// fastjet 3
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/FunctionOfPseudoJet.hh"

// ROOT
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TObjArray.h"
#include "TString.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "TBranch.h"
#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TStopwatch.h"

// TStarJetPico
#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"

#ifndef READ_IN_HH
#define READ_IN_HH

namespace io {

  //Event cuts:
  const double vertexZcut = 30;

  
    // default nEvents for all events in chain
    const int		allEvents = -1;									// selects all events in the chain
    
    // TriggerStrings
    const std::string triggerAll = "All";				// accept all events
    
    // Helper to build the TChain, used to decide which input format
    bool HasEnding (std::string const &full_string, std::string const &ending);

    // Converts TStarJetPicoVectors into PseudoJets
    void ConvertTStarJetVector( TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet> & particles, bool ClearVector = true );
  
    // Initializes the TStarJetPicoReader, so we dont have
    // To have all that code hanging around in the analysis
    // Collision Type is 'AuAu' or 'pp'
    void InitReader( TStarJetPicoReader & reader, TChain* chain, std::string collisionType, std::string triggerString, int nEvents );

    double LookupXsec( TString filename );
}



#endif
