//
//  read_in.cxx
//  
//
//  Created by Isaac Mooney on 10/25/16.
//
//

#include <stdio.h>

// The majority of the jetfinding
// And correlation code is located in
// corrFunctions.hh
#include "readfuncs.hh"

// All reader and histogram settings
// Are located in corrParameters.hh
//#include "corrParameters.hh"

// ROOT is used for histograms and
// As a base for the TStarJetPico library
// ROOT Headers
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
#include "TRandom3.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TLatex.h"

// Make use of std::vector,
// std::string, IO and algorithm
// STL Headers
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <cstring>
#include <vector>
#include <string>
#include <limits.h>
#include <unistd.h>

// Data is read in by TStarJetPico
// Library, we convert to FastJet::PseudoJet
// TStarJetPico headers
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

// The analysis is run on FastJet::PseudoJets
// We make use of the jetfinding tools
// And the convenient FastJet::Selectors
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

// Used for year 7 tracking efficiency corrections,
// if they are being used
//#include "ktTrackEff.hh"

// -------------------------
// Command line arguments: ( Defaults
// Defined for debugging in main )
// [0]: input data: can be a single .root or a .txt or .list of root files

int main (int argc, const char** argv) {
    
    //Start a timer
    TStopwatch TimeKeeper;
    TimeKeeper.Start( );
    
    // Read in command line arguments
    // ------------------------------
    // Defaults
    std::string 	executable    = "./bin/read_in";                    // placeholder
    std::string		outputDir     = "tmp/";								// directory where everything will be saved
    std::string	 	inputFile	  = "/Volumes/Promise Pegasus/STAR/DATA/AddedGeantPythia/*";			// input file: can be .root, .txt, .list
    std::string 	chainName     = "JetTreeMc";                          // tree name in input file
    
    // Now check to see if we were given modifying arguments
    switch ( argc ) {
        case 1: // Default case
	  std::cout << "Using Default Settings" << std::endl;
            break;
        case 2: { // Custom case
	    std::cout << "Using Custom Settings" << std::endl;
            std::vector<std::string> arguments( argv+1, argv+argc );
            
            // Set non-default values
            // ----------------------
            
            // output and file names
            inputFile 		= arguments[0];
	    break;
        }
        default: { // Error: invalid custom settings
	  std::cerr << "Invalid number of command line arguments" << std::endl;
            return -1;
            break;
        }
            
    }
    
    //building input
    TChain* chain = new TChain( chainName.c_str() );
    
    // Check to see if the input is a .root file or a .txt
    bool inputIsRoot = io::HasEnding( inputFile.c_str(), ".root" );
    bool inputIsTxt  = io::HasEnding( inputFile.c_str(), ".txt"  );
    bool inputIsList = io::HasEnding( inputFile.c_str(), ".list" );
    
    // If its a recognized file type, build the chain
    // If its not recognized, exit
    if ( inputIsRoot )		 	{ chain->Add( inputFile.c_str() ); }
    else if ( inputIsTxt )  { chain = TStarJetPicoUtils::BuildChainFromFileList( inputFile.c_str() ); }
    else if ( inputIsList)  { chain = TStarJetPicoUtils::BuildChainFromFileList( inputFile.c_str() ); }
    else 		    { std::cerr << "data file is not recognized type: .root or .txt only." << std::endl; return -1; }
    
    // Intialize the reader and set the chain
    // All analysis parameters are located in
    // corrParameters.hh
    // --------------------------------------
    TStarJetPicoReader reader;
    io::InitReader( reader, chain, "pp", io::triggerAll, io::allEvents );
    
    // Data classes
    TStarJetVectorContainer<TStarJetVector>* container;
    TStarJetVector* sv; // TLorentzVector* would be sufficient
    TStarJetPicoEventHeader* header;
    TStarJetPicoEvent* event;
    TClonesArray* triggerObjs;
    //TClonesArray* towers;

    // Particle container
    std::vector<fastjet::PseudoJet> particles;
    
    // Now everything is set up
    // We can start the event loop
    // First, our counters
    int nEvents = 0;
    int nHardDijets = 0;
    int nMatchedHard = 0;

    //make a histogram!
    TFile f1("histo.root","RECREATE");
    // TCanvas *c1 = new TCanvas("c1","canvas", 400, 300);
    //Event level kinematics
    TH1F *mult = new TH1F("multiplicity","multiplicity",1000,0,50);
    TH1F *zvert = new TH1F("z vertex","z vertex",1000,-50,50);
    TH1F *zvertcut = new TH1F("z vertex w cut","z vertex w/ cut",1000,-50,50);
    
    //Track level kinematics
    TH1F *pxpart = new TH1F("px particle","px particle",10000,-2,2);
    TH1F *pypart = new TH1F("py particle","py particle",10000,-2,2);
    TH1F *pzpart = new TH1F("pz particle","pz particle",10000,-2,2);
    TH1F *Epart = new TH1F("E particle","E particle",10000,0,4);
    TH1F *phipart = new TH1F("phi particle","phi particle",10000,0,6.3);
    TH1F *rappart = new TH1F("rapidity particle","rapidity particle",10000,-1,1);
    TH1F *prappart = new TH1F("pseudorapidity particle","eta particle",10000,-1,1);
    TH1F *mpart = new TH1F("mass particle","mass particle",10000,-1000,1000);
    TH1F *trackpt = new TH1F("track pt","track p_{T}",1000,0,30);
    TH1F *trackptcut = new TH1F("track pt cut","track p_{T} cut",1000,0,30);
    TH1F *tracketa = new TH1F("track eta","track #eta",1000,-2,2);
    TH1F *tracketacut = new TH1F("track eta cut","track #eta cut",1000,-2,2);
    TH1F *trackphi = new TH1F("track phi","track #phi",1000,0,6.31);
    TH1F *trackphicut = new TH1F("track phi cut","track #phi cut",1000,0,6.31);

    //Towers
    TH1F *towerpt = new TH1F("tower pt","tower p_{T}",1000,0,30);
    TH1F *towerptcut = new TH1F("tower pt cut","tower p_{T} cut",1000,0,30);
    TH1F *towereta = new TH1F("tower eta","tower #eta",1000,-2,2);
    TH1F *toweretacut = new TH1F("tower eta cut","tower #eta cut",1000,-2,2);
    TH1F *towerphi = new TH1F("tower phi","tower #phi",1000,0,6.31);
    TH1F *towerphicut = new TH1F("tower phi cut","tower #phi cut",1000,0,6.31);


    try{
        while ( reader.NextEvent() ) {
            
            // Count the event
            nEvents++;
            
            // Print out reader status every 10 seconds
            reader.PrintStatus(10);
            
            // Get the event header and event
            event = reader.GetEvent();
            header = event->GetHeader();
	    
	    // Get the output container from the reader
	    container = reader.GetOutputContainer();
	    // Convert TStarJetVector to PseudoJet
	    io::ConvertTStarJetVector( container, particles, true );
	    
	    Int_t multiplicity = header->GetGReferenceMultiplicity();
	    mult->Fill(multiplicity);

	    // Find vertex Z bin
	    double vertexZ = header->GetPrimaryVertexZ();
	    zvert->Fill(vertexZ);
	    if (vertexZ >= -30 && vertexZ <= 30) {
	      zvertcut->Fill(vertexZ);
	    }
	    //Tracks & towers
	    TStarJetVectorContainer<TStarJetVector> *tracktow = reader.GetOutputContainer();
	    for (Int_t i = 0; i < tracktow->GetEntries(); ++ i) {
	      if (tracktow->Get(i)->GetCharge() != 0) {
		Float_t ptcalc = TMath::Sqrt(tracktow->Get(i)->px()*tracktow->Get(i)->px()+tracktow->Get(i)->py()*tracktow->Get(i)->py()); 
		trackpt->Fill(ptcalc);
		tracketa->Fill(tracktow->Get(i)->pseudorapidity());
		trackphi->Fill(tracktow->Get(i)->phi());
		if (vertexZ >= -30 && vertexZ <= 30) {
		  trackptcut->Fill(ptcalc);
		  tracketacut->Fill(tracktow->Get(i)->pseudorapidity());
		  trackphicut->Fill(tracktow->Get(i)->phi());
		}     
	      }
	      else {
		Float_t towerptcalc = TMath::Sqrt(tracktow->Get(i)->px()*tracktow->Get(i)->px()+tracktow->Get(i)->py()*tracktow->Get(i)->py());
		towerpt->Fill(towerptcalc);
		towereta->Fill(tracktow->Get(i)->pseudorapidity());
		towerphi->Fill(tracktow->Get(i)->phi());
		if (vertexZ >= -30 && vertexZ <= 30) {
		  towerptcut->Fill(towerptcalc);
		  toweretacut->Fill(tracktow->Get(i)->pseudorapidity());
		  towerphicut->Fill(tracktow->Get(i)->phi());
		}
	      }
	    }
	    //TList* tracks = reader.GetListOfSelectedTracks();
	    //tracks->Print();
	    //Towers
	    //TList* towerinfo = reader.GetListOfSelectedTowers(); 
	    //towerinfo->Print();
	    //std::vector<Float_t> kinematics(8,0);
	    /*
	    Float_t px = 0, py = 0, pz = 0, E = 0, phi = 0, rap = 0, prap = 0, m = 0;
	    
	    for ( int i = 0; i < particles.size(); ++i ) {
	      //kinematics[0] = particles[i].px();
	      px = particles[i].px();py = particles[i].py();pz = particles[i].pz();E = particles[i].E();
	      phi = particles[i].phi();rap = particles[i].rap();prap = particles[i].eta();m = particles[i].m();
	      
	      pxpart->Fill(px);pypart->Fill(py);pzpart->Fill(pz);Epart->Fill(E);phipart->Fill(phi);rappart->Fill(rap);prappart->Fill(prap);mpart->Fill(m);
	    }
	    */
	    /*
	    TIter next(tracks);
	    while (TObject *obj = next()) {
	      TStarJetPicoPrimaryTrack * track = (TStarJetPicoPrimaryTrack*)obj;
	      //Float_t test = track->GetPx();
	      trackpt->Fill(track->GetPt());
	      tracketa->Fill(track->GetEta());
	      trackphi->Fill(track->GetPhi());
	      if (vertexZ >= -30 && vertexZ <= 30) {
		trackptcut->Fill(track->GetPt());
		tracketacut->Fill(track->GetEta());
		trackphicut->Fill(track->GetPhi());
	      }
	    }
	    TIter next2(towerinfo);
	    while (TObject *obj = next2()) {
	      TStarJetPicoTower * tower = (TStarJetPicoTower*)obj;
	      //Float_t test = track->GetPx();
	      towerenergy->Fill(tower->GetEnergy());
	      towereta->Fill(tower->GetEta());
	      towerphi->Fill(tower->GetPhi());
	      if (vertexZ >= -30 && vertexZ <= 30) {
		towerenergycut->Fill(tower->GetEnergy());
		toweretacut->Fill(tower->GetEta());
		towerphicut->Fill(tower->GetPhi());
	      }
	    }
	    */
        }
	std::cout << "Number of events: " << nEvents << std::endl;
	mult->Write();
	zvert->Write();
	zvertcut->Write();
	trackpt->Write();
	trackptcut->Write();
	tracketa->Write();
	tracketacut->Write();
	trackphi->Write();
	trackphicut->Write();
	towerpt->Write();
	towerptcut->Write();
	towereta->Write();
	toweretacut->Write();
	towerphi->Write();
	towerphicut->Write();
	pxpart->Write();pypart->Write();pzpart->Write();Epart->Write();phipart->Write();rappart->Write();prappart->Write();mpart->Write();
	f1.Close();
    }catch ( std::exception& e) {
            std::cerr << "Caught " << e.what() << std::endl;
            return -1;
    }
    /*
    Float_t weight=1;
    if ( InPattern.Contains("Geant") ){
      TString currentfile = reader.GetInputChain()->GetCurrentFile()->GetName();
      weight=io::LookupXsec ( currentfile );
    }
    */
    return 0;
}
