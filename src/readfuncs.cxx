#include "readfuncs.hh"

namespace io {
    
    // -------------------------
    // IO/OS Manip functionality
    // -------------------------
    
    
    // Used to understand which format of input file is being used
    // ( .root file, .txt, .list, etc )
    // ---------------------------------------------------------------------
    bool HasEnding (std::string const &full_string, std::string const &ending) {
        if (full_string.length() >= ending.length()) {
            return (0 == full_string.compare (full_string.length() - ending.length(), ending.length(), ending) );
        } else {
            return false;
        }
    }


  // Fills our working container after converting TStarJetVectors into PseudoJets
  // Also makes sure to empty the containers if they are full
  //
  // ------------------------------------------------------------------------------
  void ConvertTStarJetVector( TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet> & particles, bool ClearVector ) {
    // Empty the container
    // if called for
    if ( ClearVector )
      particles.clear();

    // Transform TStarJetVectors into (FastJet) PseudoJets
    // ---------------------------------------------------
    TStarJetVector* sv;
    for ( int i=0; i < container->GetEntries() ; ++i ){
      sv = container->Get(i);

      fastjet::PseudoJet tmpPJ = fastjet::PseudoJet( *sv );
      tmpPJ.set_user_index( sv->GetCharge() );
      particles.push_back( tmpPJ );

    }
  }
  
    // Used to initialized the reader - will set the event cuts,
    // Tower cuts, track cuts and hadronic correction
    // ---------------------------------------------------------------------
    void InitReader( TStarJetPicoReader & reader, TChain* chain, std::string collisionType, std::string triggerString, int nEvents ) {
        
        // set the chain
        reader.SetInputChain( chain );
        
        // Initialize the reader
        reader.Init( nEvents ); //runs through all events with -1


	// Tracks cuts
	TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
	towerCuts->AddBadTowers( "src/emptyList.txt" );
    }
    
}

double LookupXsec( TString filename ){
  // Some data for geant
  // -------------------
  //cross-sections for simulated GEANT data sample
  // From Renee.
  // also available via
  // http://people.physics.tamu.edu/sakuma/star/jets/c101121_event_selection/s0150_mclist_001/w\
  eb.php
    // Double_t MinbXsec=28.12;
    // Double_t Xsec[12];
    // Xsec[0]=28.11;//2
    // Xsec[1]=1.287;//3
    // Xsec[2]=0.3117;//4
    // Xsec[3]=0.1360;//5
    // Xsec[4]=0.02305;//7
    // Xsec[5]=0.005494;//9
    // Xsec[6]=0.002228;//11
    // Xsec[7]=0.0003895;//15
    // Xsec[8]=0.00001016;//25
    // Xsec[9]=0.0000005010;//35
    // Xsec[10]=0.0000000283;//45
    // Xsec[11]=0.000000001443;//55

    static const Double_t MinbXsec=28.12;
  // static const Double_t Xsec[12] = {
  //   28.11, // 2-3
  //   1.287, // 3-4
  //   0.3117, // 4-5
  //   0.1360, // 5-7
  //   0.02305, // 7-9
  //   0.005494, // 9-11
  //   0.002228, // 11-15
  //   0.0003895, // 15-25
  //   0.00001016, // 25-35
  //   0.0000005010, // 35-45
  //   0.0000000283, // 45-55
  //   0.000000001443 // 55-65
  // };

  static const Double_t Xsec[12] = {
    1.0,        // Placeholder for 2-3
    1.30E+09, // 3-4
    3.15E+08, // 4-5
    1.37E+08, // 5-7
    2.30E+07, // 7-9
    5.53E+06, // 9-11
    2.22E+06, // 11-15
    3.90E+05, // 15-25
    1.02E+04, // 25-35
    5.01E+02, // 35-45
    2.86E+01, // 45-55
    1.46E+00 // 55-65
  };

  static const Double_t Nmc[12] = {
    1, // 2-3
    672518, // 3-4
    672447, // 4-5
    393498, // 5-7
    417659, // 7-9
    412652, // 9-11
    419030, // 11-15
    396744, // 15-25
    399919, // 25-35
    119995, // 35-45
    117999, // 45-55
    119999 // 55-65
  };

  Double_t w[12];
  for ( int i=0; i<12 ; ++i ){
    w[i] = Xsec[i] / Nmc[i];
    // w[i] = Nmc[i] / Xsec[i] ;
  }

  // static const Double_t w[12] = {
  //   1, // Placeholder
  //   1.90E+03,
  //   6.30E+02,
  //   3.43E+02,
  //   5.49E+01,
  //   1.33E+01,
  //   5.30E+00,
  //   9.81E-01,
  //   2.56E-02,
  //   4.56E-03,
  //   2.43E-04,
  //   1.20E-05
  // };

  if ( filename.Contains("picoDst_3_4") ) return w[1];
  if ( filename.Contains("picoDst_4_5") ) return w[2];
  if ( filename.Contains("picoDst_5_7") ) return w[3];
  if ( filename.Contains("picoDst_7_9") ) return w[4];
  if ( filename.Contains("picoDst_9_11") ) return w[5];
  if ( filename.Contains("picoDst_11_15") ) return w[6];
  if ( filename.Contains("picoDst_15_25") ) return w[7];
  if ( filename.Contains("picoDst_25_35") ) return w[8];
  if ( filename.Contains("picoDst_35_45") ) return w[9];
  if ( filename.Contains("picoDst_45_55") ) return w[10];
  if ( filename.Contains("picoDst_55_65") ) return w[11];

  return 1;

}

    
