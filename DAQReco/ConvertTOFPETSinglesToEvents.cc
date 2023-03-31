//================================================================================================================
// Done by abdollah.mohammadi@cern.ch and Si.Xie@cern.ch
//================================================================================================================
// include std libraries
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>
#include <string.h>
#include <sstream>
// include ROOT libraries
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TFolder.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TProfile.h"
// Other headrfiles are in Tofhir Class
#include "Tofhir.h"

using namespace std;

int main(int argc, char* argv[]){

  if(argc < 1) {
    cerr << "Please give the run number " <<endl;
    return -1;
  }

  const char *Run = argv[1];
  // const char *inputFileName = argv[2];
  // const char *inputTrackerFileName = argv[3];
  // const char *outFileName   = argv[4];

  std::string runNumber(Run);

  //=====================Open TOFHIR input files=========================
  //  string singleFile= "/eos/uscms/store/group/cmstestbeam/2023_03_cmstiming_BTL/TOFHIR/RecoData/run"+runNumber+"_e.root";
  //string singleFile= "/uscms_data/d2/meridian/MTD/FNALTB_2023/BTLReco/reco/run"+runNumber+"_s.root";
  string singleFile= "/Users/meridian/scratch/FNALTB_2023/data/run"+runNumber+"_e.root";
  TFile *tofhirFile=new TFile(singleFile.c_str(),"read");
  TTree *tofhirTree = (TTree*)tofhirFile->Get("data");
  int nTriggers = tofhirTree->Draw("channelID","channelID==96","goff");
  if( tofhirTree != NULL ) cout << "\n>>> got TOFHIR tree from file " << singleFile << " with " << tofhirTree->GetEntries() << " entries and " << nTriggers << " triggers" << endl;
  else exit(-1);

  TOFHIR TOF_(tofhirTree);

  //=====================Open fast Tracker input files=========================
  string trackingFileFast = "/Users/meridian/scratch/FNALTB_2023/data/Run"+runNumber+"_CMSTiming_FastTriggerStream_converted.root";
  TFile *trackerFileFast = new TFile(trackingFileFast.c_str(),"READ");
  TTree *trackerTreeFast = (TTree*)trackerFileFast->Get("CMSTiming");
  if( trackerTreeFast != NULL ) cout << ">>> got track tree from file " << trackingFileFast << " with " << trackerTreeFast->GetEntries() << " entries" << endl;
  else exit(-1);

  TRACKER TRK_Fast(trackerTreeFast);

  if (TOF_.triggeredTofhirEv.size() == 0) {
    cout << "Warning: TOFHIR has 0 triggered events \n";
  } else {
    //=====================Files Sanity Check=========================
    // Usually the Tracker files have 20% more events than TOFHIR files
    if (trackerTreeFast->GetEntries()/TOF_.triggeredTofhirEv.size() < 1.0 || trackerTreeFast->GetEntries()/TOF_.triggeredTofhirEv.size() > 1.4 )
    {
      cout <<"Tracker has much less or much more events than TOFHIR !!!  trk="<<  trackerTreeFast->GetEntries() <<"   v.s. tof= "<<TOF_.triggeredTofhirEv.size()<< "\n";
      //exit(-1);
    }
  }


  //=====================RECREATE outputFile=========================
  std::string OutName="outFile_"+runNumber+"_match.root";
  TFile *outFile = new TFile(OutName.c_str(),"recreate");
  TTree *outTree = new TTree("data","data");
  outTree->SetAutoSave();


  if (TOF_.triggeredTofhirEv.size() > 0) {

    //=================================================================
    //Step 1             Find the begining of spill approximately ====
    //=================================================================

    TOFHIR::SpillInfo  Ev1stSpill1st = TOF_.get1stEvent1stSplill(0);
    TOFHIR::SpillInfo  EvLastSpill1st =  TOF_.getLastEvent1stSplill(Ev1stSpill1st.Index );

    TRACKER::SpillInfo  Ev1stSpill1st_trk = TRK_Fast.get1stEvent1stSplill(trackerTreeFast, 0);
    TRACKER::SpillInfo  EvLastSpill1st_trk =  TRK_Fast.getLastEvent1stSplill(trackerTreeFast, Ev1stSpill1st_trk.Index );

    //=================================================================
    //Step 2             Find the begining of spill precisely     ====
    //=================================================================

    int Matched_1st=0;
    float Itr_1st=0;
    while (!Matched_1st && Itr_1st < 100){
      Matched_1st =TOF_.find1stSpill(TRK_Fast, tofhirTree, trackerTreeFast, Itr_1st,Ev1stSpill1st.Index,Ev1stSpill1st_trk.Index);
      Itr_1st++;
    }


    //=================================================================
    //Step 3          Find the Match and fill the output file Tree ====
    //=================================================================

    TOF_.MatchAndFill(outTree, TRK_Fast,  tofhirTree, trackerTreeFast, TOF_.SPILLIndex,atoi(Run));
  }


  outFile->Write();
  outFile->Close();
  tofhirFile->Close();
  trackerFileFast->Close();
}
