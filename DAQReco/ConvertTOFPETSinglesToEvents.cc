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
  string singleFile= "/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RecoData/v1/RecoWithTracks/run"+runNumber+"_singles.root";
  TFile *tofhirFile=new TFile(singleFile.c_str(),"read");
  TTree *tofhirTree = (TTree*)tofhirFile->Get("data");
  int nTriggers = tofhirTree->Draw("channelID","channelID==255","goff");
  if( tofhirTree != NULL ) cout << "\n>>> got TOFHIR tree from file " << singleFile << " with " << tofhirTree->GetEntries() << " entries and " << nTriggers << " triggers" << endl;
  else exit(-1);

  TOFHIR TOF_(tofhirTree);

  //=====================Open fast Tracker input files=========================
  string trackingFileFast = "/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/Tracks/Run"+runNumber+"_CMSTiming_FastTriggerStream_converted.root";
  TFile *trackerFileFast = new TFile(trackingFileFast.c_str(),"READ");
  TTree *trackerTreeFast = (TTree*)trackerFileFast->Get("CMSTiming");
  if( trackerTreeFast != NULL ) cout << ">>> got track tree from file " << trackingFileFast << " with " << trackerTreeFast->GetEntries() << " entries" << endl;
  else exit(-1);

  TRACKER TRK_Fast(trackerTreeFast);

  //=====================Open Slow Tracker input files=========================
  string trackingFileSlow = "/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/Tracks/Run"+runNumber+"_CMSTiming_SlowTriggerStream_converted.root";
  TFile *trackerFileSlow = new TFile(trackingFileSlow.c_str(),"READ");
  TTree *trackerTreeSlow = (TTree*)trackerFileSlow->Get("CMSTiming");
  if( trackerTreeSlow != NULL ) cout << ">>> got track tree from file " << trackingFileSlow << " with " << trackerTreeSlow->GetEntries() << " entries" << endl;
  else exit(-1);

  SlowTrigTracker TRK_Slow(trackerTreeSlow);


  //=====================Open VME input files=========================
  //string vmeFile = "/home/daq/2020_02_cmstiming_BTL/VME/RecoData/RecoWithTracks/v2/RawDataSaver0CMSVMETiming_Run"+runNumber+"_0_Raw.root";
  //TFile *VMEInputFile = new TFile(vmeFile.c_str(),"READ");
  //TTree *VMETree = (TTree*)VMEInputFile->Get("pulse");
  //if( VMETree != NULL ) cout << ">>> got VME tree from file " << vmeFile << " with " << VMETree->GetEntries() << " entries " << endl;
  //else exit(-1);

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
  std::string OutName="outFile_"+runNumber+".root";
  TFile *outFile = new TFile(OutName.c_str(),"recreate");
  TTree *outTree = new TTree("data","data");
  outTree->SetAutoSave();


  if (TOF_.triggeredTofhirEv.size() > 0) {

    //=================================================================
    //Step 1             Find the begining of spill approximately ====
    //=================================================================

    TOFHIR::SpillInfo  Ev1stSpill1st = TOF_.get1stEvent1stSplill(0);
    TOFHIR::SpillInfo  EvLastSpill1st =  TOF_.getLastEvent1stSplill(Ev1stSpill1st.Index );
    TOFHIR::SpillInfo  Ev1stSpill2nd =  TOF_.get1stEvent2ndSplill(EvLastSpill1st.Index);
    TOFHIR::SpillInfo  EvLastSpill2nd =  TOF_.getLastEvent2ndSplill(Ev1stSpill2nd.Index);
    TOFHIR::SpillInfo  Ev1stSpill3rd =  TOF_.get1stEvent3rdSplill(EvLastSpill2nd.Index);

    TRACKER::SpillInfo  Ev1stSpill1st_trk = TRK_Fast.get1stEvent1stSplill(trackerTreeFast, 0);
    TRACKER::SpillInfo  EvLastSpill1st_trk =  TRK_Fast.getLastEvent1stSplill(trackerTreeFast, Ev1stSpill1st_trk.Index );
    TRACKER::SpillInfo  Ev1stSpill2nd_trk =  TRK_Fast.get1stEvent2ndSplill(trackerTreeFast, EvLastSpill1st_trk.Index);
    TRACKER::SpillInfo  EvLastSpill2nd_trk =  TRK_Fast.getLastEvent2ndSplill(trackerTreeFast, Ev1stSpill2nd_trk.Index);
    TRACKER::SpillInfo  Ev1stSpill3rd_trk =  TRK_Fast.get1stEvent3rdSplill(trackerTreeFast, EvLastSpill2nd_trk.Index);

    //=================================================================
    //Step 2             Find the begining of spill precisely     ====
    //=================================================================

    int Matched_1st=0;
    float Itr_1st=0;
    while (!Matched_1st && Itr_1st < 100){
      Matched_1st =TOF_.find1stSpill(TRK_Fast, tofhirTree, trackerTreeFast, Itr_1st,Ev1stSpill1st.Index,Ev1stSpill1st_trk.Index);
      Itr_1st++;
    }

    int Matched_2nd=0;
    float Itr_2nd=0;
    while (!Matched_2nd && Itr_2nd < 100){
      Matched_2nd =TOF_.find1stSpill(TRK_Fast, tofhirTree, trackerTreeFast, Itr_2nd,Ev1stSpill2nd.Index,Ev1stSpill2nd_trk.Index);
      Itr_2nd++;
    }
    if (Itr_2nd==100)  { TOF_.SPILLIndex.push_back(make_pair(-1,-1));}


    int Matched_3rd=0;
    float Itr_3rd=0;
    while (!Matched_3rd && Itr_3rd < 100){
      Matched_3rd =TOF_.find1stSpill(TRK_Fast, tofhirTree, trackerTreeFast, Itr_3rd,Ev1stSpill3rd.Index,Ev1stSpill3rd_trk.Index);
      Itr_3rd++;
    }
    if (Itr_3rd==100) { TOF_.SPILLIndex.push_back(make_pair(-1,-1));}

    //=================================================================
    //Step 3          Find the Match and fill the output file Tree ====
    //=================================================================

    TOF_.MatchAndFill(outTree, TRK_Fast, TRK_Slow, tofhirTree, trackerTreeFast, trackerTreeSlow ,TOF_.SPILLIndex,atoi(Run));
  }


  outFile->Write();
  outFile->Close();
  tofhirFile->Close();
  trackerFileFast->Close();
  trackerFileSlow->Close();
}
