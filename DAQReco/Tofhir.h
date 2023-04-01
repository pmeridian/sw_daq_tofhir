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
#include <utility>

#include "Tracker.h"

#define SEARCHWINDOWSIZE 5000 //in ps units
// from June2019 Test beam
//#define MATCH_THRESHOLD 0.000005 // in seconds
//#define ClockScaleFactor 1.04617
//Updated for Feb2020 Test beam
// #define MATCH_THRESHOLD 0.000004 // in seconds
// #define ClockScaleFactor 1.0461686
//Updated for Mar2023 Test beam
#define MATCH_THRESHOLD 0.000001 // in seconds
#define ClockScaleFactor 1.04163497

#define DEBUG_Tofhir 1
#define DEBUG_Deep_Tofhir 0

using namespace std;
//std::cout << std::setprecision(10);

class TOFHIR {
public:

  struct SpillInfo {
    int Index;
    Long64_t Time;
  };


  explicit TOFHIR (TTree *);

  const UInt_t triggerChannelID = 96;
  vector<SpillInfo> triggeredTofhirEv;

  Long64_t *  makeTimeTot(TTree * , int );

  double getScaleFactor(int);

  SpillInfo get1stEvent1stSplill(Int_t);
  SpillInfo getLastEvent1stSplill(Int_t);

  int find1stSpill(TRACKER ,TTree *, TTree *, int , int , int );
  void MatchAndFill(TTree * ,   TRACKER , TTree *, TTree *, vector< pair <int,int> > , int);

  Int_t           channelIdx[256];
  vector<Long64_t> *time;
  
  vector< pair <int,int> > SPILLIndex ;

private:


};


TOFHIR::TOFHIR(TTree * tofhirTree)
//Member initialization in constructors
{
  tofhirTree->SetBranchAddress("channelIdx",channelIdx);
  tofhirTree->SetBranchAddress("time",&time);
  // tofhirTree->SetBranchAddress("tot",&tot);
  // tofhirTree->SetBranchAddress("energy",&energy);
  // tofhirTree->SetBranchAddress("qfine",&qfine);
  // tofhirTree->SetBranchAddress("step1",&step1_);
  // tofhirTree->SetBranchAddress("step2",&step2_);


  SpillInfo test;
  for (Int_t ix=0;ix<tofhirTree->GetEntries(); ix++) {

    tofhirTree->GetEntry(ix);

    if (channelIdx[triggerChannelID]>=0) {

      test.Time= (*time)[channelIdx[triggerChannelID]];
      test.Index= ix;

      triggeredTofhirEv.push_back(test);
    }
  }
}


double TOFHIR::getScaleFactor(int run)
{
  std::cout << "==>> READING RESIDUAL SCALE FACTORS <<==" << std::endl;
  ifstream inFile("/eos/home-m/mtd/www/BTL/MTDTB_FNAL_Feb2020/triggerTimeCalibration/fitParameters.txt");
  std::string dummyLine;
  getline(inFile, dummyLine);

  int rnum;
  double intercept, slope;
  std::string chi2;
  while(true)
  {
    inFile >> rnum >> intercept >> slope >> chi2;
    if(inFile.eof())
      break;

    if(run == rnum && chi2!="inf" && atof(chi2.c_str()) <= 15)
    {
      std::cout << "SLOPE: " << slope << std::endl;
      return slope;
    }
  }
  return 0.;
}

//================================================================================================================
// First Event in 1st spill  in TOFHIR data
//================================================================================================================

TOFHIR::SpillInfo TOFHIR::get1stEvent1stSplill(Int_t startEvent){

  TOFHIR::SpillInfo spl;
  for (Int_t i=startEvent;i<triggeredTofhirEv.size(); i++) {

    Long64_t timeNextToNextEvent = triggeredTofhirEv[i+2].Time;
    Long64_t timeNextEvent = triggeredTofhirEv[i+1].Time;
    Long64_t timeThisEvent = triggeredTofhirEv[i].Time;

    //if events are within 2ms of each other then we have entered the spill
    if ( abs(timeNextEvent-timeThisEvent) *1e-12 < 2e-3 && abs(timeNextToNextEvent-timeNextEvent) * 1e-12 < 2e-2  ) {
      spl.Index = i;
      spl.Time = timeThisEvent;
      if (DEBUG_Tofhir) cout << "\nTOFHIR: first event's index and time in spill 1: " << i << " : " << timeThisEvent << "\n";
      break;
    }
  }
  return  spl;
}

//================================================================================================================
//  Last Event in 1st spill  in TOFHIR data
//================================================================================================================

TOFHIR::SpillInfo TOFHIR::getLastEvent1stSplill(Int_t startEvent){

  TOFHIR::SpillInfo spl;
  for (Int_t i=startEvent;i<triggeredTofhirEv.size(); i++) {

    Long64_t timeNextToNextEvent = triggeredTofhirEv[i+2].Time;
    Long64_t timeNextEvent = triggeredTofhirEv[i+1].Time;
    Long64_t timeThisEvent = triggeredTofhirEv[i].Time;

    //Check whether we have reached the end of spill
    if (i == triggeredTofhirEv.size()-1 ){
      cout << i<<" \t\t we reached to the end of spill \n";
      spl.Index = i;
      spl.Time = timeThisEvent;
      cout << "TOFHIR: last event's index and time in spill 1:  " << i << " : " << timeThisEvent << "\n";
      break;
    }

    //if next event is more than  10ms and the next-to-next is more than 50 ms then we have left the spill
    if ( abs(timeNextEvent-timeThisEvent) * 1e-12 > 1e-2  && abs(timeNextToNextEvent-timeNextEvent) * 1e-12 > 5e-2 ) {
      spl.Index = i;
      spl.Time = timeThisEvent;
      if (DEBUG_Tofhir) cout << "TOFHIR: last event's index and time in spill 1:  " << i << " : " << timeThisEvent << "\n";
      break;
    }
  }
  return  spl;
}


//*********************************************************************************************************************
//*********************************************************************************************************************
//                                         Finding Begining of Spill precisely
//*********************************************************************************************************************
//*********************************************************************************************************************

int TOFHIR::find1stSpill(TRACKER TRK_Fast,TTree *tofhirTree, TTree *trackerTree, int Itr, int firstIndexTOF, int firstIndexTRK){

  int n_tof=Itr%10-1;
  int n_trg=Itr/10-1;

  int SpilBegin=firstIndexTOF+n_tof;
  if (SpilBegin < 0) SpilBegin=0;
  int PreviousMatchIndex = firstIndexTRK+n_trg;
  if (PreviousMatchIndex < 0) PreviousMatchIndex=0;

  Long64_t  TofhirTimeZero=triggeredTofhirEv[SpilBegin].Time;
  trackerTree->GetEntry(PreviousMatchIndex);
  Long64_t   TrackerTimeZero = TRK_Fast.trackerEvent->bco;

  int firstOfTrack=PreviousMatchIndex;
  int totalNumEveMatched=1;

  //======================== Look for coincidences =================================================================
  for (Int_t q=SpilBegin;q<SpilBegin+50; q++) {

    if (q-SpilBegin-totalNumEveMatched > 20) {
      return 0;
    }

    double TOFHIRTriggerTimestamp = (triggeredTofhirEv[q].Time - TofhirTimeZero)*1e-12;

    for(int k=PreviousMatchIndex; k<PreviousMatchIndex+50; k++) {

      trackerTree->GetEntry( k );
      Long64_t TrackerBCOThisEvent = TRK_Fast.trackerEvent->bco;
      if (TRK_Fast.trackerEvent->bco== 4294967295) continue;
      double trackerTime = (TrackerBCOThisEvent - TrackerTimeZero) * 144e-9 * ClockScaleFactor;
      if (fabs (TOFHIRTriggerTimestamp - trackerTime ) < MATCH_THRESHOLD  ){
        totalNumEveMatched++;
        PreviousMatchIndex = k;
        break;
      }
    }
  }

  SPILLIndex.push_back(make_pair(SpilBegin,firstOfTrack));
  return 1;
}

//*********************************************************************************************************************
//*********************************************************************************************************************
//                                         Matching Function
//*********************************************************************************************************************
//*********************************************************************************************************************


void TOFHIR::MatchAndFill(TTree * outTree,   TRACKER TRK_Fast, TTree *tofhirTree, TTree *trackerTree, vector< pair <int,int> > SpillInterval, int run){

  cout<<"\nIndex of first spill in TOFHIR= "<< SpillInterval[0].first << "\t and first spill in TRACKER = "<< SpillInterval[0].second<<"\n";

  TGraph* g_temp = new TGraph();

  int SpilBegin=SpillInterval[0].first;
  int PreviousMatchIndex= SpillInterval[0].second;

  Long64_t  TofhirTimeZero=triggeredTofhirEv[SpilBegin].Time;
  trackerTree->GetEntry(PreviousMatchIndex);
  Long64_t   TrackerTimeZero = TRK_Fast.trackerEvent->bco;

  //==============get scale factor from text file============
  double ResidualClockScaleFactor = ClockScaleFactor;


  //==============Other variables================
  int totalNumEve=triggeredTofhirEv.size();
  int totalNumEveMatched=1;
  int bothTriggeredFired=0;
  int lastMatchedIndex=0;
  float TimeDiff=-10;
  int matchEff=0;

  int iev_TOFHIR;
  int iev_TRACKER;
  double t_TOFHIR;
  unsigned long t_TRACKER_bco;
  float xIntercept=-9999;
  float yIntercept=-9999;
  float xSlope=-9999;
  float ySlope=-9999;
  float x_dut=-9999;
  float y_dut=-9999;
  float chi2=-9999;
  int nplanes=-1;

  
  outTree->Branch("iev_TOFHIR",&iev_TOFHIR,"iev_TOFHIR/I");
  outTree->Branch("iev_TRACKER",&iev_TRACKER,"iev_TRACKER/I");
  outTree->Branch("t_TOFHIR",&t_TOFHIR,"t_TOFHIR/D");
  outTree->Branch("t_TRACKER_bco",&t_TRACKER_bco,"t_TRACKER_bco/l");
  outTree->Branch("xIntercept", &xIntercept, "xIntercept/F");
  outTree->Branch("yIntercept", &yIntercept, "yIntercept/F");
  outTree->Branch("xSlope", &xSlope, "xSlope/F");
  outTree->Branch("ySlope", &ySlope, "ySlope/F");
  outTree->Branch("x_dut", &x_dut, "x_dut/F");
  outTree->Branch("y_dut", &y_dut, "y_dut/F");
  outTree->Branch("chi2", &chi2, "chi2/F");
  outTree->Branch("nplanes", &nplanes, "nplanes/I");

  //    //================================================================================================================
  //    //================================================================================================================
  //    //======================== Look for coincidences =================================================================
  //    //================================================================================================================
  //    //================================================================================================================

  PreviousMatchIndex--;


  //    for (Int_t q=SpilBegin;q<triggeredTofhirEv.size(); q++) {
  int lastIndex=0;
  for (Int_t q=SpilBegin;q<triggeredTofhirEv.size(); q++) {
    iev_TOFHIR=-1;
    iev_TRACKER=-1;
    xIntercept = -9999;
    yIntercept = -9999;
    xSlope = -9999;
    ySlope = -9999;
    x_dut = -9999;
    y_dut = -9999;
    chi2 = -9999;
    nplanes = -1;

    //Fill empty informations in events without EXT trigger 
    for(Int_t tofhirIndex=lastIndex;tofhirIndex<triggeredTofhirEv[q].Index;++tofhirIndex)
	outTree->Fill();

    lastIndex=triggeredTofhirEv[q].Index+1;
    iev_TOFHIR=q;

    // tofhirTree->GetEntry(triggeredTofhirEv[q].Index);
    TimeDiff = -10;
    int NMatchedTracks = 0;
    matchEff=0;

        
    double TOFHIRTriggerTimestamp = (triggeredTofhirEv[q].Time - TofhirTimeZero)*1e-12;
    t_TOFHIR=TOFHIRTriggerTimestamp; //save time in s
    t_TRACKER_bco=-1;

    //================================================================================================================
    if (DEBUG_Deep_Tofhir) cout<<"\n\nq= "<<q<<"   time is= "<<triggeredTofhirEv[q].Time*1e-12 <<"  dif=" <<(triggeredTofhirEv[q].Time-TofhirTimeZero)*1e-12<<"\n";
					 
    for(int k=PreviousMatchIndex+1; k<trackerTree->GetEntries(); k++) {
					   
      // trackerTree->GetEntry( k + 1);
      // Long64_t TrackerBCONextEvent = TRK_Fast.trackerEvent->bco;

      trackerTree->GetEntry( k );
      Long64_t TrackerBCOThisEvent = TRK_Fast.trackerEvent->bco;

      if (TRK_Fast.trackerEvent->bco== 4294967295) continue;

      //double trackerTime = (TrackerBCOThisEvent - TrackerTimeZero) * 144e-9 * ClockScaleFactor * (1.+1.76719e-06);
      double trackerTime = (TrackerBCOThisEvent - TrackerTimeZero) * 144e-9 * ClockScaleFactor; 

      if (DEBUG_Deep_Tofhir) cout<<"\nk = "<<k<<"  time is= "<<TrackerBCOThisEvent<< " - "<<TrackerTimeZero<<" =>   dif= "<<trackerTime<< " xIn,yIn= "<< TRK_Fast.trackerEvent->xIntercept << "," << TRK_Fast.trackerEvent->yIntercept << "\n";
      if (DEBUG_Deep_Tofhir) cout<<"diff(TOFHIR-tracker) = " << fabs (TOFHIRTriggerTimestamp - trackerTime ) << "   - MATCH_THRESHOLD: " << MATCH_THRESHOLD << "\n";

      //================================================================================================================
      if (fabs ((TOFHIRTriggerTimestamp - trackerTime) ) < MATCH_THRESHOLD  ){
        if(q%10000==0  || DEBUG_Deep_Tofhir)    cout <<">>> Event with index of Tofhir= "<<q <<"  MATCHED w/ event w/ index of Tracker= " << k << " : " << TOFHIRTriggerTimestamp << "   | "<< trackerTime<<"\n";

        if (DEBUG_Tofhir)
          g_temp -> SetPoint(g_temp->GetN(),TOFHIRTriggerTimestamp,trackerTime-TOFHIRTriggerTimestamp);

        // *(1+1.76719e-06)
	
        TimeDiff= TOFHIRTriggerTimestamp - trackerTime;
        totalNumEveMatched++;
        NMatchedTracks++;
        PreviousMatchIndex = k;


        // //================================================================================================================
        // // Fill variables of the FastTriggerMode of Tracker
        // //================================================================================================================

	// //using INFO from NEXT EVENT... both 2020 & 2023 ????
        trackerTree->GetEntry( k + 1);

	iev_TRACKER=k+1;
	t_TRACKER_bco=TrackerBCOThisEvent;
        xIntercept=TRK_Fast.trackerEvent->xIntercept * 1e-3;
        yIntercept=TRK_Fast.trackerEvent->yIntercept * 1e-3;
        xSlope=TRK_Fast.trackerEvent->xSlope * 1e-3;
        ySlope=TRK_Fast.trackerEvent->ySlope * 1e-3;
        x_dut= (TRK_Fast.trackerEvent->xIntercept + TRK_Fast.trackerEvent->xSlope * 2.0e5) * 1e-3;
        y_dut= (TRK_Fast.trackerEvent->yIntercept + TRK_Fast.trackerEvent->ySlope * 2.0e5) * 1e-3;
        chi2 = TRK_Fast.trackerEvent->chi2;
        nplanes = TRK_Fast.trackerEvent->nPlanes;
	break;
      }// Check if Tofhir evet matched the trigger
      else if (TOFHIRTriggerTimestamp - trackerTime < -0.001 ) {
        break;
      }

    }

    if (DEBUG_Deep_Tofhir) std::cout << xIntercept << "," << yIntercept << std::endl;
    outTree->Fill();
  }

  iev_TOFHIR=-1;
  iev_TRACKER=-1;
  xIntercept = -9999;
  yIntercept = -9999;
  xSlope = -9999;
  ySlope = -9999;
  x_dut = -9999;
  y_dut = -9999;
  chi2 = -9999;
  nplanes = -1;
  //Fill empty info for last chunck of events
  for(Int_t tofhirIndex=lastIndex;tofhirIndex<tofhirTree->GetEntries();++tofhirIndex)
    outTree->Fill();

  if (DEBUG_Tofhir)
  {
    TFile* outFile = new TFile("temp.root","RECREATE");
    g_temp -> Write("g_temp");
    outFile -> Close();
  }

  cout<<"\n=>  Matching efficiency is "<<totalNumEveMatched <<"/" <<totalNumEve <<" = "<< 1.0*totalNumEveMatched/totalNumEve<<"\n";
  cout<<"======================================================================================================\n\n";
}
