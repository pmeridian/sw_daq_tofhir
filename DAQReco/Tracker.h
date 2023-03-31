#ifndef INCLUDE_TRACKER_H_
#define INCLUDE_TRACKER_H_

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

#define DEBUG_TRK 1

using namespace std;

//================================================================================================================
// Tracker Slow Class
//================================================================================================================

class SlowTrigTracker{

public:
  explicit SlowTrigTracker (TTree *);

  struct FTBFSlowStreamPixelEvent {
    double xSlope;
    double ySlope;
    double xIntercept;
    double yIntercept;
    double chi2;
    int trigger;
    int runNumber;
    int nPlanes;
    Long64_t timestamp;
    Long64_t bco;
  };

  FTBFSlowStreamPixelEvent *trackerEventSlow = new FTBFSlowStreamPixelEvent;

private:
};

SlowTrigTracker::SlowTrigTracker(TTree * trackerTree){

  trackerTree->SetBranchAddress("event", trackerEventSlow);
}


//================================================================================================================
// Tracker Fast Class
//================================================================================================================

class TRACKER{

public:
  explicit TRACKER (TTree *);

  struct SpillInfo {
    int Index;
    Long64_t Time;
  };

  struct FTBFFastStreamPixelEvent {
    double xSlope;
    double ySlope;
    double xIntercept;
    double yIntercept;
    double chi2;
    double xResidBack;
    double yResidBack;
    double xErrDUT;
    double yErrDUT;
    double xErr04;
    double yErr04;
    double xErr05;
    double yErr05;
    double xErrPix0;
    double yErrPix0;
    double xResid04;
    double yResid04;
    double xResid05;
    double yResid05;
    double xResidPix0;
    double yResidPix0;
    int  trigger;
    int  runNumber;
    int  nPlanes;
    int  numPixels;
    int  numBackPlanes;
    int  numTracks;
    int  numClustersPix;
    int  numClustersStripsOdd;
    int  numClustersStripsEven;
    int  numStripsWith2Clusters;
    Long64_t timestamp;
    Long64_t bco;
  };

  FTBFFastStreamPixelEvent *trackerEvent = new FTBFFastStreamPixelEvent;

  SpillInfo get1stEvent1stSplill(TTree *,Int_t);
  SpillInfo getLastEvent1stSplill(TTree *,Int_t);

private:
};

TRACKER::TRACKER(TTree * trackerTree){
  trackerTree->SetBranchAddress("event", trackerEvent);
}


//================================================================================================================
// First Event in 1st spill  in Tracker data
//================================================================================================================

TRACKER::SpillInfo TRACKER::get1stEvent1stSplill(TTree * trackerTree, Int_t startEvent){

  TRACKER::SpillInfo spl;

  for(int k=startEvent; k<trackerTree->GetEntries()-1; k++) {

    trackerTree->GetEntry( startEvent );
    Long64_t TrackerBCOStarter = trackerEvent->bco;

    trackerTree->GetEntry( k+1 );
    Long64_t TrackerBCONextEvent = trackerEvent->bco;
    if (trackerEvent->bco== 4294967295) continue;

    trackerTree->GetEntry( k );
    Long64_t TrackerBCOThisEvent = trackerEvent->bco;
    if (trackerEvent->bco== 4294967295) continue;

    if (TrackerBCONextEvent < TrackerBCOThisEvent) continue;
    if (TrackerBCOThisEvent < TrackerBCOStarter) continue;

    //              if events are within 2ms of each other then we have entered the spill
    if ( (TrackerBCONextEvent - TrackerBCOThisEvent)*144e-9 < 0.002 ) {
      spl.Time = TrackerBCOThisEvent;
      spl.Index = k;
      if (DEBUG_TRK)        cout << "\nTRACKER: first event's index and time in spill 1: " << k << " : " << TrackerBCOThisEvent << "\n";
      break;
    }
  }
  return  spl;
}


//================================================================================================================
// Last Event in 1st spill  in Tracker data
//================================================================================================================

TRACKER::SpillInfo TRACKER::getLastEvent1stSplill(TTree * trackerTree,Int_t startEvent){

  TRACKER::SpillInfo spl;


  for(int k=startEvent; k<trackerTree->GetEntries(); k++) {

    trackerTree->GetEntry( startEvent );
    Long64_t TrackerBCOStarter = trackerEvent->bco;

    trackerTree->GetEntry( k+1 );
    Long64_t TrackerBCONextEvent = trackerEvent->bco;
    if (trackerEvent->bco== 4294967295) continue;
    trackerTree->GetEntry( k );
    Long64_t TrackerBCOThisEvent = trackerEvent->bco;
    if (trackerEvent->bco== 4294967295) continue;

    if ( k<trackerTree->GetEntries()-1 && TrackerBCONextEvent < TrackerBCOThisEvent) continue;
    if ( k<trackerTree->GetEntries()-1  && TrackerBCOThisEvent < TrackerBCOStarter) continue;
    
    //Check whether we have reached the end of spill
    if (k == trackerTree->GetEntries()-1 ){
      spl.Time = TrackerBCOThisEvent;
      spl.Index = k;
      if (DEBUG_TRK)    cout << "TRACKER: last event's index and time in spill 1:  " <<  k << " : " << TrackerBCOThisEvent << "\n";
      break;
    }

    //if events are more than 0.1s of each other then we have left the spill // to get ride of some pathological fearure and events
    if ( (TrackerBCONextEvent - TrackerBCOThisEvent)*144e-9 > 0.5) {
      spl.Time = TrackerBCOThisEvent;
      spl.Index = k;
      if (DEBUG_TRK)   cout << "TRACKER: last event's index and time in spill 1:  " <<  k << " : " << TrackerBCOThisEvent << "\n";
      break;
    }
  }
  return  spl;
}

#endif
