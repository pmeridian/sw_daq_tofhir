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

#define MAX_TOFPET_CHANNEL 400
#define SEARCHWINDOWSIZE 5000 //in ps units
#define MATCH_THRESHOLD 0.000001 // in seconds
#define DEBUG 0

using namespace std;

 struct FTBFFastStreamPixelEvent {
    double xSlope;
    double ySlope;
    double xIntercept;
    double yIntercept;
    double chi2;
    int trigger;
    int runNumber;
    int nPlanes;
    int numPixels;
    int numBackPlanes;
    Long64_t timestamp;
    Long64_t bco;
  };
  
int main(int argc, char* argv[]){
   if(argc < 3) {
    cerr << "Please give 3 arguments " << "inputTOFHIRFileName " << " inputTrackerFileName " << "outputFileName" <<endl;
    return -1;
  }


  //=====================Open Tracker input files=========================
  const char *inputTrackerFileName = argv[2];
  string trackingFile = "/home/daq/2019_04_April_CMSTiming/Tracks/Run14592_CMSTiming_FastTriggerStream_converted.root";
  TFile *trackerFile = new TFile(inputTrackerFileName,"READ");
  TTree *trackerTree = (TTree*)trackerFile->Get("CMSTiming");
  FTBFFastStreamPixelEvent *trackerEvent = new FTBFFastStreamPixelEvent;
  trackerTree->SetBranchAddress("event", trackerEvent);




  //=====================Open TOFHIR input files=========================    
  const char *inputFileName = argv[1];
  const char *outFileName   = argv[3];
  TFile *AnaFile=new TFile(inputFileName,"read");
  TTree *anatree = (TTree*)AnaFile->Get("data");
  
  TFile *outFile = new TFile(outFileName,"recreate");
  TTree *outTree = new TTree("data","data");
  outTree->SetAutoSave();

  //=====================Config Variables========================= 
  const UInt_t triggerChannelID = 384;


  //=====================Initialize input tree variables========================= 
  UInt_t channelID;
  Long64_t time;
  float tot;
  float xIntercept;
  float yIntercept;
  float xSlope;
  float ySlope;
  float x_dut;
  float y_dut;
  float chi2;
  int ntracks;
  int nplanes;

  channelID=-9999;
  time=-9999;
  tot=-9999;
  xIntercept=-9999;
  yIntercept=-9999;
  xSlope=-9999;
  ySlope=-9999;
  x_dut=-9999;
  y_dut=-9999;
  chi2 = -9999;
  ntracks = -1;
  nplanes = -1;

  //==============set Branch addresses for all the input variables================  
  anatree->SetBranchAddress("channelID",&channelID);
  anatree->SetBranchAddress("time",&time);
  anatree->SetBranchAddress("tot",&tot);

  //=====================Initialize output tree variables=========================
  float step1;
  float step2;
  Long64_t chTime[400];
  float chtot[400];
  Int_t event;

 //initialize for event 1
  event=1;
  for(int k=0;k<400;k++){
    chTime[k]=-9999;
    chtot[k]=-9999;
  }
  
  //==============set Branch addresses for all the output variables================  
  outTree->Branch("event",&event,"event/I");
  outTree->Branch("step1", &step1, "step1/F");
  outTree->Branch("step2", &step2, "step2/F");
  outTree->Branch("chTime",&chTime,"chTime[400]/L");
  outTree->Branch("chtot",&chtot,"chtot[400]/F");
  outTree->Branch("xIntercept", &xIntercept, "xIntercept/F");
  outTree->Branch("yIntercept", &yIntercept, "yIntercept/F");
  outTree->Branch("xSlope", &xSlope, "xSlope/F");
  outTree->Branch("ySlope", &ySlope, "ySlope/F");
  outTree->Branch("x_dut", &x_dut, "x_dut/F");
  outTree->Branch("y_dut", &y_dut, "y_dut/F");  
  outTree->Branch("chi2", &chi2, "chi2/F");
  outTree->Branch("ntracks", &ntracks, "ntracks/I");
  outTree->Branch("nplanes", &nplanes, "nplanes/I");
		


 //========================Initialize local variables=================== 
  Int_t nentries_ana = (Int_t)anatree->GetEntries();

  //======================== Find the start of spill in TOFHIR events =======================================
  int index_firstEventInSpill = -1; 
  Long64_t TOFHIRTriggerTimeOfFirstEventInSpill = -1;
  Long64_t starttimer = -1;

  vector<Long64_t > TOFHIR_TriggerTimestamps;
  vector<Int_t> TOFHIR_TriggerIndices;
//  Abdollah: Adding this section cause to sizable delay 
  for (Int_t i=0;i<nentries_ana; i++) {
//  for (Int_t i=0;i<50000; i++) {

    anatree->GetEntry(i);    
    if (channelID == triggerChannelID) {
     TOFHIR_TriggerTimestamps.push_back(time);
      TOFHIR_TriggerIndices.push_back(i);
    }
  }



 

  for (Int_t i=0;i<TOFHIR_TriggerIndices.size()-1; i++) {
 
    Long64_t timeNextEvent = TOFHIR_TriggerTimestamps[i+1];
    Long64_t timeThisEvent = TOFHIR_TriggerTimestamps[i];
    if (i==0) starttimer = TOFHIR_TriggerTimestamps[i];



    // cout << "entry " << i << " : " << timeNextEvent << " - " << timeThisEvent << " = " << timeNextEvent-timeThisEvent << " passes? " 
    // 	 << bool( abs(timeNextEvent-timeThisEvent) < 1000000000 ) << "\n";
    cout << "entry " << i << " : " << (timeThisEvent - starttimer)*1e-12 << "\n";

    //if events are within 1ms of each other then we have entered the spill
     if ( abs(timeNextEvent-timeThisEvent) < 2000000000 ) {
       index_firstEventInSpill = i;
       TOFHIRTriggerTimeOfFirstEventInSpill = timeThisEvent;
       break;
     }
  }
  cout << "first event in spill for TOFHIR : " << index_firstEventInSpill << " : " << TOFHIRTriggerTimeOfFirstEventInSpill << "\n";

  //======================== Find the start of spill in Tracker data =======================================

  double TrackerBCOOfFirstEventInSpill = 0;
  int IndexTrackerFirstEventInSpill = -1;
  const double clockScaleFactor = 1.0420863;
  
  //Find First Event In Spill in Tracker data
  for(int k=0; k<trackerTree->GetEntries()-1; k++) {
    trackerTree->GetEntry( k+1 );
    Long64_t TrackerBCONextEvent = trackerEvent->bco;
    trackerTree->GetEntry( k );
    Long64_t TrackerBCOThisEvent = trackerEvent->bco;
    
    cout << "Tracker Trigger: " << trackerEvent->trigger << " : " << TrackerBCONextEvent << " - " 
	 << TrackerBCOThisEvent << " = " << (TrackerBCONextEvent - TrackerBCOThisEvent)*144e-9*clockScaleFactor << " | "
	 << " pass ? " << bool( (TrackerBCONextEvent - TrackerBCOThisEvent)*144e-9*clockScaleFactor < 0.002 ) << "\n";
    
    //if events are within 2ms of each other then we have entered the spill
    if ( (TrackerBCONextEvent - TrackerBCOThisEvent)*144e-9*clockScaleFactor < 0.002 ) {
      TrackerBCOOfFirstEventInSpill = TrackerBCOThisEvent;
      IndexTrackerFirstEventInSpill = k;
      cout << "Found it : " << trackerEvent->trigger << "\n";
      break;
    }            
  }

  cout << "first event in spill for tracker : " << IndexTrackerFirstEventInSpill 
       << " : " << TrackerBCOOfFirstEventInSpill*144e-9*clockScaleFactor << "\n";


  //======================== Look for coincidences =======================================
  int PreviousMatchIndex = IndexTrackerFirstEventInSpill;


  cout << "debug: " << index_firstEventInSpill << "\n";


  //Loop over the trigger signals
  for (Int_t q=index_firstEventInSpill;q<TOFHIR_TriggerIndices.size(); q++) {
    //for (Int_t q=index_firstEventInSpill;q<1000; q++) { // FIXME   added by Abdollah


    if (q%1==0) cout << "q = " << q << "\n";
    anatree->GetEntry(TOFHIR_TriggerIndices[q]);

   // cout << "nentries_ana= "<<nentries_ana<<"\n";
   // cout <<channelID <<"  "<<triggerChannelID<<"\n";

    //Do event building seeded on the triggers
    if (channelID != triggerChannelID) continue;
    
    //initialize at beginning of each event
    for(int k=0;k<400;k++){
      chTime[k]=-9999;
      chtot[k]=-9999;
    }
    
    //populate trigger data
    chTime[triggerChannelID] = time;
    chtot[triggerChannelID] = tot;
    
    //cout << "TOFHIR Event : " << q << " " << time << "\n";
    
    //find signal hits that correspond to the given trigger
    // We find that sipms are typically around 170ns EARLIER than the trigger timestamp.
    // So look within a 10ns window around that.
    Long64_t tTrigger = time;
    const int maxSearchSteps = 100;
    int j_start = TOFHIR_TriggerIndices[q]-maxSearchSteps; if (j_start < 0) j_start = 0;
    for (Int_t j=j_start;j<TOFHIR_TriggerIndices[q]+50; j++) {
      anatree->GetEntry(j);
      
      Long64_t tdiff = tTrigger - time;
      
      // cout << "channel " << channelID << " " << time << " | " << tdiff << "\n";
	
      //If channel falls within our search window, then populate the data for that channel
      if (tdiff >= 170000 - SEARCHWINDOWSIZE &&  tdiff <= 170000 + SEARCHWINDOWSIZE) {
	chTime[channelID] = time;
	chtot[channelID] = tot;
	//cout << "Coincidence: channel " << channelID << " " << time << " | " << tdiff << "\n";
      } 
      //else if (tdiff >= 170000 - 10000 &&  tdiff <= 170000 + 10000) {  
      //cout << "channel " << channelID << " " << time << " | " << tdiff << "\n";
      //}
      
    }


    //match TOFHIR trigger time to tracker timestamp  
    const double TOFHIRClockScaleFactor = 0.99609467;
    double TOFHIRTriggerTimestamp = (tTrigger - TOFHIRTriggerTimeOfFirstEventInSpill)*1e-12*TOFHIRClockScaleFactor;

    if (DEBUG==1) cout << std::setprecision(10) << "TOFHIR Time: " << TOFHIRTriggerTimestamp << "\n";

  

    //Now Match tracker event with TOFHIR event
    int NMatchedTracks = 0;
    for(int k=PreviousMatchIndex; k<trackerTree->GetEntries(); k++) {
      trackerTree->GetEntry( k );
      Long64_t TrackerBCOThisEvent = trackerEvent->bco;

      if (DEBUG==1) {
        cout << std::setprecision(10)  
	     <<"Trigger " << trackerEvent->trigger << " : " 
	     << " " << trackerEvent->bco << " "
	     << (trackerEvent->bco - TrackerBCOOfFirstEventInSpill)* 144e-9 * clockScaleFactor << "\n";      
      }
      double trackerTime = (trackerEvent->bco - TrackerBCOOfFirstEventInSpill) * 144e-9 * clockScaleFactor;
      
      //anything that came before the first trigger of the spill should be skipped
      //if (TOFHIRTriggerTimestamp - trackerTime < -0.00001) continue;
      //if (TOFHIRTriggerTimestamp - trackerTime > 0.001) continue;
 
      //now find a match with TOFHIR trigger time

      if (DEBUG==1) cout<< "Time difference is "<< TOFHIRTriggerTimestamp - trackerTime << "\t if it is smaller than "<< MATCH_THRESHOLD << " it finds the match\n";

      if (fabs (TOFHIRTriggerTimestamp - trackerTime ) < MATCH_THRESHOLD){
	if (DEBUG==1) cout << "MATCHED : " << event << " : " << TOFHIRTriggerTimestamp << "   | "<< trackerTime<<"\n";	
	NMatchedTracks++;
	PreviousMatchIndex = k;	


	//populate tracking data
	xIntercept=trackerEvent->xIntercept * 1e-3;
	yIntercept=trackerEvent->yIntercept * 1e-3;
	xSlope=trackerEvent->xSlope * 1e-3;
	ySlope=trackerEvent->ySlope * 1e-3;
	x_dut= (trackerEvent->xIntercept + trackerEvent->xSlope * 2.0e5) * 1e-3;
	y_dut= (trackerEvent->yIntercept + trackerEvent->ySlope * 2.0e5) * 1e-3;
	chi2 = trackerEvent->chi2;
	ntracks = NMatchedTracks;
	nplanes = trackerEvent->nPlanes;    

      }
// Note: I have checked and found out that the fraction of events having 2 tracks is ~ 0.5%, so it might be safe to discard them?!
//
//      cout << "Trigger: " << trackerEvent->trigger  << " : " << TOFHIRTriggerTimestamp << " " << trackerTime << " "
//	   << TOFHIRTriggerTimestamp - trackerTime << "\n";
      if (TOFHIRTriggerTimestamp - trackerTime < -0.001) {
	//cout << "BREAK\n";
     	break;
      }

    }
    
    //cout << "\n\n";
    outTree->Fill();
    event++;
             
  }
  
  //outTree->Write();
  outFile->Write();
  outFile->Close();
  AnaFile->Close();
  //trackerFile->Close();
}
