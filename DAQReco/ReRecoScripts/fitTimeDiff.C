//short list at /tmp/deguio/shortList.txt
//full list at /afs/cern.ch/work/d/deguio/MTD/TBStudies/macros/runList.txt

void fitTimeDiff(std::string fileList="/afs/cern.ch/work/d/deguio/MTD/TBStudies/macros/runList.txt")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(1);
  std::ofstream outfile("fitParameters.txt");
  outfile << "#fileName\tintercept\tslope\tchi2/ndof\n";

  std::string prefix = "/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RecoData/v1/RecoWithTracks/run";
  std::string suffix = "_events.root";

  ifstream inFile(fileList.c_str());

  std::vector<std::string> vFile;
  std::string filePath;
  while (true) {
		inFile >> filePath;
    if (inFile.eof())
    {
      break;
    }
    vFile.push_back(filePath);
  }

  TF1* myLine = new TF1("myLine","[0]+x*[1]");
  myLine->SetLineWidth(3);
  for(auto iFile = vFile.begin(); iFile!=vFile.end(); ++iFile)
  {
    std::cout << "-->> " << *iFile << std::endl;
    TFile* myfile = TFile::Open((prefix+*iFile+suffix).c_str());

    Long64_t time[256];
    TTree* ntu = (TTree*)myfile->Get("data");
    if(ntu==nullptr || ntu->GetEntries() < 1)
    {
      myfile->Close();
      continue;
    }
    ntu->SetBranchAddress("time",time);
    ntu->GetEntry();
    double initialTime = time[255]/1e12;

    TProfile* timeTrend = new TProfile("timeTrend","timeTrend", 100, initialTime, initialTime+5);
    ntu->Draw("TimeDiff:(time[255]/1e12) >> timeTrend","","goff");

    myLine->SetParameters(-1e-5, 1e-7);
    myLine->SetRange(initialTime, initialTime+1.5);
    timeTrend->Fit(myLine,"QR");

    outfile << std::fixed << std::setprecision(10) << *iFile << "\t"
            << myLine->GetParameter(0) << "\t" << myLine->GetParameter(1) << "\t" << myLine->GetChisquare()/myLine->GetNDF() << "\n";

    if(myLine->GetChisquare()/myLine->GetNDF() > 15)
    {
      TCanvas* c1 = new TCanvas();
      timeTrend->DrawClone();
      c1->Print(("~/www/mtd/triggerTimeCalibration/pics/"+*iFile+".png").c_str());
    }

    myfile->Close();
  }
  outfile.close();

  return;
}
