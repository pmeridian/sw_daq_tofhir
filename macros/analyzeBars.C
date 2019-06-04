void analyzeBars(const std::string& inputDir, const std::string& runs)
{
  //--- open files and make the tree chain
  TChain* tree = new TChain("data","data");
  
  std::stringstream ss(runs);
  std::string token;
  while( std::getline(ss,token,',') )
    {
      std::stringstream ss2(token);
      std::string token2;
      int runMin = -1;
      int runMax = -1;
      while( std::getline(ss2,token2,'-') )
	{
	  if( runMin != -1 && runMax == -1 ) runMax = atoi(token2.c_str());
	  if( runMin == -1 ) runMin = atoi(token2.c_str());
	}
      for(int run = runMin; run <= runMax; ++run)
	{
	  std::string fileName = Form("%s/*%04d*.root",inputDir.c_str(),run);
	  std::cout << ">>> Adding flle " << fileName << std::endl;
	  tree -> Add(fileName.c_str());
	}
    }
  
  
  //--- define channels
  std::vector<std::pair<std::string,std::string> > channelPairs;
  channelPairs.push_back(std::make_pair("ch0","ch14"));
  channelPairs.push_back(std::make_pair("ch1","ch15"));
  
  
  //--- define branches
  float step1, step2;
  std::vector<float>* tot = 0;
  std::vector<long long>* time = 0;
  std::map<std::string,int> channelIdx;
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1",1); tree -> SetBranchAddress("step1",&step1);
  tree -> SetBranchStatus("step2",1); tree -> SetBranchAddress("step2",&step2);
  tree -> SetBranchStatus("tot",  1); tree -> SetBranchAddress("tot",  &tot);
  tree -> SetBranchStatus("time", 1); tree -> SetBranchAddress("time", &time);
  for(auto channelPair : channelPairs)
    {
      tree -> SetBranchStatus(("%s",channelPair.first.c_str()), 1);  tree -> SetBranchAddress(("%s",channelPair.first.c_str()), &channelIdx[channelPair.first]);
      tree -> SetBranchStatus(("%s",channelPair.second.c_str()),1);  tree -> SetBranchAddress(("%s",channelPair.second.c_str()),&channelIdx[channelPair.second]);
    }
  
  
  //--- define histograms
  std::vector<std::string> stepLabels;
  
  std::map<std::string,TH1F*> h1_tot;
  std::map<std::string,TH2F*> h2_tot_corr;
  
  std::map<std::string,TH1F*> h1_timeDiff;
  
  
  
  
  //--------------------
  //--- loop over events
  int nEntries = tree->GetEntries();
  for(int entry = 0; entry < nEntries; ++entry)
    {
      tree -> GetEntry(entry);
      if( entry%1000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
      
      std::string stepLabel = Form("Vov%.1f_th%.0f",step1,step2);
      
      
      //--- create histograms, if needed
      for(auto channelPair : channelPairs)
	{
	  std::string ch1 = channelPair.first;
	  std::string ch2 = channelPair.second;
	  std::string label1 = Form("%s_%s",ch1.c_str(),stepLabel.c_str());
	  std::string label2 = Form("%s_%s",ch2.c_str(),stepLabel.c_str());
	  std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());

	  if( h1_tot[label1] == NULL )
	    {
	      h1_tot[label1] = new TH1F(Form("h1_tot_%s",label1.c_str()),"",1000,0.,1000.);
	      h1_tot[label2] = new TH1F(Form("h1_tot_%s",label2.c_str()),"",1000,0.,1000.);
	      h2_tot_corr[label12] = new TH2F(Form("h2_tot_corr_%s",label12.c_str()),"",100,0.,1000.,100,0.,1000.);
	      
	      h1_timeDiff[label12] = new TH1F(Form("h1_timeDiff_%s",label12.c_str()),"",5000,-20000.,20000.);
	      
	      stepLabels.push_back(stepLabel);
	    }
	}
      
      
      //--- fill histograms
      for(auto channelPair : channelPairs)
	{
	  std::string ch1 = channelPair.first;
	  std::string ch2 = channelPair.second;
	  std::string label1 = Form("%s_%s",ch1.c_str(),stepLabel.c_str());
	  std::string label2 = Form("%s_%s",ch2.c_str(),stepLabel.c_str());
	  std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
	  
	  if( channelIdx[ch1] >= 0 ) h1_tot[label1] -> Fill( tot->at(channelIdx[ch1])/1000. );
	  if( channelIdx[ch2] >= 0 ) h1_tot[label2] -> Fill( tot->at(channelIdx[ch2])/1000. );
	  if( channelIdx[ch1] >= 0 && channelIdx[ch2] >= 0 )
	    {
	      h2_tot_corr[label12] -> Fill( tot->at(channelIdx[ch1])/1000.,tot->at(channelIdx[ch2])/1000. );
	      
	      if( tot->at(channelIdx[ch1])/1000. > 200. && tot->at(channelIdx[ch2])/1000. > 200. )
		{
		  h1_timeDiff[label12] -> Fill( time->at(channelIdx[ch2])-time->at(channelIdx[ch1]) );
		}
	    }
	}
      
      if( channelPairs.size() == 2)
	{
	  std::string ch1 = channelPairs.at(0).first;
	  std::string ch2 = channelPairs.at(0).second;
	  std::string label1 = Form("%s_%s",ch1.c_str(),stepLabel.c_str());
	  std::string label2 = Form("%s_%s",ch2.c_str(),stepLabel.c_str());
	  
	  std::string ch3 = channelPairs.at(1).first;
	  std::string ch4 = channelPairs.at(1).second;
	  std::string label3 = Form("%s_%s",ch3.c_str(),stepLabel.c_str());
	  std::string label4 = Form("%s_%s",ch4.c_str(),stepLabel.c_str());
	  
	  std::string label1234 = Form("%s-%s_%s-%s_%s",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str(),stepLabel.c_str());
	  
	  if( channelIdx[ch1] >= 0 && channelIdx[ch2] >= 0 && 
	      channelIdx[ch3] >= 0 && channelIdx[ch4] >= 0 )
            {
	      if( tot->at(channelIdx[ch1])/1000. > 200. && tot->at(channelIdx[ch2])/1000. > 200. && 
		  tot->at(channelIdx[ch3])/1000. > 200. && tot->at(channelIdx[ch4])/1000. > 200. )
		{
		  if( h1_timeDiff[label1234] == NULL ) 
		    h1_timeDiff[label1234] = new TH1F(Form("h1_timeDiff_%s",label1234.c_str()),"",5000,-20000.,20000.);
		  
		  h1_timeDiff[label1234] -> Fill( 0.5*(time->at(channelIdx[ch3])+time->at(channelIdx[ch4])) - 0.5*(time->at(channelIdx[ch1])+time->at(channelIdx[ch2])) );
		}
	    }
	}
    }
  std::cout << std::endl;
  
  std::vector<std::string>::iterator iter;
  iter = std::unique(stepLabels.begin(),stepLabels.end());
  stepLabels.resize( std::distance(stepLabels.begin(),iter) );  
  
  
  
  
  //--------------
  //--- draw plots
  for(auto stepLabel : stepLabels)
    {
      for(auto channelPair : channelPairs)
	{
	  std::string ch1 = channelPair.first;
	  std::string ch2 = channelPair.second;
	  std::string label1 = ch1 + "_" + stepLabel;
	  std::string label2 = ch2 + "_" + stepLabel;
	  std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
	  
	  
	  TCanvas* c1 = new TCanvas(Form("c1_tot_%s",label12.c_str()),
				    Form("c1_tot_%s",label12.c_str()));
	  gPad -> SetLogy();
	  
	  h1_tot[label1] -> SetTitle(";ToT [ns];entries");
	  h1_tot[label1] -> SetLineColor(kRed);
	  h1_tot[label1] -> Draw();
	  
	  h1_tot[label2] -> SetTitle(";ToT [ns];entries");
	  h1_tot[label2] -> SetLineColor(kBlue);
	  h1_tot[label2] -> Draw("same");
	  
	  
	  TCanvas* c2 = new TCanvas(Form("c2_tot_corr_%s",label12.c_str()),Form("c2_tot_corr_%s",label12.c_str()));
	  
	  h2_tot_corr[label12] -> SetTitle(Form(";%s ToT [ns];%s ToT [ns]",ch1.c_str(),ch2.c_str()));
	  h2_tot_corr[label12] -> Draw("colz");
	  
	  
	  TCanvas* c3 = new TCanvas(Form("c3_timeDiff_%s",label12.c_str()),Form("c3_timeDiff_%s",label12.c_str()));
	  gPad -> SetLogy();
	  
	  h1_timeDiff[label12] -> SetTitle(Form(";#Deltat [ps];entries"));
	  h1_timeDiff[label12] -> Draw("");
	}
      
      if( channelPairs.size() == 2)
	{
	  std::string ch1 = channelPairs.at(0).first;
	  std::string ch2 = channelPairs.at(0).second;
	  std::string label1 = Form("%s_%s",ch1.c_str(),stepLabel.c_str());
	  std::string label2 = Form("%s_%s",ch2.c_str(),stepLabel.c_str());
	  
	  std::string ch3 = channelPairs.at(1).first;
	  std::string ch4 = channelPairs.at(1).second;
	  std::string label3 = Form("%s_%s",ch3.c_str(),stepLabel.c_str());
	  std::string label4 = Form("%s_%s",ch4.c_str(),stepLabel.c_str());
	  
	  std::string label1234 = Form("%s-%s_%s-%s_%s",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str(),stepLabel.c_str());
	  
	  
	  TCanvas* c3 = new TCanvas(Form("c3_timeDiff_%s",label1234.c_str()),Form("c3_timeDiff_%s",label1234.c_str()));
	  gPad -> SetLogy();
	  
	  h1_timeDiff[label1234] -> SetTitle(Form(";#Deltat [ps];entries"));
	  h1_timeDiff[label1234] -> Draw("");
	}
    }
}
