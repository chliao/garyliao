/* Simple xAna analysis example. */
#include <string>
#include <vector>
#include <iostream>
#include <TH1D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TStyle.h>
#include "puweicalc.h"
#include "untuplizer.h"
#include "ElectronSelections.h"
#include "MuonSelections.h"
#include "PhotonSelections.h"

void xAna(std::string inputFile_) {

  // get TTree from file ...
  TreeReader data(inputFile_.data()); // v5.3.12

  // ... or make a chain of root files with TTrees
  // const char* paths[] = {"ggtree_mc_1.root", "ggtree_mc_2.root"};
  // TreeReader data(paths, 2);
  
  // useful to determine which type of variable to use for which branches
  // data.Print();
  
  // do whathever preparations are necessary for, if MC information is present
  // pileup reweighting for MC                                                                                          
  Float_t puweiEle, puweiMu;
  PUWeightCalculator puCalcEle;
  PUWeightCalculator puCalcMu;
  if (data.HasMC()) {
    puCalcEle.Init("mcweights/mcwei_PU_RD1_ele.root");
    puCalcMu.Init("mcweights/mcwei_PU_RD1_muo.root");
  }
  
  TH1D* hM = new TH1D("hM", "Two-lepton invariant mass", 90, 50, 140);
  TH1D* hM_mu = (TH1D*)hM->Clone("hM_mu");
  TH1D* hM_el = (TH1D*)hM->Clone("hM_el");

  // event loop
  // for (Long64_t ev = 0; ev < 10000; ev++) {
  for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++) {
    // print progress
    if (ev % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());
    
    data.GetEntry(ev);

    Int_t run = data.GetInt("run");

    // IMPORTANT: branches with counters must be loaded BEFORE dependent branches
    // accessing HLT information
    //    Int_t nHLT = data.GetInt("nHLT");
    //    Int_t* HLTIndex = data.GetPtrInt("HLTIndex");

    //for (Int_t i = 0; i < nHLT; ++i) printf(" %i", HLTIndex[i]);
    //printf("\n");
    
    // IMPORTANT: branches with counters must be loaded BEFORE dependent branches
    Float_t* elePt    = data.GetPtrFloat("elePt");
    Float_t* eleEta   = data.GetPtrFloat("eleEta");
    Float_t* elePhi   = data.GetPtrFloat("elePhi");


    // electron selection
    vector<int> acc_ele1;
    vector<int> acc_ele2;
    eID2012(data, acc_ele1, 1);       // cut based eID (1: loose WP)
 //   eID2012(data, acc_ele2, 3);
    // two-lepton loop
    for (size_t i = 0; i < acc_ele1.size(); i++) {
      for (size_t j = i+1; j < acc_ele1.size(); j++) {

   // double_t loose
      //loose = acc_ele1[i]
    //double_t tight
      //tight = acc_ele2[j]
    //if (tight==loose)continue;
	// cout<<"tight_ele_index:"<<acc_ele2[j]<<" loose_ele_index:"<<acc_ele1[i]<<endl;
	// if(acc_ele1[i]==acc_ele2[j])continue;

	TLorentzVector ele1, ele2;
	ele1.SetPtEtaPhiM(elePt[acc_ele1[i]], eleEta[acc_ele1[i]], elePhi[acc_ele1[i]], 0.000511);
	ele2.SetPtEtaPhiM(elePt[acc_ele1[j]], eleEta[acc_ele1[j]], elePhi[acc_ele1[j]], 0.000511);
	
	// fill histo
	TLorentzVector Z = ele1 + ele2;
	hM_el->Fill(Z.M());
      } // for2
    } // for1

  }// event loop
  
  fprintf(stderr, "Processed all events\n");
  TFile* outFile = new TFile("zmass_emu.root","recreate");

  gStyle->SetOptFit(11111);
  hM_el->Fit("gaus");
  hM_el->Draw();
  
  hM_el->Write();
  hM_mu->Write();

  outFile->Write();

}
