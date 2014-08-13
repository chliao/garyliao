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
#include "TCanvas.h"

void xAna5(std::string inputFile_) {

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
  TH1D* hM_el_Ntp_B = (TH1D*)hM->Clone("hM_el_Ntp_B");
  TH1D* hM_el_Ntp_E = (TH1D*)hM->Clone("hM_el_Ntp_E");
  TH1D* hM_el_Ntt = (TH1D*)hM->Clone("hM_el_Ntt");
  TH1D* hM_el_Ntf_B = (TH1D*)hM->Clone("hM_el_Ntf");
  TH1D* hM_el_Ntf_E = (TH1D*)hM->Clone("hM_el_Ntf");
//TH1D* hM_mu = (TH1D*)hM->Clone("hM_mu");
  
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
    Float_t* eleSCEta = data.GetPtrFloat("eleSCEta");
    Int_t nEle        = data.GetInt("nEle");


    // electron selection
    // Ntp
    vector<int> acc_ele1;
    vector<int> acc_ele2;
    eID2012(data, acc_ele1, 1);       // cut based eID (1: loose WP)
    eID2012(data, acc_ele2, 3);
    // two-lepton loop
    for (size_t i = 0; i < acc_ele1.size(); i++) {
      for (size_t j = 0; j < acc_ele2.size(); j++) {
   	
        if(acc_ele1[i] == acc_ele2[j])continue;
        
        if (fabs(eleSCEta[acc_ele1[i]])>2.5) continue;
        if (fabs(eleSCEta[acc_ele2[j]])>2.5) continue;

        if (fabs(eleSCEta[acc_ele1[i]]) < 1.4442 && fabs(eleSCEta[acc_ele2[j]]) < 1.4442 ) {	
	TLorentzVector ele1, ele2;
	ele1.SetPtEtaPhiM(elePt[acc_ele1[i]], eleEta[acc_ele1[i]], elePhi[acc_ele1[i]], 0.000511);
	ele2.SetPtEtaPhiM(elePt[acc_ele2[j]], eleEta[acc_ele2[j]], elePhi[acc_ele2[j]], 0.000511);
    // fill histo
        TLorentzVector Z = ele1 + ele2;
        hM_el_Ntp_B->Fill(Z.M());
	}
        else if (fabs(eleSCEta[acc_ele1[i]]) > 1.566 && fabs(eleSCEta[acc_ele1[i]]) < 2.5 && fabs(eleSCEta[acc_ele2[j]]) < 1.4442 ){
        TLorentzVector ele1, ele2;
        ele1.SetPtEtaPhiM(elePt[acc_ele1[i]], eleEta[acc_ele1[i]], elePhi[acc_ele1[i]], 0.000511);
        ele2.SetPtEtaPhiM(elePt[acc_ele2[j]], eleEta[acc_ele2[j]], elePhi[acc_ele2[j]], 0.000511);     
	// fill histo
      TLorentzVector Z = ele1 + ele2;
      hM_el_Ntp_E->Fill(Z.M());
       }
        else {continue;}
      } // for2
    }// for1
    // Ntt
    vector<int> acc_ele;
    eID2012(data, acc_ele, 3);       // cut based eID (1: loose WP)
    // two-lepton loop
    for (size_t i = 0; i < acc_ele.size(); i++) {                 
      for (size_t j = i + 1; j < acc_ele.size(); j++) {

        if (fabs(eleSCEta[acc_ele[i]])>2.5) continue;
        if (fabs(eleSCEta[acc_ele[j]])>2.5) continue;
        if (fabs(eleSCEta[acc_ele[i]]) < 1.4442 && fabs(eleSCEta[acc_ele[j]]) < 1.4442 ){
        TLorentzVector ele1, ele2;
        ele1.SetPtEtaPhiM(elePt[acc_ele[i]], eleEta[acc_ele[i]], elePhi[acc_ele[i]], 0.000511);
        ele2.SetPtEtaPhiM(elePt[acc_ele[j]], eleEta[acc_ele[j]], elePhi[acc_ele[j]], 0.000511);
                        				
         // fill histo
       TLorentzVector Z = ele1 + ele2;
       hM_el_Ntt->Fill(Z.M());
        }
        else {continue;}
      }
    }
 
    // Ntf
    vector<int> fail_accepted;
        
    
     for (int k = 0; k < nEle; k++) {

      if (elePt[k] < 20.) continue;
      if (fabs(eleSCEta[k])>2.5) continue;
   
     for (size_t i = 0; i < acc_ele1.size(); i++) {

   //  cout <<"i="<< i << endl;
  //   cout <<"k="<< k << endl;
   // int tempt = acc_ele1[i];
   // cout <<"tempt:"<< tempt <<endl;
      if (k == acc_ele1[i]) continue;
  }
    fail_accepted.push_back(k);
 }
       
 
    vector<int> acc_ele3;

    eID2012(data, acc_ele3, 3);      // cut based eID (1: loose WP)  
 
    // two-lepton loop
    for (size_t i = 0; i < acc_ele3.size(); i++) {                       							                
     for (size_t j = 0; j < fail_accepted.size(); j++) {

      if (acc_ele3[i] == fail_accepted[j]) continue ;
       
      if (fabs(eleSCEta[acc_ele3[i]])>2.5) continue;
      if (fabs(eleSCEta[fail_accepted[j]])>2.5) continue;
      if (fabs(eleSCEta[acc_ele3[i]]) < 1.4442 && fabs(eleSCEta[fail_accepted[j]]) < 1.4442 ){
       TLorentzVector ele1, ele2;
       ele1.SetPtEtaPhiM(elePt[acc_ele3[i]], eleEta[acc_ele3[i]], elePhi[acc_ele3[i]], 0.000511);
       ele2.SetPtEtaPhiM(elePt[fail_accepted[j]], eleEta[fail_accepted[j]], elePhi[fail_accepted[j]], 0.000511);
        // fill histo     
      TLorentzVector Z = ele1 + ele2;
      hM_el_Ntf_B->Fill(Z.M());
          } 
    else if (fabs(eleSCEta[acc_ele3[i]]) < 1.4442 && fabs(eleSCEta[fail_accepted[j]]) > 1.566 && fabs (eleSCEta[fail_accepted[j]]) < 2.5 ){
       TLorentzVector ele1, ele2;
       ele1.SetPtEtaPhiM(elePt[acc_ele3[i]], eleEta[acc_ele3[i]], elePhi[acc_ele3[i]], 0.000511);
       ele2.SetPtEtaPhiM(elePt[fail_accepted[j]], eleEta[fail_accepted[j]], elePhi[fail_accepted[j]], 0.000511);
        // fill histo     
           TLorentzVector Z = ele1 + ele2;
           hM_el_Ntf_E->Fill(Z.M());
           }
    else {continue;} 
                   
        }
      }
                           							                                                                  						              
  }// event loop
   // numbers of entries
   double Ntt_n   = hM_el_Ntt->GetEntries();
   double Ntp_B_n = hM_el_Ntp_B->GetEntries();
   double Ntp_E_n = hM_el_Ntp_E->GetEntries();
   double Ntf_B_n = hM_el_Ntf_B->GetEntries();
   double Ntf_E_n = hM_el_Ntf_E->GetEntries();

   cout << "Ntt_n:" << Ntt_n << endl;
   cout << "Ntp_B_n:" << Ntp_B_n << endl;
   cout << "Ntp_E_n:" << Ntp_E_n << endl;
   cout << "Ntf_B_n:" << Ntf_B_n << endl;
   cout << "Ntt_E_n:" << Ntf_E_n << endl;

   double  efficiency_B = (2*Ntt_n + Ntp_B_n)/(2*Ntt_n + Ntp_B_n + Ntf_B_n);
   double  efficiency_E = (Ntp_E_n)/(Ntp_E_n + Ntf_E_n);

   cout << "efficiency_B:" << efficiency_B << endl;
   cout << "efficiency_E:" << efficiency_E << endl;

  TCanvas * Ntp_B = new TCanvas("Ntp_B","Ntp_B",500,500);
  hM_el_Ntp_B->Draw();

  TCanvas * Ntp_E = new TCanvas("Ntp_E","Ntp_E",500,500);
  hM_el_Ntp_E->Draw();
  
  TCanvas * Ntt = new TCanvas("Ntt","Ntt",500,500);
  hM_el_Ntt->Draw();

  TCanvas * Ntf_B = new TCanvas("Ntf_B","Ntf_B",500,500);
  hM_el_Ntf_B->Draw();

  TCanvas * Ntf_E = new TCanvas("Ntf_E","Ntf_E",500,500);
  hM_el_Ntf_E->Draw();

  fprintf(stderr, "Processed all events\n");
  TFile* outFile = new TFile("zmass_emu.root","recreate");

  gStyle->SetOptFit(11111);
  hM_el_Ntp_B->Fit("gaus");

  gStyle->SetOptFit(11111);
  hM_el_Ntp_E->Fit("gaus");

  gStyle->SetOptFit(11111);
  hM_el_Ntt->Fit("gaus");

  gStyle->SetOptFit(11111);
  hM_el_Ntf_B->Fit("gaus");

  gStyle->SetOptFit(11111);
  hM_el_Ntf_E->Fit("gaus");
 
 // hM_el->Draw();
  
  hM_el_Ntp_B->Write();
  hM_el_Ntp_E->Write();
  hM_el_Ntt->Write();
  hM_el_Ntf_B->Write();
  hM_el_Ntf_E->Write();

 // hM_mu->Write();

  outFile->Write();

}
