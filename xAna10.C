/* Simple xAna analysis example. */
#include <string>
#include <vector>
#include <iostream>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TStyle.h>
#include "puweicalc.h"
#include "untuplizer.h"
#include "modified_ElectronSelections.h"
#include "MuonSelections.h"
#include "PhotonSelections.h"
#include "TCanvas.h"
#include <cmath>

void xAna10(std::string inputFile_) {

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
    cout<<"pleasr cin pt_max"<<endl;
  double pt_max;
  cin >> pt_max;
  cout<<"please cin bin number"<<endl;
  int n;
  cin >> n;
  double ef_B[n];
  ef_B[n] = {0};
  double unc_B[n];
  unc_B[n]={0};
  double ef_E[n];
  ef_E[n]={0};
  double unc_E[n];
  unc_E[n]={0};

  double relative_unc_B[n];
  relative_unc_B[n]={0};
  double relative_unc_E[n];
  relative_unc_E[n]={0};
/*
  for (int c = 0; c < n; c++){
  double min = c*((pt_max-20)/n)+20;
  double max = (c+1)*((pt_max-20)/n)+20;
*/
  TH1D* hM = new TH1D("hM", "Two-lepton invariant mass", 90, 50, 140);
  TH1D* hM_el_Ntp_B = (TH1D*)hM->Clone("hM_el_Ntp_B");
  TH1D* hM_el_Ntp_E = (TH1D*)hM->Clone("hM_el_Ntp_E");
  TH1D* hM_el_Ntt_B = (TH1D*)hM->Clone("hM_el_Ntt_B");
  TH1D* hM_el_Ntt_E = (TH1D*)hM->Clone("hM_el_Ntt_E");
  TH1D* hM_el_Ntf_B = (TH1D*)hM->Clone("hM_el_Ntf_B");
  TH1D* hM_el_Ntf_E = (TH1D*)hM->Clone("hM_el_Ntf_E");

//  TH1D* hM2 = new TH1D("hM2", "The Tag and Probe Efficiency of Different Pt", pt_max-20, 20, pt_max, 100, 0, 1);
  TH1D* hM2 = new TH1D("hM2", "The Tag and Probe Efficiency of Different Pt", n, 20, pt_max);
  TH1D* hM2_el_efficiency_B = (TH1D*)hM2->Clone("hM2_el_efficiency_B");
  TH1D* hM2_el_efficiency_E = (TH1D*)hM2->Clone("hM2_el_efficiency_E");

  for (int c = 0; c < n; c++){
  double min = c*((pt_max-20)/n)+20;
  double max = (c+1)*((pt_max-20)/n)+20;

/*
  TH2D* hM2_el_efficiency_E = (TH2D*)hM2->Clone("hM2_el_efficiency_E");
TH2D* hM3 = new TH2D("hM3", "The Tag and Probe Uncertainty of Different Pt", pt_max-20, 20, pt_max, 100, 0, 1);
  TH2D* hM3_el_uncertainty_B = (TH2D*)hM3->Clone("hM3_el_uncertainty_B");
  TH2D* hM3_el_uncertainty_E = (TH2D*)hM3->Clone("hM3_el_uncertainty_E");
*/
 
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
    vector<int> acc_ele5;
    bool find = true;
    eID2012(data, acc_ele1, 1, min, max);       // cut based eID (1: loose WP)
    eID2012(data, acc_ele2, 3, 20, pt_max);
    // two-lepton loop
    for (size_t i = 0; i < acc_ele1.size(); i++) {
    for (size_t j = 0; j < acc_ele2.size(); j++) {
//	cout << "ii=" << i << endl;
//      cout << "j=" << j << endl; 
        if (acc_ele1[i] == acc_ele2[j]) {/*cout<<"break"<< endl;*/ find=false; break;}
        else {find=true;}
      
     }    
 
      if (find) { /*cout << "find=" << find << endl;
                       cout << "i=" << i << endl;
                       cout << "push_back" << endl;*/
                       acc_ele5.push_back(acc_ele1[i]);}
   }
 

/*
    vector<int> fail_tight;

     for (int l = 0; l < nEle; l++) {

      if (elePt[l] < 20.) continue;
      if (fabs(eleSCEta[l])>2.5) continue;

      for (size_t j = 0; j < acc_ele2.size(); j++) {

         if (l == acc_ele2[j]) continue;
         }
       fail_tight.push_back(l);
    }
    

     for (size_t i = 0; i < acc_ele1.size(); i++){
       for (size_t m = 0; m < fail_tight.size(); m++) {
         for (size_t j = 0; j < acc_ele2.size(); j++) {
         
         if (fabs(eleSCEta[acc_ele1[i]])>2.5) continue;
         if (fabs(eleSCEta[acc_ele2[j]])>2.5) continue;
         if (fabs(eleSCEta[fail_tight[m]])>2.5) continue;

         if (acc_ele1[i] == fail_tight[m]) {

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
         } 
             else {continue;}
      }
    }
  }
*/
        for (size_t l = 0; l < acc_ele5.size(); l++) {
        for (size_t j = 0; j < acc_ele2.size(); j++) {
     
        if (fabs(eleSCEta[acc_ele5[l]])>2.5) continue;
        if (fabs(eleSCEta[acc_ele2[j]])>2.5) continue;

        if (fabs(eleSCEta[acc_ele5[l]]) < 1.4442 && fabs(eleSCEta[acc_ele2[j]]) < 1.4442 ) {	
	TLorentzVector ele1, ele2;
	ele1.SetPtEtaPhiM(elePt[acc_ele5[l]], eleEta[acc_ele5[l]], elePhi[acc_ele5[l]], 0.000511);
	ele2.SetPtEtaPhiM(elePt[acc_ele2[j]], eleEta[acc_ele2[j]], elePhi[acc_ele2[j]], 0.000511);
    // fill histo
        TLorentzVector Z = ele1 + ele2;
        hM_el_Ntp_B->Fill(Z.M());
	}
        else if (fabs(eleSCEta[acc_ele5[l]]) > 1.566 && fabs(eleSCEta[acc_ele5[l]]) < 2.5 && fabs(eleSCEta[acc_ele2[j]]) < 1.4442 ){
        TLorentzVector ele1, ele2;
        ele1.SetPtEtaPhiM(elePt[acc_ele5[l]], eleEta[acc_ele5[l]], elePhi[acc_ele5[l]], 0.000511);
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
    vector<int> acc_ele4;
    eID2012(data, acc_ele, 3, 20, pt_max);       // cut based eID (1: loose WP)
    eID2012(data, acc_ele4, 3, min, max); 
    // two-lepton loop
    for (size_t i = 0; i < acc_ele.size(); i++) {                 
      for (size_t j = i + 1; j < acc_ele4.size(); j++) {

        if (fabs(eleSCEta[acc_ele[i]])>2.5) continue;
        if (fabs(eleSCEta[acc_ele4[j]])>2.5) continue;
        if (fabs(eleSCEta[acc_ele[i]]) < 1.4442 && fabs(eleSCEta[acc_ele4[j]]) < 1.4442 ){
        TLorentzVector ele1, ele2;
        ele1.SetPtEtaPhiM(elePt[acc_ele[i]], eleEta[acc_ele[i]], elePhi[acc_ele[i]], 0.000511);
        ele2.SetPtEtaPhiM(elePt[acc_ele4[j]], eleEta[acc_ele4[j]], elePhi[acc_ele4[j]], 0.000511);
                        				
         // fill histo
       TLorentzVector Z = ele1 + ele2;
       hM_el_Ntt_B->Fill(Z.M());
        }
      }
    }

    vector<int> acc_ele6;
    eID2012(data, acc_ele6, 3, min, max);
    for (size_t i = 0; i < acc_ele.size(); i++) {
      for (size_t j = 0 ; j < acc_ele6.size(); j++) {

        if (fabs(eleSCEta[acc_ele[i]])>2.5) continue;
        if (fabs(eleSCEta[acc_ele6[j]])>2.5) continue;
        if (acc_ele[i] == acc_ele6[j]) continue;
        if(fabs(eleSCEta[acc_ele[i]]) < 1.4442 && fabs(eleSCEta[acc_ele6[j]]) > 1.566 && fabs(eleSCEta[acc_ele6[j]]) < 2.5){
        TLorentzVector ele1, ele2;
        ele1.SetPtEtaPhiM(elePt[acc_ele[i]], eleEta[acc_ele[i]], elePhi[acc_ele[i]], 0.000511);
        ele2.SetPtEtaPhiM(elePt[acc_ele6[j]], eleEta[acc_ele6[j]], elePhi[acc_ele6[j]], 0.000511);

          // fill histo
       TLorentzVector Z = ele1 + ele2;
       hM_el_Ntt_E->Fill(Z.M());
        }
      }
    }
 
    // Ntf
    vector<int> fail_accepted;
    bool found = false ;    
    
     for (int k = 0; k < nEle; k++) {

      if (elePt[k] < 20.) continue;
      if (fabs(eleSCEta[k])>2.5) continue;
   
     for (size_t i = 0; i < acc_ele1.size(); i++) {

    // cout<<"kk="<<k<<endl;
    // cout<<"i="<<i<<endl;

      if (k == acc_ele1[i]) {/*cout<<"break"<<endl;*/ found = true; break;}
      else {found = false ;/* cout<<"found1="<<found<<endl;*/}   
   }
    // cout<<"found2="<<found<<endl;
      if (!found) { /*cout<<"k="<<k<<endl; cout<<"push_back"<<endl;*/ fail_accepted.push_back(k); }
 }
       
 
    vector<int> acc_ele3;

    eID2012(data, acc_ele3, 3, 20, pt_max);      // cut based eID (1: loose WP)  
 
    // two-lepton loop
    for (size_t i = 0; i < acc_ele3.size(); i++) {                       							                
     for (size_t j = 0; j < fail_accepted.size(); j++) {

   //   if (acc_ele3[i] == fail_accepted[j]) continue ;
      
     if (elePt[fail_accepted[j]] < min) continue;
     if (elePt[fail_accepted[j]] >= max) continue;

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
   double Ntt_B_n = hM_el_Ntt_B->GetEntries();
   double Ntt_E_n = hM_el_Ntt_E->GetEntries();
   double Ntp_B_n = hM_el_Ntp_B->GetEntries();
   double Ntp_E_n = hM_el_Ntp_E->GetEntries();
   double Ntf_B_n = hM_el_Ntf_B->GetEntries();
   double Ntf_E_n = hM_el_Ntf_E->GetEntries();
/*
   cout << "Ntt_B_n:" << Ntt_B_n << endl;
   cout << "Ntt_E_n:" << Ntt_E_n << endl;
   cout << "Ntp_B_n:" << Ntp_B_n << endl;
   cout << "Ntp_E_n:" << Ntp_E_n << endl;
   cout << "Ntf_B_n:" << Ntf_B_n << endl;
   cout << "Ntf_E_n:" << Ntf_E_n << endl;
*/
   double  efficiency_B_n = (2*Ntt_B_n + Ntp_B_n)/(2*Ntt_B_n + Ntp_B_n + Ntf_B_n);
   double  efficiency_E_n = (Ntt_E_n + Ntp_E_n)/(Ntt_E_n + Ntp_E_n + Ntf_E_n);
   double  uncertainty_B_n = sqrt((4*Ntt_B_n*Ntf_B_n*Ntf_B_n + Ntp_B_n*Ntf_B_n*Ntf_B_n + Ntf_B_n*(4*Ntt_B_n*Ntt_B_n + 4*Ntt_B_n*Ntp_B_n + Ntp_B_n*Ntp_B_n))/((2*Ntt_B_n+Ntp_B_n+ Ntf_B_n)*(2*Ntt_B_n+Ntp_B_n+ Ntf_B_n)*(2*Ntt_B_n+Ntp_B_n+ Ntf_B_n)*(2*Ntt_B_n+Ntp_B_n+ Ntf_B_n)));
   double  uncertainty_E_n = sqrt((Ntt_E_n*Ntf_E_n*Ntf_E_n + Ntp_E_n*Ntf_E_n*Ntf_E_n + Ntf_E_n*(Ntt_E_n + Ntp_E_n)*(Ntt_E_n + Ntp_E_n))/((Ntt_E_n + Ntp_E_n + Ntf_E_n)*(Ntt_E_n + Ntp_E_n + Ntf_E_n)*(Ntt_E_n + Ntp_E_n + Ntf_E_n)*(Ntt_E_n + Ntp_E_n + Ntf_E_n)));
//   double  relative_unc_B = uncertainty_B/efficiency_B;
//   double  relative_unc_E = uncertainty_E/efficiency_E;

   ef_B[c]=efficiency_B_n;
   unc_B[c]=uncertainty_B_n;
   ef_E[c]=efficiency_E_n;
   unc_E[c]=uncertainty_E_n;
//   relative_unc_B[c]= relative_unc_B;
//   relative_unc_E[c]= relative_unc_E;
   relative_unc_B[c]= unc_B[c]/ef_B[c];
   relative_unc_E[c]= unc_E[c]/ef_E[c];


  // cout << "efficiency_B:" << efficiency_B_n << endl;
  // cout << "efficiency_E:" << efficiency_E_n << endl;
   
  // cout << "uncertainty_B:" << uncertainty_B_n << endl;
  // cout << "uncertainty_E:" << uncertainty_E_n << endl;

   FILE *fout;
   fout = fopen("efficiency","w+t");
   fprintf(fout,"Ntt_B_n:%f\n", Ntt_B_n);
   fprintf(fout,"Ntt_E_n:%f\n", Ntt_E_n);
   fprintf(fout,"Ntp_B_n:%f\n", Ntp_B_n);
   fprintf(fout,"Ntp_E_n:%f\n", Ntp_E_n);
   fprintf(fout,"Ntf_B_n:%f\n", Ntf_B_n);
   fprintf(fout,"Ntf_E_n:%f\n", Ntf_E_n);
   fprintf(fout,"efficiency_B = (2*Ntt_B_n + Ntp_B_n)/(2*Ntt_B_n + Ntp_B_n + Ntf_B_n):%f\n", efficiency_B_n);
   fprintf(fout,"efficiency_E = (Ntt_E_n + Ntp_E_n)/(Ntt_E_n + Ntp_E_n + Ntf_E_n):%f\n", efficiency_E_n);
   fclose(fout);
/*
  TCanvas * Ntp_B = new TCanvas("Ntp_B","Ntp_B",500,500);
  hM_el_Ntp_B->Draw();

  TCanvas * Ntp_E = new TCanvas("Ntp_E","Ntp_E",500,500);
  hM_el_Ntp_E->Draw();
  
  TCanvas * Ntt_B = new TCanvas("Ntt_B","Ntt_B",500,500);
  hM_el_Ntt_B->Draw();

  TCanvas * Ntt_E = new TCanvas("Ntt_E","Ntt_E",500,500);
  hM_el_Ntt_E->Draw();

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
  hM_el_Ntt_B->Fit("gaus");

  gStyle->SetOptFit(11111);
  hM_el_Ntt_E->Fit("gaus");

  gStyle->SetOptFit(11111);
  hM_el_Ntf_B->Fit("gaus");

  gStyle->SetOptFit(11111);
  hM_el_Ntf_E->Fit("gaus");
*/
 TCanvas * efficiency_B = new TCanvas("efficiency_B","efficiency_B",500,500);
  hM2_el_efficiency_B->SetBinContent(c+1, ef_B[c]);
  hM2_el_efficiency_B->SetBinError(c+1, unc_B[c]);
  hM2_el_efficiency_B->Draw();

//hM2_el_efficiency_B->SetBinContent(c, ef_B[c]);
//hM2_el_efficiency_B->SetBinContent(42, ef_B[0]);
//hM2_el_efficiency_B->SetBinContent(43, ef_B[0]);
//hM2_el_efficiency_B->SetBinContent(44, ef_B[0]);
//hM2_el_efficiency_B->SetBinContent(45, ef_B[0]);

 TCanvas * efficiency_E = new TCanvas("efficiency_E","efficiency_E",500,500);
  hM2_el_efficiency_E->SetBinContent(c+1, ef_E[c]);
  hM2_el_efficiency_E->SetBinError(c+1, unc_E[c]);
  hM2_el_efficiency_E->Draw();

 // hM_el->Draw();

/*  
  hM_el_Ntp_B->Write();
  hM_el_Ntp_E->Write();
  hM_el_Ntt_B->Write();
  hM_el_Ntt_E->Write();
  hM_el_Ntf_B->Write();
  hM_el_Ntf_E->Write();
*/

 // hM_mu->Write();

 // outFile->Write();
  }//for k

  for (int k = 0; k < n; k++){
    cout<<"ef_B["<<k<<"]="<<ef_B[k]<<endl;
    cout<<"unc_B["<<k<<"]="<<unc_B[k]<<endl;
    cout<<"ef_E["<<k<<"]="<<ef_E[k]<<endl;
    cout<<"unc_E["<<k<<"]="<<unc_E[k]<<endl;
    cout<<"relative_unc_B["<<k<<"]="<<relative_unc_B[k]<<endl;
    cout<<"relative_unc_E["<<k<<"]="<<relative_unc_E[k]<<endl;
  }




  }
