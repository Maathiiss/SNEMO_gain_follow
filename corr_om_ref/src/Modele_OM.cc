#define _USE_MATH_DEFINES
#include <fstream>
#include <sstream>
#include <fstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TVector.h"
#include "TVector3.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TFractionFitter.h>
#include <TLegend.h>
#include <TParameter.h>

using namespace std;
const int gain_n_bin = 10001; // Precision au 10 000 (me
const float gain_bin_min = 0.5;
const float gain_bin_max = 1.5;
const float gain_bin_width = (gain_bin_max-gain_bin_min)/(gain_n_bin-1);

const int amplitude_n_bin = 1024;
const float amplitude_bin_min = 0e-05;
const float amplitude_bin_max = 2e5;
const float amplitude_bin_width = (amplitude_bin_max-amplitude_bin_min)/(amplitude_n_bin-1);

void MC_Simu(int run_number){

  TFile *file = new TFile(Form("entree/Modele/Modele_OM_%d.root",run_number), "READ");

  double charge;
  int om_num;

  TRandom3 rando;

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_num);
  tree->SetBranchStatus("charge_tree",1);
  tree->SetBranchAddress("charge_tree", &charge);
  
  double time = 0;
  TParameter<double> *param1 = new TParameter<double>("start_time", time);
  param1 = (TParameter<double>*)(file->Get("start_time"));
  time = param1->GetVal();


  TParameter<double> param("start_time",time);

  TFile *newfile_0 = new TFile(Form("entree/Modele/Modele_OM_712_%d.root",run_number), "RECREATE");  
  newfile_0->WriteTObject(&param);
  TFile *newfile_1 = new TFile(Form("entree/Modele/Modele_OM_713_%d.root",run_number), "RECREATE");
  //newfile_1->WriteTObject(&param);  
  TFile *newfile_2 = new TFile(Form("entree/Modele/Modele_OM_714_%d.root",run_number), "RECREATE");  
  //newfile_2->WriteTObject(&param);
  TFile *newfile_3 = new TFile(Form("entree/Modele/Modele_OM_715_%d.root",run_number), "RECREATE");  
  //newfile_3->WriteTObject(&param);
  TFile *newfile_4 = new TFile(Form("entree/Modele/Modele_OM_716_%d.root", run_number), "RECREATE");  
  //newfile_4->WriteTObject(&param);



  TH2D* Modele_712 = new TH2D("Modele_Ref_OM_712", "Modele_Ref_OM_712",gain_n_bin, gain_bin_min - (gain_bin_width/2), gain_bin_max + (gain_bin_width/2),amplitude_n_bin, amplitude_bin_min, amplitude_bin_max);
  TH2D* Modele_713 = new TH2D("Modele_Ref_OM_713", "Modele_Ref_OM_713",gain_n_bin, gain_bin_min - (gain_bin_width/2), gain_bin_max + (gain_bin_width/2),amplitude_n_bin, amplitude_bin_min, amplitude_bin_max);
  TH2D* Modele_714 = new TH2D("Modele_Ref_OM_714", "Modele_Ref_OM_714",gain_n_bin, gain_bin_min - (gain_bin_width/2), gain_bin_max + (gain_bin_width/2),amplitude_n_bin, amplitude_bin_min, amplitude_bin_max);
  TH2D* Modele_715 = new TH2D("Modele_Ref_OM_715", "Modele_Ref_OM_715",gain_n_bin, gain_bin_min - (gain_bin_width/2), gain_bin_max + (gain_bin_width/2),amplitude_n_bin, amplitude_bin_min, amplitude_bin_max);
  TH2D* Modele_716 = new TH2D("Modele_Ref_OM_716", "Modele_Ref_OM_716",gain_n_bin, gain_bin_min - (gain_bin_width/2), gain_bin_max + (gain_bin_width/2),amplitude_n_bin, amplitude_bin_min, amplitude_bin_max);

  std::cout<<tree->GetEntries()<<std::endl;  
  for (int i = 0; i < tree->GetEntries(); i++) {
    if (i%1000000 == 0) {
      std::cout << " entry = " << i << std::endl;
      }    
    double E_kolmo =0;
    tree->GetEntry(i);
    
      if (om_num ==  712 && charge > 0) {
	for(float gain = gain_bin_min; gain <= gain_bin_max; gain += gain_bin_width) {
	  E_kolmo = charge*gain;
	  Modele_712->Fill(gain, E_kolmo);
	}
      }
      else if (om_num ==  713 && charge > 0) {
        for(float gain = gain_bin_min; gain <= gain_bin_max; gain += gain_bin_width) {
          E_kolmo = charge*gain;
          Modele_713->Fill(gain, E_kolmo);
        }
      }
      else if (om_num ==  714 && charge > 0) {
        for(float gain = gain_bin_min; gain <= gain_bin_max; gain += gain_bin_width) {
          E_kolmo = charge*gain;
          Modele_714->Fill(gain, E_kolmo);
        }
      }
      else if (om_num ==  715 && charge > 0) {
        for(float gain = gain_bin_min; gain <= gain_bin_max; gain += gain_bin_width) {
          E_kolmo = charge*gain;
          Modele_715->Fill(gain, E_kolmo);
        }
      }
      else if (om_num ==  716 && charge > 0) {
        for(float gain = gain_bin_min; gain <= gain_bin_max; gain += gain_bin_width) {
          E_kolmo = charge*gain;
          Modele_716->Fill(gain, E_kolmo);
        }
      }

      
  }
  
  newfile_0->cd();
  Modele_712->Write();
  //newfile_0->WriteTObject(&param);
  newfile_0->Close();
 
  newfile_1->cd();
  Modele_713->Write();
  //newfile_1->WriteTObject(&param);
  newfile_1->Close();

  newfile_2->cd();
  Modele_714->Write();
  //newfile_2->WriteTObject(&param);
  newfile_2->Close();

  newfile_3->cd();
  Modele_715->Write();
  //newfile_3->WriteTObject(&param);
  newfile_3->Close();

  newfile_4->cd();
  //newfile_4->WriteTObject(&param);
  Modele_716->Write();
  newfile_4->Close();
  std::cout<<"fin du code";
}


int main(int argc, char const *argv[])
{
  int run_number = std::stoi(argv[1]);
  
  MC_Simu(run_number);
  return 0;
}
