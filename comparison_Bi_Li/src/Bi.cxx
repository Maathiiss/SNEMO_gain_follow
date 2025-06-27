#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstring>
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TApplication.h"
#include "TMultiGraph.h"
#include "TFeldmanCousins.h"
#include "TGaxis.h"
#include "TLeaf.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TParameter.h>
#include <TSpectrum.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TSystem.h>
//#include "/sps/nemo/scratch/granjon/bb_data/sndisplay.cc"


using namespace std; 


void fit_Bi_energy(int run_number){
  double chi2= 0.0, mean_pic= 0.0, mean_pic_error= 0.0,time_bi= 0.0,ndf= 0.0, unix_start_time =0.0;
  int om_number,nb_entries;
  TFile *file_create = new TFile(Form("sortie/Bi_result/result_run_Bi_%d.root",run_number), "RECREATE");
  TTree final_tree("final_tree","");
  final_tree.Branch("Chi2",&chi2);
  final_tree.Branch("om_number",&om_number);
  final_tree.Branch("run_number",&run_number);
  final_tree.Branch("mean_pic",&mean_pic);
  final_tree.Branch("mean_pic_error",&mean_pic_error);
  final_tree.Branch("time_bi",&time_bi);
  final_tree.Branch("nb_entries",&nb_entries);
  final_tree.Branch("ndf",&ndf);

  TFile *file = new TFile(Form("../Bi/root/file_%d.root",run_number), "READ");
  //TFile *file = new TFile(Form("../Bi/root/file_%d/extracted_data.root",run_number), "READ");
  TTree* Result_tree = (TTree*)file->Get("Event");
  vector<int>* num_om = nullptr;
  vector<int>* number_of_kinks_per_track = nullptr;
  vector<double>* charge_elec = nullptr;

  int number_of_electrons;
  gROOT->cd();
  Result_tree->SetBranchStatus("*",0);
  Result_tree->SetBranchStatus("*",0);
  Result_tree->SetBranchStatus("num_om",1);
  Result_tree->SetBranchAddress("num_om", &num_om);
  Result_tree->SetBranchStatus("number_of_kinks_per_track",1);
  Result_tree->SetBranchAddress("number_of_kinks_per_track", &number_of_kinks_per_track);
  Result_tree->SetBranchStatus("number_of_electrons",1);
  Result_tree->SetBranchAddress("number_of_electrons", &number_of_electrons);
  Result_tree->SetBranchStatus("charge_elec",1);
  Result_tree->SetBranchAddress("charge_elec", &charge_elec);
  Result_tree->SetBranchStatus("unix_start_time",1);
  Result_tree->SetBranchAddress("unix_start_time", &unix_start_time);


  std::array<TH1D*, 712> histograms;
  int nb_bins = 200;
  int max_charge = 60000;
  int min_charge = 0;
  for(int i=0; i<712; i++){
    if(i<520){
    histograms[i] = new TH1D(Form("om_%d",i),Form("om_%d",i),nb_bins,min_charge,max_charge);
    
    }
    else{
      histograms[i] = new TH1D(Form("om_%d",i),Form("om_%d",i),nb_bins/2,min_charge,max_charge);
    }
  }
  
  for(int entry=0; entry<Result_tree->GetEntries();entry++){
    Result_tree->GetEntry(entry);
    for(int k=0; k<number_of_electrons; k++){
      if(number_of_kinks_per_track->at(k)==0){
        histograms[num_om->at(k)]->Fill(charge_elec->at(k));
      }
    }
  }
  
  TF1* f_MultipleGaus = new TF1 ("f_MultipleGaus","[0]*(7.08*TMath::Gaus(x[0],[1],[2]) + 1.84*TMath::Gaus(x[0],[1]*(1047.8/975.7),[2]*sqrt((1047.8/975.7))) + 0.44*TMath::Gaus(x[0],([1]*(1059.8/975.7)),[2]*sqrt((1059.8/975.7))))", 0, 60000);
  file_create->cd();
  for(int i=0; i<712; i++){
    //fit                                                                                               
    //cout<<"hist"<<i<<endl;                                                                            
    //cout<<histograms[i][0]->GetEntries()<<endl;
    nb_entries = histograms[i]->GetEntries();
    if(histograms[i]->GetEntries()>500 && histograms[i]->GetMean()>0){
      TFitResultPtr fitResult;      
      for (int j = 0; j < 10; j++) {
        if(j!=0){
	  f_MultipleGaus->SetParLimits(2, 800, 8000);
	  //f_MultipleGaus->SetParLimits(1,histograms[i]->GetMaximumBin()*(max_charge/nb_bins)-histograms[i]->GetMaximumBin()*(max_charge/nb_bins)/5,histograms[i]->GetMaximumBin()*(max_charge/nb_bins)+histograms[i]->GetMaximumBin()*(max_charge/nb_bins)/5);
	  TSpectrum spectrum;
          if (spectrum.Search(histograms[i], 7, "", 0.11) == 2){
            Double_t *xPeaks = spectrum.GetPositionX();
            f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);
	  f_MultipleGaus->SetParLimits(1,xPeaks[0]-xPeaks[0]/5,xPeaks[0]+xPeaks[0]/5);
	  f_MultipleGaus->SetRange(f_MultipleGaus->GetParameter(1)-1.3*f_MultipleGaus->GetParameter(2), f_MultipleGaus->GetParameter(1)+3*f_MultipleGaus->GetParameter(2));
	  }
	  else{
	    Double_t *xPeaks = spectrum.GetPositionX();
            f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);
	    f_MultipleGaus->SetParLimits(1,xPeaks[0]-xPeaks[0]/5,xPeaks[0]+xPeaks[0]/5);
	    f_MultipleGaus->SetRange(f_MultipleGaus->GetParameter(1)-1.3*f_MultipleGaus->GetParameter(2), f_MultipleGaus->GetParameter(1)+3*f_MultipleGaus->GetParameter(2));
	  }
	}
	else{//if first time	  
	  TSpectrum spectrum;                                                                
	  if (spectrum.Search(histograms[i], 7, "", 0.11) == 2){
	    Double_t *xPeaks = spectrum.GetPositionX();
	    f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);		    
	  }
	  else{
	    Double_t *xPeaks = spectrum.GetPositionX();
            f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);
	    cout<<"histo "<<i<<"has " << spectrum.Search(histograms[i], 7, "", 0.11)<<" pics"<<endl;
	  }
	}
        fitResult = histograms[i]->Fit(f_MultipleGaus, "RQ0");
      }
      ndf = nb_bins*(f_MultipleGaus->GetParameter(1)+3*f_MultipleGaus->GetParameter(2)-(f_MultipleGaus->GetParameter(1)-1*f_MultipleGaus->GetParameter(2)))/max_charge; //500 bins de 0 a 4 MeV
     
      TCanvas* canvas = new TCanvas(Form("om_%d",i), Form("om_%d",i), 800, 600);
      if(fitResult == 1){
	//cout<<"om"<<i<<endl;	
        canvas->cd();
        gStyle->SetOptFit(1111);
        histograms[i]->Draw();
	//canvas->SaveAs(Form("/home/granjon/Bi/histo/histo_%d.png",i));                     
	canvas->Write();
        canvas->Close();
	om_number=i;
        chi2 = 10000;
	mean_pic = 0;
	mean_pic_error = 0;
	run_number = run_number;
        time_bi = unix_start_time;
	ndf = 1;
        final_tree.Fill();
        //cout<<"hist "<<i<<"pb fit"<<endl;                                                             
      }
      else{
	//cout<<" om bon"<<i<<endl;
        //cout<<"histo"<<i<<"chi2 "<<f_MultipleGaus->GetChisquare()<<endl;                              
	canvas->cd();
	gStyle->SetOptFit(1111);
	histograms[i]->Draw();
	f_MultipleGaus->Draw("same");
	//histograms[i]->ShowPeaks(5, "", 0.2);
	// gSystem->mkdir(Form("sortie/Comparison_Bi_Li/run_%d",run_number));
	// canvas->SaveAs(Form("sortie/Comparison_Bi_Li/run_%d/histo_%d.png",run_number,i));
	canvas->Write();
	canvas->Close();
	om_number=i;
	chi2 = f_MultipleGaus->GetChisquare();
	mean_pic = f_MultipleGaus->GetParameter(1);
	mean_pic_error = f_MultipleGaus->GetParError(1);
	run_number = run_number;
	time_bi = unix_start_time;
	final_tree.Fill();
	//histograms[i][0]->Write();                                                           

        //cout<<"hist "<<i<<"good"<<endl;                                                               
      }
    }
    else{
      TCanvas* canvas = new TCanvas(Form("om_%d",i), Form("om_%d",i), 800, 600);
      canvas->cd();
      gStyle->SetOptFit(1111);
      histograms[i]->Draw();
      //canvas->SaveAs(Form("/home/granjon/Bi/histo/run_%d/histo_%d.png",run_number,i));
      canvas->Write();
      canvas->Close();
      f_MultipleGaus->SetParameters(0,0,0);
      om_number=i;
      chi2 = 2000;
      ndf = 1;
      run_number = run_number;
      time_bi = unix_start_time;
      mean_pic = 0;
      mean_pic_error = 0;
      final_tree.Fill();      
      //cout<<"hist without entry "<<i<<endl;                                                          
    }
    	histograms[i]->Delete();
  }
  final_tree.Write();
  file_create->Close();
  file->Close();
  delete file;
  delete file_create;

}




void file_merger_Bi(std::vector<int> run){
 TFile *newfile = new TFile(Form("sortie/Bi_result/Merged_Bi_Fit_%d-%d.root", run[0],run[run.size()-1]), "RECREATE");
 double mean_pic= 0.0, mean_pic_error= 0.0, Chi2= 0.0, time_bi= 0.0, ndf= 0.0;
  int om_number, run_number, nb_entries;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("mean_pic", &mean_pic);
  Result_tree.Branch("mean_pic_error", &mean_pic_error);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("time_bi", &time_bi);
  Result_tree.Branch("nb_entries",&nb_entries);
  Result_tree.Branch("ndf",&ndf);

  for (size_t j = 0; j < run.size(); j++) {
    TFile *file1 = new TFile(Form("sortie/Bi_result/result_run_Bi_%d.root",run[j]), "READ");    
    TTree* tree = (TTree*)file1->Get("final_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("om_number",1);
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchStatus("Chi2",1);
    tree->SetBranchAddress("Chi2", &Chi2);
    tree->SetBranchStatus("mean_pic",1);
    tree->SetBranchAddress("mean_pic", &mean_pic);
    tree->SetBranchStatus("mean_pic_error",1);
    tree->SetBranchAddress("mean_pic_error", &mean_pic_error);
    tree->SetBranchStatus("run_number",1);
    tree->SetBranchAddress("run_number", &run_number);
    tree->SetBranchStatus("time_bi",1);
    tree->SetBranchAddress("time_bi", &time_bi);
    tree->SetBranchStatus("nb_entries",1);
    tree->SetBranchAddress("nb_entries", &nb_entries);
    tree->SetBranchStatus("ndf",1);
    tree->SetBranchAddress("ndf", &ndf);
    
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      Result_tree.Fill();
    }
    file1->Close();
    delete file1;
  }
  newfile->cd();
  Result_tree.Write();
  newfile->Close();
  delete newfile;
}




void Evolution_Li_SN_graph(std::vector<int> run, int n_run, int start, int stop,std::vector<int> run_bi,int nb_pics){//comparaison graph corrige et pas corrige
  //plutot que multigraph tu peux dessiner les j pics de l'histo dans le meme canvas
  //Bi part
  int n_run_bi = run_bi.size();
  double mean_pic_value[712][n_run_bi];
  double mean_pic_value_error[712][n_run_bi];
  double mean_pic_vec[712][n_run_bi];
  double Chi2_vec[712][n_run_bi];
  double ndf_vec[712][n_run_bi];
  double time_bi_vec[712][n_run_bi];
  int nb_entries_vec[712][n_run_bi];
  double mean_pic_error_vec[712][n_run_bi];
  double mean_pic_error, mean_pic, time_bi, Chi2_bi, bi_relative, bi_relative_error,mean_pic_bi, mean_pic_bi_error,Chi2, ndf_bi, ndf, time_bi_f, gain_bi, gain_bi_error;
  int run_number_bi, om_number_bi, srun_bi,nb_entries,nb_entries_bi;
  double time_vec_bi[n_run_bi];
  
  TFile *bi_file = new TFile(Form("sortie/Bi_result/Bi_final_%d-%d.root",start,stop),"RECREATE");
  TTree bi_tree("bi_tree","");
  bi_tree.Branch("om_number_bi", &om_number_bi);
  bi_tree.Branch("mean_pic_bi", &mean_pic_bi);
  bi_tree.Branch("mean_pic_bi_error", &mean_pic_bi_error);
  bi_tree.Branch("run_number_bi", &run_number_bi);
  bi_tree.Branch("gain_bi",&gain_bi);
  bi_tree.Branch("gain_bi_error",&gain_bi_error);
  bi_tree.Branch("bi_relative",&bi_relative);
  bi_tree.Branch("bi_relative_error",&bi_relative_error);
  bi_tree.Branch("Chi2_bi",&Chi2_bi);
  bi_tree.Branch("nb_entries_bi",&nb_entries_bi);
  bi_tree.Branch("ndf_bi",&ndf_bi);
  bi_tree.Branch("time_bi_f",&time_bi_f);

  TFile file_bi(Form("sortie/Bi_result/Merged_Bi_Fit_%d-%d.root",run_bi[0],run_bi[run_bi.size()-1]), "READ");
  TTree* tree_bi = (TTree*)file_bi.Get("Result_tree");
  tree_bi->SetBranchStatus("*",0);
  tree_bi->SetBranchStatus("mean_pic",1);
  tree_bi->SetBranchAddress("mean_pic", &mean_pic);
  tree_bi->SetBranchStatus("mean_pic_error",1);
  tree_bi->SetBranchAddress("mean_pic_error", &mean_pic_error);
  tree_bi->SetBranchStatus("run_number",1);
  tree_bi->SetBranchAddress("run_number", &run_number_bi);
  tree_bi->SetBranchStatus("om_number",1);
  tree_bi->SetBranchAddress("om_number", &om_number_bi);
  tree_bi->SetBranchStatus("time_bi",1);
  tree_bi->SetBranchAddress("time_bi", &time_bi);
  tree_bi->SetBranchStatus("Chi2",1);
  tree_bi->SetBranchAddress("Chi2", &Chi2);
  tree_bi->SetBranchStatus("nb_entries",1);
  tree_bi->SetBranchAddress("nb_entries", &nb_entries);
  tree_bi->SetBranchStatus("ndf",1);
  tree_bi->SetBranchAddress("ndf", &ndf);

  int scompteur_bi=0;
  srun_bi = run_bi[0];
  for (int i = 0; i < tree_bi->GetEntries(); i++) {
    tree_bi->GetEntry(i);
    if(i==0){
      time_vec_bi[0]=time_bi;
    }
    if (run_number_bi != srun_bi) {
      srun_bi = run_number_bi;
      scompteur_bi++;
      time_vec_bi[scompteur_bi]=time_bi;
    } 
      mean_pic_vec[om_number_bi][scompteur_bi] = mean_pic;
      mean_pic_error_vec[om_number_bi][scompteur_bi] = mean_pic_error;
      Chi2_vec[om_number_bi][scompteur_bi] = Chi2;      
      nb_entries_vec[om_number_bi][scompteur_bi] = nb_entries;
      ndf_vec[om_number_bi][scompteur_bi] = ndf;
      time_bi_vec[om_number_bi][scompteur_bi] = time_bi;
  }
  
  double xaxis_bi[n_run_bi];
  double xaxis_error_bi[n_run_bi];
  for(int indice = 0; indice<n_run_bi; indice++){
    xaxis_bi[indice]=time_vec_bi[indice];
    xaxis_error_bi[indice]=0;
  }
    std::array<TGraphErrors*, 712> variation_bi_vec;

    for (int i = 0; i < 712; i++) { //om num       
      double yaxis_bi[n_run_bi];
      double yaxis_error_bi[n_run_bi];
      for (int l = 0; l < n_run_bi; l++) {
	if (mean_pic_vec[i][0] > 0.1){
	  yaxis_bi[l] = mean_pic_vec[i][l]/mean_pic_vec[i][0];
	  gain_bi = mean_pic_vec[i][l]/mean_pic_vec[i][0];
	  yaxis_error_bi[l] = sqrt(pow(mean_pic_error_vec[i][l]/mean_pic_vec[i][0],2) + pow(mean_pic_vec[i][l]*mean_pic_error_vec[i][0]/pow(mean_pic_vec[i][0],2),2));
	  gain_bi_error = yaxis_error_bi[l];
	  mean_pic_value[i][l] = mean_pic_vec[i][l]/mean_pic_vec[i][0];
	  mean_pic_value_error[i][l] = sqrt(pow(mean_pic_error_vec[i][l]/mean_pic_vec[i][0],2) + pow(mean_pic_vec[i][l]*mean_pic_error_vec[i][0]/pow(mean_pic_vec[i][0],2),2));
	  om_number_bi = i;
	  mean_pic_bi = mean_pic_vec[i][l];
	  mean_pic_bi_error = mean_pic_error_vec[i][l];
	  run_number_bi = run_bi[l];
	  bi_relative = mean_pic_value[i][l];
	  bi_relative_error = mean_pic_value_error[i][l];
	  Chi2_bi = Chi2_vec[i][l];
	  ndf_bi = ndf_vec[i][l];
	  nb_entries_bi=nb_entries_vec[i][l];
	  time_bi_f = time_bi_vec[i][l];
	  //cout<<nb_entries_bi<<endl;
	  bi_tree.Fill();
	}
	else{
	  gain_bi=0;
	  gain_bi_error=0.;
	  yaxis_bi[l]=0.;
	  yaxis_error_bi[l]=0.;
	  mean_pic_value[i][l] = 0.;
	  mean_pic_value_error[i][l] = 0.;
	  om_number_bi = i;
	  mean_pic_bi = mean_pic_vec[i][l];
	  mean_pic_bi_error = 0;
	  run_number_bi = run_bi[l];
	  bi_relative = 0.;
	  bi_relative_error = 0.;
	  Chi2_bi = Chi2_vec[i][l];	  
	  nb_entries_bi=nb_entries_vec[i][l];
	  ndf_bi = ndf_vec[i][l];
	  time_bi_f = time_bi_vec[i][l];
	  bi_tree.Fill();	  
	}
      }
      variation_bi_vec[i] = new TGraphErrors(n_run_bi, xaxis_bi, yaxis_bi, xaxis_error_bi, yaxis_error_bi);
      variation_bi_vec[i]->SetLineColor(kBlack);
      variation_bi_vec[i]->SetLineWidth(3);
      variation_bi_vec[i]->SetTitle("Bi 1 MeV pic");
    }
    bi_file->cd();
    bi_tree.Write();
    bi_file->Close();
    delete bi_file;
    file_bi.Close();

  //Li part
  double amp_scor[712][nb_pics][n_run]; //om pic run
  double amp_scor_error[712][nb_pics][n_run];
  double amp_cor[712][nb_pics][n_run];
  double amp_cor_error[712][nb_pics][n_run];
  double Amplitude_error, Amplitude, time, Amplitude_corr, Amplitude_uncorr=0.0, Amplitude_corr_error, gain, gain_non_corr, time_li, gain_error;
  int run_number, om_number, pic, bundle;
  double time_vec[n_run];
  double time_li_vec[712][nb_pics][n_run];
  int bundle_ref_vec[712][nb_pics][n_run];
  int srun = start;//mettre le premier run ??
  int scompteur = 0;
  TFile *file_sortie = new TFile(Form("sortie/Comparison_Bi_Li/histo_graph_%d-%d.root",start,stop),"RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Amplitude", &Amplitude);
  Result_tree.Branch("Amplitude_uncorr", &Amplitude_uncorr);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("gain",&gain);
  Result_tree.Branch("gain_error",&gain_error);
  Result_tree.Branch("gain_non_corr",&gain_non_corr);
  Result_tree.Branch("bundle",&bundle);
  Result_tree.Branch("time_li",&time_li);

  TFile file_cor(Form("../Li/sortie/SN_Li/SN_gain_%d-%d.root",start,stop), "READ");
  TTree* tree_cor = (TTree*)file_cor.Get("Result_tree");
  tree_cor->SetBranchStatus("*",0);
  tree_cor->SetBranchStatus("Amplitude",1);
  tree_cor->SetBranchAddress("Amplitude", &Amplitude);
  tree_cor->SetBranchStatus("Amplitude_error",1);
  tree_cor->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree_cor->SetBranchStatus("Amplitude_corr",1);
  tree_cor->SetBranchAddress("Amplitude_corr", &Amplitude_corr);
  tree_cor->SetBranchStatus("Amplitude_corr_error",1);
  tree_cor->SetBranchAddress("Amplitude_corr_error", &Amplitude_corr_error);
  tree_cor->SetBranchStatus("bundle",1);
  tree_cor->SetBranchAddress("bundle", &bundle);  
  tree_cor->SetBranchStatus("run_number",1);
  tree_cor->SetBranchAddress("run_number", &run_number);
  tree_cor->SetBranchStatus("om_number",1);
  tree_cor->SetBranchAddress("om_number", &om_number);
  tree_cor->SetBranchStatus("pic",1);
  tree_cor->SetBranchAddress("pic", &pic);
  tree_cor->SetBranchStatus("time",1);
  tree_cor->SetBranchAddress("time", &time);

  srun = start;
  scompteur = 0;
    
  for (int i = 0; i < tree_cor->GetEntries(); i++) {
    tree_cor->GetEntry(i);
    if(i==0){
      time_vec[0]=time;	   
    }
    if (run_number != srun) {
      scompteur++;
      srun = run_number;
      time_vec[scompteur]=time;
    }   
      amp_cor[om_number][pic-1][scompteur] = Amplitude_corr;
      amp_cor_error[om_number][pic-1][scompteur] = Amplitude_corr_error;
      amp_scor[om_number][pic-1][scompteur] = Amplitude;
      amp_scor_error[om_number][pic-1][scompteur] = Amplitude_error;
      bundle_ref_vec[om_number][pic-1][scompteur] = bundle;
      time_li_vec[om_number][pic-1][scompteur] = time;
  }  

  double norm_amp_cor[712][nb_pics][n_run];
  double norm_amp_cor_error[712][nb_pics][n_run];
  double norm_amp_scor[712][nb_pics][n_run];
  double norm_amp_scor_error[712][nb_pics][n_run];
  double color[5] = {kBlack,kBlue,kGreen+1,kOrange-3,kRed}; //change intensities
  for (int i = 0; i < 712; i++) {//om num
    for (int j = 0; j < nb_pics; j++) {//pic num
      for (int l = 0; l < n_run; l++) {//run num
        if (amp_scor[i][j][0] > 0.1 && amp_cor[i][j][0] > 0.1 && amp_scor[i][j][l]>0.1 && amp_cor[i][j][0]>0.1) {
          norm_amp_cor[i][j][l] = amp_cor[i][j][l]/amp_cor[i][j][0]; //divide the run l by the first one to normalise
          norm_amp_cor_error[i][j][l] = sqrt(pow(amp_cor_error[i][j][l]/amp_cor[i][j][0],2) + pow(amp_cor[i][j][l]*amp_cor_error[i][j][0]/pow(amp_cor[i][j][0],2),2));
	  norm_amp_scor[i][j][l] = amp_scor[i][j][l]/amp_scor[i][j][0]; //do the same without correction
          norm_amp_scor_error[i][j][l] = sqrt(pow(amp_scor_error[i][j][l]/amp_scor[i][j][0],2) + pow(amp_scor[i][j][l]*amp_scor_error[i][j][0]/pow(amp_scor[i][j][0],2),2));
	  gain = norm_amp_cor[i][j][l];
	  gain_non_corr = norm_amp_scor[i][j][l];
	  gain_error = norm_amp_cor_error[i][j][l];
	  run_number = run[l];
	  Amplitude = amp_cor[i][j][l];
	  Amplitude_uncorr = amp_scor[i][j][l];
	  bundle = bundle_ref_vec[i][j][l];
	  time_li = time_li_vec[i][j][l];
	  om_number = i;
	  pic =j+1;
	  Result_tree.Fill();
        }
	else{
          gain = 0.;
	  gain_error=0.;
          gain_non_corr = 0.;
          run_number = run[l];
          Amplitude = amp_cor[i][j][l];
	  Amplitude_uncorr = amp_scor[i][j][l];
          bundle = bundle_ref_vec[i][j][l];
          om_number = i;
          pic =j+1;
	  time_li = time_li_vec[i][j][l];		    
          Result_tree.Fill();
	}
      }
    }
  }
  double xaxis[n_run];
  double xaxis_error[n_run];
  for(int indice = 0; indice<n_run; indice++){
    xaxis[indice]=time_vec[indice];
    xaxis_error[indice]=0;
    //std::cout<<" time_vec "<< time_vec[indice]<<std::endl;
  }
    double min_x = *std::min_element(xaxis,xaxis+n_run);
    double max_x = *std::max_element(xaxis,xaxis+n_run);
    cout<<"max = "<<max_x<<endl;
    cout<<"min = "<<min_x<<endl;    
  for (int i = 0; i < 712; i++) { //om num
    TMultiGraph *multiGraph = new TMultiGraph(Form("om_%d",i),Form("om_%d",i));
    multiGraph->GetYaxis()->SetTitle("Relative variation");
    multiGraph->GetXaxis()->SetTitle("date");
    multiGraph->GetYaxis()->SetRangeUser(0.9,1.1);
    multiGraph->GetXaxis()->SetTimeDisplay(1);    
    for (int j = 0; j < nb_pics; j++) { // pic num
      //TCanvas* c = new TCanvas("","",1200, 600);                        
      double yaxis[n_run];
      double yaxis_error[n_run];
      double syaxis[n_run];
      double syaxis_error[n_run];
      for (int l = 0; l < n_run; l++) {
        if (norm_amp_cor[i][j][l] >0.1 && norm_amp_scor[i][j][l] > 0.1) {
          yaxis[l] = norm_amp_cor[i][j][l];
          yaxis_error[l] = norm_amp_cor_error[i][j][l];	
          syaxis[l] = norm_amp_scor[i][j][l];
          syaxis_error[l] = norm_amp_scor_error[i][j][l];
        }
	else{
	  yaxis[l]=0;
	  yaxis_error[l]=0;
	  syaxis[l]=0;
	  syaxis_error[l]=0;
	}	
      }
      TGraphErrors *variation_scor = new TGraphErrors(n_run, xaxis, syaxis, xaxis_error, syaxis_error);
      TGraphErrors *variation_cor = new TGraphErrors(n_run, xaxis, yaxis, xaxis_error, yaxis_error);
      variation_scor->SetLineColor(color[j+1]);
      variation_scor->SetLineStyle(2);
      variation_scor->SetLineWidth(2);
      //variation_scor->SetTitle(Form("Intensity %d", j+1));
      variation_cor->SetLineColor(color[j+1]);
      variation_cor->SetLineWidth(2);
      variation_cor->SetTitle(Form("gain variation intensity %d", j+1));
      variation_cor->Draw("APL"/*"same"*/);
      variation_scor->SetTitle(Form("charge variation intensity %d", j+1));
      variation_scor->Draw("APL"/*"same"*/);
      multiGraph->Add(variation_cor);
      multiGraph->Add(variation_scor);
      // delete variation_scor;
      // delete variation_cor;
    }//boucle j
    multiGraph->Add(variation_bi_vec[i]);    
    file_sortie->cd();
    multiGraph->GetXaxis()->SetLimits(min_x-10000,max_x+10000);
    //cout<<"multi"<<multiGraph<<endl;
    multiGraph->Write();
    TCanvas* c = new TCanvas("","",1200, 600);
    c->cd();
    multiGraph->Draw("APL");
    c->SaveAs(Form("sortie/Comparison_Bi_Li/histo_save/om_%d.png",i));
    c->Close();
    //delete multiGraph;
    variation_bi_vec[i]->Delete();
  }//boucle i
  Result_tree.Write();
  file_sortie->Close();
  delete file_sortie;
  file_cor.Close();
}





void comparaison_Bi_Li(int start, int stop, int stop_bi){
  double gain, gain_non_corr, bi_relative, bi_relative_error, li_relative, li_relative_error, diff_li_bi,Chi2,Chi2_bi, li_relative_non_corr, diff_li_bi_non_corr, shift_Li_with_1, shift_Bi_with_1, ndf_bi, time_li, time_bi,time_bi_f, mean_pic_bi, li_amplitude, li_amplitude_uncorr=0.0,Amplitude, diff_bi_start_stop, diff_li_start_stop, Amplitude_uncorr=0.0, gain_bi=0.0, gain_bi_error=0.0, gain_error=0.0;
  int bundle, run_number, pic, om_number, om_number_bi, om_number_li, pic_li,bundle_li, run_number_bi, run_number_li, nb_entries_bi;
  bool is_before;

  TFile *file_final = new TFile(Form("sortie/Comparison_Bi_Li/final_comparison_%d-%d.root",start,stop),"RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number_bi", &om_number_bi);
  Result_tree.Branch("om_number_li", &om_number_li);
  Result_tree.Branch("pic_li", &pic_li);
  Result_tree.Branch("bundle_li", &bundle_li);
  Result_tree.Branch("run_number_li", &run_number_li);
  Result_tree.Branch("run_number_bi", &run_number_bi);
  Result_tree.Branch("bi_relative",&bi_relative);
  Result_tree.Branch("bi_relative_error",&bi_relative_error);
  Result_tree.Branch("li_relative",&li_relative);
  Result_tree.Branch("li_relative_error",&li_relative_error);
  Result_tree.Branch("li_relative_non_corr",&li_relative_non_corr);
  Result_tree.Branch("diff_li_bi",&diff_li_bi);
  Result_tree.Branch("is_before",&is_before);
  Result_tree.Branch("Chi2_bi",&Chi2_bi);
  Result_tree.Branch("nb_entries_bi",&nb_entries_bi);
  Result_tree.Branch("diff_li_bi_non_corr",&diff_li_bi_non_corr);
  Result_tree.Branch("shift_with_Li_1",&shift_Li_with_1);
  Result_tree.Branch("shift_with_Bi_1",&shift_Bi_with_1);
  Result_tree.Branch("ndf_bi",&ndf_bi);
  Result_tree.Branch("time_li",&time_li);
  Result_tree.Branch("time_bi",&time_bi);
  Result_tree.Branch("li_amplitude",&li_amplitude);
  Result_tree.Branch("li_amplitude_uncorr",&li_amplitude_uncorr);
  Result_tree.Branch("mean_pic_bi",&mean_pic_bi);
  //Result_tree.Branch("gain_bi",&gain_bi);
  Result_tree.Branch("diff_bi_start_stop",&diff_bi_start_stop);
  Result_tree.Branch("diff_li_start_stop",&diff_li_start_stop);  
  
  
  TFile file_cor(Form("sortie/Comparison_Bi_Li/histo_graph_%d-%d.root",start,stop), "READ");    
  TTree* tree_cor = (TTree*)file_cor.Get("Result_tree");
  tree_cor->SetBranchStatus("*",0);
  tree_cor->SetBranchStatus("gain",1);
  tree_cor->SetBranchAddress("gain", &gain);
  tree_cor->SetBranchStatus("gain_error",1);
  tree_cor->SetBranchAddress("gain_error", &gain_error);
  tree_cor->SetBranchStatus("gain_non_corr",1);
  tree_cor->SetBranchAddress("gain_non_corr", &gain_non_corr);
  tree_cor->SetBranchStatus("bundle",1);
  tree_cor->SetBranchAddress("bundle", &bundle);
  tree_cor->SetBranchStatus("run_number",1);
  tree_cor->SetBranchAddress("run_number", &run_number);
  tree_cor->SetBranchStatus("om_number",1);
  tree_cor->SetBranchAddress("om_number", &om_number);
  tree_cor->SetBranchStatus("pic",1);
  tree_cor->SetBranchAddress("pic", &pic);
  tree_cor->SetBranchStatus("time_li",1);
  tree_cor->SetBranchAddress("time_li", &time_li);
  tree_cor->SetBranchStatus("Amplitude",1);
  tree_cor->SetBranchAddress("Amplitude", &Amplitude);
  tree_cor->SetBranchStatus("Amplitude_uncorr",1);
  tree_cor->SetBranchAddress("Amplitude_uncorr", &Amplitude_uncorr);
 
  TFile file_bi(Form("sortie/Bi_result/Bi_final_%d-%d.root",start,stop), "READ");
  TTree* tree = (TTree*)file_bi.Get("bi_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("bi_relative",1);
  tree->SetBranchAddress("bi_relative", &bi_relative);
  tree->SetBranchStatus("bi_relative_error",1);
  tree->SetBranchAddress("bi_relative_error", &bi_relative_error);
  tree->SetBranchStatus("run_number_bi",1);
  tree->SetBranchAddress("run_number_bi", &run_number_bi);
  tree->SetBranchStatus("om_number_bi",1);
  tree->SetBranchAddress("om_number_bi", &om_number_bi);
  tree->SetBranchStatus("Chi2_bi",1);
  tree->SetBranchAddress("Chi2_bi", &Chi2);
  tree->SetBranchStatus("nb_entries_bi",1);
  tree->SetBranchAddress("nb_entries_bi", &nb_entries_bi);
  tree->SetBranchStatus("ndf_bi",1);
  tree->SetBranchAddress("ndf_bi", &ndf_bi);
  tree->SetBranchStatus("time_bi_f",1);
  tree->SetBranchAddress("time_bi_f", &time_bi_f);
  tree->SetBranchStatus("mean_pic_bi",1);
  tree->SetBranchAddress("mean_pic_bi", &mean_pic_bi);
  tree->SetBranchStatus("gain_bi",1);
  tree->SetBranchAddress("gain_bi", &gain_bi);
  tree->SetBranchStatus("gain_bi_error",1);
  tree->SetBranchAddress("gain_bi_error", &gain_bi_error);

 for (int j = 0; j < tree_cor->GetEntries(); j++) {
   tree_cor->GetEntry(j);
   diff_li_start_stop = 0;
   diff_bi_start_stop = 0;
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      om_number_li = om_number;
      om_number_bi = om_number_bi;
      if(om_number_li==om_number_bi){
	pic_li = pic;
	bundle_li=bundle;
	run_number_li = run_number;
	run_number_bi = run_number_bi;
	bi_relative = bi_relative;
	li_relative = gain;
	//cout<<"run num bi stop "<<run_number_bi<<" "<<stop_bi<<"run num li "<<run_number_li<<" "<<stop<<endl;
	if(run_number_bi==stop_bi && run_number_li==stop){
	  diff_li_start_stop = li_relative;
	}
        if(run_number_bi==stop_bi && run_number_li==stop){
	  diff_bi_start_stop = bi_relative;
	}
	shift_Li_with_1 = fabs(1-gain);
	shift_Bi_with_1 = fabs(1-bi_relative);
	li_relative_non_corr = gain_non_corr;
	Chi2_bi = Chi2;
	nb_entries_bi = nb_entries_bi;
	ndf_bi = ndf_bi;
	time_li = time_li;
	time_bi = time_bi_f;
	mean_pic_bi = mean_pic_bi;
	gain_bi=gain_bi;
	bi_relative_error = gain_bi_error;
	li_amplitude = Amplitude;
	li_relative_error = gain_error;
	li_amplitude_uncorr = Amplitude_uncorr;
	if(bi_relative!=0){
	diff_li_bi = 100*(li_relative-bi_relative)/bi_relative;	
	diff_li_bi_non_corr = 100*(li_relative_non_corr-bi_relative)/bi_relative;
	}
	else{
	  diff_li_bi=0;
	  diff_li_bi_non_corr = 0;
	}
	if(run_number_bi==run_number+2){
          is_before = true;
          Result_tree.Fill();
        }	
	else if(run_number_bi==run_number-6 || (run_number_bi== 1716 && run_number==1723)){
	  is_before = false;
	  Result_tree.Fill();
	}
      }
    }
 }
 
 file_final->cd();
 Result_tree.Write();
 file_final->Close();
 delete file_final;
 file_cor.Close();
 file_bi.Close();
}
 




void each_OM_mean(int start, int stop){
  double diff_li_bi = 0.0, Chi2_bi = 0.0, li_relative = 0.0, li_relative_error=0.0, bi_relative = 0.0, bi_relative_error=0.0, ndf_bi = 0.0, Mean = 0.0, RMS = 0.0, variation_li = 0.0, Mean_li = 0.0, variation_bi = 0.0,Mean_bi = 0.0, diff_li_start_stop = 0.0, diff_bi_start_stop = 0.0, min_somme_li_bi_save = 0.0, min_somme_li_bi_save_test = 0.0, li_amplitude = 0.0, li_amplitude_uncorr=0.0, mean_pic_bi=0.0, best_charge_diff = 0.0, charge_diff = 0.0, bi_amplitude = 0.0, li_amplitude1 = 0.0, best_li_amplitude1 = 0.0, min_charge_diff=0.0, max_li_value = 0.0, min_li_value = 0.0, Mean_derivative_li = 0.0, variation_derivative_li=0.0, time_li=0.0, time_bi=0.0;
  int pic_li, nb_entries_bi, om_number_li, pic_li_after, om_number_li_after, bundle_li, run_number_li, run_number_bi, om_number_bi, pic_associated, bundle_li1;
  bool li_drop;

  //sncalo = new sndisplay::calorimeter ("sndiplay_dis",1);
  TFile* file = TFile::Open(Form("sortie/Comparison_Bi_Li/final_comparison_%d-%d.root",start,stop), "READ");
  TTree* Result_tree = (TTree*)file->Get("Result_tree");
  Result_tree->SetBranchStatus("*",0);
  Result_tree->SetBranchStatus("diff_li_bi",1);
  Result_tree->SetBranchAddress("diff_li_bi",&diff_li_bi);
  Result_tree->SetBranchStatus("pic_li",1);
  Result_tree->SetBranchAddress("pic_li",&pic_li);
  Result_tree->SetBranchStatus("Chi2_bi",1);
  Result_tree->SetBranchAddress("Chi2_bi",&Chi2_bi);
  Result_tree->SetBranchStatus("nb_entries_bi",1);
  Result_tree->SetBranchAddress("nb_entries_bi",&nb_entries_bi);
  Result_tree->SetBranchStatus("ndf_bi",1);
  Result_tree->SetBranchAddress("ndf_bi",&ndf_bi);
  Result_tree->SetBranchStatus("nb_entries_bi",1);
  Result_tree->SetBranchAddress("nb_entries_bi",&nb_entries_bi);
  Result_tree->SetBranchStatus("li_relative",1);
  Result_tree->SetBranchAddress("li_relative",&li_relative);
  Result_tree->SetBranchStatus("li_relative_error",1);
  Result_tree->SetBranchAddress("li_relative_error",&li_relative_error);
  Result_tree->SetBranchStatus("bi_relative",1);
  Result_tree->SetBranchAddress("bi_relative",&bi_relative);
  Result_tree->SetBranchStatus("bi_relative_error",1);
  Result_tree->SetBranchAddress("bi_relative_error",&bi_relative_error);
  Result_tree->SetBranchStatus("time_bi",1);
  Result_tree->SetBranchAddress("time_bi",&time_bi);
  Result_tree->SetBranchStatus("om_number_li",1);
  Result_tree->SetBranchAddress("om_number_li",&om_number_li);
  Result_tree->SetBranchStatus("om_number_bi",1);
  Result_tree->SetBranchAddress("om_number_bi",&om_number_bi);
  Result_tree->SetBranchStatus("run_number_li",1);
  Result_tree->SetBranchAddress("run_number_li",&run_number_li);
  Result_tree->SetBranchStatus("run_number_bi",1);
  Result_tree->SetBranchAddress("run_number_bi",&run_number_bi);
  Result_tree->SetBranchStatus("bundle_li",1);
  Result_tree->SetBranchAddress("bundle_li",&bundle_li);
  Result_tree->SetBranchStatus("diff_li_start_stop",1);
  Result_tree->SetBranchAddress("diff_li_start_stop",&diff_li_start_stop);
  Result_tree->SetBranchStatus("diff_bi_start_stop",1);
  Result_tree->SetBranchAddress("diff_bi_start_stop",&diff_bi_start_stop);
  Result_tree->SetBranchStatus("li_amplitude",1);
  Result_tree->SetBranchAddress("li_amplitude",&li_amplitude);
  Result_tree->SetBranchStatus("li_amplitude_uncorr",1);
  Result_tree->SetBranchAddress("li_amplitude_uncorr",&li_amplitude_uncorr);
  Result_tree->SetBranchStatus("time_li",1);
  Result_tree->SetBranchAddress("time_li",&time_li);
  Result_tree->SetBranchStatus("mean_pic_bi",1);
  Result_tree->SetBranchAddress("mean_pic_bi",&mean_pic_bi);

  
  TFile *file_create = new TFile(Form("sortie/Comparison_Bi_Li/each_OM_save_%d-%d.root",start,stop), "RECREATE");
  TTree final_tree("Result_tree","");
  final_tree.Branch("Mean",&Mean);
  final_tree.Branch("RMS",&RMS);
  final_tree.Branch("pic_li",&pic_li_after);
  final_tree.Branch("om_number_li",&om_number_li_after);
  final_tree.Branch("bundle_li",&bundle_li1);
  final_tree.Branch("RMS_li",&variation_li);
  final_tree.Branch("Mean_derivative_li",&Mean_derivative_li);
  final_tree.Branch("variation_derivative_li",&variation_derivative_li);
  final_tree.Branch("max_li_value",&max_li_value);
  final_tree.Branch("min_li_value",&min_li_value);
  final_tree.Branch("Mean_li",&Mean_li);
  final_tree.Branch("diff_li_start_stop",&diff_li_start_stop);
  final_tree.Branch("li_drop",&li_drop);
  final_tree.Branch("variation_bi",&variation_bi);
  final_tree.Branch("Mean_bi",&Mean_bi);
  final_tree.Branch("diff_bi_start_stop",&diff_bi_start_stop);
  final_tree.Branch("pic_associated",&pic_associated); //which intensity is associated with 
  final_tree.Branch("min_somme_li_bi_save",&min_somme_li_bi_save);
  final_tree.Branch("min_somme_li_bi_save_test",&min_somme_li_bi_save_test);
  final_tree.Branch("best_charge_diff",&best_charge_diff); // charge diff with better result in relative part (li-bi compare to 1 diff)
  final_tree.Branch("charge_diff",&charge_diff);
  final_tree.Branch("min_charge_diff",&min_charge_diff); //min charge diff in terms of absolute comparison between Li and Bi charge valeur (ex : 20000-40000)
  final_tree.Branch("li_amplitude",&li_amplitude1);
  final_tree.Branch("li_amplitude_uncorr",&li_amplitude_uncorr);
  final_tree.Branch("best_li_amplitude",&best_li_amplitude1); //best in terms of relative results
  final_tree.Branch("bi_amplitude",&bi_amplitude); 
  
  std::array<std::array<TH1D*,4>, 712> histograms; //712 OM et 4 pics
  std::array<std::array<TH1D*,4>, 712> histograms_variation; //712 OM et 4 pics
  std::array<std::array<TH1D*,4>, 712> derivate_li; //712 OM et 4 pics
  std::array<std::array<std::vector<double>,4>, 712> RMS_Mean_Li_per_OM; //712 OM et 4 pics
  std::array<std::array<std::vector<double>,4>, 712> time_li_per_OM; //712 OM et 4 pics
  std::array<std::vector<double>, 712> time_bi_per_OM;
  std::array<std::array<std::vector<double>,4>, 712> error_li_relative; //712 OM et 4 pics
  std::array<std::vector<double>, 712> bi_relative_vec_draw;
  std::array<std::vector<double>, 712> bi_relative_vec_draw_error;
  
  std::array<TH1D*,712> histograms_bi;
  double diff_li_start_stop_vec[712][4] = {};
  double somme_diff_li_bi_vec[712][4] = {};
  double somme_charge_li_bi_vec[712][4] = {};
  double diff_bi_start_stop_vec[712] = {};
  double bundle_li_vec[712] = {};
  double li_amplitude_vec[712][4] = {};
  double li_amplitude_uncorr_vec[712][4] = {};
  double bi_amplitude_vec[712] = {};
  vector<vector<double>> bi_amplitude_vec_size(712);
  double best_pic[712] = {};
  double best_charge_pic[712] = {};
  
  for(int i=0;i<712;i++){
    histograms_bi[i] = new TH1D(Form("om_%d",i),Form("om_%d",i),200,-100,100);
    for (size_t j = 0; j <4; j++){
      histograms[i][j] = new TH1D(Form("om_%d_pic_%lu",i,j+1),Form("om_%d_pic_%lu",i,j+1),200,-100,100);
      histograms_variation[i][j] = new TH1D(Form("om_var_%d_pic_%lu",i,j+1),Form("om_var_%d_pic_%lu",i,j+1),2000,-100,100);      
    }
  }
  
  for (int j = 0; j < Result_tree->GetEntries(); j++) {
    Result_tree->GetEntry(j);
    // if(diff_li_start_stop!=0){
    //   diff_li_start_stop_vec[om_number_li][pic_li-1] = diff_li_start_stop;
    // }
    // if(diff_bi_start_stop!=0){
    //   diff_bi_start_stop_vec[om_number_bi] = diff_bi_start_stop;
    // }
    if(/*Chi2_bi!=2000 && nb_entries_bi>50 && Chi2_bi/ndf_bi<3 && */bi_relative>0.1 && pic_li==3){
      //we choose 1 intensity to have only 1 time Bi points
      bi_relative_vec_draw[om_number_bi].push_back(bi_relative);
      bi_relative_vec_draw_error[om_number_bi].push_back(bi_relative_error);
      time_bi_per_OM[om_number_bi].push_back(time_bi-25 * 365.25 * 86400);
    }
    if(li_relative!=0){
      histograms_variation[om_number_li][pic_li-1]->Fill(li_relative);
      RMS_Mean_Li_per_OM[om_number_li][pic_li-1].push_back(li_relative);
      time_li_per_OM[om_number_li][pic_li-1].push_back(time_li-25 * 365.25 * 86400);
      error_li_relative[om_number_li][pic_li-1].push_back(li_relative_error);
      li_amplitude_vec[om_number_li][pic_li-1] += li_amplitude;
      li_amplitude_uncorr_vec[om_number_li][pic_li-1] += li_amplitude_uncorr;
      bundle_li_vec[om_number_li] = bundle_li;
      if(/*Chi2_bi!=2000 && nb_entries_bi>50 && Chi2_bi/ndf_bi<3 && */bi_relative>0.1){
	histograms_bi[om_number_bi]->Fill(bi_relative);
	histograms[om_number_li][pic_li-1]->Fill(diff_li_bi);
	//histograms[om_number_li][pic_li-1]->Fill(abs((li_amplitude-mean_pic_bi)/mean_pic_bi));
	somme_diff_li_bi_vec[om_number_li][pic_li-1]+=abs(diff_li_bi)/*diff_li_bi*diff_li_bi/(100*100)*/;
	somme_charge_li_bi_vec[om_number_li][pic_li-1]+=abs((li_amplitude-mean_pic_bi)/mean_pic_bi);	
	bi_amplitude_vec[om_number_bi] += mean_pic_bi;
	bi_amplitude_vec_size[om_number_bi].push_back(mean_pic_bi);
	//cout<<"om num bi "<<om_number_bi<<" value "<<mean_pic_bi<<endl;
	if(run_number_bi==1415 && run_number_li==1417){
	  diff_li_start_stop_vec[om_number_li][pic_li-1] = li_relative;
	  diff_bi_start_stop_vec[om_number_bi]=bi_relative;
	}
      }
    }
  }

  std::vector<std::vector<bool>> li_drop_vec(712, std::vector<bool>(4, false));
  double color[4] = {kBlue,kGreen+1,kOrange-3,kRed};
  std::array<std::array<std::vector<double>,4>, 712> li_derivative_vec; 
  for(int i=0; i<712; i++){
    TMultiGraph *mg = new TMultiGraph();
    bundle_li1 = bundle_li_vec[i];
    if(histograms_bi[i]->GetEntries()>0){
      bi_amplitude = bi_amplitude_vec[i]/histograms_bi[i]->GetEntries();
    }
    else{
      bi_amplitude = 0;
    }
    charge_diff = 0;
    best_charge_diff = 0;
    om_number_li_after = i;
    variation_bi = histograms_bi[i]->GetRMS();
    Mean_bi = histograms_bi[i]->GetMean();
    diff_bi_start_stop = diff_bi_start_stop_vec[i];
    pic_associated = -1;
    double min_somme_li_bi = 1000;
    double min_charge_diff_temp = 1000;
    int j_store=0;
    for (size_t j = 0; j <4; j++){
      derivate_li[i][j] = new TH1D(Form("om_derivative_%d_pic_%lu",i,j+1),Form("om_derivative_%d_pic_%lu",i,j+1),1000,0,0);
      if(min_somme_li_bi>somme_diff_li_bi_vec[i][j]/histograms[i][j]->GetEntries() && somme_diff_li_bi_vec[i][j]>0.0001 && somme_diff_li_bi_vec[i][j]<10e3 && histograms[i][j]->GetEntries() > 0){
	min_somme_li_bi = somme_diff_li_bi_vec[i][j]/histograms[i][j]->GetEntries();
	pic_associated = j+1;
	best_charge_diff = somme_charge_li_bi_vec[i][j]/histograms[i][j]->GetEntries();
	best_li_amplitude1 = li_amplitude_vec[i][j]/histograms_variation[i][j]->GetEntries();
      }
      if(min_charge_diff_temp>somme_charge_li_bi_vec[i][j]/histograms[i][j]->GetEntries()){
	min_charge_diff_temp = somme_charge_li_bi_vec[i][j]/histograms[i][j]->GetEntries();
	j_store = j+1;
      }      
    }
    for (size_t j = 0; j <4; j++){
      li_amplitude1 = li_amplitude_vec[i][j]/histograms_variation[i][j]->GetEntries();
      li_amplitude_uncorr = li_amplitude_uncorr_vec[i][j]/histograms_variation[i][j]->GetEntries();
      pic_li_after = j+1;
      diff_li_start_stop = diff_li_start_stop_vec[i][j];
      Mean = histograms[i][j]->GetMean();
      RMS = histograms[i][j]->GetRMS();	
      // variation_li = histograms_variation[i][j]->GetRMS();
      // Mean_li = histograms_variation[i][j]->GetMean();
      //RMS_Mean_Li_per_OM contains every li_relative
      li_drop=false;
      // cout<<" i ="<<i<<" j= "<<j<<endl;
      // cout<<" size = "<<RMS_Mean_Li_per_OM[i][j].size()<<endl;
      for(int k=0; k<RMS_Mean_Li_per_OM[i][j].size(); k++){
      	if(k>0){
      	  //double drop_li = abs(RMS_Mean_Li_per_OM[i][j][k]-RMS_Mean_Li_per_OM[i][j][k-1])/(time_li_per_OM[i][j][k]-time_li_per_OM[i][j][k-1]);
      	  double drop_li = abs(RMS_Mean_Li_per_OM[i][j][k]-RMS_Mean_Li_per_OM[i][j][k-1]);
      	  li_derivative_vec[i][j].push_back(drop_li);
      	  derivate_li[i][j]->Fill(drop_li);
	  // if(i==177 && j==3){
	  //   cout<<"drop = "<<drop_li<<endl;
	  // }
      	  if(drop_li>0.02){
	    // if(i==177 && j==3){
	    //   cout<<"DROP = "<<drop_li<<endl;
	    // }
      	    li_drop=true;
      	    li_drop_vec[i][j]=true;
      	  }
      	}
      }
      Mean_derivative_li = TMath::Mean(li_derivative_vec[i][j].size(), &li_derivative_vec[i][j][0]);
      variation_derivative_li = TMath::StdDev(li_derivative_vec[i][j].size(), &li_derivative_vec[i][j][0]);
      
      Mean_li = TMath::Mean(RMS_Mean_Li_per_OM[i][j].size(), &RMS_Mean_Li_per_OM[i][j][0]);
      variation_li = TMath::StdDev(RMS_Mean_Li_per_OM[i][j].size(), &RMS_Mean_Li_per_OM[i][j][0]);
      if(RMS_Mean_Li_per_OM[i][j].size()){
      max_li_value = *std::max_element(RMS_Mean_Li_per_OM[i][j].begin(),RMS_Mean_Li_per_OM[i][j].end());
      min_li_value = *std::min_element(RMS_Mean_Li_per_OM[i][j].begin(),RMS_Mean_Li_per_OM[i][j].end());
      }
      if(somme_charge_li_bi_vec[i][j]>0.0001 && histograms[i][j]->GetEntries() > 0){
	charge_diff = somme_charge_li_bi_vec[i][j]/histograms[i][j]->GetEntries();
	min_charge_diff = min_charge_diff_temp;
	//sncalo->setcontent(i,j_store);
	best_charge_pic[i] = j_store;
      }
      else{
	charge_diff=-10.0;
	min_charge_diff=0;
      }
      if(j+1 == pic_associated){
	min_somme_li_bi_save = somme_diff_li_bi_vec[i][j]/histograms[i][j]->GetEntries();
	best_pic[i] = pic_associated;
	best_charge_diff = best_charge_diff;
      }		  
      else{
	min_somme_li_bi_save = 0;
      }
      min_somme_li_bi_save_test = somme_diff_li_bi_vec[i][j]/histograms[i][j]->GetEntries();
    
      final_tree.Fill();
      file_create->cd();      
      derivate_li[i][j]->Write();
      std::vector<double> x(RMS_Mean_Li_per_OM[i][j].size(), 0.0);
      TGraphErrors *g1 = new TGraphErrors(RMS_Mean_Li_per_OM[i][j].size(), time_li_per_OM[i][j].data() ,RMS_Mean_Li_per_OM[i][j].data(),x.data(),error_li_relative[i][j].data());
      g1->SetLineColor(color[j]);
      if(j+1 == j_store){
	g1->SetLineStyle(1);
	g1->SetLineWidth(2);
      }
      else{
	//g1->SetLineStyle(2);
      }
      //g1->SetLineWidth(2);
      g1->SetMarkerColor(color[j]);
      mg->Add(g1);
    }
    std::vector<double> x_bi(bi_relative_vec_draw[i].size(), 0.0);
    TGraphErrors *g2 = new TGraphErrors(bi_relative_vec_draw[i].size(), time_bi_per_OM[i].data() ,bi_relative_vec_draw[i].data(),x_bi.data(),bi_relative_vec_draw_error[i].data());
    g2->SetLineColor(kBlack);
    g2->SetLineStyle(1);
    g2->SetLineWidth(2);
    mg->Add(g2);
    TAxis *xAxis = mg->GetXaxis();
    xAxis->SetTimeDisplay(1);
    xAxis->SetTimeFormat("%Y-%m-%d");
    xAxis->SetTitle("time (date)");
    mg->GetYaxis()->SetTitle("Relative variation");
    TCanvas *c1 = new TCanvas("c1", "Canvas Example", 800, 600);
    c1->cd();
    mg->Draw("APL");
    c1->SaveAs(Form("sortie/Comparison_Bi_Li/new_every_OM/om_%d.root",i));
    c1->Close();

    //    sncalo->settext(i, Form("%d", i));
  
  }
  file_create->cd();
  final_tree.Write();


  


  bool drop, has_Bi, has_Li, is_stable, is_perfect, decrease, li_pic_drop, best_pic_li;
  int om_num, pic;
  double Bi_and_Li_close=0.0;
  TFile *file_tree = new TFile("sortie/Comparison_Bi_Li/summary_all_bool.root", "RECREATE");
  TTree end_tree("Result_tree","");
  end_tree.Branch("om_number",&om_num);
  end_tree.Branch("pic_li",&pic);
  end_tree.Branch("best_pic_li",&best_pic_li);
  end_tree.Branch("has_Bi",&has_Bi);
  end_tree.Branch("has_Li",&has_Li);
  end_tree.Branch("Bi_and_Li_close",&Bi_and_Li_close);
  end_tree.Branch("li_pic_drop",&li_pic_drop);
  end_tree.Branch("drop",&drop);
  end_tree.Branch("is_stable",&is_stable);
  end_tree.Branch("is_perfect",&is_perfect);
  end_tree.Branch("decrease",&decrease);
  
  
  for(int i=0; i<520; i++){    
    om_num=i;
    has_Bi=true;
    is_perfect=false;
    int compteur_li_drop=0;
    Bi_and_Li_close=0.0;
    if(bi_amplitude_vec_size[i].size()==0){
      has_Bi=false;
    }
    for(int j=0; j<4; j++){
      //cout<<" i = "<<i<< " j = "<<j<<endl;
      pic=j+1;
      has_Li=true;
      drop=false;
      li_pic_drop=false;
      best_pic_li=false;
      is_stable=true;
      decrease=false;
      if(j+1==best_charge_pic[i]){
	best_pic_li=true;
      }
      if(RMS_Mean_Li_per_OM[i][j].size()==0){
  	has_Li=false;
	end_tree.Fill();	      
  	continue;
      }
      //if(somme_diff_li_bi_vec[i][j]/histograms[i][j]->GetEntries()>4){
      //we need less than 4% of differences
      if(li_drop_vec[i][j]==true){
	li_pic_drop=true;
	compteur_li_drop++;
	//cout<<"OM = "<<i<<endl;
	if(compteur_li_drop==4){
	  //cout<<"OM found = "<<i<<endl;
	  drop=true;
	}
      }
      
      if(/*pic==4 &&*/ has_Bi && has_Li){
	//cout<<" li "<<diff_li_start_stop_vec[i][j]<<" bi "<<diff_bi_start_stop_vec[i]<<endl;
  	Bi_and_Li_close = 100*abs(diff_li_start_stop_vec[i][j]-diff_bi_start_stop_vec[i])/diff_bi_start_stop_vec[i];       
  	double Mean_derivative_li_f = TMath::Mean(li_derivative_vec[i][j].size(), &li_derivative_vec[i][j][0]);
  	double variation_derivative_li_f = TMath::StdDev(li_derivative_vec[i][j].size(), &li_derivative_vec[i][j][0]);

  	if(Mean_derivative_li_f>0.010 && variation_derivative_li_f>0.010){
  	  is_stable=false;
  	}
  	else{//if stable
	  if(li_pic_drop==false){
	    if(diff_li_start_stop_vec[i][j]>0.98){
	      is_perfect=true;
	    }
	    else{
	      decrease=true;
	    }
	  }	
	}
      }
      end_tree.Fill();
    }
  }
  file_tree->cd();
  end_tree.Write();
  file_tree->Close();
  

  double best_diff_li_bi=0.0, best_diff_li_bi_charge = 0.0, RMS_li=0.0, RMS_bi=0.0, best_li_amplitude = 0.0,best_li_amplitude_charge =0.0, bi_amplitude1 = 0.0; 
  int best_pic_li_end;
  TFile *file_result = new TFile("sortie/Comparison_Bi_Li/final_file_best_%.root", "RECREATE");
  TTree result_tree("Result_tree","");
  result_tree.Branch("li_relative",&li_relative);
  result_tree.Branch("bundle_li",&bundle_li);
  result_tree.Branch("run_number_li",&run_number_li);
  result_tree.Branch("Chi2_bi",&Chi2_bi);
  result_tree.Branch("nb_entries_bi",&nb_entries_bi);
  result_tree.Branch("ndf_bi",&ndf_bi);
  result_tree.Branch("bi_relative",&bi_relative);
  result_tree.Branch("om_number_li",&om_number_li);
  result_tree.Branch("pic_li",&pic_li);
  result_tree.Branch("diff_li_bi",&diff_li_bi);
  result_tree.Branch("best_diff_li_bi",&best_diff_li_bi);
  result_tree.Branch("best_diff_li_bi_charge",&best_diff_li_bi_charge);
  result_tree.Branch("RMS_li",&RMS_li);
  result_tree.Branch("RMS_bi",&RMS_bi);
  result_tree.Branch("mean_pic_bi",&mean_pic_bi);
  result_tree.Branch("li_amplitude",&li_amplitude);
  result_tree.Branch("best_li_amplitude",&best_li_amplitude);
  result_tree.Branch("best_li_amplitude_charge",&best_li_amplitude_charge);
  result_tree.Branch("best_pic_li_end",&best_pic_li_end);


  for (int j = 0; j < Result_tree->GetEntries(); j++) {
    Result_tree->GetEntry(j);
    if(li_relative!=0 /*&& li_relative!=1*/){
      if(!((bundle_li ==4 && run_number_li>1300 && run_number_li<1365)|| (bundle_li ==9 && run_number_li>1300 && run_number_li<1365))){
        if(Chi2_bi!=2000 && nb_entries_bi>50 && Chi2_bi/ndf_bi<3 && /*bi_relative!=1 &&*/ bi_relative!=0&& diff_li_bi!=1){
  	  if(best_pic[om_number_li]==pic_li && diff_li_bi!=0){
	    RMS_li = histograms_variation[om_number_li][pic_li-1]->GetRMS();	    
	    RMS_bi = histograms_bi[om_number_li]->GetRMS();	    
	    best_diff_li_bi = diff_li_bi;
	    best_li_amplitude = li_amplitude;
	    bi_amplitude1 = mean_pic_bi;
	    best_pic_li_end=pic_li;
  	  }
	  else{
	    best_diff_li_bi = 0.0;
	    //best_li_amplitude = 0.0;	    	    
	  }
	  if(best_charge_pic[om_number_li]==pic_li){
	    best_diff_li_bi_charge = diff_li_bi;
	    best_li_amplitude_charge = li_amplitude;
	  }
	  else{
	    best_diff_li_bi_charge = 0.0;
	  }
        }
      }
    }
    result_tree.Fill();
  }
  for(int i=0;i<712;i++){
    delete histograms_bi[i];
    for(int j=0; j<4; j++){
      delete histograms[i][j];
      delete histograms_variation[i][j];
      delete derivate_li[i][j];
    }
  }
  
  file_result->cd();
  result_tree.Write();
  file_result->Close();
  file_create->Close();
  delete file_create;
  
  delete file;
  delete file_result;
  
  //sncalo->draw1();
  // sncalo->setrange(0,4);
  // for(int i =0; i<712; i++){sncalo->settext(i, Form("%d", i));}
  //sncalo->draw1();
  
  //sncalo->canvas->SaveAs("/home/granjon/Bi/sortie/Comparison_Bi_Li/each_OM/calo_display_li_bi.png");
 }






int main(int argc, char const *argv[]){
  int nb_intensities = 8;
  std::vector<int> ref_run_number, ref_time, energy_run_number/*, run_number*/;
  string file, ref_correction, calo_correction;
  int n_files_LED = std::atoi(argv[1]);
  std::vector<int> LED_files;
  for (int i = 2; i < 2 + n_files_LED; ++i) {
    LED_files.push_back(std::atoi(argv[i]));
  }
  std::vector<int> vecteur_Bi_vector;
  for (int i = 2+n_files_LED; i < argc; ++i) {
    vecteur_Bi_vector.push_back(std::stoi(argv[i]));
  }

  //Passage au calo 
  std::cout << "" << '\n';
  std::cout << "START OF THE CALORIMETER FIT" << '\n';
  std::cout << "" << '\n';
  
  for(int run_value : vecteur_Bi_vector){
    fit_Bi_energy(run_value);
    std::cout<<"fit_BI ok "<<run_value<<std::endl;     
  }
  
  file_merger_Bi(vecteur_Bi_vector);
  cout<<"merger Bi ok"<<endl;
  cout<<"file "<<vecteur_Bi_vector[0]<<" - "<<vecteur_Bi_vector[vecteur_Bi_vector.size()-1] <<" created "<<endl;
  
  Evolution_Li_SN_graph(LED_files,LED_files.size(),LED_files[0],LED_files[LED_files.size()-1],vecteur_Bi_vector, nb_intensities/2);
  std::cout<<"evolution Li ok"<<std::endl;
  comparaison_Bi_Li(LED_files[0],LED_files[LED_files.size()-1], vecteur_Bi_vector[vecteur_Bi_vector.size()-1]);

  cout<<"comparison ok"<<endl;
  each_OM_mean(LED_files[0],LED_files[LED_files.size()-1]);
  std::cout<<"each OM ok"<<endl;
  std::cout<<"end"<<std::endl;
  
  return 0;
}


