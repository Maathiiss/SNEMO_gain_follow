// Standard library:
#include <iostream>
#include <fstream>
#include <exception>
#include <cstdlib>
#include <numeric> // for mean, sigma
#include <algorithm>
#include <string>

#include "TParameter.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

#include <snfee/snfee.h>
#include <snfee/io/multifile_data_reader.h>
#include <snfee/data/raw_trigger_data.h>

// TO ADD
#include <snfee/data/calo_waveform_drawer.h>
#include <snfee/data/calo_waveform_data.h>
#include <snfee/algo/calo_waveform_analysis.h>

// This project:
#include <sncabling/sncabling.h>
#include <sncabling/om_id.h>
#include <sncabling/calo_hv_id.h>
#include <sncabling/calo_hv_cabling.h>
#include <sncabling/label.h>

#include <sncabling/service.h>
#include <sncabling/calo_signal_cabling.h>

using namespace std;


bool debug = true;


bool usage(){

  std::clog<<std::endl;
  std::clog<<"+--------------------------------------------------+"<<std::endl;
  std::clog<<"| SuperNEMO calorimeter commissioning tutorial lv2 |"<<std::endl;
  std::clog<<"+--------------------------------------------------+"<<std::endl;

  std::clog<<"How to : "<<std::endl;
  std::clog<<" "<<std::endl;
  std::clog<<"List of option : "<<std::endl;
  std::clog<<"    |->   -i or --input_file "<<std::endl;
  std::clog<<std::endl;
  std::clog<<std::endl;

  std::clog<<"Example : "<<std::endl;
  std::clog<<"snfee-tuto_lv1 -i toto.data.gz"<<std::endl;

  return 1;
  
}

// Main program

int main(int argc, char **argv)
{
  sncabling::initialize();  
  int error_code = EXIT_SUCCESS;
  
  try {
    
    if (argc < 2){
      usage();
      exit;
    }

    int n_event = -1;
    std::string filename="default.dat";
    std::string output_filename="default.root";
    std::string crate_file="histo";
    int run_number = -1;

    for(int i = 0; i<argc; i++){
      if(std::string(argv[i]) =="--input_file" || std::string(argv[i]) =="-i" ){filename = argv[i+1];}
      if(std::string(argv[i]) =="--output_file" || std::string(argv[i]) =="-o" ){output_filename = argv[i+1];}
      if(std::string(argv[i]) =="--event_number" || std::string(argv[i]) =="-n"){n_event = atoi(argv[i+1]);} 
      if(std::string(argv[i]) =="--run_number" || std::string(argv[i]) =="-r"){run_number = atoi(argv[i+1]);}
    }

    double start_time = 0;
    std::ifstream  timefile(Form("/sps/nemo/snemo/snemo_data/raw_data/CBD/run-%d/snemo_crate-1_run-%d.log", run_number, run_number));

    if(timefile.is_open())        // Get starting time
      {
	string ligne; 
	int compteur = 0;
	while(getline(timefile, ligne)) 
	  {
	    if (compteur == 3){
	      start_time = atof(&ligne[20]);
	    }
	    compteur++;
	  }
      } 
    else
      {
	cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
      }
    
    TParameter<double> param("start_time", start_time);



    TFile file(Form("%s", output_filename.c_str()),"RECREATE");
    std::ifstream  charge("/sps/nemo/scratch/aguerre/sw/rtd_analysis/ReadRTD/gain_energy_692.txt");
  
    int ch_lt;
    uint16_t lt_trigger_counter = 0;
    uint32_t lt_trigger_time = 0;
    ULong64_t tdc = 0;
    int om_num = 0; 
    TH2F histo_pm_amplitude("histo_pm_amplitude", "amplitude", 717, 0, 717, 1024, 0, 1250);
    TH2F histo_pm_charge("histo_pm_charge", "charge", 717, 0, 717, 1024, 0, 200000);
    TH2F histo_pm_charge_50("histo_pm_charge_50", "charge", 717, 0, 717, 1024, 0, 50000);
    TH2F histo_pm_charge_energie("histo_pm_charge_energie", "charge", 717, 0, 717, 1024, 0, 10);
    TTree *Result_tree;
    Result_tree= new TTree("Result_tree","");
    Result_tree->Branch("om_number", &om_num);
    Result_tree->Branch("lt_trigger_counter", &lt_trigger_counter);
    Result_tree->Branch("lt_trigger_time", &lt_trigger_time);
    Result_tree->Branch("tdc", &tdc);
    Result_tree->Branch("ch_lt", &ch_lt);

    double charge_valeur_fit[717];
    memset(charge_valeur_fit, -1, 717*sizeof(float));   

    int charge_om_num;

    while (charge >> charge_om_num)
      {
	charge >> charge_valeur_fit[charge_om_num];
      }





    // Cf tutorial lv0
    sncabling::service snCabling;
    snCabling.initialize_simple();
    
    snfee::io::multifile_data_reader::config_type reader_cfg;
    reader_cfg.filenames.push_back(filename);
    
    snfee::io::multifile_data_reader rtd_source(reader_cfg);
    snfee::data::raw_trigger_data rtd; 


    const sncabling::calo_signal_cabling & caloSignalCabling
      = snCabling.get_calo_signal_cabling();
    
    
    bool display_wf    = true;
    /// Waveform drawer:
    std::unique_ptr<snfee::data::calo_waveform_drawer> calo_drawer;
    calo_drawer.reset(new snfee::data::calo_waveform_drawer(0,
							    snfee::model::feb_constants::SAMLONG_DEFAULT_TDC_LSB_NS
							    * snfee::model::feb_constants::SAMLONG_MAX_NUMBER_OF_SAMPLES,
							    snfee::model::feb_constants::SAMLONG_ADC_MIN_VOLTAGE_MV,
							    snfee::model::feb_constants::SAMLONG_ADC_MAX_VOLTAGE_MV / 5,
							    false,
							    true));
    
    
    uint64_t first_tdc = 0;
    uint64_t last_tdc  = 0;

    std::size_t rtd_counter = 0;
    //Main loop
    while (rtd_source.has_record_tag() ) {
      rtd_counter++;
      rtd_source.load(rtd);
      int32_t trigger_id = rtd.get_trigger_id();
      int32_t run_id     = rtd.get_run_id();

      //if(rtd_counter < 10)std::clog<<"In Run : "<<run_id<<" Trigger # "<<trigger_id <<std::endl;
      
      std::size_t calo_counter = 0;
      // Loop on calo hit
      for (const auto & p_calo_hit : rtd.get_calo_hits()) {
	const snfee::data::calo_hit_record & calo_hit = *p_calo_hit;
	tdc             = calo_hit.get_tdc();        // TDC timestamp (48 bits)
	int32_t  crate_num       = calo_hit.get_crate_num();  // Crate number (0,1,2)
	int32_t  board_num       = calo_hit.get_board_num();  // Board number (0-19)
	int32_t  chip_num        = calo_hit.get_chip_num();   // Chip number (0-7)
        bool     has_waveforms   = calo_hit.has_waveforms(); // Default: true
      

	if(rtd_counter == 1) first_tdc = tdc;
	else last_tdc = tdc;
	
	// if(rtd_counter < 10 ){
	// std::clog <<"--------------------"<<std::endl;
	// std::clog<<"   |-> tdc          : "<< tdc <<std::endl;
	// std::clog<<"   |-> tdc(ns)      : "<<  snfee::model::feb_constants::SAMLONG_DEFAULT_TDC_LSB_NS*tdc <<std::endl;
	//}



	
	// loop on channels
	for (int ichannel = 0; ichannel < snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS; ichannel++) {
	  const snfee::data::calo_hit_record::channel_data_record & ch_data = calo_hit.get_channel_data(ichannel);
	  ch_lt           = ch_data.is_lt();            // Low threshold flag
	  bool    ch_ht           = ch_data.is_ht();            // High threshold flag
	  int32_t ch_baseline     = ch_data.get_baseline();     // Computed baseline       (LSB: ADC unit/16)
	  int32_t ch_peak         = ch_data.get_peak();         // Computed peak amplitude (LSB: ADC unit/8)
	  int32_t ch_charge       = ch_data.get_charge();       // Computed charge
          lt_trigger_counter = calo_hit.get_lt_trigger_counter(ichannel);  // nombre de coups pendant le calo hit (pour calcul temps mort)
          lt_trigger_time = 1.25*calo_hit.get_lt_time_counter(ichannel);    // duree du coup (pour le temps mort)

	  snfee::data::channel_id ch_id (crate_num, board_num, snfee::model::feb_constants::SAMLONG_NUMBER_OF_CHANNELS * chip_num + ichannel);
	  sncabling::calo_signal_id readout_id;
	  DT_THROW_IF(!ch_id.export_to(readout_id), std::logic_error, "Cannot export calorimeter channel ID to SNCabling readout channel ID!");
	 

	  int wall_num = -1;
	  int side_num = -1;
	  int column_num = -1;
	  int row_num = -1;
	  om_num = -1;

	  
	  if (caloSignalCabling.has_channel(readout_id)) {
	    const sncabling::om_id & calo_id = caloSignalCabling.get_om(readout_id);
	    
            if (calo_id.is_main()) {
              side_num = calo_id.get_side();
              column_num = calo_id.get_column();
              row_num = calo_id.get_row();
              om_num = side_num*20*13 + column_num*13 + row_num;
            }
            else if (calo_id.is_xwall()) {
              side_num = calo_id.get_side();
              wall_num = calo_id.get_wall();
              column_num = calo_id.get_column();
              row_num = calo_id.get_row();
              om_num = 520 + side_num*64 +  wall_num*32 + column_num*16 + row_num;
            }
            else if (calo_id.is_gveto()) {
              side_num = calo_id.get_side();
              wall_num = calo_id.get_wall();
              column_num = calo_id.get_column();
              om_num = 520 + 128 + side_num*32 + wall_num*16 + column_num;	    }
	  } // if (has channel)
	    
	  if (crate_num == 1){
	    if ((board_num > - 1) && (board_num < 5)){
	      if (chip_num*2+ichannel == 14){
		om_num = 712 + board_num;
	      }
	    }
	  } 
	  if (om_num != -1){
	    
	    histo_pm_amplitude.Fill(om_num, -(ch_peak)*2500/(4096.0*8) );
	    histo_pm_charge.Fill(om_num, -ch_charge);
	    histo_pm_charge_50.Fill(om_num, -ch_charge);
	      
	    if (charge_valeur_fit[om_num] >0){
	      histo_pm_charge_energie.Fill(om_num, (-ch_charge)*charge_valeur_fit[om_num] );
	    }
	    Result_tree->Fill();	  
	  }
	  

	  // if (caloSignalCabling.has_channel(readout_id)) {
	   
	  //   if(rtd_counter < 10 ){
	  //     std::clog <<"ch_ht : "<<ch_ht<<std::endl;
	  //     std::clog <<"readout channel: "<<readout_id.to_label()<<std::endl;
	  //     std::clog <<"   |-> crate   : "<<crate_num<<std::endl;
	  //     std::clog <<"   |-> board   : "<<board_num<<std::endl;
	  //     std::clog <<"   |-> channel : "<<chip_num*2+ichannel<<std::endl;
	  //     std::clog <<"om id : "<<calo_id.to_label()<<std::endl;
	  //     std::clog <<"   |-> side   : "<<calo_id.get_side()<<std::endl;
	  //     std::clog <<"   |-> column : "<<calo_id.get_column()<<std::endl;
	  //     std::clog <<"   |-> row    : "<<calo_id.get_row()<<std::endl;
	  //     std::clog <<"Peak        : "<<ch_peak<<std::endl;
	  //     std::clog <<"Peak(mV)    : "<<ch_peak * snfee::model::feb_constants::SAMLONG_ADC_VOLTAGE_LSB_MV / 8<<std::endl;
	  //     std::clog <<"Charge  : "<<ch_charge<<std::endl;
	  //     std::clog <<"--------------------"<<std::endl;
	  //   }
	  // }

	  // // Declare the waveform array for this SAMLONG channel:
          // std::vector<uint16_t> ch_waveform;
	  // if (has_waveforms)
	  //   {
	  //   ch_waveform.reserve(snfee::model::feb_constants::SAMLONG_MAX_NUMBER_OF_SAMPLES);//1024
	  //   const sncabling::om_id & calo_id = caloSignalCabling.get_om(readout_id);
	  //   for (int isample = 0; isample <  1024; isample++)
	  //     {
	  // 	int16_t adc = calo_hit.get_waveforms().get_adc(isample, ichannel); // 0-4095 (ADC)


	  // 	    ch_waveform.push_back(adc);
		  
	  //     }
	  //   }// can play with WF
	    

	  // //channel to look after 
	  // sncabling::calo_signal_id control_readout_id(sncabling::CALOSIGNAL_CHANNEL,
	  // 				       0, 
	  // 				       5, 
	  // 				       1);

	 
	  // if (display_wf == true){
	  //   // Working waveform data structure:
	  //   snfee::data::calo_waveform_info waveform_info;
	  //   if (calo_drawer)
	  //     {
	  //     if (waveform_info.waveform.size() == 0) {
	  // 	// If waveform samples are not set yet, fetch them from the RTD
	  // 	// and draw the waveform shape:
	  // 	int16_t adc_zero  = snfee::model::feb_constants::SAMLONG_ADC_ZERO;
	  // 	snfee::algo::calo_waveform_analysis::populate_waveform(ch_waveform,
	  // 							       waveform_info.waveform,
	  // 							       snfee::model::feb_constants::SAMLONG_DEFAULT_TDC_LSB_NS,
	  // 							       snfee::model::feb_constants::SAMLONG_ADC_ZERO,
	  // 							       snfee::model::feb_constants::SAMLONG_ADC_VOLTAGE_LSB_MV);
	  //     }
	      
	  //     std::string counter_str = std::to_string(rtd_counter);
	  //     std::string const tmp_title { "event # " + counter_str + " ch : " + readout_id.to_label().to_string()};

	  //     if(rtd_counter == 1 || readout_id == control_readout_id)

	      
	  // 	{
	  // 	  //calo_drawer->draw(waveform_info, "", tmp_title);
				  
	  // 	}
	  //     //calo_drawer->draw(waveform_info, "", "TOTO");
	  //     //calo_drawer->draw(waveform_info, "", readout_id.to_label().to_string());
	  //     //
	  //     }
	  // }
	  
	}//end of channels
      }//end of calohit

      if (rtd_counter == n_event)break ;
  


    }//end of file


    uint64_t run_duration=snfee::model::feb_constants::SAMLONG_DEFAULT_TDC_LSB_NS*(last_tdc-first_tdc);
    
    std::clog<<"====== Summary ========="<<std::endl;
    std::clog<<std::endl;
    std::clog<<std::endl;

    std::clog<<"Run duration : "<<run_duration<<" ns"<<std::endl;
    std::clog<<"Nb of entry  : "<<rtd_counter<<std::endl;
    histo_pm_amplitude.Write();
    histo_pm_charge.Write();
    histo_pm_charge_50.Write();
    histo_pm_charge_energie.Write();
    Result_tree->Write();
    file.WriteTObject(&param);
    file.Close();

  } catch (std::exception & error) {
    std::cerr << "[error] " << error.what() << std::endl;
    error_code = EXIT_FAILURE;
  }
  sncabling::terminate();
  return error_code;
}
