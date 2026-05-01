#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPad.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TParameter.h>
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
// #include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <stdio.h>
#include <string.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <bits/stdc++.h>
#include "Strlen.h"
#include "TAttMarker.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TVirtualX.h"
#include "TVirtualPadEditor.h"
#include "TColor.h"
#include "Rtypes.h"

using std::cout; 
using std::endl;
using namespace std;

/** @brief Properties recorded for each detected pulse in a waveform */
struct pulse_temp {
  double start;  /* Start time of pulse (10% peak) in waveform (ns) */
  double end;    /* End time of pulse (reach baseline) in waveform (ns) */
  double peak;   /* Max amplitude of pulse (photo-electrons) */
  double energy; /* Energy (integral) of pulse (photo-electrons) */
};

/** @brief Properties recorded for each detected pulse in a waveform */
struct pulse {
  int file_num;             /* Entry # from inputFile TTree T */
  int entry;                /* Entry # from inputFile TTree T */
  long long unix_time;      /* Unix starttime value from inputFile TTree T */
  long long start;          /* Universal start time of pulse (10% peak) in waveform (ns) */
  long long end;            /* Universal end time of pulse (reach baseline) in waveform (ns) */
  double peak;              /* Max amplitude of pulse (photo-electrons) */
  double energy;            /* Energy (integral) of pulse (photo-electrons) */
  int number;               /* Number of channels in which we see pulse (photo-electrons) */
  bool single;              /* Is pulse (photo-electrons) timing consistent across all channels */
  bool beam;                /* Tracks whether beam is on or off */
  int trigger;              /* Tracks whether trigger is external (2) or internal (16) */
  int length;               /* Length of waveform in number of bins */
  double side_vp_energy;    /* Energy (integral) of pulse (photo-electrons) in SIDE veto panels */
  double top_vp_energy;     /* Energy (integral) of pulse (photo-electrons) in TOP veto panel */
  double all_vp_energy;     /* Energy (integral) of pulse (photo-electrons) in ALL veto panels */
  long long last_cosmics_time; /* Time value of most recent detected cosmic particle event */
  long long last_det_time;  /* Time value of most recent detected detector event */
  long long wf_time;        /* nsTime of event */
  long long vp_start_16;    /* Veto panel 16 start time */
  long long vp_start_17;    /* Veto panel 17 start time */
  long long vp_start_18;    /* Veto panel 18 start time */
  long long vp_start_19;    /* Veto panel 19 start time */
  long long vp_start_20;    /* Veto panel 20 start time */
  long long vp_start_21;    /* Veto panel 21 start time */
  long long vp_start_22;    /* Veto panel 22 start time */
  long long vp_start_23;    /* Veto panel 23 start time */
  long long vp_start_top;   /* Top veto panel start time */
  long long vp_int_16;      /* Veto panel 16 integral value */
  long long vp_int_17;      /* Veto panel 17 integral value */
  long long vp_int_18;      /* Veto panel 18 integral value */
  long long vp_int_19;      /* Veto panel 19 integral value */
  long long vp_int_20;      /* Veto panel 20 integral value */
  long long vp_int_21;      /* Veto panel 21 integral value */
  long long vp_int_22;      /* Veto panel 22 integral value */
  long long vp_int_23;      /* Veto panel 23 integral value */
  long long vp_int_top;     /* Top veto panel integral value */
  bool problem;             /* Tracks whether there is a potential problem with the event */
  long long ev61_start_time;/* Start time of Event 61 pulse */
  bool ev61_significant;    /* Flag to keep track of whether Event 61 event is also potentially significant */
  int selection_cut;        /* Tacks which selection cut(s) a given event satisfies */
};

/** @brief Constants for pulse and pulse-edge detection */
const int PULSE_THRESHOLD = 30;     /* Pulse detected if read above this value */
const int BS_UNCERTAINTY = 5;       /* Baseline uncertainty */
const int PULSE_PE_THRESHOLD = 40;  /* Large pulse detected if read above this value (given in units of photoelectrons) */
const int EV61_THRESHOLD = 1200;    /* SNS Beam assumed to be on if Event 61 read above this value (given in units of ADC) */

/** @brief Maximum number of waveforms to process from input root file */
const int MAX_NUM_ENTRIES = 1900000;        // run4144: 2700000 ; run4176: 2500000 ; run4193: 1900000

/** @brief These values should change for each PMT */
const char* inputFilePath = "";
const char* inputFileName = "processed";
const char* outputFilePath = "";                          // What does "" mean?
const char* outputFileName = "7894EventAnalysis";         // Why is output file going into ~ directory?
const char* outputStatsName = "PMTAnalysisStats";

std::vector<double> amplitudeToPE{ 26,28,30,29,27,12,30,28,30,30,24,30 };                // 12/3: run4176

std::vector<double> MichelSiPMNoiseAvg1{ 187.887,188.1116667,314.3093333,195.7781667,190.6345,202.5948333,209.708,201.6861667,243.2 };                          // Jul 17, 2023 - Aug 29, 2023

std::vector<double> MichelSiPMNoiseAvg2{ 185.7410476,174.2948095,461.303619,192.1831429,190.1772857,189.2594762,187.5979524,195.2732857,282.2847619 };          // Sep 3, 2023 - Jan 30, 2024

std::vector<double> MichelSiPMNoiseAvg3{ 220.497,202.8251905,545.8977619,220.8942857,220.4852381,222.7695238,220.9086667,231.9756667,300.433381 };              // Feb 9, 2024 - Jul 21, 2024

std::vector<double> MichelSiPMNoiseAvg4{ 222.2114286,213.4123571,231.4738571,222.3042143,220.602,234.8156429,237.2396429,240.2786429,282.2042143 };             // Jul 31, 2024 - Nov 23, 2024

std::vector<double> MichelSiPMNoiseAvg5{ 214.9131176,219.8434706,227.0522353,211.6268824,215.7442353,222.5197059,217.8758824,230.3478235,258.0880588 };         // Jan 12, 2025 - May 24, 2025

std::vector<double> MichelSiPMNoiseAvg6{ 202.2175,214.0718,190.558875,202.427125,209.80825,212.1195,195.1185,225.44875,227.461 };                               // May 24, 2025 - Aug 23, 2025

std::vector<double> MinMuonSiPMCut1{ 400,500,800,1000,300,400,350,250,400 };            // Jul 17, 2023 - Jan 30, 2024

std::vector<double> MinMuonSiPMCut2{ 600,700,900,1100,400,500,450,350,500 };            // Feb 9, 2024 - Dec 15, 2024

std::vector<double> MinMuonSiPMCut3{ 500,600,900,1100,400,400,350,300,500 };            // Jan 12, 2025 - Aug 24, 2025

template<typename T>
double getAverage(std::vector<T> const& v) {

    if (v.empty()) {

        return 0;

    }

    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

template<typename T>
int mostFrequent(std::vector<T> const& v) {

    if (v.empty()) {

        return 0;

    }

    int maxcount = 0;
    int element_having_max_freq;
    for (int i = 0; i < v.size(); i++) {
        int count = 0;
        for (int j = 0; j < v.size(); j++) {
            if (v[i] == v[j]) {
                count++;
            }
        }

        if (count > maxcount) {
            maxcount = count;
            element_having_max_freq = v[i];
        }
    }

    if (maxcount > 1) {

        return element_having_max_freq;

    }

    else {

        return std::accumulate(v.begin(), v.end(), 0.0) / v.size();

    }
    
}

template<typename T>
int mostFrequentIndex(std::vector<T> const& v) {

    if (v.empty()) {

        return -1;

    }

    int maxcount = 0;
    int index_having_max_freq;
    for (int i = 0; i < v.size(); i++) {
        int count = 0;
        for (int j = 0; j < v.size(); j++) {
            if (v[i] == v[j]) {
                count++;
            }
        }

        if (count > maxcount) {
            maxcount = count;
            index_having_max_freq = i;
        }
    }

    if (maxcount > 1) {

        return index_having_max_freq;

    }

    else {

        return -2;

    }
    
}

template<typename T>
T variance(const std::vector<T>& vec) {

    const size_t sz = vec.size();

    if (vec.empty()) {
        return 0.0;
    } else if (sz == 1) {
        return 0.0;
    }

    // Calculate the mean
    const T mean = std::accumulate(vec.begin(), vec.end(), 0.0) / sz;

    // Now calculate the variance
    auto variance_func = [&mean, &sz](T accumulator, const T& val) {
        return accumulator + ((val - mean) * (val - mean) / (sz - 1));      // No - 1?
    };

    return std::accumulate(vec.begin(), vec.end(), 0.0, variance_func);
}

template<typename T>
double rmsValue(std::vector<T>& vec)
{
    double square = 0.0;
    double mean = 0.0, root = 0.0;
 
    // Calculate square.
    for (int i = 0; i < vec.size(); i++) {
        square += pow(vec[i], 2);
    }
 
    // Calculate Mean.
    mean = (square / vec.size());
 
    // Calculate Root.
    root = sqrt(mean);
 
    return root;
}

int main(int argc, char *argv[])
{
    int run {0};
    int last_run {0};
    if (argc == 3) {
        sscanf(argv[1], "%i", &run);
        sscanf(argv[2], "%i", &last_run);
    } else if(argc == 2){ 
        sscanf(argv[1], "%i", &run);
        last_run = run;
    } else {
        cout <<"Usase: "<< argv[0] <<" [run]" << " [last run] (optional)" << endl;
        return -1; 
    }   
    cout << "\n" << "Run: " << run << " Last Run: " << last_run << endl;

    int ADCSIZE = 45;
    TH1D *h_wf = new TH1D("h_wf", "Waveform", ADCSIZE, 0, ADCSIZE);

    TH1D* h_intev61 = new TH1D("h_intev61", "Distribution of Event 61 Integral Values", 2200, -200, 2000);

    TH1D* h_intlowLED = new TH1D("h_intlowLED", "Distribution of Low Light LED Integral Values", 150, -500, 1000);

    TH1D* h_inthighLED = new TH1D("h_inthighLED", "Distribution of High Light LED Integral Values", 200, 0, 2000);            //, 11000, -1000, 10000);

    TH1D* h_intminLED = new TH1D("h_intminLED", "Distribution of Minimum Bias LED Integral Values", 150, -500, 1000);

    // Create a text string, which is used to output the text file
    
    std::vector<int> runlist;

    string intval1;

    // Read from the text file

    // ifstream ReadRunListFile("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/event_rate_over_time_runlist.txt");
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/yearoneevents/yearoneevents_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/eventrunlist/eventrunlist_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/sum2024runlist/sum2024runlist_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/sum2024goldenrunlist/sum2024goldenrunlist_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/2024allgoldenruns/2024allgoldenruns_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/weekly2024goldenruns/weekly2024goldenruns_%i.txt", run));
    ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/allweeklyruns/allweeklyruns_%i.txt", run));
    // 16151-16174 = 1721442002-1721524801 = 7/19/24-7/20/24 = file ~48 // 16351-16374 = 1722090001-1722172801 = 7/27/24-7/28/24 = file ~49
    // 16647-16671 = 1723155601-1723242001 = 8/08/24-8/09/24 = file 50  // 16855-16878 = 1723870861-1723953661 = 8/17/24-8/18/24 = file 51
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/F24-S25_GR_Runs/F24-S25_GR_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/F24-S25_BO_Runs/F24-S25_BO_%i.txt", run));

    // Use a while loop together with the getline() function to read the file line by line

    int iLine = 0;

    while (getline(ReadRunListFile, intval1)) {

        // Output the text from the file

        runlist.push_back(stod(intval1));

        iLine += 1;

        if (iLine == 24) {break;}
        
    }
    
    // Close the file

    ReadRunListFile.close();

    ReadRunListFile.clear();

    // Create a text string, which is used to output the text file

    int HLmeancounter = 0;

    double sumHLmean = 0.0;

    string intval2;

    // Read from the text file

    ifstream ReadHLMFile("Vertical_Muon_Peak_Center_1-124_Run_Vals.txt");           // hl_mean_val_all_1.txt          // Vertical_Muon_Peak_Center_1-124_Val_Only.txt

    // Use a while loop together with the getline() function to read the file line by line

    while (getline(ReadHLMFile, intval2)) {

        if (stod(intval2) > 0.0) {

            HLmeancounter += 1;

            sumHLmean += stod(intval2);

        }
        
    }

    double defaultHLmean = 1813.94;         // sumHLmean / HLmeancounter;

    double currentHLmean = defaultHLmean;

    // Close the file

    ReadHLMFile.close();

    ReadHLMFile.clear();

    /* Create vectors for plotting histograms */

    std::vector<double> intlowLED;

    std::vector<double> inthighLED;

    std::vector<double> intminLED;

    int ev_61_dt_max = pow(10, 5);

    int muon_dt_max = 16 * pow(10, 3);

    int muon_dt_extra = 500 * pow(10, 3);

    bool cent_event_found = false;

    bool take_sample = false;

    bool bt_LL_n_HL = false;

    int num_events_low_dt = 0;

    int num_events_low_dt_cent_spike = 0;

    int num_events_low_dt_no_spike_neg = 0;

    int num_events_low_dt_no_spike_pos = 0;

    int last_cosmics_counter = 0;

    int last_det_counter = 0;

    int last_event_counter = 0;

    int num_large_found_1 = 0;
    
    int num_large_found_2 = 0;

    int num_large_found_3 = 0;

    int num_large_found_4 = 0;

    int num_large_found_5 = 0;

    int num_large_found_6 = 0;

    int num_large_found_7 = 0;

    int num_large_found_8 = 0;

    int num_large_found_9 = 0;

    int num_large_found_10 = 0;

    int num_veto_found_1 = 0;

    int num_veto_found_2 = 0;

    int num_veto_found_3 = 0;

    int num_veto_found_4 = 0;

    int num_veto_found_5 = 0;

    int num_veto_found_6 = 0;

    int num_veto_found_7 = 0;

    int num_veto_found_8 = 0;

    int num_veto_found_9 = 0;

    int num_veto_found_10 = 0;

    // for(int iRun = run; iRun <= last_run ; iRun++) {
    
    for (size_t iRun = 0; iRun < runlist.size(); iRun++) {

        // int run_iterable = iRun;

        int run_iterable = runlist[iRun];

        // Below system makes code work for old and new SiPM/Event61 channel maps (won't work for runs 10367 - 10678)

        // int data_num = 9; int vers_num = 4;

        int data_num = 41; int vers_num = 5;

        if (run_iterable <= 14976) {data_num = 9; vers_num = 4;}

        else if (run_iterable > 14976 && run_iterable < 15696) {data_num = 41; vers_num = 4;}

        else if (run_iterable >= 15696) {data_num = 41; vers_num = 5;}

        int echo_cut = 1000;

        if (run_iterable < 20000) {echo_cut = 1000;}

        else if (run_iterable >= 20000) {echo_cut = 2000;}

        int vp18intcut = 400;

        if (run_iterable < 16431) {vp18intcut = 400;}

        else if (run_iterable >= 16431) {vp18intcut = 150;}

        bool old_channel_map = false; int Event61Chan = 15; int SiPM1Chan = 16; int SiPM2Chan = 17; int SiPM3Chan = 18; int SiPM4Chan = 19; int SiPM5Chan = 20; int SiPM6Chan = 21; int SiPM7Chan = 22; int SiPM8Chan = 23; int SiPM9Chan = 24; int SiPM10Chan = 25;

        if (vers_num < 5) {old_channel_map = true; Event61Chan = 15; SiPM1Chan = 16; SiPM2Chan = 17; SiPM3Chan = 18; SiPM4Chan = 19; SiPM5Chan = 20; SiPM6Chan = 21; SiPM7Chan = 22; SiPM8Chan = 23; SiPM9Chan = 24; SiPM10Chan = 25;}

        else if (vers_num == 5) {old_channel_map = false; SiPM1Chan = 12; SiPM2Chan = 13; SiPM3Chan = 14; SiPM4Chan = 15; SiPM5Chan = 16; SiPM6Chan = 17; SiPM7Chan = 18; SiPM8Chan = 19; SiPM9Chan = 20; SiPM10Chan = 21; Event61Chan = 22;}

        TFile *f;
        //your root file location here
        // if(gSystem->AccessPathName(Form("/data%i/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v%i.root", data_num, run_iterable, vers_num))){
        if (gSystem->AccessPathName(Form("/data%i/coherent/data/d2o/processedData/run%i_processed_v%i.root", data_num, run_iterable, vers_num))) {
            // if(gSystem->AccessPathName(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/run%i_processed_v%i.root", run_iterable, vers_num))){
            cout << "Could not open file " << run_iterable << endl;
            continue;           // return -1; 
        } else{
            // f = new TFile(Form("/data%i/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v%i.root", data_num, run_iterable, vers_num));
            f = new TFile(Form("/data%i/coherent/data/d2o/processedData/run%i_processed_v%i.root", data_num, run_iterable, vers_num));
            // f = new TFile(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/run%i_processed_v%i.root", run_iterable, vers_num));
        }

        auto tsstart = (TParameter<Long64_t> *) f->Get("starttime");
        
        Long64_t run_starttime = tsstart->GetVal();

        // Create a text string, which is used to output the text file

        std::vector<double> num_vec;

        num_vec.clear();

        string int_val_1;

        string num_str;

        int_val_1.clear();

        num_str.clear();

        // Read from the text file

        ifstream ReadHLMFile("Vertical_Muon_Peak_Center_1-124_Run_Nums.txt");           // hl_mean_val_look_up_1.txt          // Vertical_Muon_Peak_Center_1-124_No_HL_Mean_Correction.txt

        // Use a while loop together with the getline() function to read the file line by line

        while (getline(ReadHLMFile, int_val_1)) {

            num_vec.clear();

            num_str.clear();

            int_val_1.push_back('\t');

            for (char ch : int_val_1) {
                
                if (ch != '\t') {num_str.push_back(ch);}
                
                else if (ch == '\t') {num_vec.push_back(stod(num_str)); num_str.clear();}
            
            }

            if (num_vec[0] == run_iterable && num_vec[1] > 0.0) {

                // Output the text from the file

                currentHLmean = num_vec[1];

            }
            
        }

        // Close the file

        ReadHLMFile.close();

        ReadHLMFile.clear();
        
        // Create a text string, which is used to output the text file

        std::vector<double> integralToPE;

        std::vector<double> y_val;

        std::vector<double> b_val;

        y_val.clear();

        b_val.clear();

        string int_val_2;

        string y_str;

        string b_str;

        int_val_2.clear();

        y_str.clear();

        b_str.clear();

        integralToPE.clear();

        // Read from the text file

        // ifstream ReadSPIFile("Single_Phe_Integral_Values_Run7894.txt");
        // ifstream ReadSPIFile(Form("Single_Phe_Integral_Values_Run%i.txt", run_iterable));
        ifstream ReadSPIFile("sphe_int_conv_fits.txt");

        // Use a while loop together with the getline() function to read the file line by line

        iLine = 0;

        while (getline(ReadSPIFile, int_val_2)) {

            iLine += 1;

            // Output the text from the file

            // integralToPE.push_back(stod(int_val_2));

            if (run_starttime < 1706473201) {int_val_2.push_back('\t'); if (iLine == 1) {for (char ch : int_val_2) {if (ch != '\t') {y_str.push_back(ch);} else if (ch == '\t') {y_val.push_back(stod(y_str)); y_str.clear();}}} else if (iLine == 2) {for (char ch : int_val_2) {if (ch != '\t') {b_str.push_back(ch);} else if (ch == '\t') {b_val.push_back(stod(b_str)); b_str.clear();}}}}

            else if (run_starttime >= 1706473201 && run_starttime < 1721578801) {int_val_2.push_back('\t'); if (iLine == 3) {for (char ch : int_val_2) {if (ch != '\t') {y_str.push_back(ch);} else if (ch == '\t') {y_val.push_back(stod(y_str)); y_str.clear();}}} else if (iLine == 4) {for (char ch : int_val_2) {if (ch != '\t') {b_str.push_back(ch);} else if (ch == '\t') {b_val.push_back(stod(b_str)); b_str.clear();}}}}

            else if (run_starttime >= 1721578801 && run_starttime < 1723645201) {int_val_2.push_back('\t'); if (iLine == 5) {for (char ch : int_val_2) {if (ch != '\t') {y_str.push_back(ch);} else if (ch == '\t') {y_val.push_back(stod(y_str)); y_str.clear();}}} else if (iLine == 6) {for (char ch : int_val_2) {if (ch != '\t') {b_str.push_back(ch);} else if (ch == '\t') {b_val.push_back(stod(b_str)); b_str.clear();}}}}

            else if (run_starttime >= 1723645201 && run_starttime < 1730937661) {int_val_2.push_back('\t'); if (iLine == 7) {for (char ch : int_val_2) {if (ch != '\t') {y_str.push_back(ch);} else if (ch == '\t') {y_val.push_back(stod(y_str)); y_str.clear();}}} else if (iLine == 8) {for (char ch : int_val_2) {if (ch != '\t') {b_str.push_back(ch);} else if (ch == '\t') {b_val.push_back(stod(b_str)); b_str.clear();}}}}

            else if (run_starttime >= 1730937661 && run_starttime < 1736388061) {int_val_2.push_back('\t'); if (iLine == 9) {for (char ch : int_val_2) {if (ch != '\t') {y_str.push_back(ch);} else if (ch == '\t') {y_val.push_back(stod(y_str)); y_str.clear();}}} else if (iLine == 10) {for (char ch : int_val_2) {if (ch != '\t') {b_str.push_back(ch);} else if (ch == '\t') {b_val.push_back(stod(b_str)); b_str.clear();}}}}

            else if (run_starttime >= 1736388061 && run_starttime < 1740117662) {int_val_2.push_back('\t'); if (iLine == 11) {for (char ch : int_val_2) {if (ch != '\t') {y_str.push_back(ch);} else if (ch == '\t') {y_val.push_back(stod(y_str)); y_str.clear();}}} else if (iLine == 12) {for (char ch : int_val_2) {if (ch != '\t') {b_str.push_back(ch);} else if (ch == '\t') {b_val.push_back(stod(b_str)); b_str.clear();}}}}

            else if (run_starttime >= 1740117662 && run_starttime < 1748473261) {int_val_2.push_back('\t'); if (iLine == 13) {for (char ch : int_val_2) {if (ch != '\t') {y_str.push_back(ch);} else if (ch == '\t') {y_val.push_back(stod(y_str)); y_str.clear();}}} else if (iLine == 14) {for (char ch : int_val_2) {if (ch != '\t') {b_str.push_back(ch);} else if (ch == '\t') {b_val.push_back(stod(b_str)); b_str.clear();}}}}

            else if (run_starttime >= 1748473261 && run_starttime < 1750147261) {int_val_2.push_back('\t'); if (iLine == 15) {for (char ch : int_val_2) {if (ch != '\t') {y_str.push_back(ch);} else if (ch == '\t') {y_val.push_back(stod(y_str)); y_str.clear();}}} else if (iLine == 16) {for (char ch : int_val_2) {if (ch != '\t') {b_str.push_back(ch);} else if (ch == '\t') {b_val.push_back(stod(b_str)); b_str.clear();}}}}

            else if (run_starttime >= 1750147261 && run_starttime < 1752739262) {int_val_2.push_back('\t'); if (iLine == 17) {for (char ch : int_val_2) {if (ch != '\t') {y_str.push_back(ch);} else if (ch == '\t') {y_val.push_back(stod(y_str)); y_str.clear();}}} else if (iLine == 18) {for (char ch : int_val_2) {if (ch != '\t') {b_str.push_back(ch);} else if (ch == '\t') {b_val.push_back(stod(b_str)); b_str.clear();}}}}

            else if (run_starttime >= 1752739262) {int_val_2.push_back('\t'); if (iLine == 19) {for (char ch : int_val_2) {if (ch != '\t') {y_str.push_back(ch);} else if (ch == '\t') {y_val.push_back(stod(y_str)); y_str.clear();}}} else if (iLine == 20) {for (char ch : int_val_2) {if (ch != '\t') {b_str.push_back(ch);} else if (ch == '\t') {b_val.push_back(stod(b_str)); b_str.clear();}}}}
            
        }

        // Close the file

        ReadSPIFile.close();

        ReadSPIFile.clear();

        // Select appropriate SiPM Noise values

        std::vector<double> vp_cut_vals;

        if (run_starttime <= 1693282801) {vp_cut_vals = MichelSiPMNoiseAvg1;}

        else if (run_starttime > 1693282801 && run_starttime <= 1706628002) {vp_cut_vals = MichelSiPMNoiseAvg2;}

        else if (run_starttime > 1706628002 && run_starttime <= 1721521201) {vp_cut_vals = MichelSiPMNoiseAvg3;}

        else if (run_starttime > 1721521201 && run_starttime <= 1734260462) {vp_cut_vals = MichelSiPMNoiseAvg4;}

        else if (run_starttime > 1734260462 && run_starttime <= 1748775661) {vp_cut_vals = MichelSiPMNoiseAvg5;}

        else if (run_starttime > 1748775661) {vp_cut_vals = MichelSiPMNoiseAvg6;}

        // Select appropriate Minimum Muon values

        std::vector<double> min_muon_vals;

        if (run_starttime <= 1706628002) {min_muon_vals = MinMuonSiPMCut1;}

        else if (run_starttime > 1706628002 && run_starttime <= 1734260462) {min_muon_vals = MinMuonSiPMCut2;}

        else if (run_starttime > 1734260462) {min_muon_vals = MinMuonSiPMCut3;}

        // TTree *t = (TTree *) f->Get("wf");
        TTree* t = (TTree*)f->Get("tree");

        // Declaration of leaf types
        Int_t           eventID;
        Int_t           nSamples[32];
        ///        UInt_t          adcTime;
        ///        Int_t           adcSize[32];
        Short_t         adcVal[32][45];
        ///        Int_t           trigPattern;
        ///        ULong64_t       timeStamp_extTrig;
        Double_t        baselineMean[32];
        Double_t        baselineRMS[32];
        ///        Int_t           peakBin[32];
        Double_t        pulseH[32];
        Int_t           peakPosition[32];
        Double_t        area[32];
        Long64_t        nsTime;
        Int_t           triggerBits;

        // List of branches
        TBranch        *b_eventID;
        TBranch        *b_nSamples;
        ///        TBranch        *b_adcTime;
        ///        TBranch        *b_adcSize;
        TBranch        *b_adcVal;
        ///        TBranch        *b_trigPattern;
        ///        TBranch        *b_timeStamp_extTrig;
        TBranch        *b_baselineMean;
        TBranch        *b_baselineRMS;
        ///        TBranch        *b_peakBin;
        TBranch        *b_pulseH;
        TBranch        *b_peakPosition;
        TBranch        *b_area;
        TBranch        *b_nsTime;
        TBranch        *b_triggerBits;

        t->SetBranchAddress("eventID", &eventID, &b_eventID);
        t->SetBranchAddress("nSamples", &nSamples, &b_nSamples);
        ///        t->SetBranchAddress("adcTime", &adcTime, &b_adcTime);
        ///        t->SetBranchAddress("adcSize", adcSize, &b_adcSize);
        t->SetBranchAddress("adcVal", adcVal, &b_adcVal);
        ///        t->SetBranchAddress("trigPattern", &trigPattern, &b_trigPattern);
        ///        t->SetBranchAddress("timeStamp_extTrig", &timeStamp_extTrig, &b_timeStamp_extTrig);
        t->SetBranchAddress("baselineMean", baselineMean, &b_baselineMean);
        t->SetBranchAddress("baselineRMS", baselineRMS, &b_baselineRMS);
        ///        t->SetBranchAddress("peakBin", peakBin, &b_peakBin);
        t->SetBranchAddress("pulseH", pulseH, &b_pulseH);
        t->SetBranchAddress("peakPosition", &peakPosition, &b_peakPosition);
        t->SetBranchAddress("area", area, &b_area);
        t->SetBranchAddress("nsTime", &nsTime, &b_nsTime);
        t->SetBranchAddress("triggerBits", &triggerBits, &b_triggerBits);

        // Create output root file of processed waveform data
        TString outputFile;
        if (strcmp(outputFilePath, "") == 0) {
            outputFile.Form("%s.root", Form("ev61_analysis_week%i_run%i", run, run_iterable));
        } else {
            outputFile.Form("%s/%s.root", outputFilePath, Form("ev61_analysis_week%i_run%i", run, run_iterable));
        }
        TFile *fileOut = new TFile(outputFile, "RECREATE");
        TTree *eventTree = new TTree("eventTree", "eventTree");
        eventTree->SetDirectory(fileOut);

        // Data to keep track of for each eventTree entry
        Int_t *br_fnum;     /* File # from inputFile TTree T */
        Int_t *br_entry;    /* Entry # from inputFile TTree T */
        Long64_t *br_utime; /* Unix starttime value from inputFile TTree T */
        Long64_t *br_st;    /* Universal start time (10% peak) of event peak (nanoseconds) */
        Long64_t *br_ed;    /* Universal end time (10% peak) of event peak (nanoseconds) */
        Double_t *br_pk;    /* Sum of peaks (max amplitude) of all event peaks across all PMT channels (photo-electrons) */
        Double_t *br_ey;    /* Sum of energies of all event peaks across all PMT channels (photo-electrons) */
        Int_t *br_nr;       /* Multiplicity of event (number of PMT channels in which a peak was seen) */
        Bool_t *br_se;      /* Is first pulse (photo-electrons) timing consistent across all channels */
        Bool_t *br_bm;      /* Tracks whether beam is on or off for first pulse */
        Int_t *br_tr;       /* Type of event (triggerBits value) */
        Int_t *br_lh;       /* Size of waveform in bins */
        Double_t *br_svpe;  /* Sum of all side SiPM energies (photo-electrons) */
        Double_t *br_tvpe;  /* Sum of all top SiPM energies (photo-electrons) */
        Double_t *br_avpe;  /* Sum of all SiPM energies (photo-electrons) */
        Long64_t *br_lct;   /* Universal time of most recent muon event before current event (nanoseconds) */
        Long64_t *br_ldt;   /* Universal time of most recent detector event before current event (nanoseconds) */
        Bool_t *br_issue;   /* Flag to keep track of unusual eventTree entries */
        Long64_t *br_wftime;    /* nsTime of event */
        Long64_t *br_vptime16;  /* Veto panel 16 start time */
        Long64_t *br_vptime17;  /* Veto panel 17 start time */
        Long64_t *br_vptime18;  /* Veto panel 18 start time */
        Long64_t *br_vptime19;  /* Veto panel 19 start time */
        Long64_t *br_vptime20;  /* Veto panel 20 start time */
        Long64_t *br_vptime21;  /* Veto panel 21 start time */
        Long64_t *br_vptime22;  /* Veto panel 22 start time */
        Long64_t *br_vptime23;  /* Veto panel 23 start time */
        Long64_t *br_vptimetop; /* Top Veto panel start time */
        Long64_t *br_vpint16;   /* Veto panel 16 integral value */
        Long64_t *br_vpint17;   /* Veto panel 17 integral value */
        Long64_t *br_vpint18;   /* Veto panel 18 integral value */
        Long64_t *br_vpint19;   /* Veto panel 19 integral value */
        Long64_t *br_vpint20;   /* Veto panel 20 integral value */
        Long64_t *br_vpint21;   /* Veto panel 21 integral value */
        Long64_t *br_vpint22;   /* Veto panel 22 integral value */
        Long64_t *br_vpint23;   /* Veto panel 23 integral value */
        Long64_t *br_vpinttop;  /* Top veto panel integral value */
        Long64_t *br_e61st;     /* Start time of Event 61 pulse */
        Bool_t *br_e61sig;      /* Flag to keep track of whether Event 61 event is also potentially significant */
        Int_t *br_selcut;       /* Tacks which selection cut(s) a given event satisfies */

        eventTree->Branch("fnum", &br_fnum, "fnum/I");
        eventTree->Branch("entry", &br_entry, "entry/I");
        eventTree->Branch("utime", &br_utime, "utime/L");
        eventTree->Branch("st", &br_st, "st/L");
        eventTree->Branch("ed", &br_ed, "ed/L");
        eventTree->Branch("pk", &br_pk, "pk/d");
        eventTree->Branch("ey", &br_ey, "ey/d");
        eventTree->Branch("nr", &br_nr, "nr/I");
        eventTree->Branch("se", &br_se, "se/O");
        eventTree->Branch("bm", &br_bm, "bm/O");
        eventTree->Branch("tr", &br_tr, "tr/I");
        eventTree->Branch("lh", &br_lh, "lh/I");
        eventTree->Branch("svpe", &br_svpe, "svpe/d");
        eventTree->Branch("tvpe", &br_tvpe, "tvpe/d");
        eventTree->Branch("avpe", &br_avpe, "avpe/d");
        eventTree->Branch("lct", &br_lct, "lct/L");
        eventTree->Branch("ldt", &br_ldt, "ldt/L");
        eventTree->Branch("issue", &br_issue, "issue/O");
        eventTree->Branch("wftime", &br_wftime, "wftime/L");
        eventTree->Branch("vptime16", &br_vptime16, "vptime16/L");
        eventTree->Branch("vptime17", &br_vptime17, "vptime17/L");
        eventTree->Branch("vptime18", &br_vptime18, "vptime18/L");
        eventTree->Branch("vptime19", &br_vptime19, "vptime19/L");
        eventTree->Branch("vptime20", &br_vptime20, "vptime20/L");
        eventTree->Branch("vptime21", &br_vptime21, "vptime21/L");
        eventTree->Branch("vptime22", &br_vptime22, "vptime22/L");
        eventTree->Branch("vptime23", &br_vptime23, "vptime23/L");
        eventTree->Branch("vptimetop", &br_vptimetop, "vptimetop/L");
        eventTree->Branch("vpint16", &br_vpint16, "vpint16/L");
        eventTree->Branch("vpint17", &br_vpint17, "vpint17/L");
        eventTree->Branch("vpint18", &br_vpint18, "vpint18/L");
        eventTree->Branch("vpint19", &br_vpint19, "vpint19/L");
        eventTree->Branch("vpint20", &br_vpint20, "vpint20/L");
        eventTree->Branch("vpint21", &br_vpint21, "vpint21/L");
        eventTree->Branch("vpint22", &br_vpint22, "vpint22/L");
        eventTree->Branch("vpint23", &br_vpint23, "vpint23/L");
        eventTree->Branch("vpinttop", &br_vpinttop, "vpinttop/L");
        eventTree->Branch("e61st", &br_e61st, "e61st/L");
        eventTree->Branch("e61sig", &br_e61sig, "e61sig/O");
        eventTree->Branch("selcut", &br_selcut, "selcut/I");

        // Get statistics for up to 1 million entries from fileIn TTree T
        int numEntries = std::min((int)t->GetEntries(), MAX_NUM_ENTRIES);

        std::vector<double> chan_lengths;

        std::vector<double> peak_pos_RMS;

        std::vector<double> chan_start_no_outliers;

        std::vector<struct pulse> pulses;

        long long last_event_time = 0.0;

        long long LL_time = 0.0;

        long long HL_time = 0.0;

        bool record_pulses = false;

        bool record_ext_trig = false;

        double avg_peak_pos_RMS = 0.0;

        double var_peak_pos_RMS = 0.0;

        int ch15_on = 0;

        double last_take_sample_time = 0.0;
    
        // Can use either "t->GetEntries()" or "numEntries" or number
        for (int iEnt = 0; iEnt < t->GetEntries(); iEnt++) {

            // cout << "\n" << "Processing event " << iEnt + 1 << " of " << t->GetEntries() << endl;

            Long64_t tentry = t->LoadTree(iEnt);

            b_eventID->GetEntry(tentry);
            b_nSamples->GetEntry(tentry);
            // b_adcTime->GetEntry(tentry);
            // b_adcSize->GetEntry(tentry);
            b_adcVal->GetEntry(tentry);
            // b_trigPattern->GetEntry(tentry);
            // b_timeStamp_extTrig->GetEntry(tentry);
            b_baselineMean->GetEntry(tentry);
            b_baselineRMS->GetEntry(tentry);
            // b_peakBin->GetEntry(tentry);
            b_pulseH->GetEntry(tentry);
            b_peakPosition->GetEntry(tentry);
            b_area->GetEntry(tentry);
            b_nsTime->GetEntry(tentry);
            b_triggerBits->GetEntry(tentry);

            // Create vector to hold info of all pulses detected
            // std::vector<struct pulse> pulses;
            std::vector<struct pulse_temp> pulses_temp;

            // Digitizer cannot read above this so flag entries w/ p1,p2 > maxPeak
            double maxPeak = 15776 / amplitudeToPE[0];

            // Create variables to hold info of each current pulse
            bool onPulse = false, onLastPulseTail = false;
            int thresholdBin = 0, peakBin = 0;
            double tailWindow = 250;
            double peak = 0.;
            double pulseEnergy = 0.;
            double AllPulseEnergy = 0.;
            double Ev61Energy = 0.;
            double Ev61Peak = 0.;
            double chan_15_start = 0.;

            std::vector<double> all_chan_start;
            std::vector<double> all_chan_start_adj;
            std::vector<double> all_chan_end;
            std::vector<double> all_chan_peak;
            std::vector<double> all_chan_energy;
            std::vector<double> peak_energy;
            std::vector<double> ten_chan_peak;
            std::vector<double> all_chan_peakbin;
            std::vector<double> ten_chan_peakbin;
            std::vector<double> side_veto_panel_energy;
            std::vector<double> top_veto_panel_energy;
            std::vector<double> all_chan_start_pmt;
            std::vector<double> all_chan_start_vp;
            std::vector<double> all_chan_start_pmt_adj;
            std::vector<double> all_chan_start_vp_adj;

            std::vector<int> event_bin_val;

            bool all_chan_beam = false;

            bool top_vp_event = false;

            bool pulse_at_end = false;

            int pulse_at_end_count = 0;

            int num_chan = 0;

            int num_chan_over_1phe = 0;

            int num_chan_simple = 0;

            double num_chan_vp = 0;

            double sphe_int = 0;

            double p_int_16 = 0.0;
            double p_int_17 = 0.0;
            double p_int_18 = 0.0;
            double p_int_19 = 0.0;
            double p_int_20 = 0.0;
            double p_int_21 = 0.0;
            double p_int_22 = 0.0;
            double p_int_23 = 0.0;
            double p_int_24 = 0.0;
            double p_int_25 = 0.0;

            double p_start_16 = 0.0;
            double p_start_17 = 0.0;
            double p_start_18 = 0.0;
            double p_start_19 = 0.0;
            double p_start_20 = 0.0;
            double p_start_21 = 0.0;
            double p_start_22 = 0.0;
            double p_start_23 = 0.0;
            double p_start_24 = 0.0;
            double p_start_25 = 0.0;

            chan_lengths.clear();

            for (int iChan = 0; iChan < 26; iChan++) {

                if (old_channel_map) {if (iChan == 12 || iChan == 13 || iChan == 14) {all_chan_start.push_back(1000); continue;}}

                else if (!old_channel_map) {if (iChan == 23 || iChan == 24 || iChan == 25) {all_chan_start.push_back(1000); continue;}}

                for (int i = 0; i < ADCSIZE; i++) {h_wf->SetBinContent(i + 1, adcVal[iChan][i] - baselineMean[iChan]); event_bin_val.push_back(adcVal[iChan][i] - baselineMean[iChan]);}

                if (iChan <= 11) {sphe_int = run_starttime * y_val[iChan] + b_val[iChan];}

                // h_wf->GetXaxis()->SetTitle("Time (1 bin = 16 ns)");
                // h_wf->GetYaxis()->SetTitle("Counts");
                // h_wf->Draw();
                // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/wf_image/wf_run_%i_event_%i_ch_%i.png", run_iterable, iEnt + 1, iChan));
                // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/wf_image/wf_run_%i_event_%i_ch_%i.png", run_iterable, iEnt + 1, iChan));

                if (iChan == Event61Chan) {          // Move this to line 641?

                    Ev61Energy = 0.0;

                    for (int iBin = 1; iBin <= h_wf->GetNbinsX(); iBin++) {

                        double iBinContent = h_wf->GetBinContent(iBin);

                        Ev61Energy += iBinContent;

                    }

                    if (Ev61Energy > EV61_THRESHOLD) {

                        all_chan_beam = true;

                        ch15_on += 1;

                    }

                }

                pulses_temp.clear();

                onLastPulseTail = false;

                bool peak_in_any_chan = false;

                bool peak_in_chan = false;

                bool peak_in_chan_simple = false;

                for (int iBin = 1; iBin <= h_wf->GetNbinsX(); iBin++) {

                    double iBinContent = h_wf->GetBinContent(iBin);
                    // cout << iBinContent << "\t";

                    if (iChan <= 11 && !peak_in_chan_simple && iBinContent >= PULSE_THRESHOLD) {

                        peak_in_chan_simple = true;

                        num_chan_simple += 1;

                    }
                    /*
                    if (iChan >= SiPM1Chan && !peak_in_chan_veto && iBinContent >= PULSE_THRESHOLD) {

                        peak_in_chan_veto = true;

                        num_chan_veto += 1;

                    }
                    */
                    if (iBin >= 20) {

                        AllPulseEnergy += iBinContent;

                    }

                    if (pulses_temp.size() > 0) {    // Could get rid of this check b/c event window is so small

                        onLastPulseTail = iBin * 16 - pulses_temp.back().start < tailWindow;

                    }

                    // Find pulse
                    if (!onPulse && !onLastPulseTail && iBinContent >= PULSE_THRESHOLD && iBin >= 20) {

                        onPulse = true;
                        thresholdBin = iBin;
                        peakBin = iBin;
                        peak = iBinContent;
                        pulseEnergy += iBinContent;

                    }

                    // Pulse found. Find pulse duration & energy
                    else if (onPulse) {

                        // Check if pulse is at end of waveform
                        if (iChan <= 11 && iBin == ADCSIZE && iBinContent > 100) {pulse_at_end_count += 1; if (pulse_at_end_count >= 10) {pulse_at_end = true;}}            // cout << "\n" << "pulse_at_end_count = " << pulse_at_end_count << ", Run, Event, triggerBits = " << run_iterable << " " << iEnt + 1 << " " << triggerBits << endl;}}

                        // if (iChan <= 11 && run_iterable == 7894 && iEnt + 1 == 389830) {cout << "\n" << "Run 7894, Event 389830: pulse_at_end_count = " << pulse_at_end_count << "; iChan = " << iChan << "; iBin = " << iBin << "; iBinContent = " << iBinContent << endl;}

                        // Accumlate energy of pulse after threshold bin
                        pulseEnergy += iBinContent;

                        // Update pulse's peak
                        if (peak < iBinContent) {
                            peak = iBinContent;
                            peakBin = iBin;
                        }

                        // Search for end of pulse (falls below noiselevel)
                        // Assumes no pulse pileup
                        if (iBinContent < BS_UNCERTAINTY || (iBin == ADCSIZE && iBinContent < peak * 0.6)) {

                            // cout << "\n" << "Found a pulse in Event " << iEnt + 1 << " Channel " << iChan << endl;

                            // Create pulse info
                            struct pulse_temp p;        // Clear these struct variables?
                            p.start = thresholdBin * 16.;
                            p.peak = peak / amplitudeToPE[iChan];
                            // Record end of pulse
                            p.end = iBin * 16.;
                            // Accumlate energy of pulse before threshold bin
                            for (int j = peakBin - 1; BS_UNCERTAINTY < h_wf->GetBinContent(j); j--)
                            {
                                if (j < thresholdBin) {     // Technically, this if statement is a bit time inefficient
                                    pulseEnergy += h_wf->GetBinContent(j);

                                    // Record start of pulse (10% of peak)
                                    if (h_wf->GetBinContent(j) > peak * 0.1) {
                                        p.start = j * 16.;
                                    }

                                }

                                if (j == 0) {break;}
                                
                            }

                            if (iChan <= 11 && peak > 100) {peak_in_any_chan = true; all_chan_start.push_back(p.start);}

                            if (iChan >= SiPM1Chan && iChan <= SiPM8Chan && peak > 50) {peak_in_any_chan = true; all_chan_start.push_back(p.start);}

                            if ((iChan == SiPM9Chan || iChan == SiPM10Chan) && peak > 30) {peak_in_any_chan = true; all_chan_start.push_back(p.start);}

                            if (iChan == Event61Chan && peak > 100) {peak_in_any_chan = true; all_chan_start.push_back(p.start);}

                            if (iChan == Event61Chan) {Ev61Peak = peak;}

                            if (iChan <= 11) {

                                // Record energy of pulse
                                p.energy = pulseEnergy / sphe_int;          // integralToPE[iChan];

                                all_chan_start_pmt.push_back(p.start);
                                all_chan_end.push_back(p.end);
                                peak_energy.push_back(p.energy);
                                // all_chan_peak.push_back(p.peak);

                                if (!peak_in_chan) {

                                    num_chan += 1;

                                    all_chan_peak.push_back(p.peak);

                                    ten_chan_peak.push_back(p.peak);

                                    all_chan_peakbin.push_back(peakBin);

                                    ten_chan_peakbin.push_back(peakBin);

                                    // all_chan_energy.push_back(p.energy);

                                    peak_in_chan = true;

                                }

                                // Add pulse info to vector of pulses
                                pulses_temp.push_back(p);

                            }

                            if (iChan == Event61Chan) {

                                chan_15_start = p.start;

                            }

                            if (iChan >= SiPM1Chan && iChan <= SiPM10Chan) {

                                all_chan_start_vp.push_back(p.start);

                                if (iChan == SiPM1Chan) {p_start_16 = p.start;}

                                else if (iChan == SiPM2Chan) {p_start_17 = p.start;}

                                else if (iChan == SiPM3Chan) {p_start_18 = p.start;}

                                else if (iChan == SiPM4Chan) {p_start_19 = p.start;}

                                else if (iChan == SiPM5Chan) {p_start_20 = p.start;}

                                else if (iChan == SiPM6Chan) {p_start_21 = p.start;}

                                else if (iChan == SiPM7Chan) {p_start_22 = p.start;}

                                else if (iChan == SiPM8Chan) {p_start_23 = p.start;}

                                else if (iChan == SiPM9Chan) {p_start_24 = p.start;}

                                else if (iChan == SiPM10Chan) {p_start_25 = p.start;}

                            }

                            // Clear current pulse variables to look for new pulse
                            peak = 0.;
                            peakBin = 0;
                            pulseEnergy = 0.;
                            thresholdBin = 0;
                            onPulse = false;
                        }

                        if (iBin == ADCSIZE) {

                            if (iChan >= SiPM1Chan && iChan <= SiPM10Chan) {

                                if (iChan == SiPM1Chan) {p_start_16 = thresholdBin * 16.;}

                                else if (iChan == SiPM2Chan) {p_start_17 = thresholdBin * 16.;}

                                else if (iChan == SiPM3Chan) {p_start_18 = thresholdBin * 16.;}

                                else if (iChan == SiPM4Chan) {p_start_19 = thresholdBin * 16.;}

                                else if (iChan == SiPM5Chan) {p_start_20 = thresholdBin * 16.;}

                                else if (iChan == SiPM6Chan) {p_start_21 = thresholdBin * 16.;}

                                else if (iChan == SiPM7Chan) {p_start_22 = thresholdBin * 16.;}

                                else if (iChan == SiPM8Chan) {p_start_23 = thresholdBin * 16.;}

                                else if (iChan == SiPM9Chan) {p_start_24 = thresholdBin * 16.;}

                                else if (iChan == SiPM10Chan) {p_start_25 = thresholdBin * 16.;}

                            }

                            // Clear current pulse variables to look for new pulse
                            peak = 0.;
                            peakBin = 0;
                            pulseEnergy = 0.;
                            thresholdBin = 0;
                            onPulse = false;

                        }

                    }

                } // Waveform loop

                if (iChan <= 11) {

                    all_chan_energy.push_back(AllPulseEnergy / sphe_int);           // integralToPE[iChan]

                    if (AllPulseEnergy / sphe_int > 1) {            // integralToPE[iChan]

                        num_chan_over_1phe += 1;

                    }

                    if (!peak_in_chan) {

                        ten_chan_peak.push_back(0);

                        ten_chan_peakbin.push_back(0);

                    }

                }

                if (iChan >= SiPM1Chan && iChan <= SiPM8Chan) {

                    if (AllPulseEnergy > 200) {

                        num_chan_vp += 1;

                    }

                    if (iChan == SiPM1Chan) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_16 = AllPulseEnergy;

                    }

                    else if (iChan == SiPM2Chan) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_17 = AllPulseEnergy;

                    }

                    else if (iChan == SiPM3Chan) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_18 = AllPulseEnergy;

                    }

                    else if (iChan == SiPM4Chan) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_19 = AllPulseEnergy;

                    }

                    else if (iChan == SiPM5Chan) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_20 = AllPulseEnergy;

                    }

                    else if (iChan == SiPM6Chan) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_21 = AllPulseEnergy;

                    }

                    else if (iChan == SiPM7Chan) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_22 = AllPulseEnergy;

                    }

                    else if (iChan == SiPM8Chan) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_23 = AllPulseEnergy;

                    }

                }

                if (iChan == SiPM9Chan || iChan == SiPM10Chan) {

                    if (!top_vp_event && AllPulseEnergy > 200) {

                        num_chan_vp += 1;

                        top_vp_event = true;
                        
                    }

                    if (iChan == SiPM9Chan) {

                        top_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_24 = 1.07809 * AllPulseEnergy;

                    }

                    else if (iChan == SiPM10Chan) {

                        top_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_25 = AllPulseEnergy;

                    }

                }

                if (!peak_in_any_chan) {all_chan_start.push_back(1000); if (iChan == Event61Chan) {Ev61Peak = 0.0;}}

                AllPulseEnergy = 0.;

                h_wf->Reset();

                chan_lengths.push_back(nSamples[iChan]);

            }   // Channel loop

            double first_peak_time = 0.0;

            if (std::accumulate(all_chan_start.begin(), all_chan_start.end(), 0.0) != 26000) {first_peak_time = *min_element(all_chan_start.begin(), all_chan_start.end());}

            long long pulse_start_time = mostFrequent(all_chan_start_pmt);

            long long pulse_end_time = mostFrequent(all_chan_end);

            long long vp_start_time = nsTime + mostFrequent(all_chan_start_vp);

            double chan_15_start_adj = chan_15_start;

            if (first_peak_time == 0.0) {pulse_start_time = 0.0; pulse_end_time = 0.0; vp_start_time = nsTime; chan_15_start_adj = 0.0;}

            else {

                for (size_t iPeak = 0; iPeak < all_chan_start.size(); iPeak++) {

                    if (iPeak <= 11 && all_chan_start[iPeak] != 1000) {all_chan_start_pmt_adj.push_back(all_chan_start[iPeak] - first_peak_time);}

                    else if (iPeak >= SiPM1Chan && iPeak <= SiPM10Chan && all_chan_start[iPeak] != 1000) {all_chan_start_vp_adj.push_back(all_chan_start[iPeak] - first_peak_time);}

                    else if (iPeak == Event61Chan && all_chan_start[iPeak] != 1000) {chan_15_start_adj = all_chan_start[iPeak] - first_peak_time;}

                    else if (iPeak == Event61Chan && all_chan_start[iPeak] == 1000) {chan_15_start_adj = 0.0;}

                }

                pulse_start_time = mostFrequent(all_chan_start_pmt_adj);

                vp_start_time = nsTime + mostFrequent(all_chan_start_vp_adj);

            }

            if (num_chan_over_1phe >= 10) {int pulse_start_index = mostFrequentIndex(all_chan_start_pmt);}

            // if (pulse_start_time == 0) {continue;}

            struct pulse avg_pulse;

            avg_pulse.file_num = run_iterable;
            avg_pulse.entry = iEnt;
            avg_pulse.unix_time = run_starttime;
            avg_pulse.start = nsTime + pulse_start_time;                                                   // Need a way to weed out secondary pulses - use weighted average?
            avg_pulse.end = nsTime + pulse_end_time;
            avg_pulse.energy = std::accumulate(all_chan_energy.begin(), all_chan_energy.end(), 0.0) / (currentHLmean / defaultHLmean);           // Integrals still sum over secondary peaks discarded in .start and .end values
            avg_pulse.peak = std::accumulate(all_chan_peak.begin(), all_chan_peak.end(), 0.0);
            avg_pulse.number = num_chan_over_1phe;
            avg_pulse.beam = all_chan_beam;
            avg_pulse.trigger = triggerBits;
            avg_pulse.length = getAverage(chan_lengths);
            avg_pulse.side_vp_energy = std::accumulate(side_veto_panel_energy.begin(), side_veto_panel_energy.end(), 0.0);
            avg_pulse.top_vp_energy = std::accumulate(top_veto_panel_energy.begin(), top_veto_panel_energy.end(), 0.0);
            avg_pulse.all_vp_energy = avg_pulse.side_vp_energy + avg_pulse.top_vp_energy;
            avg_pulse.wf_time = nsTime;
            avg_pulse.vp_start_16 = p_start_16;
            avg_pulse.vp_start_17 = p_start_17;
            avg_pulse.vp_start_18 = p_start_18;
            avg_pulse.vp_start_19 = p_start_19;
            avg_pulse.vp_start_20 = p_start_20;
            avg_pulse.vp_start_21 = p_start_21;
            avg_pulse.vp_start_22 = p_start_22;
            avg_pulse.vp_start_23 = p_start_23;
            avg_pulse.vp_start_top = 0.5 * (p_start_24 + p_start_25);
            avg_pulse.vp_int_16 = p_int_16;
            avg_pulse.vp_int_17 = p_int_17;
            avg_pulse.vp_int_18 = p_int_18;
            avg_pulse.vp_int_19 = p_int_19;
            avg_pulse.vp_int_20 = p_int_20;
            avg_pulse.vp_int_21 = p_int_21;
            avg_pulse.vp_int_22 = p_int_22;
            avg_pulse.vp_int_23 = p_int_23;
            avg_pulse.vp_int_top = p_int_24 + p_int_25;
            avg_pulse.problem = false;
            avg_pulse.ev61_start_time = nsTime + chan_15_start;

            for (int iPeak = 0; iPeak < all_chan_start.size(); iPeak++) {

                if (all_chan_start[iPeak] < (pulse_start_time + 10 * 16) && all_chan_start[iPeak] > (pulse_start_time - 10 * 16)) {chan_start_no_outliers.push_back(all_chan_start[iPeak]);}

                else {

                    // Check if pulse outside this range occurs more than once

                    for (int jPeak = 0; jPeak < all_chan_start.size(); jPeak++) {

                        if (iPeak == jPeak) {continue;}

                        else if (all_chan_start[iPeak] < (all_chan_start[jPeak] + 1 * 16) && all_chan_start[iPeak] > (all_chan_start[jPeak] - 1 * 16)) {chan_start_no_outliers.push_back(all_chan_start[iPeak]);}

                    }

                }
                
            }

            double var_val = variance(chan_start_no_outliers);

            if (var_val < 5 * 16) {avg_pulse.single = true;}            // Pulses largely not consistent, need to weed out outlier and secondary pulses

            else {avg_pulse.single = false;}

            // Is this event a muon? (Looking at events with peak at the very end of the waveform)

            // if (pulse_at_end && (p_int_16 > 750 / 2 || p_int_17 > 950 / 2 || p_int_18 > 1200 / 2 || p_int_19 > 1375 / 2 || p_int_20 > 525 / 2 || p_int_21 > 700 / 2 || p_int_22 > 700 / 2 || p_int_23 > 500 / 2 || avg_pulse.top_vp_energy > 450 / 2)) {         // * (currentHLmean / defaultHLmean)

            if (pulse_at_end && (p_int_16 > min_muon_vals[0] / 2 || p_int_17 > min_muon_vals[1] / 2 || p_int_18 > min_muon_vals[2] / 2 || p_int_19 > min_muon_vals[3] / 2 || p_int_20 > min_muon_vals[4] / 2 || p_int_21 > min_muon_vals[5] / 2 || p_int_22 > min_muon_vals[6] / 2 || p_int_23 > min_muon_vals[7] / 2 || avg_pulse.vp_int_top > min_muon_vals[8] / 2)) {
                
                avg_pulse.last_cosmics_time = vp_start_time;

                last_cosmics_counter += 1;

            }

            // Is this event a muon? (Looking at "normal" events)

            // else if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {            // * (currentHLmean / defaultHLmean)

            else if (p_int_16 > min_muon_vals[0] || p_int_17 > min_muon_vals[1] || p_int_18 > min_muon_vals[2] || p_int_19 > min_muon_vals[3] || p_int_20 > min_muon_vals[4] || p_int_21 > min_muon_vals[5] || p_int_22 > min_muon_vals[6] || p_int_23 > min_muon_vals[7] || avg_pulse.vp_int_top > min_muon_vals[8]) {

                avg_pulse.last_cosmics_time = vp_start_time;

                last_cosmics_counter += 1;

            }

            // Record every high light LED event, enable recording of vector pulses

            if (avg_pulse.trigger == 8) {

                if (avg_pulse.energy < 0 || avg_pulse.peak >= maxPeak || avg_pulse.length != 45 || !avg_pulse.single || avg_pulse.trigger != 8) {avg_pulse.problem = true;}

                pulses.push_back(avg_pulse);
                
                record_pulses = true;
                
            }
            
            // Record every low light LED event, disable recording of vector pulses
            
            if (avg_pulse.trigger == 16) {
                
                if (avg_pulse.energy < 0 || avg_pulse.peak >= maxPeak || avg_pulse.length != 45 || !avg_pulse.single || avg_pulse.trigger != 16) {avg_pulse.problem = true;}

                pulses.push_back(avg_pulse);

                if (pulses.size() < 3 || (pulses.size() == 3 && pulses[0].trigger == 8 && (pulses[1].trigger == 1 || pulses[1].trigger == 33) && pulses[2].trigger == 16)) {}

                else {

                    for (size_t iVec = 0; iVec < pulses.size(); iVec++) {

                        Int_t fnum = pulses[iVec].file_num;
                        Int_t entry = pulses[iVec].entry;
                        Long64_t utime = pulses[iVec].unix_time;
                        Long64_t st = pulses[iVec].start;
                        Long64_t ed = pulses[iVec].end;
                        Double_t pk = pulses[iVec].peak;
                        Double_t ey = pulses[iVec].energy;
                        Int_t nr = pulses[iVec].number;
                        Bool_t se = pulses[iVec].single;
                        Bool_t bm = pulses[iVec].beam;
                        Int_t tr = pulses[iVec].trigger;
                        Int_t lh = pulses[iVec].length;
                        Double_t svpe = pulses[iVec].side_vp_energy;
                        Double_t tvpe = pulses[iVec].top_vp_energy;
                        Double_t avpe = pulses[iVec].all_vp_energy;
                        Long64_t lct = pulses[iVec].last_cosmics_time;
                        Long64_t ldt = pulses[iVec].last_det_time;
                        Long64_t wftime = pulses[iVec].wf_time;
                        Long64_t vptime16 = pulses[iVec].vp_start_16;
                        Long64_t vptime17 = pulses[iVec].vp_start_17;
                        Long64_t vptime18 = pulses[iVec].vp_start_18;
                        Long64_t vptime19 = pulses[iVec].vp_start_19;
                        Long64_t vptime20 = pulses[iVec].vp_start_20;
                        Long64_t vptime21 = pulses[iVec].vp_start_21;
                        Long64_t vptime22 = pulses[iVec].vp_start_22;
                        Long64_t vptime23 = pulses[iVec].vp_start_23;
                        Long64_t vptimetop = pulses[iVec].vp_start_top;
                        Long64_t vpint16 = pulses[iVec].vp_int_16;
                        Long64_t vpint17 = pulses[iVec].vp_int_17;
                        Long64_t vpint18 = pulses[iVec].vp_int_18;
                        Long64_t vpint19 = pulses[iVec].vp_int_19;
                        Long64_t vpint20 = pulses[iVec].vp_int_20;
                        Long64_t vpint21 = pulses[iVec].vp_int_21;
                        Long64_t vpint22 = pulses[iVec].vp_int_22;
                        Long64_t vpint23 = pulses[iVec].vp_int_23;
                        Long64_t vpinttop = pulses[iVec].vp_int_top;
                        Bool_t issue = pulses[iVec].problem;
                        Long64_t e61st = pulses[iVec].ev61_start_time;
                        Bool_t e61sig = pulses[iVec].ev61_significant;
                        Int_t selcut = pulses[iVec].selection_cut;

                        eventTree->SetBranchAddress("fnum", &fnum);
                        eventTree->SetBranchAddress("entry", &entry);
                        eventTree->SetBranchAddress("utime", &utime);
                        eventTree->SetBranchAddress("st", &st);
                        eventTree->SetBranchAddress("ed", &ed);
                        eventTree->SetBranchAddress("pk", &pk);
                        eventTree->SetBranchAddress("ey", &ey);
                        eventTree->SetBranchAddress("nr", &nr);
                        eventTree->SetBranchAddress("se", &se);
                        eventTree->SetBranchAddress("bm", &bm);
                        eventTree->SetBranchAddress("tr", &tr);
                        eventTree->SetBranchAddress("lh", &lh);
                        eventTree->SetBranchAddress("svpe", &svpe);
                        eventTree->SetBranchAddress("tvpe", &tvpe);
                        eventTree->SetBranchAddress("avpe", &avpe);
                        eventTree->SetBranchAddress("lct", &lct);
                        eventTree->SetBranchAddress("ldt", &ldt);
                        eventTree->SetBranchAddress("issue", &issue);
                        eventTree->SetBranchAddress("wftime", &wftime);
                        eventTree->SetBranchAddress("vptime16", &vptime16);
                        eventTree->SetBranchAddress("vptime17", &vptime17);
                        eventTree->SetBranchAddress("vptime18", &vptime18);
                        eventTree->SetBranchAddress("vptime19", &vptime19);
                        eventTree->SetBranchAddress("vptime20", &vptime20);
                        eventTree->SetBranchAddress("vptime21", &vptime21);
                        eventTree->SetBranchAddress("vptime22", &vptime22);
                        eventTree->SetBranchAddress("vptime23", &vptime23);
                        eventTree->SetBranchAddress("vptimetop", &vptimetop);
                        eventTree->SetBranchAddress("vpint16", &vpint16);
                        eventTree->SetBranchAddress("vpint17", &vpint17);
                        eventTree->SetBranchAddress("vpint18", &vpint18);
                        eventTree->SetBranchAddress("vpint19", &vpint19);
                        eventTree->SetBranchAddress("vpint20", &vpint20);
                        eventTree->SetBranchAddress("vpint21", &vpint21);
                        eventTree->SetBranchAddress("vpint22", &vpint22);
                        eventTree->SetBranchAddress("vpint23", &vpint23);
                        eventTree->SetBranchAddress("vpinttop", &vpinttop);
                        eventTree->SetBranchAddress("e61st", &e61st);
                        eventTree->SetBranchAddress("e61sig", &e61sig);
                        eventTree->SetBranchAddress("selcut", &selcut);

                        eventTree->Fill();

                    }

                }
                
                record_pulses = false; pulses.clear();
                
            }

            // Is a given Event 61 event potentially a physical signal?

            avg_pulse.ev61_significant = false;

            avg_pulse.selection_cut = 0;

            if (record_pulses && avg_pulse.beam && avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {

                if (avg_pulse.number >= 10 && avg_pulse.energy >= 50 && avg_pulse.start - last_event_time > echo_cut) {

                    avg_pulse.selection_cut = 1;

                    avg_pulse.ev61_significant = true;

                    // if (p_int_16 <= 750 && p_int_17 <= 950 && p_int_18 <= 1200 && p_int_19 <= 1375 && p_int_20 <= 525 && p_int_21 <= 700 && p_int_22 <= 700 && p_int_23 <= 500 && avg_pulse.top_vp_energy <= 450) {

                    if (p_int_16 <= min_muon_vals[0] && p_int_17 <= min_muon_vals[1] && p_int_18 <= min_muon_vals[2] && p_int_19 <= min_muon_vals[3] && p_int_20 <= min_muon_vals[4] && p_int_21 <= min_muon_vals[5] && p_int_22 <= min_muon_vals[6] && p_int_23 <= min_muon_vals[7] && avg_pulse.vp_int_top <= min_muon_vals[8]) {

                        avg_pulse.selection_cut = 2;
                            
                        if (avg_pulse.start - avg_pulse.last_det_time > muon_dt_max) {
                                
                            avg_pulse.selection_cut = 3;

                            if (avg_pulse.start - avg_pulse.last_cosmics_time > muon_dt_max) {

                                avg_pulse.selection_cut = 4;

                                if (p_int_16 <= vp_cut_vals[0] && p_int_17 <= vp_cut_vals[1] && p_int_18 <= vp_cut_vals[2] && p_int_19 <= vp_cut_vals[3] && p_int_20 <= vp_cut_vals[4] && p_int_21 <= vp_cut_vals[5] && p_int_22 <= vp_cut_vals[6] && p_int_23 <= vp_cut_vals[7] && avg_pulse.vp_int_top <= vp_cut_vals[8]) {

                                    avg_pulse.selection_cut = 5;

                                    if (avg_pulse.energy <= 500) {

                                        avg_pulse.selection_cut = 6;

                                        if (avg_pulse.start - avg_pulse.last_cosmics_time > muon_dt_extra) {

                                            avg_pulse.selection_cut = 7;

                                        }

                                    }

                                }

                            }
                            
                        }

                    }

                }
            
            }

            // Record every Event 61 event

            if (record_pulses && avg_pulse.beam && (avg_pulse.trigger == 1 || avg_pulse.trigger == 3 || avg_pulse.trigger == 33 || avg_pulse.trigger == 35)) {
                
                if (avg_pulse.energy < 0 || avg_pulse.peak >= maxPeak || avg_pulse.length != 45 || !avg_pulse.single || (avg_pulse.trigger != 1 && avg_pulse.trigger != 3 && avg_pulse.trigger != 33 && avg_pulse.trigger != 35)) {avg_pulse.problem = true;}

                pulses.push_back(avg_pulse);

            }

            // Record internally triggered detector events, all if statements are just cuts on which events I record

            avg_pulse.selection_cut = 0;

            if (record_pulses && !avg_pulse.beam && avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {
                
                if (avg_pulse.number >= 10 && avg_pulse.energy >= 50 && avg_pulse.start - last_event_time > echo_cut) {

                    avg_pulse.selection_cut = 1;

                    // if (p_int_16 <= 750 && p_int_17 <= 950 && p_int_18 <= 1200 && p_int_19 <= 1375 && p_int_20 <= 525 && p_int_21 <= 700 && p_int_22 <= 700 && p_int_23 <= 500 && avg_pulse.top_vp_energy <= 450) {

                    if (p_int_16 <= min_muon_vals[0] && p_int_17 <= min_muon_vals[1] && p_int_18 <= min_muon_vals[2] && p_int_19 <= min_muon_vals[3] && p_int_20 <= min_muon_vals[4] && p_int_21 <= min_muon_vals[5] && p_int_22 <= min_muon_vals[6] && p_int_23 <= min_muon_vals[7] && avg_pulse.vp_int_top <= min_muon_vals[8]) {

                        avg_pulse.selection_cut = 2;
                            
                        if (avg_pulse.start - avg_pulse.last_det_time > muon_dt_max) {

                            avg_pulse.selection_cut = 3;

                            if (avg_pulse.start - avg_pulse.last_cosmics_time > muon_dt_max) {

                                avg_pulse.selection_cut = 4;

                                if (p_int_16 <= vp_cut_vals[0] && p_int_17 <= vp_cut_vals[1] && p_int_18 <= vp_cut_vals[2] && p_int_19 <= vp_cut_vals[3] && p_int_20 <= vp_cut_vals[4] && p_int_21 <= vp_cut_vals[5] && p_int_22 <= vp_cut_vals[6] && p_int_23 <= vp_cut_vals[7] && avg_pulse.vp_int_top <= vp_cut_vals[8]) {

                                    avg_pulse.selection_cut = 5;

                                    if (avg_pulse.energy <= 500) {

                                        avg_pulse.selection_cut = 6;

                                        if (avg_pulse.start - avg_pulse.last_cosmics_time > muon_dt_extra) {

                                            avg_pulse.selection_cut = 7;

                                        }

                                    }

                                }

                            }

                        }

                    }

                }

            }

            if (record_pulses && !avg_pulse.beam && avg_pulse.selection_cut >= 5) {

                if (avg_pulse.energy < 0 || avg_pulse.peak >= maxPeak || avg_pulse.length != 45 || !avg_pulse.single || (avg_pulse.trigger != 2 && avg_pulse.trigger != 34)) {avg_pulse.problem = true;}

                pulses.push_back(avg_pulse);
                
            }

            // Is this event a possible Michel event?

            if (avg_pulse.number >= 10 && avg_pulse.energy >= 50 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {avg_pulse.last_det_time = avg_pulse.start; last_det_counter += 1;}         // * (currentHLmean / defaultHLmean)

            // Are large energy detector events also being detected in veto panels?

            if (avg_pulse.number >= 10 && avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {

                if (avg_pulse.energy > 0 && avg_pulse.energy <= 200) {num_large_found_1 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_1 += 1;}}            // * (currentHLmean / defaultHLmean)
                
                if (avg_pulse.energy > 200 && avg_pulse.energy <= 400) {num_large_found_2 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_2 += 1;}}

                if (avg_pulse.energy > 400 && avg_pulse.energy <= 600) {num_large_found_3 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_3 += 1;}}

                if (avg_pulse.energy > 600 && avg_pulse.energy <= 800) {num_large_found_4 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_4 += 1;}}

                if (avg_pulse.energy > 800 && avg_pulse.energy <= 1000) {num_large_found_5 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_5 += 1;}}

                if (avg_pulse.energy > 1000 && avg_pulse.energy <= 1200) {num_large_found_6 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_6 += 1;}}
            
                if (avg_pulse.energy > 1200 && avg_pulse.energy <= 1400) {num_large_found_7 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_7 += 1;}}
            
                if (avg_pulse.energy > 1400 && avg_pulse.energy <= 1600) {num_large_found_8 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_8 += 1;}}
            
                if (avg_pulse.energy > 1600 && avg_pulse.energy <= 1800) {num_large_found_9 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_9 += 1;}}
            
                if (avg_pulse.energy > 1800 && avg_pulse.energy <= 2000) {num_large_found_10 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_10 += 1;}}
            
            }

            last_event_time = avg_pulse.start; last_event_counter += 1;

            if (avg_pulse.start - avg_pulse.last_det_time > muon_dt_max && avg_pulse.start - avg_pulse.last_cosmics_time > muon_dt_max) {
                
                h_intev61->Fill(Ev61Energy);

                if (avg_pulse.trigger == 16) {h_intlowLED->Fill(avg_pulse.energy); intlowLED.push_back(avg_pulse.energy);}

                else if (avg_pulse.trigger == 8) {h_inthighLED->Fill(avg_pulse.energy); inthighLED.push_back(avg_pulse.energy);}

                else if (avg_pulse.trigger == 4) {h_intminLED->Fill(avg_pulse.energy); intminLED.push_back(avg_pulse.energy);}
                
            }

            // Reset variables at end of loop

            p_int_16 = 0.0;
            p_int_17 = 0.0;
            p_int_18 = 0.0;
            p_int_19 = 0.0;
            p_int_20 = 0.0;
            p_int_21 = 0.0;
            p_int_22 = 0.0;
            p_int_23 = 0.0;
            p_int_24 = 0.0;
            p_int_25 = 0.0;
            p_start_16 = 0.0;
            p_start_17 = 0.0;
            p_start_18 = 0.0;
            p_start_19 = 0.0;
            p_start_20 = 0.0;
            p_start_21 = 0.0;
            p_start_22 = 0.0;
            p_start_23 = 0.0;
            p_start_24 = 0.0;
            p_start_25 = 0.0;
            num_chan = 0;
            num_chan_vp = 0;
            num_chan_simple = 0;
            num_chan_over_1phe = 0;

            all_chan_beam = false;
            top_vp_event = false;

            all_chan_start.clear();
            all_chan_start_adj.clear();
            all_chan_end.clear();
            all_chan_peak.clear();
            all_chan_energy.clear();
            peak_energy.clear();
            ten_chan_peak.clear();
            all_chan_peakbin.clear();
            ten_chan_peakbin.clear();
            chan_start_no_outliers.clear();
            side_veto_panel_energy.clear();
            top_veto_panel_energy.clear();
            all_chan_start_pmt.clear();
            all_chan_start_vp.clear();
            all_chan_start_pmt_adj.clear();
            all_chan_start_vp_adj.clear();
            chan_lengths.clear();

            event_bin_val.clear();

        } // Event loop

        eventTree->Write();

        fileOut->Close();

        ch15_on = 0;

        integralToPE.clear();

        avg_peak_pos_RMS = getAverage(peak_pos_RMS);

        var_peak_pos_RMS = variance(peak_pos_RMS);

        peak_pos_RMS.clear();

        // cout << "\n" << "Completed run " << run_iterable << endl;

        f->Close();

    } // Run loop

    TCanvas* c_intev61 = new TCanvas("c_intev61", "Distribution of Event 61 Integral Values", 1200, 700);
    c_intev61->SetLogy();
    c_intev61->cd();
    h_intev61->GetXaxis()->SetTitle("Integral (ADC)");
    h_intev61->GetYaxis()->SetTitle("Counts");
    h_intev61->Draw();
    c_intev61->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_intev61_%i.png", run));
    c_intev61->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_intev61.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_intev61.png"));
    h_intev61->Reset();

    TCanvas* c_intlowLED = new TCanvas("c_intlowLED", "Low Light LED Integral Values", 1200, 700);
    c_intlowLED->SetLogy();
    c_intlowLED->cd();
    h_intlowLED->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intlowLED->GetYaxis()->SetTitle("Counts");
    h_intlowLED->Draw();
    TLatex *lat_intlowLED = new TLatex(0.75, 0.7, Form("RMS = %.2f", rmsValue(intlowLED)));
    lat_intlowLED->SetNDC();
    lat_intlowLED->SetTextColor(1);
    lat_intlowLED->SetTextSize(0.035);
    lat_intlowLED->Draw();
    c_intlowLED->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_intlowLED_%i.png", run));
    c_intlowLED->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_intlowLED.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intlowLED.png"));
    h_intlowLED->Reset();

    TCanvas* c_inthighLED = new TCanvas("c_inthighLED", "High Light LED Integral Values", 1200, 700);
    c_inthighLED->SetLogy();
    c_inthighLED->cd();
    h_inthighLED->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_inthighLED->GetYaxis()->SetTitle("Counts");
    h_inthighLED->Draw();
    TLatex *lat_inthighLED = new TLatex(0.75, 0.7, Form("RMS = %.2f", rmsValue(inthighLED)));
    lat_inthighLED->SetNDC();
    lat_inthighLED->SetTextColor(1);
    lat_inthighLED->SetTextSize(0.035);
    lat_inthighLED->Draw();
    c_inthighLED->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_inthighLED_%i.png", run));
    c_inthighLED->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_inthighLED.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_inthighLED.png"));
    h_inthighLED->Reset();

    TCanvas* c_intminLED = new TCanvas("c_intminLED", "Min Bias LED Integral Values", 1200, 700);
    c_intminLED->SetLogy();
    c_intminLED->cd();
    h_intminLED->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intminLED->GetYaxis()->SetTitle("Counts");
    h_intminLED->Draw();
    TLatex *lat_intminLED = new TLatex(0.75, 0.7, Form("RMS = %.2f", rmsValue(intminLED)));
    lat_intminLED->SetNDC();
    lat_intminLED->SetTextColor(1);
    lat_intminLED->SetTextSize(0.035);
    lat_intminLED->Draw();
    c_intminLED->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_intminLED_%i.png", run));
    c_intminLED->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_intminLED.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intminLED.png"));
    h_intminLED->Reset();
    
    cout << "\n" << "Number of cosmic events = " << last_cosmics_counter << endl;

    cout << "\n" << "Number of potential Michel events = " << last_det_counter << endl;

    cout << "\n" << "Number of all events = " << last_event_counter << endl;

    double veto_det_ratio_1 = num_veto_found_1 * pow(num_large_found_1, -1);

    double veto_det_ratio_2 = num_veto_found_2 * pow(num_large_found_2, -1);

    double veto_det_ratio_3 = num_veto_found_3 * pow(num_large_found_3, -1);

    double veto_det_ratio_4 = num_veto_found_4 * pow(num_large_found_4, -1);

    double veto_det_ratio_5 = num_veto_found_5 * pow(num_large_found_5, -1);

    double veto_det_ratio_6 = num_veto_found_6 * pow(num_large_found_6, -1);

    double veto_det_ratio_7 = num_veto_found_7 * pow(num_large_found_7, -1);

    double veto_det_ratio_8 = num_veto_found_8 * pow(num_large_found_8, -1);

    double veto_det_ratio_9 = num_veto_found_9 * pow(num_large_found_9, -1);

    double veto_det_ratio_10 = num_veto_found_10 * pow(num_large_found_10, -1);

    cout << "\n" << "num_veto_found_1 / num_large_found_1 = " << num_veto_found_1 << " / " << num_large_found_1 << " = " << veto_det_ratio_1 << endl;
    
    cout << "\n" << "num_veto_found_2 / num_large_found_2 = " << num_veto_found_2 << " / " << num_large_found_2 << " = " << veto_det_ratio_2 << endl;

    cout << "\n" << "num_veto_found_3 / num_large_found_3 = " << num_veto_found_3 << " / " << num_large_found_3 << " = " << veto_det_ratio_3 << endl;

    cout << "\n" << "num_veto_found_4 / num_large_found_4 = " << num_veto_found_4 << " / " << num_large_found_4 << " = " << veto_det_ratio_4 << endl;

    cout << "\n" << "num_veto_found_5 / num_large_found_5 = " << num_veto_found_5 << " / " << num_large_found_5 << " = " << veto_det_ratio_5 << endl;

    cout << "\n" << "num_veto_found_6 / num_large_found_6 = " << num_veto_found_6 << " / " << num_large_found_6 << " = " << veto_det_ratio_6 << endl;

    cout << "\n" << "num_veto_found_7 / num_large_found_7 = " << num_veto_found_7 << " / " << num_large_found_7 << " = " << veto_det_ratio_7 << endl;

    cout << "\n" << "num_veto_found_8 / num_large_found_8 = " << num_veto_found_8 << " / " << num_large_found_8 << " = " << veto_det_ratio_8 << endl;

    cout << "\n" << "num_veto_found_9 / num_large_found_9 = " << num_veto_found_9 << " / " << num_large_found_9 << " = " << veto_det_ratio_9 << endl;

    cout << "\n" << "num_veto_found_10 / num_large_found_10 = " << num_veto_found_10 << " / " << num_large_found_10 << " = " << veto_det_ratio_10 << endl;

    cout << "\n" << "End of code." << "\n" << endl;

    return 0;

}
