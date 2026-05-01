#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPad.h>
#include <iostream>
#include <fstream>
#include <string>
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
  double start;   /* Universal start time of pulse (10% peak) in waveform (ns) */
  double end;     /* Universal end time of pulse (reach baseline) in waveform (ns) */
  double peak;    /* Max amplitude of pulse (photo-electrons) */
  double energy;  /* Energy (integral) of pulse (photo-electrons) */
  double number;  /* Number of channels in which we see pulse (photo-electrons) */
  bool single;    /* Is pulse (photo-electrons) timing consistent across all channels */
  bool beam;      /* Tracks whether beam is on or off */
  double trigger; /* Tracks whether trigger is external (2) or internal (16) */
  double length;  /* Length of waveform in number of bins */
  double side_vp_energy; /* Energy (integral) of pulse (photo-electrons) in SIDE veto panels */
  double top_vp_energy; /* Energy (integral) of pulse (photo-electrons) in TOP veto panel */
  long long last_muon_time; /* Time value of most recent detected muon event */
  long long last_lead_time; /* Time value of most recent detected muon stopping in lead event */
  long long last_det_time; /* Time value of most recent detected detector event */
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
const char* outputFilePath = "";                            // What does "" mean?
const char* outputFileName = "MichelEvents";         // Why is output file going into ~ directory?
const char* outputStatsName = "PMTAnalysisStats";
// double integralToPE = 163.43;
// double amplitudeToPE = 60.33;

// std::vector<double> amplitudeToPE{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15.9,16.4,15.65,15.9,16.3,0,10.0,15.3,14.9,0,14.75,14,14.45,15.35 };                    // 10/5: run2757

// std::vector<double> amplitudeToPE{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15.4,16.5,15.95,15.3,15.7,11.0,14.95,14.8,16.1,16.0,16.15,16.25 };                        // 11/7: run3529

// std::vector<double> amplitudeToPE{ 28,27,30,29,26,11,29,28,30,29,24,29 };             // 12/2: run4144

std::vector<double> amplitudeToPE{ 26,28,30,29,27,12,30,28,30,30,24,30 };                // 12/3: run4176

// std::vector<double> amplitudeToPE{ 25,27,30,29,25,10,27,27,30,30,23,30 };             // 12/4: run4193

// std::vector<double> integralToPETemp{ 102.068,104.926,103.562,103.482,106.928,69.5778,101.677,100.208,100.963,120.495,102.519,104.107 };                     // 01/23: run2001

std::vector<double> integralToPETemp{ 103.536,102.226,108.643,101.837,107.954,109.53,111.911,102.516,104.326,90.6458,105.732,105.475 };                         // runs 25120 - 25143

std::vector<double> MichelSiPMNoiseAvg1{ 187.887,188.1116667,314.3093333,195.7781667,190.6345,202.5948333,209.708,201.6861667,243.2 };                          // Jul 17, 2023 - Aug 29, 2023

std::vector<double> MichelSiPMNoiseAvg2{ 185.7410476,174.2948095,461.303619,192.1831429,190.1772857,189.2594762,187.5979524,195.2732857,282.2847619 };          // Sep 3, 2023 - Jan 30, 2024

std::vector<double> MichelSiPMNoiseAvg3{ 220.497,202.8251905,545.8977619,220.8942857,220.4852381,222.7695238,220.9086667,231.9756667,300.433381 };              // Feb 9, 2024 - Jul 21, 2024

std::vector<double> MichelSiPMNoiseAvg4{ 222.2114286,213.4123571,231.4738571,222.3042143,220.602,234.8156429,237.2396429,240.2786429,282.2042143 };             // Jul 31, 2024 - Nov 23, 2024

std::vector<double> MichelSiPMNoiseAvg5{ 214.9131176,219.8434706,227.0522353,211.6268824,215.7442353,222.5197059,217.8758824,230.3478235,258.0880588 };         // Jan 12, 2025 - May 24, 2025

std::vector<double> MichelSiPMNoiseAvg6{ 202.2175,214.0718,190.558875,202.427125,209.80825,212.1195,195.1185,225.44875,227.461 };                               // May 24, 2025 - Aug 23, 2025

// 251.0594118,246.5616875,243.1353529,238.1902941,234.0304118,222.7274706,213.5342941,240.3322941,354.5087647

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
    cout << "\n" << "Run: " << run << " Last Run: " << last_run << "\n" << endl;

    int ADCSIZE = 45;
    TH1D *h_wf = new TH1D("h_wf", "Waveform", ADCSIZE, 0, ADCSIZE);

    int write_out_counter = 0;

    int run_counter = 0;

    int hour12_run = 0;

    int HL_dt_17ms_counter = 0;

    Long64_t hour12_st = 0.0;

    bool int_trig = false;

    bool min_bias_trig = false;

    bool hi_li_trig = false;

    bool low_li_trig = false;

    bool michel_event = false;

    int cosmic_events = 0;

    int hist_events = 0;

    double tot_michel_energy = 0.0;

    /* Create vectors for plotting histograms */

    std::vector<double> pulse_time_sep;

    std::vector<double> p1_trig_type;

    std::vector<double> p2_trig_type;

    std::vector<int> trig_type;

    std::vector<double> peak_height;

    std::vector<double> peak_time;

    std::vector<int> multiplicity;

    std::vector<int> multiplicity_ext;

    std::vector<int> multiplicity_int;

    std::vector<int> multiplicity_data;

    std::vector<double> pulse_time_l1s;

    std::vector<double> event_time_l1s;

    std::vector<int> event_num;

    std::vector<double> all_nsTime;

    std::vector<double> no_muon_cut_time_count;

    std::vector<double> simult_muon_cut_time_count;

    std::vector<double> intall;

    std::vector<double> intlowLED;

    std::vector<double> inthighLED;

    std::vector<double> intminLED;

    std::vector<double> intallint;

    std::vector<double> intallmultcut;

    std::vector<double> intallmultmuoncut;

    int ev_61_dt_max = pow(10, 4); // 2 * pow(10, 7);

    int muon_dt_max = 16 * pow(10, 3);

    int no_muons_time_cut = 5000;

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

    /* Initialize histograms */                 // NUMBER OF BINS SHOULD BE A MULTIPLE OF THE SMALLEST AXIS UNIT!!! (OR LENGTH OF AXIS DIVIDED BY SMALLEST AXIS UNIT?)

    TH1D* h_amp = new TH1D("h_amp", "Distribution of All Internally Triggered Event Amplitude Values", 2000, 0, 2000);

    TH1D* h_int = new TH1D("h_int", "Distribution of All Internally Triggered Event Integral Values", 5500, -500, 5000);

    TH2D* h_dten = new TH2D("h_dten", "Delta-T vs Event Number", 1000, 0, muon_dt_max, 1000, 0, MAX_NUM_ENTRIES);

    TH1D* h_tB = new TH1D("h_tB", "TriggerBits Distribution", 40, 0, 40);

    TH1D* h_blM_pmt = new TH1D("h_blM_pmt", "PMT BaselineMean Distribution", 400, -100, 300);

    TH1D* h_blM_vp = new TH1D("h_blM_vp", "SiPM BaselineMean Distribution", 400, -100, 300);

    TH1D* h_blM_ev61 = new TH1D("h_blM_ev61", "Event 61 BaselineMean Distribution", 400, -100, 300);

    TH1D* h_dt = new TH1D("h_dt", "Distribution of Time Separations Between Peaks", 1000, 0, muon_dt_max);

    TH1D* h_tc = new TH1D("h_tc", "Distribution of Time Separations Between Peaks (time_count)", 1000, 0, muon_dt_max);

    TH1D* h_time_bt_ev61 = new TH1D("h_time_bt_ev61", "Distribution of Time Separations Between Event 61 Peaks", 1000, 0, 6 * pow(10, 7));

    TH1D* h_time_bt_minLED = new TH1D("h_time_bt_minLED", "Distribution of Time Separations Between Min Bias LED Peaks", 1000, 0, 6 * pow(10, 7));

    TH1D* h_time_bt_lowLED = new TH1D("h_time_bt_lowLED", "Distribution of Time Separations Between Low Light LED Peaks", 1000, 0, 6 * pow(10, 7));

    TH1D* h_time_bt_highLED = new TH1D("h_time_bt_highLED", "Distribution of Time Separations Between High Light LED Peaks", 1000, 0, 6 * pow(10, 7));

    TH1D* h_dtint = new TH1D("h_dtint", "Distribution of Time Separations Between Internally Triggered Peaks", 1000, 0, muon_dt_max);

    TH1D* h_dtmin = new TH1D("h_dtmin", "Distribution of Time Separations Between Minimum Bias Triggered Peaks", 1000, 0, muon_dt_max);

    TH1D* h_dthlt = new TH1D("h_dthlt", "Distribution of Time Separations Between High Light Triggered Peaks", 1000, 0, muon_dt_max);

    TH1D* h_dtllt = new TH1D("h_dtllt", "Distribution of Time Separations Between Low Light Triggered Peaks", 1000, 0, muon_dt_max);

    TH1D* h_intall = new TH1D("h_intall", "Distribution of All Event Integral Values", 5500, -500, 5000);

    TH1D* h_intlowLED = new TH1D("h_intlowLED", "Distribution of Low Light LED Integral Values", 150, -500, 1000);

    TH1D* h_inthighLED = new TH1D("h_inthighLED", "Distribution of High Light LED Integral Values", 400, -2000, 2000);            //, 11000, -1000, 10000);

    TH1D* h_intminLED = new TH1D("h_intminLED", "Distribution of Minimum Bias LED Integral Values", 150, -500, 1000);

    TH1D* h_MB_HL_dt = new TH1D("h_MB_HL_dt", "Delta-t Between High Light Events and Most Recent Min Bias Event", 2250, 0, 36 * pow(10, 6));

    TH1D* h_intallint = new TH1D("h_intallint", "Distribution of All Internally Triggered Event Integral Values", 5500, -500, 5000);

    TH1D* h_intallmultcut = new TH1D("h_intallmultcut", "Distribution of All Internally Triggered Event Integral Values With Multiplicity 10+", 5500, -500, 5000);

    TH1D* h_intallmultmuoncut = new TH1D("h_intallmultmuoncut", "Distribution of All Internally Triggered Event Integral Values With Multiplicity 10+ And No Simultaneous Muons", 5500, -500, 5000);

    TH1D* h_ampdiff = new TH1D("h_ampdiff", "Distribution of Amplitude Differences Between Peak Pairs", 100, 0, 2000);

    TH2D* h_ampboth = new TH2D("h_ampboth", "Distribution of Amplitude Values", 100, 0, 2000, 100, 0, 2000);

    TH2D* h_intboth = new TH2D("h_intboth", "Distribution of Integral Values", 3000, 0, 3000, 3000, 0, 3000);

    TH1D* h_ampP1 = new TH1D("h_ampP1", "Distribution of First Peak Amplitude Values", 100, 0, 2000);

    TH1D* h_ampP2 = new TH1D("h_ampP2", "Distribution of Second Peak Amplitude Values", 100, 0, 2000);

    TH1D* h_ampOn = new TH1D("h_ampOn", "Distribution of Peak Amplitude Values for SNS Beam On vs Off", 100, 0, 200);

    TH1D* h_ampOff = new TH1D("h_ampOff", "Distribution of Peak Amplitude Values for SNS Beam On vs Off", 100, 0, 200);

    TH1D* h_ml = new TH1D("h_ml", "Internally Triggered Event Multiplicity", 17, -2, 15);

    TH2D* h_ml_v_int = new TH2D("h_ml_v_int", "Internally Triggered Event Multiplicity vs Integral Values", 17, -2, 15, 5500, -500, 5000);

    TH1D* h_pst = new TH1D("h_pst", "Distribution of Peak Start Times in Each Channel", 45, 0, 720);

    TH1D* h_mfpst = new TH1D("h_mfpst", "Distribution of Most Frequent Peak Start Times in Each Channel", 45, 0, 720);

    TH2D* h_pst_v_pint = new TH2D("h_pst_v_pint", "Distribution of Peak Start Times in Each Channel vs Peak Integral Values", 45, 0, 720, 100, 0, 1000);

    TH2D* h_mfpst_v_mfpint = new TH2D("h_mfpst_v_mfpint", "Distribution of Most Frequent Peak Start Times vs Most Frequent Peak Integral Values", 45, 0, 720, 100, 0, 1000);

    TH1D* h_var = new TH1D("h_var", "Distribution of Variance in Peak Start Times", 48, 0, 480);

    TH1D* h_spheLL = new TH1D("h_spheLL", "Distribution of Single Channel Low Light Integral Values", 1000, -500, 500);

    TH1D* h_spheMB = new TH1D("h_spheMB", "Distribution of Single Channel Minimum Bias Integral Values", 1000, -500, 500);

    TH1D* h_intmuons = new TH1D("h_intmuons", "Distribution of Muon Integral Values", 200, 0, 2000);

    TH1D* h_intmichels = new TH1D("h_intmichels", "Distribution of Michel Electron Integral Values", 80, 0, 800);

    TH2D* h_dt_v_intmichels = new TH2D("h_dt_v_intmichels", "Delta-T vs Michel PMT Integral Value", 1000, 0, muon_dt_max, 80, 0, 800);

    TH1D* h_intcosmics = new TH1D("h_intcosmics", "Distribution of Cosmic Integral Values", 400, 0, 4000);

    TH1D* h_intev61 = new TH1D("h_intev61", "Distribution of Event 61 Integral Values", 2200, -200, 2000);

    TH1D* h_intneu = new TH1D("h_intneu", "Distribution of Beam On Event Integral Values", 1650, -150, 1500);

    TH1D* h_dt_61 = new TH1D("h_dt_61", "Distribution of Time Separations Between Ev61 and Detector Peaks", 100, - ev_61_dt_max, ev_61_dt_max);

    TH1D* h_no_prior_muons_dt = new TH1D("h_no_prior_muons_dt", "Distribution of Time Separations Between Muons and Non-Muons", 100, 0, no_muons_time_cut + 5000);

    TH1D* h_int_61 = new TH1D("h_int_61", "Distribution of Detector Integral Values After Ev61 Events", 400, 0, 4000);

    TH1D* h_dt_61_vc = new TH1D("h_dt_61_vc", "Distribution of Time Separations Between Ev61 and Detector Peaks", 100, - ev_61_dt_max, ev_61_dt_max);

    TH1D* h_int_61_vc = new TH1D("h_int_61_vc", "Distribution of Detector Integral Values After Ev61 Events", 400, 0, 4000);

    TH1D* h_dt_61_smc = new TH1D("h_dt_61_smc", "Distribution of Time Separations Between Ev61 and Veto Panel Peaks", 100, - ev_61_dt_max, ev_61_dt_max);

    TH1D* h_int_61_smc = new TH1D("h_int_61_smc", "Distribution of Veto Panel Integral Values After Ev61 Events", 400, 0, 4000);

    TH1D* h_dt_61_pmc = new TH1D("h_dt_61_pmc", "Distribution of Time Separations Between Ev61 and Detector Peaks, Prior Muon Cut", 100, - ev_61_dt_max, ev_61_dt_max);

    TH1D* h_int_61_pmc = new TH1D("h_int_61_pmc", "Distribution of Detector Integral Values After Ev61 Events, Prior Muon Cut", 100, 0, 1000);

    TH1D* h_dt_61_spmc = new TH1D("h_dt_61_spmc", "Distribution of Time Separations Between Ev61 and Detector Peaks, Prior + Simult. Muon Cuts", 100, - ev_61_dt_max, ev_61_dt_max);

    TH1D* h_int_61_spmc = new TH1D("h_int_61_spmc", "Distribution of Detector Integral Values After Ev61 Events, Prior + Simult. Muon Cuts", 100, 0, 1000);

    TH1D* h_int_61_spmc_inner = new TH1D("h_int_61_spmc_inner", "Distribution of Detector Integral Values After Ev61 Events, Prior + Simult. Muon Cuts", 100, 0, 1000);

    TH1D* h_int_61_spmc_outer = new TH1D("h_int_61_spmc_outer", "Distribution of Detector Integral Values After Ev61 Events, Prior + Simult. Muon Cuts", 100, 0, 1000);

    TH2D* h_dt_v_int_61 = new TH2D("h_dt_v_int_61", "Time Separations between Ev61 and Detector Peaks vs Detector Integral Values After Ev61 Events", 100, - ev_61_dt_max, ev_61_dt_max, 100, -200, 4000);

    TH2D* h_dt_v_int_61_smc = new TH2D("h_dt_v_int_61_smc", "Time Separations between Ev61 and Detector Peaks vs Detector Integral Values After Ev61 Events, Simult. Muon Cut", 100, - ev_61_dt_max, ev_61_dt_max, 100, 0, 2000);

    TH2D* h_dt_v_int_61_pmc = new TH2D("h_dt_v_int_61_pmc", "Time Separations between Ev61 and Detector Peaks vs Detector Integral Values After Ev61 Events, Prior Muon Cut", 100, - ev_61_dt_max, ev_61_dt_max, 100, 0, 2000);

    TH2D* h_dt_v_int_61_spmc = new TH2D("h_dt_v_int_61_spmc", "Time Separations between Ev61 and Detector Peaks vs Detector Integral Values After Ev61 Events, Prior + Simult. Muon Cuts", 100, - ev_61_dt_max, ev_61_dt_max, 100, 0, 2000);

    TH1D* h_cent_spke_tB_smc = new TH1D("h_cent_spke_tB_smc", "Distribution of Detector Peak TriggerBits Values for Central Spike of Ev61-Det dt Plot, Simult. Muon Cut", 40, 0, 40);

    TH1D* h_cent_spke_tB_pmc = new TH1D("h_cent_spke_tB_pmc", "Distribution of Detector Peak TriggerBits Values for Central Spike of Ev61-Det dt Plot, Prior Muon Cut", 40, 0, 40);

    TH1D* h_cent_spke_tB_spmc = new TH1D("h_cent_spke_tB_spmc", "Distribution of Detector Peak TriggerBits Values for Central Spike of Ev61-Det dt Plot, Prior + Simult. Muon Cuts", 40, 0, 40);

    TH1D* h_tvp_int = new TH1D("h_tvp_int", "Distribution of Top Veto Panel Integral Values", 3000, 0, 3000);           // After Ev61 Events

    TH1D* h_neg_dt_ev61_peak_pos = new TH1D("h_neg_dt_ev61_peak_pos", "Distribution of Event 61 Peak Positions for Negative Event 61 - Detector Peak Time Separations", 45, 0, 720);

    TH1D* h_34detint = new TH1D("h_34detint", "Distribution of All triggerBits = 34 Detector Integral Values", 200, 0, 2000);

    TH1D* h_34vpint = new TH1D("h_34vpint", "Distribution of All triggerBits = 34 Veto Panel Integral Values", 200, 0, 2000);

    TH1D* h_slice_vpintavg = new TH1D("h_slice_vpintavg", "Distribution of Averaged Veto Panel Integral Values from Background Between 200-400 Ph.e.", 100, 0, 100);

    TH1D* h_dt_61_16 = new TH1D("h_dt_61_16", "Delta-t Between Event 61 and SiPM 1 Peaks", 17000, 0, ev_61_dt_max);                                              // 1 bin = 1,000 ns

    TH1D* h_dt_61_17 = new TH1D("h_dt_61_17", "Delta-t Between Event 61 and SiPM 2 Peaks", 17000, 0, ev_61_dt_max);                                              // 1 bin = 1,000 ns

    TH1D* h_dt_61_18 = new TH1D("h_dt_61_18", "Delta-t Between Event 61 and SiPM 3 Peaks", 17000, 0, ev_61_dt_max);                                              // 1 bin = 1,000 ns

    TH1D* h_dt_61_19 = new TH1D("h_dt_61_19", "Delta-t Between Event 61 and SiPM 4 Peaks", 17000, 0, ev_61_dt_max);                                              // 1 bin = 1,000 ns

    TH1D* h_dt_61_20 = new TH1D("h_dt_61_20", "Delta-t Between Event 61 and SiPM 5 Peaks", 17000, 0, ev_61_dt_max);                                              // 1 bin = 1,000 ns

    TH1D* h_dt_61_21 = new TH1D("h_dt_61_21", "Delta-t Between Event 61 and SiPM 6 Peaks", 17000, 0, ev_61_dt_max);                                              // 1 bin = 1,000 ns

    TH1D* h_dt_61_22 = new TH1D("h_dt_61_22", "Delta-t Between Event 61 and SiPM 7 Peaks", 17000, 0, ev_61_dt_max);                                              // 1 bin = 1,000 ns

    TH1D* h_dt_61_23 = new TH1D("h_dt_61_23", "Delta-t Between Event 61 and SiPM 8 Peaks", 17000, 0, ev_61_dt_max);                                              // 1 bin = 1,000 ns

    TH1D* h_dt_61_top = new TH1D("h_dt_61_top", "Delta-t Between Event 61 and SiPM 9 & 10 Peaks", 17000, 0, ev_61_dt_max);                                       // 1 bin = 1,000 ns

    TH1D* h_vpint_16 = new TH1D("h_vpint_16", "Low dt SiPM 1 Integral Values", 150, -400, 2000);

    TH1D* h_vpint_17 = new TH1D("h_vpint_17", "Low dt SiPM 2 Integral Values", 150, -400, 2000);

    TH1D* h_vpint_18 = new TH1D("h_vpint_18", "Low dt SiPM 3 Integral Values", 150, -400, 2000);

    TH1D* h_vpint_19 = new TH1D("h_vpint_19", "Low dt SiPM 4 Integral Values", 150, -400, 2000);

    TH1D* h_vpint_20 = new TH1D("h_vpint_20", "Low dt SiPM 5 Integral Values", 150, -400, 2000);

    TH1D* h_vpint_21 = new TH1D("h_vpint_21", "Low dt SiPM 6 Integral Values", 150, -400, 2000);

    TH1D* h_vpint_22 = new TH1D("h_vpint_22", "Low dt SiPM 7 Integral Values", 150, -400, 2000);

    TH1D* h_vpint_23 = new TH1D("h_vpint_23", "Low dt SiPM 8 Integral Values", 150, -400, 2000);

    TH1D* h_vpint_top = new TH1D("h_vpint_top", "Low dt SiPM 9 & 10 Integral Values", 150, -400, 2000);

    TH2D* h_detint_vs_vpint_16 = new TH2D("h_detint_vs_vpint_16", "SiPM 1 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_17 = new TH2D("h_detint_vs_vpint_17", "SiPM 2 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_18 = new TH2D("h_detint_vs_vpint_18", "SiPM 3 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_19 = new TH2D("h_detint_vs_vpint_19", "SiPM 4 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_20 = new TH2D("h_detint_vs_vpint_20", "SiPM 5 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_21 = new TH2D("h_detint_vs_vpint_21", "SiPM 6 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_22 = new TH2D("h_detint_vs_vpint_22", "SiPM 7 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_23 = new TH2D("h_detint_vs_vpint_23", "SiPM 8 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_top = new TH2D("h_detint_vs_vpint_top", "SiPM 9 & 10 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH1D* h_intmichels_PMT0 = new TH1D("h_intmichels_PMT0", "PMT0 Michel Electron Integral Values", 100, 0, 100);

    TH1D* h_intmichels_PMT1 = new TH1D("h_intmichels_PMT1", "PMT1 Michel Electron Integral Values", 100, 0, 100);

    TH1D* h_intmichels_PMT2 = new TH1D("h_intmichels_PMT2", "PMT2 Michel Electron Integral Values", 100, 0, 100);

    TH1D* h_intmichels_PMT3 = new TH1D("h_intmichels_PMT3", "PMT3 Michel Electron Integral Values", 100, 0, 100);

    TH1D* h_intmichels_PMT4 = new TH1D("h_intmichels_PMT4", "PMT4 Michel Electron Integral Values", 100, 0, 100);

    TH1D* h_intmichels_PMT5 = new TH1D("h_intmichels_PMT5", "PMT5 Michel Electron Integral Values", 100, 0, 100);

    TH1D* h_intmichels_PMT6 = new TH1D("h_intmichels_PMT6", "PMT6 Michel Electron Integral Values", 100, 0, 100);

    TH1D* h_intmichels_PMT7 = new TH1D("h_intmichels_PMT7", "PMT7 Michel Electron Integral Values", 100, 0, 100);

    TH1D* h_intmichels_PMT8 = new TH1D("h_intmichels_PMT8", "PMT8 Michel Electron Integral Values", 100, 0, 100);

    TH1D* h_intmichels_PMT9 = new TH1D("h_intmichels_PMT9", "PMT9 Michel Electron Integral Values", 100, 0, 100);

    TH1D* h_intmichels_PMT10 = new TH1D("h_intmichels_PMT10", "PMT10 Michel Electron Integral Values", 100, 0, 100);

    TH1D* h_intmichels_PMT11 = new TH1D("h_intmichels_PMT11", "PMT11 Michel Electron Integral Values", 100, 0, 100);

    TH1D* h_muon_pmt_sipm_dt = new TH1D("h_muon_pmt_sipm_dt", "Time Difference Between PMT and SiPM Peaks in Muon Events", 45, -360, 360);

    TH1D* h_only_one_sipm = new TH1D("h_only_one_sipm", "Integral Distribution of SiPM8 with No Events Anywhere Else", 200, 0, 2000);

    std::vector<std::vector<double>> amp_on_off_val;

    std::vector<double> amp_val;

    // Create a text string, which is used to output the text file
    
    std::vector<int> runlist;

    string intval1;

    // Read from the text file

    // ifstream ReadRunListFile("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/bef_aft_aug16.txt");
    // ifstream ReadRunListFile("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/event_rate_over_time_runlist.txt");
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/yearoneevents/yearoneevents_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/sum2024goldenrunlist/sum2024goldenrunlist_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/2024allgoldenruns/2024allgoldenruns_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/weekly2024goldenruns/weekly2024goldenruns_%i.txt", run));
    ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/allweeklyruns/allweeklyruns_%i.txt", run));           // 1-46: 2023-2024; 47-62: Fall 2024; 63-82: Spring 2025; 83-94: Summer 2025 (Beam Off)
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

    string int_val_1;

    // Read from the text file

    ifstream ReadHLMFile("Vertical_Muon_Peak_Center_1-124_Run_Vals.txt");          // Vertical_Muon_Peak_Center_1-124_Val_Only.txt         // hl_mean_val_all_1.txt          // hl_mean_val_all_temp.txt

    // Use a while loop together with the getline() function to read the file line by line

    while (getline(ReadHLMFile, int_val_1)) {

        if (stod(int_val_1) > 0.0) {

            HLmeancounter += 1;

            sumHLmean += stod(int_val_1);

        }
        
    }

    double defaultHLmean = 1813.94;           // = 1813.94; (12/1/24 Vert Muon Mean Value)         // = 34.5323; (11/26/24 HL Mean Value)         // sumHLmean / HLmeancounter;           // sumHLmean / HLmeancounter = 1783.874 (Averaged value)

    double currentHLmean = defaultHLmean;

    // cout << "\n" << "defaultHLmean = " << defaultHLmean << endl;

    // Close the file

    ReadHLMFile.close();

    ReadHLMFile.clear();

    // Create and open a text file
    ofstream MIFile(Form("Michel_Integrals_Run%i.txt", run));
    ofstream MERFile(Form("Michel_Event_Rate_Run%i.txt", run));
    ofstream MBVFile(Form("Michel_Bin_Vals_Run%i.txt", run));
    ofstream uBVFile(Form("Muon_Bin_Vals_Run%i.txt", run));
    ofstream HLBVFile(Form("High_Light_Bin_Vals_Run%i.txt", run));
    ofstream LLBVFile(Form("Low_Light_Bin_Vals_Run%i.txt", run));
    ofstream MBBVFile(Form("Min_Bias_Bin_Vals_Run%i.txt", run));
    ofstream MIAFile(Form("Michel_Integral_Amplitude_Run%i.txt", run));
    ofstream MIPFile(Form("Michel_Integral_Peak_Center_Run%i.txt", run));

    // for(int iRun = run; iRun <= last_run; iRun++) {
    
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

        bool old_channel_map = false; int Event61Chan = 15; int SiPM1Chan = 16; int SiPM2Chan = 17; int SiPM3Chan = 18; int SiPM4Chan = 19; int SiPM5Chan = 20; int SiPM6Chan = 21; int SiPM7Chan = 22; int SiPM8Chan = 23; int SiPM9Chan = 24; int SiPM10Chan = 25;

        if (vers_num < 5) {old_channel_map = true; Event61Chan = 15; SiPM1Chan = 16; SiPM2Chan = 17; SiPM3Chan = 18; SiPM4Chan = 19; SiPM5Chan = 20; SiPM6Chan = 21; SiPM7Chan = 22; SiPM8Chan = 23; SiPM9Chan = 24; SiPM10Chan = 25;}

        else if (vers_num == 5) {old_channel_map = false; SiPM1Chan = 12; SiPM2Chan = 13; SiPM3Chan = 14; SiPM4Chan = 15; SiPM5Chan = 16; SiPM6Chan = 17; SiPM7Chan = 18; SiPM8Chan = 19; SiPM9Chan = 20; SiPM10Chan = 21; Event61Chan = 22;}

        TFile *f;
        //your root file location here
        //  if(gSystem->AccessPathName(Form("/data%i/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v%i.root", data_num, run_iterable, vers_num))){
        if (gSystem->AccessPathName(Form("/data%i/coherent/data/d2o/processedData/run%i_processed_v%i.root", data_num, run_iterable, vers_num))) {
        //  if(gSystem->AccessPathName(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/run%i_processed_v%i.root", run_iterable, vers_num))){
            cout << "Could not open file " << run_iterable << endl;
            continue;           // return -1;
        } else{
        //     f = new TFile(Form("/data%i/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v%i.root", data_num, run_iterable, vers_num));
            f = new TFile(Form("/data%i/coherent/data/d2o/processedData/run%i_processed_v%i.root", data_num, run_iterable, vers_num));
        //     f = new TFile(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/run%i_processed_v%i.root", run_iterable, vers_num));
        }

        auto tsstart = (TParameter<Long64_t> *) f->Get("starttime");
        
        Long64_t run_starttime = tsstart->GetVal();

        // Create a text string, which is used to output the text file

        std::vector<double> num_vec;

        num_vec.clear();

        string intval2;

        string num_str;

        intval2.clear();

        num_str.clear();

        // Read from the text file

        ifstream ReadHLMFile("Vertical_Muon_Peak_Center_1-124_Run_Nums.txt");          // Vertical_Muon_Peak_Center_1-124_No_HL_Mean_Correction.txt          // hl_mean_val_look_up_1.txt          // hl_mean_val_look_up_temp.txt

        // Use a while loop together with the getline() function to read the file line by line

        while (getline(ReadHLMFile, intval2)) {

            num_vec.clear();

            num_str.clear();

            intval2.push_back('\t');

            for (char ch : intval2) {
                
                if (ch != '\t') {num_str.push_back(ch);}
                
                else if (ch == '\t') {num_vec.push_back(stod(num_str)); num_str.clear();}
            
            }

            // cout << "\n" << "num_vec[0] = " << num_vec[0] << " num_vec[1] =  " << num_vec[1] << " num_vec.size() =  " << num_vec.size() << endl;

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
	      outputFile.Form("%s.root", Form("MichelEvents_Run%i", run_iterable));
	    } else {
	      outputFile.Form("%s/%s.root", outputFilePath, Form("MichelEvents_Run%i", run_iterable));
	    }
	    TFile *fileOut = new TFile(outputFile, "RECREATE");
	    TTree *michelTree = new TTree("michelTree", "michelTree");
	    michelTree->SetDirectory(fileOut);

	    // Data to keep track of for each michelTree entry
	    int *br_entry;  /* Entry # from inputFile TTree T */
	    double *br_e1;  /* Energy of first pulse (photo-electrons) */
	    double *br_p1;  /* Peak (max amplitude) of first pulse (photo-electrons) */
	    double *br_t1;  /* Universal start time (10% peak) of first pulse (nanoseconds) */
	    double *br_d1;  /* Duration of first pulse (nanoseconds) */
	    double *br_e2;  /* Energy of second pulse (photo-electrons) */
	    double *br_p2;  /* Peak (max amplitude) of second pulse (photo-electrons) */
	    double *br_t2;  /* Universal start time (10% peak) of second pulse (nanoseconds) */
	    double *br_d2;  /* Duration of second pulse (nanoseconds) */
	    double *br_dt;  /* Time separation between pulse onsets (nanoseconds) */
	    bool *br_issue; /* Flag to keep track of unusual michelTree entries */
        double* br_n1;  /* Number of channels in which we see first pulse (photo-electrons) */
        bool* br_s1;    /* Is first pulse (photo-electrons) timing consistent across all channels */
        bool* br_b1;    /* Tracks whether beam is on or off for first pulse */
        double* br_tr1; /* Tracks whether trigger is external (2) or internal (16) for first pulse */
        double* br_l1;  /* Length of first waveform in number of bins */
        double* br_n2;  /* Number of channels in which we see second pulse (photo-electrons) */
        bool* br_s2;    /* Is second pulse (photo-electrons) timing consistent across all channels */
        bool* br_b2;    /* Tracks whether beam is on or off for second pulse */
        double* br_tr2; /* Tracks whether trigger is external (2) or internal (16) for second pulse */
        double* br_l2;  /* Length of second waveform in number of bins */
        // double* br_c1;  /* PMT integral values of all cosmic events */

	    michelTree->Branch("entry", &br_entry, "entry/I");
	    michelTree->Branch("e1", &br_e1, "e1/d");
	    michelTree->Branch("p1", &br_p1, "p1/d");
	    michelTree->Branch("t1", &br_t1, "t1/d");
	    michelTree->Branch("d1", &br_d1, "d1/d");
	    michelTree->Branch("e2", &br_e2, "e2/d");
	    michelTree->Branch("p2", &br_p2, "p2/d");
	    michelTree->Branch("t2", &br_t2, "t2/d");
	    michelTree->Branch("d2", &br_d2, "d2/d");
	    michelTree->Branch("dt", &br_dt, "dt/d");
	    michelTree->Branch("issue", &br_issue, "issue/O");
        michelTree->Branch("n1", &br_n1, "n1/d");
        michelTree->Branch("s1", &br_s1, "s1/O");
        michelTree->Branch("b1", &br_b1, "b1/O");
        michelTree->Branch("tr1", &br_tr1, "tr1/d");
        michelTree->Branch("l1", &br_l1, "l1/d");
        michelTree->Branch("n2", &br_n2, "n2/d");
        michelTree->Branch("s2", &br_s2, "s2/O");
        michelTree->Branch("b2", &br_b2, "b2/O");
        michelTree->Branch("tr2", &br_tr2, "tr2/d");
        michelTree->Branch("l2", &br_l2, "l2/d");
        // michelTree->Branch("c1", &br_c1, "c1/d");

        // Setup a vector to contain the waveform info
        std::map<int, int> numPulses;

        // Create vector to hold info of peak times of pulses
        std::vector<struct pulse> mult_pulses;
        // std::vector<vector<struct pulse>> pulse_groups;
        double time_count = 0.0;
        bool WithinPulse = false;
        double first_pulse_start_time;

        double time_count_61 = 0.0;
        bool WithinPulse_61 = false;
        double first_pulse_start_time_61;

        double time_count_61_muons = 0.0;
        bool WithinPulse_61_muons = false;
        double first_pulse_start_time_61_muons;

        double time_count_61_rev = 0.0;
        bool WithinPulse_61_rev = false;
        double first_pulse_start_time_61_rev;
        double first_pulse_int_61_rev;

        double time_count_61_muons_rev = 0.0;
        bool WithinPulse_61_muons_rev = false;
        double first_pulse_start_time_61_muons_rev;

        double time_count_61_vc = 0.0;
        bool WithinPulse_61_vc = false;
        double first_pulse_start_time_61_vc;

        double time_count_61_vc_muons = 0.0;
        bool WithinPulse_61_vc_muons = false;
        double first_pulse_start_time_61_vc_muons;

        double time_count_61_vc_rev = 0.0;
        bool WithinPulse_61_vc_rev = false;
        double first_pulse_start_time_61_vc_rev;
        double first_pulse_int_61_vc_rev;

        double time_count_61_vc_muons_rev = 0.0;
        bool WithinPulse_61_vc_muons_rev = false;
        double first_pulse_start_time_61_vc_muons_rev;

        double time_count_61_smc = 0.0;
        bool WithinPulse_61_smc = false;
        double first_pulse_start_time_61_smc;

        double time_count_61_smc_muons = 0.0;
        bool WithinPulse_61_smc_muons = false;
        double first_pulse_start_time_61_smc_muons;

        double time_count_61_smc_rev = 0.0;
        bool WithinPulse_61_smc_rev = false;
        double first_pulse_start_time_61_smc_rev;
        double first_pulse_int_61_smc_rev;

        double time_count_61_smc_muons_rev = 0.0;
        bool WithinPulse_61_smc_muons_rev = false;
        double first_pulse_start_time_61_smc_muons_rev;

        double time_count_61_spmc = 0.0;
        bool WithinPulse_61_spmc = false;
        double first_pulse_start_time_61_spmc;

        double time_count_61_spmc_muons = 0.0;
        bool WithinPulse_61_spmc_muons = false;
        double first_pulse_start_time_61_spmc_muons;

        double time_count_61_spmc_rev = 0.0;
        bool WithinPulse_61_spmc_rev = false;
        double first_pulse_start_time_61_spmc_rev;
        double first_pulse_int_61_spmc_rev;

        double time_count_61_spmc_muons_rev = 0.0;
        bool WithinPulse_61_spmc_muons_rev = false;
        double first_pulse_start_time_61_spmc_muons_rev;

        double time_count_61_pmc = 0.0;
        bool WithinPulse_61_pmc = false;
        double first_pulse_start_time_61_pmc;

        double time_count_61_pmc_muons = 0.0;
        bool WithinPulse_61_pmc_muons = false;
        double first_pulse_start_time_61_pmc_muons;

        double time_count_61_pmc_rev = 0.0;
        bool WithinPulse_61_pmc_rev = false;
        double first_pulse_start_time_61_pmc_rev;
        double first_pulse_int_61_pmc_rev;

        double time_count_61_pmc_muons_rev = 0.0;
        bool WithinPulse_61_pmc_muons_rev = false;
        double first_pulse_start_time_61_pmc_muons_rev;

        bool MichelVetoActivity = true;

        // Get statistics for up to 1 million entries from fileIn TTree T
        int numEntries = std::min((int)t->GetEntries(), MAX_NUM_ENTRIES);

        std::vector<double> chan_lengths;

        std::vector<double> peak_pos_RMS;

        std::vector<double> chan_start_no_outliers;

        std::vector<int> u_event_bin_val;

        long long last_event_time = 0.0;

        long long last_int_trig_event_time = 0.0;

        double avg_peak_pos_RMS = 0.0;

        double var_peak_pos_RMS = 0.0;

        int ch15_on = 0;

        double ev61_time = 0.0;
        double minLED_time = 0.0;
        double lowLED_time = 0.0;
        double highLED_time = 0.0;
        double exttrig_time = 0.0;

        int counter = 0;

        long long prev_event_time = 0.0;

        int prev_event_tB = 0.0;
    
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

            all_nsTime.push_back(nsTime);

            // if (triggerBits >= 32) {triggerBits = triggerBits - 32;}         // && (triggerBits == 3 || triggerBits == 33 || triggerBits == 35)

            // if (counter < 1000) {counter += 1; cout << "\n" << "Run, Event Number, triggerBits =  " << run_iterable << " " << iEnt << " " << triggerBits << endl;}

            h_tB->Fill(triggerBits);

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

            double p_int_16 = 0;
            double p_int_17 = 0;
            double p_int_18 = 0;
            double p_int_19 = 0;
            double p_int_20 = 0;
            double p_int_21 = 0;
            double p_int_22 = 0;
            double p_int_23 = 0;
            double p_int_24 = 0;
            double p_int_25 = 0;

            double p_start_16 = 0;
            double p_start_17 = 0;
            double p_start_18 = 0;
            double p_start_19 = 0;
            double p_start_20 = 0;
            double p_start_21 = 0;
            double p_start_22 = 0;
            double p_start_23 = 0;
            double p_start_24 = 0;
            double p_start_25 = 0;

            chan_lengths.clear();

            for (int iChan = 0; iChan < 26; iChan++) {

                if (old_channel_map) {if (iChan == 12 || iChan == 13 || iChan == 14) {all_chan_start.push_back(1000); continue;}}

                else if (!old_channel_map) {if (iChan == 23 || iChan == 24 || iChan == 25) {all_chan_start.push_back(1000); continue;}}

                for (int i = 0; i < ADCSIZE; i++) {h_wf->SetBinContent(i + 1, adcVal[iChan][i] - baselineMean[iChan]); event_bin_val.push_back(adcVal[iChan][i] - baselineMean[iChan]);}

                if (iChan <= 11) {sphe_int = run_starttime * y_val[iChan] + b_val[iChan];}          // run_starttime * y_val[iChan] + b_val[iChan];}          // integralToPETemp[iChan];}

                if (iChan <= 11) {h_blM_pmt->Fill(baselineMean[iChan]);}

                else if (iChan >= SiPM1Chan && iChan <= SiPM10Chan) {h_blM_vp->Fill(baselineMean[iChan]);}

                else if (iChan == Event61Chan) {h_blM_ev61->Fill(baselineMean[iChan]);}

                // h_wf->GetXaxis()->SetTitle("Time (1 bin = 16 ns)");
                // h_wf->GetYaxis()->SetTitle("Counts");
                // h_wf->Draw();
                // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/wf_image/wf_run_%i_event_%i_ch_%i.png", run_iterable, iEnt + 1, iChan));
                // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/wf_image/wf_run_%i_event_%i_ch_%i.png", run_iterable, iEnt + 1, iChan));

                if (iChan == Event61Chan) {

                    Ev61Energy = 0.0;

                    for (int iBin = 1; iBin <= h_wf->GetNbinsX(); iBin++) {

                        double iBinContent = h_wf->GetBinContent(iBin);

                        Ev61Energy += iBinContent;

                    }

                    if (Ev61Energy > EV61_THRESHOLD) {

                        all_chan_beam = true;

                        ch15_on += 1;

                    }

                    if (triggerBits == 1 || triggerBits == 3 || triggerBits == 33 || triggerBits == 35) {h_intev61->Fill(Ev61Energy);}
                    
                    /*
                    if (Ev61Energy > 1650 && Ev61Energy < 1750) {

                        h_wf->GetXaxis()->SetTitle("Time (1 bin = 16 ns)");
                        h_wf->GetYaxis()->SetTitle("Counts");
                        h_wf->Draw();
                        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/wf_image/wf_run_%i_event_%i_ch_%i.png", run_iterable, iEnt + 1, iChan));
                        gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/wf_image/wf_run_%i_event_%i_ch_%i.png", run_iterable, iEnt + 1, iChan));

                    }
                    */

                    //if (all_chan_beam == true) {           // avg_pulse.beam == true && (avg_pulse.trigger >= 32 && avg_pulse.trigger <= 35) && avg_pulse.number >= 10

                        //h_intneu->Fill(std::accumulate(all_chan_peak.begin(), all_chan_peak.end(), 0.0));

                    //}

                    // continue;

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
                                p.energy = pulseEnergy / sphe_int;

                                all_chan_start_pmt.push_back(p.start);      // This is triggering for multiple bins above threshold? Even within single peak?
                                all_chan_end.push_back(p.end);
                                peak_energy.push_back(p.energy);
                                // all_chan_peak.push_back(p.peak);

                                // cout << "\n" << "Channel " << iChan << " pulse start time: " << p.start << endl;

                                // cout << "\n" << "Channel " << iChan << " pulse peak bin: " << peakBin << endl;

                                // cout << "\n" << "Channel " << iChan << " pulse height: " << peak << " = " << p.peak << " Ph.e." << endl;

                                h_pst->Fill(p.start);                       // Secondary pulses in each channel are afterpulsing? Ignore? Check if secondary pulse appears in multiple channels?

                                h_pst_v_pint->Fill(p.start, p.energy);

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

                    // if (iChan == 10 && triggerBits == 8) {

                        // h_inthighLED->Fill(AllPulseEnergy / (sphe_int * (currentHLmean / defaultHLmean)));

                        // inthighLED.push_back(AllPulseEnergy / (sphe_int * (currentHLmean / defaultHLmean)));

                    // }

                    if (iChan == 0 && triggerBits == 4) {

                        h_spheLL->Fill(AllPulseEnergy);

                    }

                    if (iChan == 0 && triggerBits == 16) {

                        h_spheMB->Fill(AllPulseEnergy);

                    }

                    all_chan_energy.push_back(AllPulseEnergy / sphe_int);

                    if (AllPulseEnergy / sphe_int > 1) {

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

            if (triggerBits != 0 && triggerBits != 1 && triggerBits != 4 && triggerBits != 8 && triggerBits != 16) {h_mfpst->Fill(pulse_start_time);}

            if (num_chan_over_1phe >= 10) {

                int pulse_start_index = mostFrequentIndex(all_chan_start_pmt);

                if (pulse_start_index >= 0) {h_mfpst_v_mfpint->Fill(pulse_start_time, peak_energy[pulse_start_index]);}

                else if (pulse_start_index == -1) {h_mfpst_v_mfpint->Fill(pulse_start_time, 0);}

                else if (pulse_start_index == -2) {h_mfpst_v_mfpint->Fill(pulse_start_time, (std::accumulate(peak_energy.begin(), peak_energy.end(), 0.0) / peak_energy.size()));}

            }

            // if (pulse_start_time == 0) {continue;}

            struct pulse avg_pulse;

            avg_pulse.start = nsTime + pulse_start_time;                                                        // Need a way to weed out secondary pulses - use weighted average?
            avg_pulse.end = nsTime + pulse_end_time;
            avg_pulse.energy = std::accumulate(all_chan_energy.begin(), all_chan_energy.end(), 0.0) / (currentHLmean / defaultHLmean);           // Integrals still sum over secondary peaks discarded in .start and .end values 
            avg_pulse.peak = std::accumulate(all_chan_peak.begin(), all_chan_peak.end(), 0.0);
            avg_pulse.number = num_chan_over_1phe;
            avg_pulse.beam = all_chan_beam;
            avg_pulse.trigger = triggerBits;
            avg_pulse.length = getAverage(chan_lengths);
            avg_pulse.side_vp_energy = std::accumulate(side_veto_panel_energy.begin(), side_veto_panel_energy.end(), 0.0);
            avg_pulse.top_vp_energy = std::accumulate(top_veto_panel_energy.begin(), top_veto_panel_energy.end(), 0.0);
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

            if (avg_pulse.trigger == 34) {h_34detint->Fill(avg_pulse.energy); h_34vpint->Fill(avg_pulse.side_vp_energy + avg_pulse.top_vp_energy);}

            // if (avg_pulse.beam == true) {cout << "\n" << "Run Number = " << run_iterable << ", Event Number = " << iEnt + 1 << ", triggerBits = " << avg_pulse.trigger << endl;}

            // cout << "\n" << "Run Number = " << run_iterable << ", Event Number = " << iEnt + 1 << ", triggerBits = " << avg_pulse.trigger << endl;

            // if (avg_pulse.beam == true && avg_pulse.number >= 10) {h_intneu->Fill(avg_pulse.energy);}           // && (avg_pulse.trigger >= 32 && avg_pulse.trigger <= 35)

            if ((avg_pulse.trigger >= 32 && avg_pulse.trigger <= 35) && avg_pulse.number >= 10) {amp_val.push_back(avg_pulse.peak);}

            if (avg_pulse.trigger >= 32 && avg_pulse.trigger <= 35) {

                h_int->Fill(avg_pulse.energy);

                h_amp->Fill(avg_pulse.peak);

            }

            if (all_chan_peakbin.size() > 0) {

                peak_pos_RMS.push_back(rmsValue(all_chan_peakbin));

            }

            for (int iPeak = 0; iPeak < all_chan_start.size(); iPeak++) {

                if (all_chan_start[iPeak] < (pulse_start_time + 10 * 16) && all_chan_start[iPeak] > (pulse_start_time - 10 * 16)) {

                    chan_start_no_outliers.push_back(all_chan_start[iPeak]);

                }

                else {

                    // Check if pulse outside this range occurs more than once

                    for (int jPeak = 0; jPeak < all_chan_start.size(); jPeak++) {

                        if (iPeak == jPeak) {

                            continue;

                        }

                        else if (all_chan_start[iPeak] < (all_chan_start[jPeak] + 1 * 16) && all_chan_start[iPeak] > (all_chan_start[jPeak] - 1 * 16)) {

                            chan_start_no_outliers.push_back(all_chan_start[iPeak]);

                        }

                    }

                }
            }

            double var_val = variance(chan_start_no_outliers);

            h_var->Fill(var_val);

            if (var_val < 5 * 16) {            // Pulses largely not consistent, need to weed out outlier and secondary pulses

                avg_pulse.single = true;

            }

            else {

                avg_pulse.single = false;

            }

            // cout << "\n" << "Pulse consistent across all channels? " << avg_pulse.single << endl;

            if (avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {

                h_ml->Fill(num_chan_over_1phe);

                h_ml_v_int->Fill(num_chan_over_1phe, avg_pulse.energy);

            }

            if (avg_pulse.beam) {

                // cout << "\n" << "Event 61 event" << endl;

                if (ev61_time == 0.0) {

                    ev61_time = nsTime + chan_15_start;

                }

                else {

                    // if ((nsTime + chan_15_start) - ev61_time < 16000000) {cout << "\n" << "Low Event 61 dt value: Run & Event Number = " << run_iterable << " " << iEnt + 1 << endl;}
                
                    h_time_bt_ev61->Fill((nsTime + chan_15_start) - ev61_time);

                    ev61_time = nsTime + chan_15_start;

                }

            }

            h_intall->Fill(avg_pulse.energy);

            intall.push_back(avg_pulse.energy);

            if (avg_pulse.trigger == 16) {

                // cout << "\n" << "Min bias LED event (tB == 4)" << endl;

                h_intlowLED->Fill(avg_pulse.energy);

                intlowLED.push_back(avg_pulse.energy);

                if (lowLED_time == 0.0) {

                    lowLED_time = avg_pulse.start;
                    
                }

                else {

                    h_time_bt_lowLED->Fill(avg_pulse.start - lowLED_time);

                    lowLED_time = avg_pulse.start;

                }

                if (avg_pulse.energy > 450 && avg_pulse.energy < 575) {h_MB_HL_dt->Fill(lowLED_time - exttrig_time);}

                // if (avg_pulse.energy > 475 && avg_pulse.energy < 525 && counter < 1000) {counter += 1; cout << "\n" << "LOW LIGHT PLOT: Run, Event Number, Integral, triggerBits, Previous triggerBits, dt =  " << run_iterable << " " << iEnt << " " << avg_pulse.energy << " " << triggerBits << " " << prev_event_tB << " " << avg_pulse.start - prev_event_time << endl;}

                // if (avg_pulse.energy > -100 && avg_pulse.energy < -50 && counter < 1000) {counter += 1; cout << "\n" << "LOW LIGHT PLOT: Run, Event Number, Integral, triggerBits, Previous triggerBits, dt =  " << run_iterable << " " << iEnt << " " << avg_pulse.energy << " " << triggerBits << " " << prev_event_tB << " " << avg_pulse.start - prev_event_time << endl;}

                // if (counter < 15000) {counter += 1; for (int iVec = 0; iVec < event_bin_val.size(); iVec++) {LLBVFile << event_bin_val[iVec] << "\t";} LLBVFile << "\n";}

                exttrig_time = avg_pulse.start;

            }

            else if (avg_pulse.trigger == 8) {

                // cout << "\n" << "High light LED event (tB == 8)" << endl;

                // run_iterable >= 16000 && run_iterable <= 16200           // run_iterable >= 16300 && run_iterable <= 16500

                h_inthighLED->Fill(avg_pulse.energy);

                inthighLED.push_back(avg_pulse.energy);

                if (highLED_time == 0.0) {

                    highLED_time = avg_pulse.start;

                }

                else {

                    h_time_bt_highLED->Fill(avg_pulse.start - highLED_time);

                    if (avg_pulse.start - highLED_time < 18000000) {HL_dt_17ms_counter += 1;}

                    highLED_time = avg_pulse.start;

                }

                // h_MB_HL_dt->Fill(highLED_time - minLED_time);

                // if (avg_pulse.energy < 300 && counter < 1000) {counter += 1; cout << "\n" << "HIGH LIGHT PLOT: Run, Event Number, Integral, triggerBits, Previous triggerBits, dt =  " << run_iterable << " " << iEnt << " " << avg_pulse.energy << " " << triggerBits << " " << prev_event_tB << " " << avg_pulse.start - prev_event_time << endl;}

                // if (counter < 15000) {counter += 1; for (int iVec = 0; iVec < event_bin_val.size(); iVec++) {HLBVFile << event_bin_val[iVec] << "\t";} HLBVFile << "\n";}

                exttrig_time = avg_pulse.start;

            }

            else if (avg_pulse.trigger == 4) {

                // cout << "\n" << "Low light LED event (tB == 16)" << endl;

                h_intminLED->Fill(avg_pulse.energy);

                intminLED.push_back(avg_pulse.energy);

                if (minLED_time == 0.0) {

                    minLED_time = avg_pulse.start;

                }

                else {

                    h_time_bt_minLED->Fill(avg_pulse.start - minLED_time);

                    minLED_time = avg_pulse.start;

                }

                // if (counter < 15000) {counter += 1; for (int iVec = 0; iVec < event_bin_val.size(); iVec++) {MBBVFile << event_bin_val[iVec] << "\t";} MBBVFile << "\n";}

                exttrig_time = avg_pulse.start;

            }

            else if (avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {

                h_intallint->Fill(avg_pulse.energy);

                intallint.push_back(avg_pulse.energy);

                if (avg_pulse.number >= 10) {

                    // h_intallmultcut->Fill(avg_pulse.energy);

                    intallmultcut.push_back(avg_pulse.energy);

                    // if (p_int_16 <= 750 && p_int_17 <= 950 && p_int_18 <= 1200 && p_int_19 <= 1375 && p_int_20 <= 525 && p_int_21 <= 700 && p_int_22 <= 700 && p_int_23 <= 500 && avg_pulse.top_vp_energy <= 450) {

                    if (p_int_16 <= min_muon_vals[0] && p_int_17 <= min_muon_vals[1] && p_int_18 <= min_muon_vals[2] && p_int_19 <= min_muon_vals[3] && p_int_20 <= min_muon_vals[4] && p_int_21 <= min_muon_vals[5] && p_int_22 <= min_muon_vals[6] && p_int_23 <= min_muon_vals[7] && avg_pulse.vp_int_top <= min_muon_vals[8]) {

                        // h_intallmultmuoncut->Fill(avg_pulse.energy);

                        intallmultmuoncut.push_back(avg_pulse.energy);

                    }

                }

            }

            if (avg_pulse.trigger == 4 || avg_pulse.trigger == 8 || avg_pulse.trigger == 16 || avg_pulse.trigger == 1 || avg_pulse.trigger == 3 || avg_pulse.trigger == 33 || avg_pulse.trigger == 35) {

                prev_event_tB = avg_pulse.trigger;

                prev_event_time = avg_pulse.start;

            }

            // Events in SiPM 6 or 8 but in nothing else

            if (avg_pulse.energy < 50 && p_int_16 < 50 && p_int_17 < 50 && p_int_18 < 50 && p_int_19 < 50 && p_int_20 < 50 && p_int_21 < 50 && p_int_22 < 50 && p_int_23 > 50 && avg_pulse.vp_int_top < 50) {h_only_one_sipm ->Fill(p_int_23);}

            // Search for cosmic events

            // Use these cuts:

            // (p_int_16 < 750 && p_int_17 < 950 && p_int_18 < 1200 && p_int_19 < 1375 && p_int_20 < 525 && p_int_21 < 700 && p_int_22 < 700 && p_int_23 < 500) && avg_pulse.vp_int_top > 450

            if (avg_pulse.number >= 10 && avg_pulse.energy >= 50 && avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16 && avg_pulse.start - last_event_time > echo_cut) {

                // if (p_int_16 >= 100 || p_int_17 >= 100 || p_int_18 >= 100 || p_int_19 >= 100 || p_int_20 >= 100 || p_int_21 >= 100 || p_int_22 >= 100 || p_int_23 >= 100 || avg_pulse.top_vp_energy >= 100) {

                // if ((p_int_16 < 750 && p_int_17 < 950 && p_int_18 < 1200 && p_int_19 < 1375 && p_int_20 < 525 && p_int_21 < 700 && p_int_22 < 700 && p_int_23 < 500) && avg_pulse.vp_int_top > 450) {

                // if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.vp_int_top > 450) {

                if (p_int_16 > min_muon_vals[0] || p_int_17 > min_muon_vals[1] || p_int_18 > min_muon_vals[2] || p_int_19 > min_muon_vals[3] || p_int_20 > min_muon_vals[4] || p_int_21 > min_muon_vals[5] || p_int_22 > min_muon_vals[6] || p_int_23 > min_muon_vals[7] || avg_pulse.vp_int_top > min_muon_vals[8]) {

                    if (avg_pulse.start - avg_pulse.last_muon_time > muon_dt_max && avg_pulse.start - avg_pulse.last_lead_time > muon_dt_max && avg_pulse.start - avg_pulse.last_det_time > muon_dt_max) {

                        cosmic_events += 1;

                        h_intcosmics->Fill(avg_pulse.energy);

                        // double c1 = avg_pulse.energy;

                        // michelTree->SetBranchAddress("c1", &c1);

                        // michelTree->Fill();

                        h_muon_pmt_sipm_dt->Fill(pulse_start_time - (vp_start_time - nsTime));

                        // cout << "\n" << "Run Number, Event Number, triggerBits, pulse_start_time, vp_start_time = " << run_iterable << " " << iEnt << " " << avg_pulse.trigger << " " << pulse_start_time << " " << (vp_start_time - nsTime) << endl;

                    }
                    
                }
            
            }

            // Is this event a muon stopping in lead? (Looking at events with peak at the very end of the waveform)
                
            // if (pulse_at_end && avg_pulse.energy < (20 / 2) && (p_int_16 > 750 / 2 || p_int_17 > 950 / 2 || p_int_18 > 1200 / 2 || p_int_19 > 1375 / 2 || p_int_20 > 525 / 2 || p_int_21 > 700 / 2 || p_int_22 > 700 / 2 || p_int_23 > 500 / 2 || avg_pulse.vp_int_top > 450 / 2)) {

            if (pulse_at_end && avg_pulse.energy < (20 / 2) && (p_int_16 > min_muon_vals[0] / 2 || p_int_17 > min_muon_vals[1] / 2 || p_int_18 > min_muon_vals[2] / 2 || p_int_19 > min_muon_vals[3] / 2 || p_int_20 > min_muon_vals[4] / 2 || p_int_21 > min_muon_vals[5] / 2 || p_int_22 > min_muon_vals[6] / 2 || p_int_23 > min_muon_vals[7] / 2 || avg_pulse.vp_int_top > min_muon_vals[8] / 2)) {
                
                avg_pulse.last_lead_time = vp_start_time;

            }
            
            // Is this event a muon stopping in lead? (Looking at "normal" events)

            // else if (avg_pulse.energy < 20 && (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.vp_int_top > 450)) {

            else if (avg_pulse.energy < 20 && (p_int_16 > min_muon_vals[0] || p_int_17 > min_muon_vals[1] || p_int_18 > min_muon_vals[2] || p_int_19 > min_muon_vals[3] || p_int_20 > min_muon_vals[4] || p_int_21 > min_muon_vals[5] || p_int_22 > min_muon_vals[6] || p_int_23 > min_muon_vals[7] || avg_pulse.vp_int_top > min_muon_vals[8])) {
                
                avg_pulse.last_lead_time = vp_start_time;

            }

            // Is this event a muon? (Looking at events with peak at the very end of the waveform)

            // if (pulse_at_end && avg_pulse.energy > (20 / 2) && (p_int_16 > 750 / 2 || p_int_17 > 950 / 2 || p_int_18 > 1200 / 2 || p_int_19 > 1375 / 2 || p_int_20 > 525 / 2 || p_int_21 > 700 / 2 || p_int_22 > 700 / 2 || p_int_23 > 500 / 2 || avg_pulse.vp_int_top > 450 / 2)) {

            if (pulse_at_end && avg_pulse.energy > (20 / 2) && (p_int_16 > min_muon_vals[0] / 2 || p_int_17 > min_muon_vals[1] / 2 || p_int_18 > min_muon_vals[2] / 2 || p_int_19 > min_muon_vals[3] / 2 || p_int_20 > min_muon_vals[4] / 2 || p_int_21 > min_muon_vals[5] / 2 || p_int_22 > min_muon_vals[6] / 2 || p_int_23 > min_muon_vals[7] / 2 || avg_pulse.vp_int_top > min_muon_vals[8] / 2)) {

                avg_pulse.last_muon_time = vp_start_time;

            }

            // Is this event a muon? (Looking at "normal" events)

            // else if (avg_pulse.energy > 20 && (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.vp_int_top > 450)) {

            else if (avg_pulse.energy > 20 && (p_int_16 > min_muon_vals[0] || p_int_17 > min_muon_vals[1] || p_int_18 > min_muon_vals[2] || p_int_19 > min_muon_vals[3] || p_int_20 > min_muon_vals[4] || p_int_21 > min_muon_vals[5] || p_int_22 > min_muon_vals[6] || p_int_23 > min_muon_vals[7] || avg_pulse.vp_int_top > min_muon_vals[8])) {

                avg_pulse.last_muon_time = vp_start_time;

            }

            // Is there a triple coincidence between veto panels 6 and 8 and at least 10 PMTs?

            if (p_int_21 > min_muon_vals[5] && p_int_23 > min_muon_vals[7] && (p_int_16 < min_muon_vals[0] && p_int_17 < min_muon_vals[1] && p_int_18 < min_muon_vals[2] && p_int_19 < min_muon_vals[3] && p_int_20 < min_muon_vals[4] && p_int_22 < min_muon_vals[6] && avg_pulse.vp_int_top < min_muon_vals[8])) {

                if ((avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) && avg_pulse.number >= 10 && avg_pulse.energy >= 50 && avg_pulse.start - last_event_time > echo_cut) {
                
                    // cout << "\n" << "Run Number, Event Number, triggerBits = " << run_iterable << " " << iEnt << " " << avg_pulse.trigger << endl;

                }
            
            }

            /*
            Conditions I can use to cut down on which peaks I am looking at:
                1. First pulse threshold
                2. Delta-t magnitude
                3. Event trigger type
                4. Number of channels peak appeared in
                5. Variance of peak position (after secondary peaks discarded)
            */

            // Find average pulse start time by finding pulse times that are repeated at least X number of times

            if (time_count == 0.0 && (avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) && avg_pulse.number >= 10 && avg_pulse.energy >= 50 && avg_pulse.start - last_event_time > echo_cut) {         // * (currentHLmean / defaultHLmean)            // avg_pulse.peak > PULSE_PE_THRESHOLD && 

                // if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.vp_int_top > 450) {

                // if (p_int_16 > vp_cut_vals[0] || p_int_17 > vp_cut_vals[1] || p_int_18 > vp_cut_vals[2] || p_int_19 > vp_cut_vals[3] || p_int_20 > vp_cut_vals[4] || p_int_21 > vp_cut_vals[5] || p_int_22 > vp_cut_vals[6] || p_int_23 > vp_cut_vals[7] || avg_pulse.vp_int_top > vp_cut_vals[8]) {

                if (p_int_16 > min_muon_vals[0] || p_int_17 > min_muon_vals[1] || p_int_18 > min_muon_vals[2] || p_int_19 > min_muon_vals[3] || p_int_20 > min_muon_vals[4] || p_int_21 > min_muon_vals[5] || p_int_22 > min_muon_vals[6] || p_int_23 > min_muon_vals[7] || avg_pulse.vp_int_top > min_muon_vals[8]) {
                    
                    WithinPulse = true;

                    first_pulse_start_time = avg_pulse.start;

                    MichelVetoActivity = true;

                }

            }

            if (WithinPulse && (avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) && avg_pulse.number >= 10 && avg_pulse.energy >= 50 && avg_pulse.start - last_event_time > echo_cut) {           // * (currentHLmean / defaultHLmean)                // Conditions on if statement mean not all counted pulses get stored in mult_pulses!!!
                
                // if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.vp_int_top > 450) {

                // if (p_int_16 > vp_cut_vals[0] || p_int_17 > vp_cut_vals[1] || p_int_18 > vp_cut_vals[2] || p_int_19 > vp_cut_vals[3] || p_int_20 > vp_cut_vals[4] || p_int_21 > vp_cut_vals[5] || p_int_22 > vp_cut_vals[6] || p_int_23 > vp_cut_vals[7] || avg_pulse.vp_int_top > vp_cut_vals[8]) {

                if (p_int_16 > min_muon_vals[0] || p_int_17 > min_muon_vals[1] || p_int_18 > min_muon_vals[2] || p_int_19 > min_muon_vals[3] || p_int_20 > min_muon_vals[4] || p_int_21 > min_muon_vals[5] || p_int_22 > min_muon_vals[6] || p_int_23 > min_muon_vals[7] || avg_pulse.vp_int_top > min_muon_vals[8]) {

                    u_event_bin_val.clear();

                    u_event_bin_val = event_bin_val;
                
                    if (time_count == 0.0) {
                    
                        mult_pulses.push_back(avg_pulse);

                        time_count = avg_pulse.start - first_pulse_start_time;

                        if (time_count == 0.0) {time_count = 0.5;}
                
                    }

                    else {

                        mult_pulses.clear();

                        mult_pulses.push_back(avg_pulse);

                        first_pulse_start_time = avg_pulse.start;

                        time_count = 0.5;

                    }

                }

                // else if (p_int_16 <= 750 && p_int_17 <= 950 && p_int_18 <= 1200 && p_int_19 <= 1375 && p_int_20 <= 525 && p_int_21 <= 700 && p_int_22 <= 700 && p_int_23 <= 500 && avg_pulse.vp_int_top <= 450) {

                // else if (p_int_16 <= vp_cut_vals[0] && p_int_17 <= vp_cut_vals[1] && p_int_18 <= vp_cut_vals[2] && p_int_19 <= vp_cut_vals[3] && p_int_20 <= vp_cut_vals[4] && p_int_21 <= vp_cut_vals[5] && p_int_22 <= vp_cut_vals[6] && p_int_23 <= vp_cut_vals[7] && avg_pulse.vp_int_top <= vp_cut_vals[8]) {

                else if (p_int_16 <= min_muon_vals[0] && p_int_17 <= min_muon_vals[1] && p_int_18 <= min_muon_vals[2] && p_int_19 <= min_muon_vals[3] && p_int_20 <= min_muon_vals[4] && p_int_21 <= min_muon_vals[5] && p_int_22 <= min_muon_vals[6] && p_int_23 <= min_muon_vals[7] && avg_pulse.vp_int_top <= min_muon_vals[8]) {

                    mult_pulses.push_back(avg_pulse);

                    time_count = avg_pulse.start - first_pulse_start_time;

                }

                if (time_count == 0.0) {time_count = 0.5;}
                
                if (run == 7894 && time_count > 1 && p_int_16 <= 100 && p_int_17 <= 100 && p_int_18 <= 200 && p_int_19 <= 100 && p_int_20 <= 100 && p_int_21 <= 100 && p_int_22 <= 100 && p_int_23 <= 100 && avg_pulse.vp_int_top <= 150) {MichelVetoActivity = false;}
                
                else if (run == 11241 && time_count > 1 && p_int_16 <= 100 && p_int_17 <= 100 && p_int_18 <= 300 && p_int_19 <= 100 && p_int_20 <= 100 && p_int_21 <= 100 && p_int_22 <= 100 && p_int_23 <= 100 && avg_pulse.vp_int_top <= 200) {MichelVetoActivity = false;}
                
                else if (run == 12636 && time_count > 1 && p_int_16 <= 150 && p_int_17 <= 150 && p_int_18 <= 400 && p_int_19 <= 150 && p_int_20 <= 150 && p_int_21 <= 150 && p_int_22 <= 150 && p_int_23 <= 150 && avg_pulse.vp_int_top <= 250) {MichelVetoActivity = false;}

                // else if (run == 15312 && time_count > 1 && p_int_16 <= 300 && p_int_17 <= 300 && p_int_18 <= 600 && p_int_19 <= 300 && p_int_20 <= 300 && p_int_21 <= 300 && p_int_22 <= 300 && p_int_23 <= 300 && avg_pulse.vp_int_top <= 300) {MichelVetoActivity = false;}

                else {MichelVetoActivity = true;}

            }

            if (time_count > 2000 && time_count <= muon_dt_max) {          // && MichelVetoActivity == false

                // if (counter < 1000) {counter += 1; cout << "\n" << "Run Number = " << run_iterable << ", Event Number = " << iEnt << ", triggerBits = " << avg_pulse.trigger << ", dt = " << time_count << endl;}

                // Fill histograms here?

                // run_iterable >= 16000 && run_iterable <= 16200           // run_iterable >= 16300 && run_iterable <= 16500

                hist_events += 1;

                tot_michel_energy += mult_pulses.back().energy;

                h_dt->Fill(time_count);
                            
                h_intmuons->Fill(mult_pulses.front().energy);
                            
                h_intmichels->Fill(mult_pulses.back().energy);

                h_dt_v_intmichels->Fill(time_count, mult_pulses.back().energy);

                if (mult_pulses.back().energy > 200 && mult_pulses.back().energy < 400) {h_slice_vpintavg->Fill((mult_pulses.back().side_vp_energy + mult_pulses.back().top_vp_energy) / 9);}

                // Write to the file

                for (int iVec = 0; iVec < event_bin_val.size(); iVec++) {MBVFile << event_bin_val[iVec] << "\t";} MBVFile << avg_pulse.start - last_int_trig_event_time << "\t" << "\n";

                for (int iVec = 0; iVec < u_event_bin_val.size(); iVec++) {uBVFile << u_event_bin_val[iVec] << "\t";} uBVFile << avg_pulse.start - last_int_trig_event_time << "\t" << "\n";

                // if (mult_pulses.back().vp_int_18 > 350 && mult_pulses.back().vp_int_18 < 450) {cout << "\n" << "350 - 450 ADC: " << run_iterable << " " << iEnt << endl;}

                // if (mult_pulses.back().vp_int_18 > -450 && mult_pulses.back().vp_int_18 < -350) {cout << "\n" << "-450 - -350 ADC: " << run_iterable << " " << iEnt << endl;}

                // if (mult_pulses.back().vp_int_18 > 1400 && mult_pulses.back().vp_int_18 < 1500) {cout << "\n" << "1200 - 1400 ADC: " << run_iterable << " " << iEnt << endl;}
                            
                h_intmichels_PMT0->Fill(all_chan_energy[0]);

                h_intmichels_PMT1->Fill(all_chan_energy[1]);

                h_intmichels_PMT2->Fill(all_chan_energy[2]);

                h_intmichels_PMT3->Fill(all_chan_energy[3]);

                h_intmichels_PMT4->Fill(all_chan_energy[4]);

                h_intmichels_PMT5->Fill(all_chan_energy[5]);

                h_intmichels_PMT6->Fill(all_chan_energy[6]);

                h_intmichels_PMT7->Fill(all_chan_energy[7]);

                h_intmichels_PMT8->Fill(all_chan_energy[8]);

                h_intmichels_PMT9->Fill(all_chan_energy[9]);

                h_intmichels_PMT10->Fill(all_chan_energy[10]);

                h_intmichels_PMT11->Fill(all_chan_energy[11]);

                h_dt_61_16->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_16 - mult_pulses.front().wf_time); h_vpint_16->Fill(mult_pulses.back().vp_int_16);

                h_dt_61_17->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_17 - mult_pulses.front().wf_time); h_vpint_17->Fill(mult_pulses.back().vp_int_17);

                h_dt_61_18->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_18 - mult_pulses.front().wf_time); h_vpint_18->Fill(mult_pulses.back().vp_int_18);

                h_dt_61_19->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_19 - mult_pulses.front().wf_time); h_vpint_19->Fill(mult_pulses.back().vp_int_19);

                h_dt_61_20->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_20 - mult_pulses.front().wf_time); h_vpint_20->Fill(mult_pulses.back().vp_int_20);

                h_dt_61_21->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_21 - mult_pulses.front().wf_time); h_vpint_21->Fill(mult_pulses.back().vp_int_21);

                h_dt_61_22->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_22 - mult_pulses.front().wf_time); h_vpint_22->Fill(mult_pulses.back().vp_int_22);

                h_dt_61_23->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_23 - mult_pulses.front().wf_time); h_vpint_23->Fill(mult_pulses.back().vp_int_23);

                h_dt_61_top->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_top - mult_pulses.front().wf_time); h_vpint_top->Fill(mult_pulses.back().vp_int_top);

                // if (run == 7894 && mult_pulses.back().vp_int_16 > 100) {h_vpint_16->Fill(mult_pulses.back().vp_int_16);} else if (run == 11241 && mult_pulses.back().vp_int_16 > 100) {h_vpint_16->Fill(mult_pulses.back().vp_int_16);} else if (run == 12636 && mult_pulses.back().vp_int_16 > 150) {h_vpint_16->Fill(mult_pulses.back().vp_int_16);} else if (run == 15312 && mult_pulses.back().vp_int_16 > 300) {h_vpint_16->Fill(mult_pulses.back().vp_int_16);} else if (run == 16100 && mult_pulses.back().vp_int_16 > 150) {h_vpint_16->Fill(mult_pulses.back().vp_int_16);}

                // if (run == 7894 && mult_pulses.back().vp_int_17 > 100) {h_vpint_17->Fill(mult_pulses.back().vp_int_17);} else if (run == 11241 && mult_pulses.back().vp_int_17 > 100) {h_vpint_17->Fill(mult_pulses.back().vp_int_17);} else if (run == 12636 && mult_pulses.back().vp_int_17 > 150) {h_vpint_17->Fill(mult_pulses.back().vp_int_17);} else if (run == 15312 && mult_pulses.back().vp_int_17 > 300) {h_vpint_17->Fill(mult_pulses.back().vp_int_17);} else if (run == 16100 && mult_pulses.back().vp_int_17 > 150) {h_vpint_17->Fill(mult_pulses.back().vp_int_17);}

                // if (run == 7894 && mult_pulses.back().vp_int_18 > 200) {h_vpint_18->Fill(mult_pulses.back().vp_int_18);} else if (run == 11241 && mult_pulses.back().vp_int_18 > 300) {h_vpint_18->Fill(mult_pulses.back().vp_int_18);} else if (run == 12636 && mult_pulses.back().vp_int_18 > 400) {h_vpint_18->Fill(mult_pulses.back().vp_int_18);} else if (run == 15312 && mult_pulses.back().vp_int_18 > 600) {h_vpint_18->Fill(mult_pulses.back().vp_int_18);} else if (run == 16100 && mult_pulses.back().vp_int_18 > 400) {h_vpint_18->Fill(mult_pulses.back().vp_int_18);}

                // if (run == 7894 && mult_pulses.back().vp_int_19 > 100) {h_vpint_19->Fill(mult_pulses.back().vp_int_19);} else if (run == 11241 && mult_pulses.back().vp_int_19 > 100) {h_vpint_19->Fill(mult_pulses.back().vp_int_19);} else if (run == 12636 && mult_pulses.back().vp_int_19 > 150) {h_vpint_19->Fill(mult_pulses.back().vp_int_19);} else if (run == 15312 && mult_pulses.back().vp_int_19 > 300) {h_vpint_19->Fill(mult_pulses.back().vp_int_19);} else if (run == 16100 && mult_pulses.back().vp_int_19 > 150) {h_vpint_19->Fill(mult_pulses.back().vp_int_19);}

                // if (run == 7894 && mult_pulses.back().vp_int_20 > 100) {h_vpint_20->Fill(mult_pulses.back().vp_int_20);} else if (run == 11241 && mult_pulses.back().vp_int_20 > 100) {h_vpint_20->Fill(mult_pulses.back().vp_int_20);} else if (run == 12636 && mult_pulses.back().vp_int_20 > 150) {h_vpint_20->Fill(mult_pulses.back().vp_int_20);} else if (run == 15312 && mult_pulses.back().vp_int_20 > 300) {h_vpint_20->Fill(mult_pulses.back().vp_int_20);} else if (run == 16100 && mult_pulses.back().vp_int_20 > 150) {h_vpint_20->Fill(mult_pulses.back().vp_int_20);}

                // if (run == 7894 && mult_pulses.back().vp_int_21 > 100) {h_vpint_21->Fill(mult_pulses.back().vp_int_21);} else if (run == 11241 && mult_pulses.back().vp_int_21 > 100) {h_vpint_21->Fill(mult_pulses.back().vp_int_21);} else if (run == 12636 && mult_pulses.back().vp_int_21 > 150) {h_vpint_21->Fill(mult_pulses.back().vp_int_21);} else if (run == 15312 && mult_pulses.back().vp_int_21 > 300) {h_vpint_21->Fill(mult_pulses.back().vp_int_21);} else if (run == 16100 && mult_pulses.back().vp_int_21 > 150) {h_vpint_21->Fill(mult_pulses.back().vp_int_21);}

                // if (run == 7894 && mult_pulses.back().vp_int_22 > 100) {h_vpint_22->Fill(mult_pulses.back().vp_int_22);} else if (run == 11241 && mult_pulses.back().vp_int_22 > 100) {h_vpint_22->Fill(mult_pulses.back().vp_int_22);} else if (run == 12636 && mult_pulses.back().vp_int_22 > 150) {h_vpint_22->Fill(mult_pulses.back().vp_int_22);} else if (run == 15312 && mult_pulses.back().vp_int_22 > 300) {h_vpint_22->Fill(mult_pulses.back().vp_int_22);} else if (run == 16100 && mult_pulses.back().vp_int_22 > 150) {h_vpint_22->Fill(mult_pulses.back().vp_int_22);}

                // if (run == 7894 && mult_pulses.back().vp_int_23 > 100) {h_vpint_23->Fill(mult_pulses.back().vp_int_23);} else if (run == 11241 && mult_pulses.back().vp_int_23 > 100) {h_vpint_23->Fill(mult_pulses.back().vp_int_23);} else if (run == 12636 && mult_pulses.back().vp_int_23 > 150) {h_vpint_23->Fill(mult_pulses.back().vp_int_23);} else if (run == 15312 && mult_pulses.back().vp_int_23 > 300) {h_vpint_23->Fill(mult_pulses.back().vp_int_23);} else if (run == 16100 && mult_pulses.back().vp_int_23 > 150) {h_vpint_23->Fill(mult_pulses.back().vp_int_23);}

                // if (run == 7894 && mult_pulses.back().vp_int_top > 150) {h_vpint_top->Fill(mult_pulses.back().vp_int_top);} else if (run == 11241 && mult_pulses.back().vp_int_top > 200) {h_vpint_top->Fill(mult_pulses.back().vp_int_top);} else if (run == 12636 && mult_pulses.back().vp_int_top > 250) {h_vpint_top->Fill(mult_pulses.back().vp_int_top);} else if (run == 15312 && mult_pulses.back().vp_int_top > 300) {h_vpint_top->Fill(mult_pulses.back().vp_int_top);} else if (run == 16100 && mult_pulses.back().vp_int_top > 250) {h_vpint_top->Fill(mult_pulses.back().vp_int_top);}

                h_detint_vs_vpint_16->Fill(mult_pulses.back().energy, mult_pulses.back().vp_int_16);

                h_detint_vs_vpint_17->Fill(mult_pulses.back().energy, mult_pulses.back().vp_int_17);

                h_detint_vs_vpint_18->Fill(mult_pulses.back().energy, mult_pulses.back().vp_int_18);

                h_detint_vs_vpint_19->Fill(mult_pulses.back().energy, mult_pulses.back().vp_int_19);

                h_detint_vs_vpint_20->Fill(mult_pulses.back().energy, mult_pulses.back().vp_int_20);

                h_detint_vs_vpint_21->Fill(mult_pulses.back().energy, mult_pulses.back().vp_int_21);

                h_detint_vs_vpint_22->Fill(mult_pulses.back().energy, mult_pulses.back().vp_int_22);

                h_detint_vs_vpint_23->Fill(mult_pulses.back().energy, mult_pulses.back().vp_int_23);

                h_detint_vs_vpint_top->Fill(mult_pulses.back().energy, mult_pulses.back().vp_int_top);

                michel_event = true;

                // cout << "\n" << "mult_pulses.size() = " << mult_pulses.size() << endl;

                // cout << "\n" << "eventID of first peak: " << eventID_fp << endl;

                // cout << "\n" << "eventID of second peak: " << eventID << endl;

                // cout << "\n" << "Pair of peaks found. Time difference (t_c) = " << time_count << endl;

                // We know the first pulse is interesting; 
                //   don't care about dt between secondary peaks
                if (mult_pulses.size() > 1) {

                    int entry = iEnt;
                    bool issue = false;
                    for (int j = 1; j < mult_pulses.size(); j++) {

                        // cout << "\n" << "mult_pulses.front().energy: " << mult_pulses.front().energy << endl;

                        // cout << "\n" << "mult_pulses[j].energy: " << mult_pulses[j].energy << endl;

                        double e1 = mult_pulses.front().energy;
                        double p1 = mult_pulses.front().peak;
                        double t1 = mult_pulses.front().start;
                        double d1 = mult_pulses.front().end - t1;
                        double e2 = mult_pulses[j].energy;
                        double p2 = mult_pulses[j].peak;
                        double t2 = mult_pulses[j].start;
                        double d2 = mult_pulses[j].end - t2;
                        double dt = t2 - t1;
                        // Potential issues, any xi < 0 or e2 > e1          (Include new variables in detecting issues)
                        if (dt < 0 || e1 < 0 || e2 < 0 || e1 < e2 ||
                            p1 > maxPeak || p2 > maxPeak) {issue = true;}
                        double n1 = mult_pulses.front().number;
                        bool s1 = mult_pulses.front().single;
                        bool b1 = mult_pulses.front().beam;
                        double tr1 = mult_pulses.front().trigger;
                        double l1 = mult_pulses.front().length;
                        double n2 = mult_pulses[j].number;
                        bool s2 = mult_pulses[j].single;
                        bool b2 = mult_pulses[j].beam;
                        double tr2 = mult_pulses[j].trigger;
                        double l2 = mult_pulses[j].length;

                        michelTree->SetBranchAddress("entry", &entry);
                        michelTree->SetBranchAddress("e1", &e1);
                        michelTree->SetBranchAddress("p1", &p1);
                        michelTree->SetBranchAddress("t1", &t1);
                        michelTree->SetBranchAddress("d1", &d1);
                        michelTree->SetBranchAddress("e2", &e2);
                        michelTree->SetBranchAddress("p2", &p2);
                        michelTree->SetBranchAddress("t2", &t2);
                        michelTree->SetBranchAddress("d2", &d2);
                        michelTree->SetBranchAddress("dt", &dt);
                        michelTree->SetBranchAddress("issue", &issue);
                        michelTree->SetBranchAddress("n1", &n1);
                        michelTree->SetBranchAddress("s1", &s1);
                        michelTree->SetBranchAddress("b1", &b1);
                        michelTree->SetBranchAddress("tr1", &tr1);
                        michelTree->SetBranchAddress("l1", &l1);
                        michelTree->SetBranchAddress("n2", &n2);
                        michelTree->SetBranchAddress("s2", &s2);
                        michelTree->SetBranchAddress("b2", &b2);
                        michelTree->SetBranchAddress("tr2", &tr2);
                        michelTree->SetBranchAddress("l2", &l2);

                        michelTree->Fill();

                        write_out_counter += 1;

                        // cout << "\n" << "Fill Michel Tree " << write_out_counter << " Time(s) " << endl;

                        // cout << "\n" << "Pair of peaks written to ROOT tree. Time difference (dt) = " << dt << endl;

                        /* Fill histograms */

                        if (e2 > 0) {                       // e2 > 50

                            // h_dt->Fill(dt);
                            
                            // h_intmuons->Fill(e1);
                            
                            // h_intmichels->Fill(e2);

                            // Old histograms

                            h_dten->Fill(dt, iEnt + 1);     // dt vs time_count?

                            h_tc->Fill(time_count);

                            h_ampdiff->Fill(p2 - p1);

                            h_ampboth->Fill(p1, p2);

                            h_ampP1->Fill(p1);

                            h_ampP2->Fill(p2);

                            h_intboth->Fill(e1, e2);

                        }

                    }
                }

                // pulse_groups.push_back(mult_pulses); // pulse_groups is not the right vector type?

                mult_pulses.clear();

                WithinPulse = false;

                time_count = 0.0;

                MichelVetoActivity = true;

            }
            
            if (time_count > muon_dt_max) {            // time_count > 10000

                // if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.vp_int_top > 450) {

                // if (p_int_16 > vp_cut_vals[0] || p_int_17 > vp_cut_vals[1] || p_int_18 > vp_cut_vals[2] || p_int_19 > vp_cut_vals[3] || p_int_20 > vp_cut_vals[4] || p_int_21 > vp_cut_vals[5] || p_int_22 > vp_cut_vals[6] || p_int_23 > vp_cut_vals[7] || avg_pulse.vp_int_top > vp_cut_vals[8]) {

                if (p_int_16 > min_muon_vals[0] || p_int_17 > min_muon_vals[1] || p_int_18 > min_muon_vals[2] || p_int_19 > min_muon_vals[3] || p_int_20 > min_muon_vals[4] || p_int_21 > min_muon_vals[5] || p_int_22 > min_muon_vals[6] || p_int_23 > min_muon_vals[7] || avg_pulse.vp_int_top > min_muon_vals[8]) {

                    u_event_bin_val.clear();

                    u_event_bin_val = event_bin_val;

                    mult_pulses.clear();

                    mult_pulses.push_back(avg_pulse);

                    time_count = 0.5;

                    WithinPulse = true;

                    first_pulse_start_time = avg_pulse.start;

                    MichelVetoActivity = true;

                    // eventID_fp = eventID;
                
                }

                else {

                    mult_pulses.clear();

                    WithinPulse = false;

                    time_count = 0.0;

                    MichelVetoActivity = true;

                    // eventID_fp = eventID;

                }
                
            }

            // Are large energy detector events also being detected in veto panels?

            if (!michel_event && avg_pulse.number >= 10 && avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {

                if (avg_pulse.energy > 0 && avg_pulse.energy <= 200) {num_large_found_1 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_1 += 1;}}
                
                if (avg_pulse.energy > 200 && avg_pulse.energy <= 400) {num_large_found_2 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_2 += 1;}}

                if (avg_pulse.energy > 400 && avg_pulse.energy <= 600) {num_large_found_3 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_3 += 1;}}

                if (avg_pulse.energy > 600 && avg_pulse.energy <= 800) {num_large_found_4 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_4 += 1;}}

                if (avg_pulse.energy > 800 && avg_pulse.energy <= 1000) {num_large_found_5 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_5 += 1;}}

                if (avg_pulse.energy > 1000) {num_large_found_6 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_6 += 1;}}
            
                // if (avg_pulse.energy > 1200 && avg_pulse.energy <= 1400) {num_large_found_7 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_7 += 1;}}
            
                // if (avg_pulse.energy > 1400 && avg_pulse.energy <= 1600) {num_large_found_8 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_8 += 1;}}
            
                // if (avg_pulse.energy > 1600 && avg_pulse.energy <= 1800) {num_large_found_9 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_9 += 1;}}
            
                // if (avg_pulse.energy > 1800 && avg_pulse.energy <= 2000) {num_large_found_10 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_10 += 1;}}
            
            }

            if (avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {last_int_trig_event_time = avg_pulse.start;}

            last_event_time = avg_pulse.start;

            // Is this event an internally triggered event?

            if (avg_pulse.number >= 10 && avg_pulse.energy >= 50 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {avg_pulse.last_det_time = avg_pulse.start;}

            // Reset variables at end of loop
            
            michel_event = false;

            num_chan = 0;
            num_chan_over_1phe = 0;
            num_chan_simple = 0;

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

            bool no_muon_cut_event_found = false;

            bool no_muon_cut_event_found_rev = false;

            num_chan_vp = 0;

            p_int_16 = 0;
            p_int_17 = 0;
            p_int_18 = 0;
            p_int_19 = 0;
            p_int_20 = 0;
            p_int_21 = 0;
            p_int_22 = 0;
            p_int_23 = 0;
            p_int_24 = 0;
            p_int_25 = 0;
            p_start_16 = 0;
            p_start_17 = 0;
            p_start_18 = 0;
            p_start_19 = 0;
            p_start_20 = 0;
            p_start_21 = 0;
            p_start_22 = 0;
            p_start_23 = 0;
            p_start_24 = 0;
            p_start_25 = 0;

            // Ev61Energy = 0.;

        } // Event loop

        // cout << "\n" << "Total number of Michel Tree Fills: " << write_out_counter << "\n" << endl;

	    michelTree->Write();

        fileOut->Close();

        // michelTree->Clear();

        // cout << "\n" << "We wrote to file. " << endl;

        // cout << "\n" << "Variance of nsTime = " << variance(all_nsTime)  << endl;

        // cout << "\n" << "RMS of nsTime = " << rmsValue(all_nsTime)  << endl;

        // cout << "\n" << "Number of events in which Channel 15 is 'on' = " << ch15_on << endl;

        ch15_on = 0;

        all_nsTime.clear();

        amp_on_off_val.push_back(amp_val);

        amp_val.clear();

        integralToPE.clear();

        vp_cut_vals.clear();

        u_event_bin_val.clear();

        avg_peak_pos_RMS = getAverage(peak_pos_RMS);

        var_peak_pos_RMS = variance(peak_pos_RMS);

        // cout << "\n" << "Average peak position RMS within each detector event = " << avg_peak_pos_RMS << endl;

        // cout << "\n" << "Variance of average peak position RMS within each detector event = " << var_peak_pos_RMS << "\n" << endl;

        peak_pos_RMS.clear();

        run_counter += 1;

        if (run_counter == 12) {hour12_st = run_starttime; hour12_run = run_iterable;}

        else if (run_counter == 24) {run_counter = 0;}

        f->Close();

    } // Run loop

    cout << "Number of HL events separated by less than 18ms, Run " << run << " = " << HL_dt_17ms_counter << "\n" << endl;

    cout << "Number of cosmic events, Run " << run << " = " << cosmic_events << "\n" << endl;          // runlist[0]

    cout << "Number of Michel electron events, Run " << run << " = " << hist_events << "\n" << endl;           // runlist[0]

    cout << "Average Michel electron energy, Run " << run << " = " << tot_michel_energy * pow(hist_events, -1) << "\n" << "\n" << endl;            // runlist[0]

    double avg_mich_en = tot_michel_energy * pow(hist_events, -1);

    // Write to the file
    MIFile << hour12_st << "\t" << avg_mich_en << "\n";

    MERFile << hour12_st << "\t" << hist_events << "\n";

    cosmic_events = 0;

    hist_events = 0;

    tot_michel_energy = 0.0;

    // Close the file
    MIFile.close();

    MIFile.clear();

    MERFile.close();

    MERFile.clear();

    MBVFile.close();

    MBVFile.clear();

    uBVFile.close();

    uBVFile.clear();

    HLBVFile.close();

    HLBVFile.clear();

    LLBVFile.close();

    LLBVFile.clear();

    MBBVFile.close();

    MBBVFile.clear();

    /* Fill & plot histograms */
    
    if (run != 11241 && run != 12636) {

        TCanvas* c_dt = new TCanvas("c_dt", "Delta-T", 900, 700);
        // c_dt->SetLogy();
        c_dt->cd();

        TF1* func1 = new TF1("func1", "[0]+[1]*exp(-x/[2])", 2100, 12100);
        double Par1[] = {100, 15000, 2100};
        func1->SetParameters(Par1);
        func1->SetParName(0, "Baseline");
        func1->SetParName(1, "Amplitude");
        func1->SetParName(2, "Time Constant");
        h_dt->Fit("func1", "R");
        Double_t chi2_1 = func1->GetChisquare();

        cout << "\n" << "Chi-square of e^x/a fit: " << chi2_1 << "\n" << endl;

        // TF1* func2 = new TF1("func2", "[0]+[1]*exp(-x/2200)+[2]*exp(-x/1800)", 2100, 12100);
        // func2->SetParNames("Baseline", "Amplitude_u+", "Amplitude_u-");
        // func2->SetParameters(100, 15000, 15000);
        // h_dt->Fit("func2", "R");
        // Double_t chi2_2 = func2->GetChisquare();

        // cout << "\n" << "Chi-square of e^x/2200 + e^x/1800 fit: " << chi2_2 << "\n" << endl;

        /*TF1* func3 = new TF1("func3", "[0]+[1]*exp(-x/[2])+[3]*exp(-x/[4])", 2100, 12100);

        double Par3[] = {100, 15000, 2200, 15000, 1800};
        func3->SetParameters(Par3);
        func3->SetParName(0, "Baseline");
        func3->SetParName(1, "Amplitude_u+");
        func3->SetParName(2, "Time Constant_u+");
        func3->SetParName(3, "Amplitude_u-");
        func3->SetParName(4, "Time Constant_u-");
        h_dt->Fit("func3", "R");
        Double_t chi2_3 = func3->GetChisquare();

        cout << "\n" << "Chi-square of e^x/a + e^x/b fit: " << chi2_3 << endl;*/

        h_dt->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt->GetYaxis()->SetTitle("Counts");
        h_dt->GetYaxis()->SetRange(0, 200);
        h_dt->Draw();
        func1->Draw("same");
        func1->SetLineColor(2);
        // func2->Draw("same");
        // func2->SetLineColor(3);
        // func3->Draw("same");
        // func3->SetLineColor(3);
        // TLegend *leg_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        // leg_dt->AddEntry(func1, "e^x/a fit", "l");
        // leg_dt->AddEntry(func2, "e^x/2200 + e^x/1800 fit", "l");
        // leg_dt->AddEntry(func3, "e^x/a + e^x/b fit", "l");
        // leg_dt->Draw();

        Double_t param10 = func1->GetParameter(0);
        Double_t param11 = func1->GetParameter(1);
        Double_t param12 = func1->GetParameter(2);

        Double_t error10 = func1->GetParError(0);
        Double_t error11 = func1->GetParError(1);
        Double_t error12 = func1->GetParError(2);

        TLatex *lat_dt_1 = new TLatex(0.3, 0.7, Form("Baseline = %.2f #pm %.2f", param10, error10));
        lat_dt_1->SetNDC();
        lat_dt_1->SetTextColor(1);
        lat_dt_1->SetTextSize(0.04);
        lat_dt_1->Draw();

        TLatex *lat_dt_2 = new TLatex(0.3, 0.65, Form("Amplitude = %.2f #pm %.2f", param11, error11));
        lat_dt_2->SetNDC();
        lat_dt_2->SetTextColor(1);
        lat_dt_2->SetTextSize(0.04);
        lat_dt_2->Draw();

        TLatex *lat_dt_3 = new TLatex(0.3, 0.6, Form("Time Constant (ns) = %.2f #pm %.2f", param12, error12));
        lat_dt_3->SetNDC();
        lat_dt_3->SetTextColor(1);
        lat_dt_3->SetTextSize(0.04);
        lat_dt_3->Draw();
        
        /*TLatex *lat_dt_9 = new TLatex(0.3, 0.7, Form("e^x/a chi-square = %.2f", chi2_1));
        lat_dt_9->SetNDC();
        lat_dt_9->SetTextColor(1);
        lat_dt_9->SetTextSize(0.02);
        lat_dt_9->Draw();

        Double_t param30 = func3->GetParameter(0);
        Double_t param31 = func3->GetParameter(1);
        Double_t param32 = func3->GetParameter(2);
        Double_t param33 = func3->GetParameter(3);
        Double_t param34 = func3->GetParameter(4);

        TLatex *lat_dt_4 = new TLatex(0.3, 0.65, Form("e^x/a + e^x/b baseline = %.2f", param30));
        lat_dt_4->SetNDC();
        lat_dt_4->SetTextColor(1);
        lat_dt_4->SetTextSize(0.02);
        lat_dt_4->Draw();

        TLatex *lat_dt_5 = new TLatex(0.3, 0.6, Form("e^x/a + e^x/b u+ amplitude = %.2f", param31));
        lat_dt_5->SetNDC();
        lat_dt_5->SetTextColor(1);
        lat_dt_5->SetTextSize(0.02);
        lat_dt_5->Draw();

        TLatex *lat_dt_6 = new TLatex(0.3, 0.55, Form("e^x/a + e^x/b u+ time constant (ns) = %.2f", param32));
        lat_dt_6->SetNDC();
        lat_dt_6->SetTextColor(1);
        lat_dt_6->SetTextSize(0.02);
        lat_dt_6->Draw();

        TLatex *lat_dt_7 = new TLatex(0.3, 0.5, Form("e^x/a + e^x/b u- amplitude = %.2f", param33));
        lat_dt_7->SetNDC();
        lat_dt_7->SetTextColor(1);
        lat_dt_7->SetTextSize(0.02);
        lat_dt_7->Draw();

        TLatex *lat_dt_8 = new TLatex(0.3, 0.45, Form("e^x/a + e^x/b u- time constant (ns) = %.2f", param34));
        lat_dt_8->SetNDC();
        lat_dt_8->SetTextColor(1);
        lat_dt_8->SetTextSize(0.02);
        lat_dt_8->Draw();

        TLatex *lat_dt_10 = new TLatex(0.3, 0.4, Form("e^x/a + e^x/b chi-square = %.2f", chi2_3));
        lat_dt_10->SetNDC();
        lat_dt_10->SetTextColor(1);
        lat_dt_10->SetTextSize(0.02);
        lat_dt_10->Draw();*/

        c_dt->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_%i.png", run));
        c_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_dt.png"));
        h_dt->Reset();

        h_intmuons->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmuons->GetYaxis()->SetTitle("Counts");
        h_intmuons->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmuons_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmuons.png"));
        h_intmuons->Reset();

        TCanvas* c_intmichels = new TCanvas("c_intmichels", "Michel Electron Integral Values", 1200, 700);
        // c_intmichels->SetLogy();
        c_intmichels->cd();

        TF1* func_intmichels = new TF1("func_intmichels", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))", 60, 280);
        double par_michels[] = {500, 4000, 150, 100};
        func_intmichels->SetParameters(par_michels);
        func_intmichels->SetParName(0, "Baseline");
        func_intmichels->SetParName(1, "Amplitude");
        func_intmichels->SetParName(2, "PeakCenter");
        func_intmichels->SetParName(3, "StandardDeviation");
        h_intmichels->Fit("func_intmichels", "R");
        double_t am_michels = func_intmichels->GetParameter(1);
        double_t pc_michels = func_intmichels->GetParameter(2);

        h_intmichels->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels->GetYaxis()->SetTitle("Counts");
        h_intmichels->Draw();
        TLatex *lat_intmichels_1 = new TLatex(0.6, 0.7, Form("Average Michel energy = %.2f", avg_mich_en));
        lat_intmichels_1->SetNDC();
        lat_intmichels_1->SetTextColor(1);
        lat_intmichels_1->SetTextSize(0.04);
        lat_intmichels_1->Draw();
        c_intmichels->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_%i.png", run));
        // c_intmichels->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_intmichels.png");
        c_intmichels->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels.png"));
        // h_intmichels->Print("all");
        h_intmichels->Reset();

        MIAFile << hour12_st << "\t" << am_michels << endl;
        MIPFile << hour12_st << "\t" << pc_michels << endl;

        TCanvas* c_dt_v_intmichels = new TCanvas("c_dt_v_intmichels", "Event 61 Detector dt vs Integral Value", 1200, 700);
        // c_dt_v_intmichels->SetLogy();
        c_dt_v_intmichels->cd();
        h_dt_v_intmichels->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_intmichels->GetYaxis()->SetTitle("Integral (Ph.e.)");
        h_dt_v_intmichels->SetMarkerStyle(7);
        h_dt_v_intmichels->Draw();
        c_dt_v_intmichels->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_v_intmichels_%i.png", run));
        // c_dt_v_intmichels->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_v_intmichels.png");
        c_dt_v_intmichels->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_intmichels.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_dt_v_intmichels.png"));
        h_dt_v_intmichels->Reset();

        h_intcosmics->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intcosmics->GetYaxis()->SetTitle("Counts");
        h_intcosmics->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intcosmics_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intcosmics.png"));
        // h_intcosmics->Print("all");
        h_intcosmics->Reset();

        TCanvas* c_muon_pmt_sipm_dt = new TCanvas("c_muon_pmt_sipm_dt", "Muon PMT vs SiPM dt", 1200, 700);
        // c_muon_pmt_sipm_dt->SetLogy();
        c_muon_pmt_sipm_dt->cd();

        TF1* func_pmt_sipm_dt = new TF1("func_pmt_sipm_dt", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -40, 60);
        double par_pmt_sipm_dt[] = {1000, 6000000, 10, 20};
        func_pmt_sipm_dt->SetParameters(par_pmt_sipm_dt);
        func_pmt_sipm_dt->SetParName(0, "Baseline");
        func_pmt_sipm_dt->SetParName(1, "Amplitude");
        func_pmt_sipm_dt->SetParName(2, "PeakCenter");
        func_pmt_sipm_dt->SetParName(3, "StandardDeviation");
        h_muon_pmt_sipm_dt->Fit("func_pmt_sipm_dt", "R");
        double_t dt_center = func_pmt_sipm_dt->GetParameter(2);
        double_t dt_sigma = func_pmt_sipm_dt->GetParameter(3);

        h_muon_pmt_sipm_dt->GetXaxis()->SetTitle("Delta-T (ns)");
        h_muon_pmt_sipm_dt->GetYaxis()->SetTitle("Counts");
        // h_muon_pmt_sipm_dt->GetYaxis()->SetRange(0, 50);
        h_muon_pmt_sipm_dt->Draw();
        TLatex *lat_pmt_sipm_dt = new TLatex(0.6, 0.7, Form("Peak Center #pm Sigma = %.2f #pm %.2f", dt_center, dt_sigma));
        lat_pmt_sipm_dt->SetNDC();
        lat_pmt_sipm_dt->SetTextColor(1);
        lat_pmt_sipm_dt->SetTextSize(0.04);
        lat_pmt_sipm_dt->Draw();
        func_pmt_sipm_dt->Draw("same");
        func_pmt_sipm_dt->SetLineColor(2);
        
        c_muon_pmt_sipm_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_muon_pmt_sipm_dt.png");
        // c_muon_pmt_sipm_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_muon_pmt_sipm_dt.png");
        c_muon_pmt_sipm_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_muon_pmt_sipm_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_muon_pmt_sipm_dt.png"));
        h_muon_pmt_sipm_dt->Reset();

        TCanvas* c_only_one_sipm = new TCanvas("c_only_one_sipm", "Only One SiPM Integrals", 1200, 700);
        c_only_one_sipm->SetLogy();
        c_only_one_sipm->cd();
        h_only_one_sipm->GetXaxis()->SetTitle("Integral (ADC)");
        h_only_one_sipm->GetYaxis()->SetTitle("Counts");
        // h_only_one_sipm->GetYaxis()->SetRange(0, 50);
        h_only_one_sipm->Draw();
        c_only_one_sipm->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_only_one_sipm.png");
        // c_only_one_sipm->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_only_one_sipm.png");
        c_only_one_sipm->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_only_one_sipm.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_only_one_sipm.png"));
        h_only_one_sipm->Reset();
        
        TCanvas* c_34detint = new TCanvas("c_34detint", "Distribution of All tB34 Det. Integral Values", 1200, 700);
        c_34detint->SetLogy();
        c_34detint->cd();
        h_34detint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_34detint->GetYaxis()->SetTitle("Counts");
        h_34detint->Draw();
        c_34detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_34detint.png");
        // c_34detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_34detint.png");
        c_34detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_34detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_34detint.png"));
        h_34detint->Reset();
        // delete TCanvas c_34detint;           // Every time you have "new", need "delete" after

        TCanvas* c_34vpint = new TCanvas("c_34vpint", "Distribution of All tB34 VP Integral Values", 1200, 700);
        c_34vpint->SetLogy();
        c_34vpint->cd();
        h_34vpint->GetXaxis()->SetTitle("Integral (ADC)");
        h_34vpint->GetYaxis()->SetTitle("Counts");
        h_34vpint->Draw();
        c_34vpint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_34vpint.png");
        // c_34vpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_34vpint.png");
        c_34vpint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_34vpint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_34vpint.png"));
        h_34vpint->Reset();

        TCanvas* c_slice_vpintavg = new TCanvas("c_slice_vpintavg", "Distribution of Averaged Veto Panel Integral Values", 1200, 700);
        c_slice_vpintavg->SetLogy();
        c_slice_vpintavg->cd();
        h_slice_vpintavg->GetXaxis()->SetTitle("Integral (ADC)");
        h_slice_vpintavg->GetYaxis()->SetTitle("Counts");
        h_slice_vpintavg->Draw();
        c_slice_vpintavg->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_slice_vpintavg.png");
        // c_slice_vpintavg->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_slice_vpintavg.png");
        c_slice_vpintavg->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_slice_vpintavg.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_slice_vpintavg.png"));
        h_slice_vpintavg->Reset();

        TCanvas* c_dt_61_16 = new TCanvas("c_dt_61_16", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_16->SetLogy();
        c_dt_61_16->cd();
        h_dt_61_16->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_16->GetYaxis()->SetTitle("Counts");
        h_dt_61_16->Draw();
        c_dt_61_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_16.png");
        // c_dt_61_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_16.png");
        c_dt_61_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_16.png"));
        h_dt_61_16->Reset();

        TCanvas* c_dt_61_17 = new TCanvas("c_dt_61_17", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_17->SetLogy();
        c_dt_61_17->cd();
        h_dt_61_17->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_17->GetYaxis()->SetTitle("Counts");
        h_dt_61_17->Draw();
        c_dt_61_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_17.png");
        // c_dt_61_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_17.png");
        c_dt_61_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_17.png"));
        h_dt_61_17->Reset();

        TCanvas* c_dt_61_18 = new TCanvas("c_dt_61_18", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_18->SetLogy();
        c_dt_61_18->cd();
        h_dt_61_18->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_18->GetYaxis()->SetTitle("Counts");
        h_dt_61_18->Draw();
        c_dt_61_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_18.png");
        // c_dt_61_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_18.png");
        c_dt_61_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_18.png"));
        h_dt_61_18->Reset();

        TCanvas* c_dt_61_19 = new TCanvas("c_dt_61_19", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_19->SetLogy();
        c_dt_61_19->cd();
        h_dt_61_19->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_19->GetYaxis()->SetTitle("Counts");
        h_dt_61_19->Draw();
        c_dt_61_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_19.png");
        // c_dt_61_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_19.png");
        c_dt_61_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_19.png"));
        h_dt_61_19->Reset();

        TCanvas* c_dt_61_20 = new TCanvas("c_dt_61_20", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_20->SetLogy();
        c_dt_61_20->cd();
        h_dt_61_20->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_20->GetYaxis()->SetTitle("Counts");
        h_dt_61_20->Draw();
        c_dt_61_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_20.png");
        // c_dt_61_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_20.png");
        c_dt_61_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_20.png"));
        h_dt_61_20->Reset();

        TCanvas* c_dt_61_21 = new TCanvas("c_dt_61_21", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_21->SetLogy();
        c_dt_61_21->cd();
        h_dt_61_21->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_21->GetYaxis()->SetTitle("Counts");
        h_dt_61_21->Draw();
        c_dt_61_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_21.png");
        // c_dt_61_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_21.png");
        c_dt_61_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_21.png"));
        h_dt_61_21->Reset();

        TCanvas* c_dt_61_22 = new TCanvas("c_dt_61_22", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_22->SetLogy();
        c_dt_61_22->cd();
        h_dt_61_22->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_22->GetYaxis()->SetTitle("Counts");
        h_dt_61_22->Draw();
        c_dt_61_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_22.png");
        // c_dt_61_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_22.png");
        c_dt_61_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_22.png"));
        h_dt_61_22->Reset();

        TCanvas* c_dt_61_23 = new TCanvas("c_dt_61_23", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_23->SetLogy();
        c_dt_61_23->cd();
        h_dt_61_23->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_23->GetYaxis()->SetTitle("Counts");
        h_dt_61_23->Draw();
        c_dt_61_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_23.png");
        // c_dt_61_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_23.png");
        c_dt_61_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_23.png"));
        h_dt_61_23->Reset();

        TCanvas* c_dt_61_top = new TCanvas("c_dt_61_top", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_top->SetLogy();
        c_dt_61_top->cd();
        h_dt_61_top->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_top->GetYaxis()->SetTitle("Counts");
        h_dt_61_top->Draw();
        c_dt_61_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_top.png");
        // c_dt_61_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_top.png");
        c_dt_61_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_top.png"));
        h_dt_61_top->Reset();

        TCanvas* c_vpint_16 = new TCanvas("c_vpint_16", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_16->SetLogy();
        c_vpint_16->cd();
        h_vpint_16->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_16->GetYaxis()->SetTitle("Counts");
        h_vpint_16->Draw();
        c_vpint_16->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_16_%i.png", run));
        // c_vpint_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_16.png");
        c_vpint_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_16.png"));
        h_vpint_16->Reset();

        TCanvas* c_vpint_17 = new TCanvas("c_vpint_17", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_17->SetLogy();
        c_vpint_17->cd();
        h_vpint_17->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_17->GetYaxis()->SetTitle("Counts");
        h_vpint_17->Draw();
        c_vpint_17->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_17_%i.png", run));
        // c_vpint_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_17.png");
        c_vpint_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_17.png"));
        h_vpint_17->Reset();

        TCanvas* c_vpint_18 = new TCanvas("c_vpint_18", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_18->SetLogy();
        c_vpint_18->cd();
        h_vpint_18->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_18->GetYaxis()->SetTitle("Counts");
        h_vpint_18->Draw();
        c_vpint_18->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_18_%i.png", run));
        // c_vpint_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_18.png");
        c_vpint_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_18.png"));
        h_vpint_18->Reset();

        TCanvas* c_vpint_19 = new TCanvas("c_vpint_19", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_19->SetLogy();
        c_vpint_19->cd();
        h_vpint_19->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_19->GetYaxis()->SetTitle("Counts");
        h_vpint_19->Draw();
        c_vpint_19->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_19_%i.png", run));
        // c_vpint_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_19.png");
        c_vpint_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_19.png"));
        h_vpint_19->Reset();

        TCanvas* c_vpint_20 = new TCanvas("c_vpint_20", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_20->SetLogy();
        c_vpint_20->cd();
        h_vpint_20->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_20->GetYaxis()->SetTitle("Counts");
        h_vpint_20->Draw();
        c_vpint_20->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_20_%i.png", run));
        // c_vpint_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_20.png");
        c_vpint_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_20.png"));
        h_vpint_20->Reset();

        TCanvas* c_vpint_21 = new TCanvas("c_vpint_21", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_21->SetLogy();
        c_vpint_21->cd();
        h_vpint_21->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_21->GetYaxis()->SetTitle("Counts");
        h_vpint_21->Draw();
        c_vpint_21->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_21_%i.png", run));
        // c_vpint_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_21.png");
        c_vpint_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_21.png"));
        h_vpint_21->Reset();

        TCanvas* c_vpint_22 = new TCanvas("c_vpint_22", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_22->SetLogy();
        c_vpint_22->cd();
        h_vpint_22->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_22->GetYaxis()->SetTitle("Counts");
        h_vpint_22->Draw();
        c_vpint_22->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_22_%i.png", run));
        // c_vpint_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_22.png");
        c_vpint_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_22.png"));
        h_vpint_22->Reset();

        TCanvas* c_vpint_23 = new TCanvas("c_vpint_23", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_23->SetLogy();
        c_vpint_23->cd();
        h_vpint_23->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_23->GetYaxis()->SetTitle("Counts");
        h_vpint_23->Draw();
        c_vpint_23->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_23_%i.png", run));
        // c_vpint_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_23.png");
        c_vpint_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_23.png"));
        h_vpint_23->Reset();

        TCanvas* c_vpint_top = new TCanvas("c_vpint_top", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_top->SetLogy();
        c_vpint_top->cd();
        h_vpint_top->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_top->GetYaxis()->SetTitle("Counts");
        h_vpint_top->Draw();
        c_vpint_top->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_top_%i.png", run));
        // c_vpint_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_top.png");
        c_vpint_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_top.png"));
        h_vpint_top->Reset();

        TCanvas* c_detint_vs_vpint_16 = new TCanvas("c_detint_vs_vpint_16", "SiPM 1 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_16->SetLogy();
        // c_detint_vs_vpint_16->SetLogx();
        c_detint_vs_vpint_16->cd();
        h_detint_vs_vpint_16->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_16->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_16->Draw();
        c_detint_vs_vpint_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_16.png");
        // c_detint_vs_vpint_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_16.png");
        c_detint_vs_vpint_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_16.png"));
        h_detint_vs_vpint_16->Reset();

        TCanvas* c_detint_vs_vpint_17 = new TCanvas("c_detint_vs_vpint_17", "SiPM 2 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_17->SetLogy();
        // c_detint_vs_vpint_17->SetLogx();
        c_detint_vs_vpint_17->cd();
        h_detint_vs_vpint_17->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_17->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_17->Draw();
        c_detint_vs_vpint_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_17.png");
        // c_detint_vs_vpint_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_17.png");
        c_detint_vs_vpint_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_17.png"));
        h_detint_vs_vpint_17->Reset();

        TCanvas* c_detint_vs_vpint_18 = new TCanvas("c_detint_vs_vpint_18", "SiPM 3 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_18->SetLogy();
        // c_detint_vs_vpint_18->SetLogx();
        c_detint_vs_vpint_18->cd();
        h_detint_vs_vpint_18->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_18->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_18->Draw();
        c_detint_vs_vpint_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_18.png");
        // c_detint_vs_vpint_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_18.png");
        c_detint_vs_vpint_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_18.png"));
        h_detint_vs_vpint_18->Reset();

        TCanvas* c_detint_vs_vpint_19 = new TCanvas("c_detint_vs_vpint_19", "SiPM 4 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_19->SetLogy();
        // c_detint_vs_vpint_19->SetLogx();
        c_detint_vs_vpint_19->cd();
        h_detint_vs_vpint_19->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_19->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_19->Draw();
        c_detint_vs_vpint_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_19.png");
        // c_detint_vs_vpint_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_19.png");
        c_detint_vs_vpint_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_19.png"));
        h_detint_vs_vpint_19->Reset();

        TCanvas* c_detint_vs_vpint_20 = new TCanvas("c_detint_vs_vpint_20", "SiPM 5 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_20->SetLogy();
        // c_detint_vs_vpint_20->SetLogx();
        c_detint_vs_vpint_20->cd();
        h_detint_vs_vpint_20->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_20->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_20->Draw();
        c_detint_vs_vpint_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_20.png");
        // c_detint_vs_vpint_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_20.png");
        c_detint_vs_vpint_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_20.png"));
        h_detint_vs_vpint_20->Reset();

        TCanvas* c_detint_vs_vpint_21 = new TCanvas("c_detint_vs_vpint_21", "SiPM 6 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_21->SetLogy();
        // c_detint_vs_vpint_21->SetLogx();
        c_detint_vs_vpint_21->cd();
        h_detint_vs_vpint_21->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_21->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_21->Draw();
        c_detint_vs_vpint_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_21.png");
        // c_detint_vs_vpint_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_21.png");
        c_detint_vs_vpint_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_21.png"));
        h_detint_vs_vpint_21->Reset();

        TCanvas* c_detint_vs_vpint_22 = new TCanvas("c_detint_vs_vpint_22", "SiPM 7 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_22->SetLogy();
        // c_detint_vs_vpint_22->SetLogx();
        c_detint_vs_vpint_22->cd();
        h_detint_vs_vpint_22->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_22->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_22->Draw();
        c_detint_vs_vpint_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_22.png");
        // c_detint_vs_vpint_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_22.png");
        c_detint_vs_vpint_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_22.png"));
        h_detint_vs_vpint_22->Reset();

        TCanvas* c_detint_vs_vpint_23 = new TCanvas("c_detint_vs_vpint_23", "SiPM 8 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_23->SetLogy();
        // c_detint_vs_vpint_23->SetLogx();
        c_detint_vs_vpint_23->cd();
        h_detint_vs_vpint_23->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_23->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_23->Draw();
        c_detint_vs_vpint_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_23.png");
        // c_detint_vs_vpint_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_23.png");
        c_detint_vs_vpint_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_23.png"));
        h_detint_vs_vpint_23->Reset();

        TCanvas* c_detint_vs_vpint_top = new TCanvas("c_detint_vs_vpint_top", "Top SiPM & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_top->SetLogy();
        // c_detint_vs_vpint_top->SetLogx();
        c_detint_vs_vpint_top->cd();
        h_detint_vs_vpint_top->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_top->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_top->Draw();
        c_detint_vs_vpint_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_top.png");
        // c_detint_vs_vpint_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_detint_vs_vpint_top.png");
        c_detint_vs_vpint_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_top.png"));
        h_detint_vs_vpint_top->Reset();

        TCanvas* c_tB = new TCanvas("c_tB", "TriggerBits Distribution", 1200, 700);
        c_tB->SetLogy();
        c_tB->cd();
        h_tB->GetXaxis()->SetTitle("triggerBits Value");
        h_tB->GetYaxis()->SetTitle("Counts");
        h_tB->Draw();
        c_tB->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_tB_%i.png", run));
        c_tB->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_tB.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_tB.png"));
        h_tB->Reset();

        TCanvas* c_blM_pmt = new TCanvas("c_blM_pmt", "PMT baselineMean Distribution", 1200, 700);
        c_blM_pmt->SetLogy();
        c_blM_pmt->cd();
        h_blM_pmt->GetXaxis()->SetTitle("baselineMean Value");
        h_blM_pmt->GetYaxis()->SetTitle("Counts");
        h_blM_pmt->Draw();
        c_blM_pmt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_blM_pmt.png");
        c_blM_pmt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_blM_pmt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_blM_pmt.png"));
        h_blM_pmt->Reset();

        TCanvas* c_blM_vp = new TCanvas("c_blM_vp", "SiPM baselineMean Distribution", 1200, 700);
        c_blM_vp->SetLogy();
        c_blM_vp->cd();
        h_blM_vp->GetXaxis()->SetTitle("baselineMean Value");
        h_blM_vp->GetYaxis()->SetTitle("Counts");
        h_blM_vp->Draw();
        c_blM_vp->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_blM_vp.png");
        c_blM_vp->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_blM_vp.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_blM_vp.png"));
        h_blM_vp->Reset();

        TCanvas* c_blM_ev61 = new TCanvas("c_blM_ev61", "Event 61 baselineMean Distribution", 1200, 700);
        c_blM_ev61->SetLogy();
        c_blM_ev61->cd();
        h_blM_ev61->GetXaxis()->SetTitle("baselineMean Value");
        h_blM_ev61->GetYaxis()->SetTitle("Counts");
        h_blM_ev61->Draw();
        c_blM_ev61->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_blM_ev61.png");
        c_blM_ev61->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_blM_ev61.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_blM_ev61.png"));
        h_blM_ev61->Reset();

        TCanvas* c_intev61 = new TCanvas("c_intev61", "Distribution of Event 61 Integral Values", 1200, 700);
        c_intev61->SetLogy();
        c_intev61->cd();
        h_intev61->GetXaxis()->SetTitle("Integral (ADC)");
        h_intev61->GetYaxis()->SetTitle("Counts");
        h_intev61->Draw();
        c_intev61->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intev61_%i.png", run));
        c_intev61->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_intev61.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intev61.png"));
        h_intev61->Reset();

        h_intmichels_PMT0->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels_PMT0->GetYaxis()->SetTitle("Counts");
        h_intmichels_PMT0->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_PMT0_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels_PMT0.png"));
        h_intmichels_PMT0->Reset();

        h_intmichels_PMT1->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels_PMT1->GetYaxis()->SetTitle("Counts");
        h_intmichels_PMT1->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_PMT1_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels_PMT1.png"));
        h_intmichels_PMT1->Reset();

        h_intmichels_PMT2->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels_PMT2->GetYaxis()->SetTitle("Counts");
        h_intmichels_PMT2->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_PMT2_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels_PMT2.png"));
        h_intmichels_PMT2->Reset();

        h_intmichels_PMT3->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels_PMT3->GetYaxis()->SetTitle("Counts");
        h_intmichels_PMT3->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_PMT3_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels_PMT3.png"));
        h_intmichels_PMT3->Reset();

        h_intmichels_PMT4->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels_PMT4->GetYaxis()->SetTitle("Counts");
        h_intmichels_PMT4->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_PMT4_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels_PMT4.png"));
        h_intmichels_PMT4->Reset();

        h_intmichels_PMT5->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels_PMT5->GetYaxis()->SetTitle("Counts");
        h_intmichels_PMT5->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_PMT5_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels_PMT5.png"));
        h_intmichels_PMT5->Reset();

        h_intmichels_PMT6->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels_PMT6->GetYaxis()->SetTitle("Counts");
        h_intmichels_PMT6->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_PMT6_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels_PMT6.png"));
        h_intmichels_PMT6->Reset();

        h_intmichels_PMT7->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels_PMT7->GetYaxis()->SetTitle("Counts");
        h_intmichels_PMT7->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_PMT7_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels_PMT7.png"));
        h_intmichels_PMT7->Reset();

        h_intmichels_PMT8->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels_PMT8->GetYaxis()->SetTitle("Counts");
        h_intmichels_PMT8->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_PMT8_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels_PMT8.png"));
        h_intmichels_PMT8->Reset();

        h_intmichels_PMT9->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels_PMT9->GetYaxis()->SetTitle("Counts");
        h_intmichels_PMT9->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_PMT9_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels_PMT9.png"));
        h_intmichels_PMT9->Reset();

        h_intmichels_PMT10->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels_PMT10->GetYaxis()->SetTitle("Counts");
        h_intmichels_PMT10->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_PMT10_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels_PMT10.png"));
        h_intmichels_PMT10->Reset();

        h_intmichels_PMT11->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels_PMT11->GetYaxis()->SetTitle("Counts");
        h_intmichels_PMT11->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intmichels_PMT11_%i.png", run));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intmichels_PMT11.png"));
        h_intmichels_PMT11->Reset();

    }

    if (run == 11241) {

        TCanvas* c_dt = new TCanvas("c_dt", "Delta-T", 900, 700);
        // c_dt->SetLogy();
        c_dt->cd();

        TF1* func1 = new TF1("func1", "[0]+[1]*exp(-x/[2])", 2100, 12100);
        func1->SetParNames("Baseline", "Amplitude", "Time Constant");
        func1->SetParameters(100, 15000, 2100);
        h_dt->Fit("func1", "R");
        Double_t chi2_1 = func1->GetChisquare();

        cout << "\n" << "Chi-square of e^x/a fit: " << chi2_1 << "\n" << endl;

        // TF1* func2 = new TF1("func2", "[0]+[1]*exp(-x/2200)+[2]*exp(-x/1800)", 2100, 12100);
        // func2->SetParNames("Baseline", "Amplitude_u+", "Amplitude_u-");
        // func2->SetParameters(100, 15000, 15000);
        // h_dt->Fit("func2", "R");
        // Double_t chi2_2 = func2->GetChisquare();

        // cout << "\n" << "Chi-square of e^x/2200 + e^x/1800 fit: " << chi2_2 << "\n" << endl;

        TF1* func3 = new TF1("func3", "[0]+[1]*exp(-x/[2])+[3]*exp(-x/[4])", 2100, 12100);
        func3->SetParNames("Baseline", "Amplitude_u+", "Time Constant_u+", "Amplitude_u-", "Time Constant_u-");
        func3->SetParameters(100, 15000, 2200, 15000, 1800);
        h_dt->Fit("func3", "R");
        Double_t chi2_3 = func3->GetChisquare();

        cout << "\n" << "Chi-square of e^x/a + e^x/b fit: " << chi2_3 << endl;

        h_dt->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt->GetYaxis()->SetTitle("Counts");
        h_dt->GetYaxis()->SetRange(0, 200);
        h_dt->Draw();
        func1->Draw("same");
        func1->SetLineColor(2);
        // func2->Draw("same");
        // func2->SetLineColor(3);
        func3->Draw("same");
        func3->SetLineColor(3);
        TLegend *leg_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_dt->AddEntry(func1, "e^x/a fit", "l");
        // leg_dt->AddEntry(func2, "e^x/2200 + e^x/1800 fit", "l");
        leg_dt->AddEntry(func3, "e^x/a + e^x/b fit", "l");
        leg_dt->Draw();
        c_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_dt.png");
        c_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots11241/h_dt.png"));
        h_dt->Reset();

        h_intmuons->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmuons->GetYaxis()->SetTitle("Counts");
        h_intmuons->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_intmuons.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots11241/h_intmuons.png"));
        h_intmuons->Reset();

        h_intmichels->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels->GetYaxis()->SetTitle("Counts");
        h_intmichels->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_intmichels.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots11241/h_intmichels.png"));
        // h_intmichels->Print("all");
        h_intmichels->Reset();

        h_intcosmics->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intcosmics->GetYaxis()->SetTitle("Counts");
        h_intcosmics->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_intcosmics.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots11241/h_intcosmics.png"));
        // h_intcosmics->Print("all");
        h_intcosmics->Reset();

        TCanvas* c_34detint = new TCanvas("c_34detint", "Distribution of All tB34 Det. Integral Values", 1200, 700);
        c_34detint->SetLogy();
        c_34detint->cd();
        h_34detint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_34detint->GetYaxis()->SetTitle("Counts");
        h_34detint->Draw();
        c_34detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_34detint.png");
        // c_34detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_34detint.png");
        c_34detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_34detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_34detint.png"));
        h_34detint->Reset();

        TCanvas* c_34vpint = new TCanvas("c_34vpint", "Distribution of All tB34 VP Integral Values", 1200, 700);
        c_34vpint->SetLogy();
        c_34vpint->cd();
        h_34vpint->GetXaxis()->SetTitle("Integral (ADC)");
        h_34vpint->GetYaxis()->SetTitle("Counts");
        h_34vpint->Draw();
        c_34vpint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_34vpint.png");
        // c_34vpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_34vpint.png");
        c_34vpint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_34vpint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_34vpint.png"));
        h_34vpint->Reset();

        TCanvas* c_slice_vpintavg = new TCanvas("c_slice_vpintavg", "Distribution of Averaged Veto Panel Integral Values", 1200, 700);
        c_slice_vpintavg->SetLogy();
        c_slice_vpintavg->cd();
        h_slice_vpintavg->GetXaxis()->SetTitle("Integral (ADC)");
        h_slice_vpintavg->GetYaxis()->SetTitle("Counts");
        h_slice_vpintavg->Draw();
        c_slice_vpintavg->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_slice_vpintavg.png");
        // c_slice_vpintavg->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_slice_vpintavg.png");
        c_slice_vpintavg->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_slice_vpintavg.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_slice_vpintavg.png"));
        h_slice_vpintavg->Reset();

        TCanvas* c_dt_61_16 = new TCanvas("c_dt_61_16", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_16->SetLogy();
        c_dt_61_16->cd();
        h_dt_61_16->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_16->GetYaxis()->SetTitle("Counts");
        h_dt_61_16->Draw();
        c_dt_61_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_dt_61_16.png");
        // c_dt_61_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_dt_61_16.png");
        c_dt_61_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_16.png"));
        h_dt_61_16->Reset();

        TCanvas* c_dt_61_17 = new TCanvas("c_dt_61_17", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_17->SetLogy();
        c_dt_61_17->cd();
        h_dt_61_17->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_17->GetYaxis()->SetTitle("Counts");
        h_dt_61_17->Draw();
        c_dt_61_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_dt_61_17.png");
        // c_dt_61_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_dt_61_17.png");
        c_dt_61_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_17.png"));
        h_dt_61_17->Reset();

        TCanvas* c_dt_61_18 = new TCanvas("c_dt_61_18", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_18->SetLogy();
        c_dt_61_18->cd();
        h_dt_61_18->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_18->GetYaxis()->SetTitle("Counts");
        h_dt_61_18->Draw();
        c_dt_61_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_dt_61_18.png");
        // c_dt_61_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_dt_61_18.png");
        c_dt_61_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_18.png"));
        h_dt_61_18->Reset();

        TCanvas* c_dt_61_19 = new TCanvas("c_dt_61_19", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_19->SetLogy();
        c_dt_61_19->cd();
        h_dt_61_19->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_19->GetYaxis()->SetTitle("Counts");
        h_dt_61_19->Draw();
        c_dt_61_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_dt_61_19.png");
        // c_dt_61_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_dt_61_19.png");
        c_dt_61_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_19.png"));
        h_dt_61_19->Reset();

        TCanvas* c_dt_61_20 = new TCanvas("c_dt_61_20", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_20->SetLogy();
        c_dt_61_20->cd();
        h_dt_61_20->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_20->GetYaxis()->SetTitle("Counts");
        h_dt_61_20->Draw();
        c_dt_61_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_dt_61_20.png");
        // c_dt_61_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_dt_61_20.png");
        c_dt_61_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_20.png"));
        h_dt_61_20->Reset();

        TCanvas* c_dt_61_21 = new TCanvas("c_dt_61_21", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_21->SetLogy();
        c_dt_61_21->cd();
        h_dt_61_21->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_21->GetYaxis()->SetTitle("Counts");
        h_dt_61_21->Draw();
        c_dt_61_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_dt_61_21.png");
        // c_dt_61_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_dt_61_21.png");
        c_dt_61_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_21.png"));
        h_dt_61_21->Reset();

        TCanvas* c_dt_61_22 = new TCanvas("c_dt_61_22", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_22->SetLogy();
        c_dt_61_22->cd();
        h_dt_61_22->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_22->GetYaxis()->SetTitle("Counts");
        h_dt_61_22->Draw();
        c_dt_61_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_dt_61_22.png");
        // c_dt_61_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_dt_61_22.png");
        c_dt_61_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_22.png"));
        h_dt_61_22->Reset();

        TCanvas* c_dt_61_23 = new TCanvas("c_dt_61_23", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_23->SetLogy();
        c_dt_61_23->cd();
        h_dt_61_23->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_23->GetYaxis()->SetTitle("Counts");
        h_dt_61_23->Draw();
        c_dt_61_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_dt_61_23.png");
        // c_dt_61_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_dt_61_23.png");
        c_dt_61_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_23.png"));
        h_dt_61_23->Reset();

        TCanvas* c_dt_61_top = new TCanvas("c_dt_61_top", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_top->SetLogy();
        c_dt_61_top->cd();
        h_dt_61_top->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_top->GetYaxis()->SetTitle("Counts");
        h_dt_61_top->Draw();
        c_dt_61_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_dt_61_top.png");
        // c_dt_61_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_dt_61_top.png");
        c_dt_61_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_top.png"));
        h_dt_61_top->Reset();
        
        TCanvas* c_vpint_16 = new TCanvas("c_vpint_16", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_16->SetLogy();
        c_vpint_16->cd();
        h_vpint_16->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_16->GetYaxis()->SetTitle("Counts");
        h_vpint_16->Draw();
        c_vpint_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_vpint_16.png");
        // c_vpint_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_vpint_16.png");
        c_vpint_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_16.png"));
        h_vpint_16->Reset();

        TCanvas* c_vpint_17 = new TCanvas("c_vpint_17", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_17->SetLogy();
        c_vpint_17->cd();
        h_vpint_17->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_17->GetYaxis()->SetTitle("Counts");
        h_vpint_17->Draw();
        c_vpint_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_vpint_17.png");
        // c_vpint_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_vpint_17.png");
        c_vpint_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_17.png"));
        h_vpint_17->Reset();

        TCanvas* c_vpint_18 = new TCanvas("c_vpint_18", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_18->SetLogy();
        c_vpint_18->cd();
        h_vpint_18->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_18->GetYaxis()->SetTitle("Counts");
        h_vpint_18->Draw();
        c_vpint_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_vpint_18.png");
        // c_vpint_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_vpint_18.png");
        c_vpint_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_18.png"));
        h_vpint_18->Reset();

        TCanvas* c_vpint_19 = new TCanvas("c_vpint_19", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_19->SetLogy();
        c_vpint_19->cd();
        h_vpint_19->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_19->GetYaxis()->SetTitle("Counts");
        h_vpint_19->Draw();
        c_vpint_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_vpint_19.png");
        // c_vpint_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_vpint_19.png");
        c_vpint_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_19.png"));
        h_vpint_19->Reset();

        TCanvas* c_vpint_20 = new TCanvas("c_vpint_20", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_20->SetLogy();
        c_vpint_20->cd();
        h_vpint_20->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_20->GetYaxis()->SetTitle("Counts");
        h_vpint_20->Draw();
        c_vpint_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_vpint_20.png");
        // c_vpint_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_vpint_20.png");
        c_vpint_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_20.png"));
        h_vpint_20->Reset();

        TCanvas* c_vpint_21 = new TCanvas("c_vpint_21", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_21->SetLogy();
        c_vpint_21->cd();
        h_vpint_21->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_21->GetYaxis()->SetTitle("Counts");
        h_vpint_21->Draw();
        c_vpint_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_vpint_21.png");
        // c_vpint_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_vpint_21.png");
        c_vpint_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_21.png"));
        h_vpint_21->Reset();

        TCanvas* c_vpint_22 = new TCanvas("c_vpint_22", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_22->SetLogy();
        c_vpint_22->cd();
        h_vpint_22->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_22->GetYaxis()->SetTitle("Counts");
        h_vpint_22->Draw();
        c_vpint_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_vpint_22.png");
        // c_vpint_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_vpint_22.png");
        c_vpint_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_22.png"));
        h_vpint_22->Reset();

        TCanvas* c_vpint_23 = new TCanvas("c_vpint_23", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_23->SetLogy();
        c_vpint_23->cd();
        h_vpint_23->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_23->GetYaxis()->SetTitle("Counts");
        h_vpint_23->Draw();
        c_vpint_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_vpint_23.png");
        // c_vpint_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_vpint_23.png");
        c_vpint_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_23.png"));
        h_vpint_23->Reset();

        TCanvas* c_vpint_top = new TCanvas("c_vpint_top", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_top->SetLogy();
        c_vpint_top->cd();
        h_vpint_top->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_top->GetYaxis()->SetTitle("Counts");
        h_vpint_top->Draw();
        c_vpint_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_vpint_top.png");
        // c_vpint_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_vpint_top.png");
        c_vpint_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_top.png"));
        h_vpint_top->Reset();

        TCanvas* c_detint_vs_vpint_16 = new TCanvas("c_detint_vs_vpint_16", "SiPM 1 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_16->SetLogy();
        // c_detint_vs_vpint_16->SetLogx();
        c_detint_vs_vpint_16->cd();
        h_detint_vs_vpint_16->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_16->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_16->Draw();
        c_detint_vs_vpint_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_16.png");
        // c_detint_vs_vpint_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_16.png");
        c_detint_vs_vpint_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_16.png"));
        h_detint_vs_vpint_16->Reset();

        TCanvas* c_detint_vs_vpint_17 = new TCanvas("c_detint_vs_vpint_17", "SiPM 2 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_17->SetLogy();
        // c_detint_vs_vpint_17->SetLogx();
        c_detint_vs_vpint_17->cd();
        h_detint_vs_vpint_17->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_17->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_17->Draw();
        c_detint_vs_vpint_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_17.png");
        // c_detint_vs_vpint_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_17.png");
        c_detint_vs_vpint_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_17.png"));
        h_detint_vs_vpint_17->Reset();

        TCanvas* c_detint_vs_vpint_18 = new TCanvas("c_detint_vs_vpint_18", "SiPM 3 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_18->SetLogy();
        // c_detint_vs_vpint_18->SetLogx();
        c_detint_vs_vpint_18->cd();
        h_detint_vs_vpint_18->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_18->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_18->Draw();
        c_detint_vs_vpint_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_18.png");
        // c_detint_vs_vpint_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_18.png");
        c_detint_vs_vpint_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_18.png"));
        h_detint_vs_vpint_18->Reset();

        TCanvas* c_detint_vs_vpint_19 = new TCanvas("c_detint_vs_vpint_19", "SiPM 4 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_19->SetLogy();
        // c_detint_vs_vpint_19->SetLogx();
        c_detint_vs_vpint_19->cd();
        h_detint_vs_vpint_19->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_19->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_19->Draw();
        c_detint_vs_vpint_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_19.png");
        // c_detint_vs_vpint_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_19.png");
        c_detint_vs_vpint_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_19.png"));
        h_detint_vs_vpint_19->Reset();

        TCanvas* c_detint_vs_vpint_20 = new TCanvas("c_detint_vs_vpint_20", "SiPM 5 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_20->SetLogy();
        // c_detint_vs_vpint_20->SetLogx();
        c_detint_vs_vpint_20->cd();
        h_detint_vs_vpint_20->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_20->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_20->Draw();
        c_detint_vs_vpint_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_20.png");
        // c_detint_vs_vpint_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_20.png");
        c_detint_vs_vpint_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_20.png"));
        h_detint_vs_vpint_20->Reset();

        TCanvas* c_detint_vs_vpint_21 = new TCanvas("c_detint_vs_vpint_21", "SiPM 6 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_21->SetLogy();
        // c_detint_vs_vpint_21->SetLogx();
        c_detint_vs_vpint_21->cd();
        h_detint_vs_vpint_21->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_21->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_21->Draw();
        c_detint_vs_vpint_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_21.png");
        // c_detint_vs_vpint_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_21.png");
        c_detint_vs_vpint_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_21.png"));
        h_detint_vs_vpint_21->Reset();

        TCanvas* c_detint_vs_vpint_22 = new TCanvas("c_detint_vs_vpint_22", "SiPM 7 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_22->SetLogy();
        // c_detint_vs_vpint_22->SetLogx();
        c_detint_vs_vpint_22->cd();
        h_detint_vs_vpint_22->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_22->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_22->Draw();
        c_detint_vs_vpint_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_22.png");
        // c_detint_vs_vpint_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_22.png");
        c_detint_vs_vpint_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_22.png"));
        h_detint_vs_vpint_22->Reset();

        TCanvas* c_detint_vs_vpint_23 = new TCanvas("c_detint_vs_vpint_23", "SiPM 8 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_23->SetLogy();
        // c_detint_vs_vpint_23->SetLogx();
        c_detint_vs_vpint_23->cd();
        h_detint_vs_vpint_23->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_23->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_23->Draw();
        c_detint_vs_vpint_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_23.png");
        // c_detint_vs_vpint_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_23.png");
        c_detint_vs_vpint_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_23.png"));
        h_detint_vs_vpint_23->Reset();

        TCanvas* c_detint_vs_vpint_top = new TCanvas("c_detint_vs_vpint_top", "Top SiPM & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_top->SetLogy();
        // c_detint_vs_vpint_top->SetLogx();
        c_detint_vs_vpint_top->cd();
        h_detint_vs_vpint_top->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_top->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_top->Draw();
        c_detint_vs_vpint_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_top.png");
        // c_detint_vs_vpint_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11241/h_detint_vs_vpint_top.png");
        c_detint_vs_vpint_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_top.png"));
        h_detint_vs_vpint_top->Reset();

        TCanvas* c_tB = new TCanvas("c_tB", "TriggerBits Distribution", 1200, 700);
        c_tB->SetLogy();
        c_tB->cd();
        h_tB->GetXaxis()->SetTitle("triggerBits Value");
        h_tB->GetYaxis()->SetTitle("Counts");
        h_tB->Draw();
        c_tB->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_tB.png");
        c_tB->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_tB.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots11241/h_tB.png"));
        h_tB->Reset();

        TCanvas* c_blM_pmt = new TCanvas("c_blM_pmt", "PMT baselineMean Distribution", 1200, 700);
        c_blM_pmt->SetLogy();
        c_blM_pmt->cd();
        h_blM_pmt->GetXaxis()->SetTitle("baselineMean Value");
        h_blM_pmt->GetYaxis()->SetTitle("Counts");
        h_blM_pmt->Draw();
        c_blM_pmt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_blM_pmt.png");
        c_blM_pmt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_blM_pmt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots11241/h_blM_pmt.png"));
        h_blM_pmt->Reset();

        TCanvas* c_blM_vp = new TCanvas("c_blM_vp", "SiPM baselineMean Distribution", 1200, 700);
        c_blM_vp->SetLogy();
        c_blM_vp->cd();
        h_blM_vp->GetXaxis()->SetTitle("baselineMean Value");
        h_blM_vp->GetYaxis()->SetTitle("Counts");
        h_blM_vp->Draw();
        c_blM_vp->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_blM_vp.png");
        c_blM_vp->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_blM_vp.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots11241/h_blM_vp.png"));
        h_blM_vp->Reset();

        TCanvas* c_blM_ev61 = new TCanvas("c_blM_ev61", "Event 61 baselineMean Distribution", 1200, 700);
        c_blM_ev61->SetLogy();
        c_blM_ev61->cd();
        h_blM_ev61->GetXaxis()->SetTitle("baselineMean Value");
        h_blM_ev61->GetYaxis()->SetTitle("Counts");
        h_blM_ev61->Draw();
        c_blM_ev61->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11241/h_blM_ev61.png");
        c_blM_ev61->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_blM_ev61.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots11241/h_blM_ev61.png"));
        h_blM_ev61->Reset();

    }

    if (run == 12636) {

        TCanvas* c_dt = new TCanvas("c_dt", "Delta-T", 900, 700);
        // c_dt->SetLogy();
        c_dt->cd();

        TF1* func1 = new TF1("func1", "[0]+[1]*exp(-x/[2])", 2100, 12100);
        func1->SetParNames("Baseline", "Amplitude", "Time Constant");
        func1->SetParameters(100, 15000, 2100);
        h_dt->Fit("func1", "R");
        Double_t chi2_1 = func1->GetChisquare();

        cout << "\n" << "Chi-square of e^x/a fit: " << chi2_1 << "\n" << endl;

        // TF1* func2 = new TF1("func2", "[0]+[1]*exp(-x/2200)+[2]*exp(-x/1800)", 2100, 12100);
        // func2->SetParNames("Baseline", "Amplitude_u+", "Amplitude_u-");
        // func2->SetParameters(100, 15000, 15000);
        // h_dt->Fit("func2", "R");
        // Double_t chi2_2 = func2->GetChisquare();

        // cout << "\n" << "Chi-square of e^x/2200 + e^x/1800 fit: " << chi2_2 << "\n" << endl;

        TF1* func3 = new TF1("func3", "[0]+[1]*exp(-x/[2])+[3]*exp(-x/[4])", 2100, 12100);
        func3->SetParNames("Baseline", "Amplitude_u+", "Time Constant_u+", "Amplitude_u-", "Time Constant_u-");
        func3->SetParameters(100, 15000, 2200, 15000, 1800);
        h_dt->Fit("func3", "R");
        Double_t chi2_3 = func3->GetChisquare();

        cout << "\n" << "Chi-square of e^x/a + e^x/b fit: " << chi2_3 << endl;

        h_dt->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt->GetYaxis()->SetTitle("Counts");
        h_dt->GetYaxis()->SetRange(0, 200);
        h_dt->Draw();
        func1->Draw("same");
        func1->SetLineColor(2);
        // func2->Draw("same");
        // func2->SetLineColor(3);
        func3->Draw("same");
        func3->SetLineColor(3);
        TLegend *leg_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_dt->AddEntry(func1, "e^x/a fit", "l");
        // leg_dt->AddEntry(func2, "e^x/2200 + e^x/1800 fit", "l");
        leg_dt->AddEntry(func3, "e^x/a + e^x/b fit", "l");
        leg_dt->Draw();
        c_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt.png");
        c_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots12636/h_dt.png"));
        h_dt->Reset();

        h_intmuons->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmuons->GetYaxis()->SetTitle("Counts");
        h_intmuons->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_intmuons.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots12636/h_intmuons.png"));
        h_intmuons->Reset();

        h_intmichels->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intmichels->GetYaxis()->SetTitle("Counts");
        h_intmichels->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_intmichels.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots12636/h_intmichels.png"));
        // h_intmichels->Print("all");
        h_intmichels->Reset();

        h_intcosmics->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_intcosmics->GetYaxis()->SetTitle("Counts");
        h_intcosmics->Draw();
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_intcosmics.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots12636/h_intcosmics.png"));
        // h_intcosmics->Print("all");
        h_intcosmics->Reset();

        TCanvas* c_34detint = new TCanvas("c_34detint", "Distribution of All tB34 Det. Integral Values", 1200, 700);
        c_34detint->SetLogy();
        c_34detint->cd();
        h_34detint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_34detint->GetYaxis()->SetTitle("Counts");
        h_34detint->Draw();
        c_34detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_34detint.png");
        // c_34detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_34detint.png");
        c_34detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_34detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_34detint.png"));
        h_34detint->Reset();

        TCanvas* c_34vpint = new TCanvas("c_34vpint", "Distribution of All tB34 VP Integral Values", 1200, 700);
        c_34vpint->SetLogy();
        c_34vpint->cd();
        h_34vpint->GetXaxis()->SetTitle("Integral (ADC)");
        h_34vpint->GetYaxis()->SetTitle("Counts");
        h_34vpint->Draw();
        c_34vpint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_34vpint.png");
        // c_34vpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_34vpint.png");
        c_34vpint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_34vpint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_34vpint.png"));
        h_34vpint->Reset();

        TCanvas* c_slice_vpintavg = new TCanvas("c_slice_vpintavg", "Distribution of Averaged Veto Panel Integral Values", 1200, 700);
        c_slice_vpintavg->SetLogy();
        c_slice_vpintavg->cd();
        h_slice_vpintavg->GetXaxis()->SetTitle("Integral (ADC)");
        h_slice_vpintavg->GetYaxis()->SetTitle("Counts");
        h_slice_vpintavg->Draw();
        c_slice_vpintavg->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_slice_vpintavg.png");
        // c_slice_vpintavg->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_slice_vpintavg.png");
        c_slice_vpintavg->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_slice_vpintavg.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_slice_vpintavg.png"));
        h_slice_vpintavg->Reset();

        TCanvas* c_dt_61_16 = new TCanvas("c_dt_61_16", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_16->SetLogy();
        c_dt_61_16->cd();
        h_dt_61_16->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_16->GetYaxis()->SetTitle("Counts");
        h_dt_61_16->Draw();
        c_dt_61_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_16.png");
        // c_dt_61_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_dt_61_16.png");
        c_dt_61_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_16.png"));
        h_dt_61_16->Reset();

        TCanvas* c_dt_61_17 = new TCanvas("c_dt_61_17", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_17->SetLogy();
        c_dt_61_17->cd();
        h_dt_61_17->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_17->GetYaxis()->SetTitle("Counts");
        h_dt_61_17->Draw();
        c_dt_61_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_17.png");
        // c_dt_61_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_dt_61_17.png");
        c_dt_61_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_17.png"));
        h_dt_61_17->Reset();

        TCanvas* c_dt_61_18 = new TCanvas("c_dt_61_18", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_18->SetLogy();
        c_dt_61_18->cd();
        h_dt_61_18->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_18->GetYaxis()->SetTitle("Counts");
        h_dt_61_18->Draw();
        c_dt_61_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_18.png");
        // c_dt_61_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_dt_61_18.png");
        c_dt_61_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_18.png"));
        h_dt_61_18->Reset();

        TCanvas* c_dt_61_19 = new TCanvas("c_dt_61_19", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_19->SetLogy();
        c_dt_61_19->cd();
        h_dt_61_19->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_19->GetYaxis()->SetTitle("Counts");
        h_dt_61_19->Draw();
        c_dt_61_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_19.png");
        // c_dt_61_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_dt_61_19.png");
        c_dt_61_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_19.png"));
        h_dt_61_19->Reset();

        TCanvas* c_dt_61_20 = new TCanvas("c_dt_61_20", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_20->SetLogy();
        c_dt_61_20->cd();
        h_dt_61_20->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_20->GetYaxis()->SetTitle("Counts");
        h_dt_61_20->Draw();
        c_dt_61_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_20.png");
        // c_dt_61_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_dt_61_20.png");
        c_dt_61_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_20.png"));
        h_dt_61_20->Reset();

        TCanvas* c_dt_61_21 = new TCanvas("c_dt_61_21", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_21->SetLogy();
        c_dt_61_21->cd();
        h_dt_61_21->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_21->GetYaxis()->SetTitle("Counts");
        h_dt_61_21->Draw();
        c_dt_61_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_21.png");
        // c_dt_61_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_dt_61_21.png");
        c_dt_61_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_21.png"));
        h_dt_61_21->Reset();

        TCanvas* c_dt_61_22 = new TCanvas("c_dt_61_22", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_22->SetLogy();
        c_dt_61_22->cd();
        h_dt_61_22->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_22->GetYaxis()->SetTitle("Counts");
        h_dt_61_22->Draw();
        c_dt_61_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_22.png");
        // c_dt_61_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_dt_61_22.png");
        c_dt_61_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_22.png"));
        h_dt_61_22->Reset();

        TCanvas* c_dt_61_23 = new TCanvas("c_dt_61_23", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_23->SetLogy();
        c_dt_61_23->cd();
        h_dt_61_23->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_23->GetYaxis()->SetTitle("Counts");
        h_dt_61_23->Draw();
        c_dt_61_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_23.png");
        // c_dt_61_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_dt_61_23.png");
        c_dt_61_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_23.png"));
        h_dt_61_23->Reset();

        TCanvas* c_dt_61_top = new TCanvas("c_dt_61_top", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_top->SetLogy();
        c_dt_61_top->cd();
        h_dt_61_top->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_top->GetYaxis()->SetTitle("Counts");
        h_dt_61_top->Draw();
        c_dt_61_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_top.png");
        // c_dt_61_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_dt_61_top.png");
        c_dt_61_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_top.png"));
        h_dt_61_top->Reset();

        TCanvas* c_vpint_16 = new TCanvas("c_vpint_16", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_16->SetLogy();
        c_vpint_16->cd();
        h_vpint_16->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_16->GetYaxis()->SetTitle("Counts");
        h_vpint_16->Draw();
        c_vpint_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_vpint_16.png");
        // c_vpint_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_vpint_16.png");
        c_vpint_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_16.png"));
        h_vpint_16->Reset();

        TCanvas* c_vpint_17 = new TCanvas("c_vpint_17", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_17->SetLogy();
        c_vpint_17->cd();
        h_vpint_17->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_17->GetYaxis()->SetTitle("Counts");
        h_vpint_17->Draw();
        c_vpint_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_vpint_17.png");
        // c_vpint_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_vpint_17.png");
        c_vpint_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_17.png"));
        h_vpint_17->Reset();

        TCanvas* c_vpint_18 = new TCanvas("c_vpint_18", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_18->SetLogy();
        c_vpint_18->cd();
        h_vpint_18->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_18->GetYaxis()->SetTitle("Counts");
        h_vpint_18->Draw();
        c_vpint_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_vpint_18.png");
        // c_vpint_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_vpint_18.png");
        c_vpint_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_18.png"));
        h_vpint_18->Reset();

        TCanvas* c_vpint_19 = new TCanvas("c_vpint_19", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_19->SetLogy();
        c_vpint_19->cd();
        h_vpint_19->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_19->GetYaxis()->SetTitle("Counts");
        h_vpint_19->Draw();
        c_vpint_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_vpint_19.png");
        // c_vpint_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_vpint_19.png");
        c_vpint_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_19.png"));
        h_vpint_19->Reset();

        TCanvas* c_vpint_20 = new TCanvas("c_vpint_20", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_20->SetLogy();
        c_vpint_20->cd();
        h_vpint_20->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_20->GetYaxis()->SetTitle("Counts");
        h_vpint_20->Draw();
        c_vpint_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_vpint_20.png");
        // c_vpint_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_vpint_20.png");
        c_vpint_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_20.png"));
        h_vpint_20->Reset();

        TCanvas* c_vpint_21 = new TCanvas("c_vpint_21", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_21->SetLogy();
        c_vpint_21->cd();
        h_vpint_21->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_21->GetYaxis()->SetTitle("Counts");
        h_vpint_21->Draw();
        c_vpint_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_vpint_21.png");
        // c_vpint_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_vpint_21.png");
        c_vpint_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_21.png"));
        h_vpint_21->Reset();

        TCanvas* c_vpint_22 = new TCanvas("c_vpint_22", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_22->SetLogy();
        c_vpint_22->cd();
        h_vpint_22->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_22->GetYaxis()->SetTitle("Counts");
        h_vpint_22->Draw();
        c_vpint_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_vpint_22.png");
        // c_vpint_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_vpint_22.png");
        c_vpint_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_22.png"));
        h_vpint_22->Reset();

        TCanvas* c_vpint_23 = new TCanvas("c_vpint_23", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_23->SetLogy();
        c_vpint_23->cd();
        h_vpint_23->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_23->GetYaxis()->SetTitle("Counts");
        h_vpint_23->Draw();
        c_vpint_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_vpint_23.png");
        // c_vpint_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_vpint_23.png");
        c_vpint_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_23.png"));
        h_vpint_23->Reset();

        TCanvas* c_vpint_top = new TCanvas("c_vpint_top", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_top->SetLogy();
        c_vpint_top->cd();
        h_vpint_top->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_top->GetYaxis()->SetTitle("Counts");
        h_vpint_top->Draw();
        c_vpint_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_vpint_top.png");
        // c_vpint_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_vpint_top.png");
        c_vpint_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_top.png"));
        h_vpint_top->Reset();

        TCanvas* c_detint_vs_vpint_16 = new TCanvas("c_detint_vs_vpint_16", "SiPM 1 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_16->SetLogy();
        // c_detint_vs_vpint_16->SetLogx();
        c_detint_vs_vpint_16->cd();
        h_detint_vs_vpint_16->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_16->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_16->Draw();
        c_detint_vs_vpint_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_16.png");
        // c_detint_vs_vpint_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_16.png");
        c_detint_vs_vpint_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_16.png"));
        h_detint_vs_vpint_16->Reset();

        TCanvas* c_detint_vs_vpint_17 = new TCanvas("c_detint_vs_vpint_17", "SiPM 2 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_17->SetLogy();
        // c_detint_vs_vpint_17->SetLogx();
        c_detint_vs_vpint_17->cd();
        h_detint_vs_vpint_17->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_17->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_17->Draw();
        c_detint_vs_vpint_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_17.png");
        // c_detint_vs_vpint_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_17.png");
        c_detint_vs_vpint_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_17.png"));
        h_detint_vs_vpint_17->Reset();

        TCanvas* c_detint_vs_vpint_18 = new TCanvas("c_detint_vs_vpint_18", "SiPM 3 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_18->SetLogy();
        // c_detint_vs_vpint_18->SetLogx();
        c_detint_vs_vpint_18->cd();
        h_detint_vs_vpint_18->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_18->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_18->Draw();
        c_detint_vs_vpint_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_18.png");
        // c_detint_vs_vpint_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_18.png");
        c_detint_vs_vpint_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_18.png"));
        h_detint_vs_vpint_18->Reset();

        TCanvas* c_detint_vs_vpint_19 = new TCanvas("c_detint_vs_vpint_19", "SiPM 4 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_19->SetLogy();
        // c_detint_vs_vpint_19->SetLogx();
        c_detint_vs_vpint_19->cd();
        h_detint_vs_vpint_19->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_19->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_19->Draw();
        c_detint_vs_vpint_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_19.png");
        // c_detint_vs_vpint_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_19.png");
        c_detint_vs_vpint_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_19.png"));
        h_detint_vs_vpint_19->Reset();

        TCanvas* c_detint_vs_vpint_20 = new TCanvas("c_detint_vs_vpint_20", "SiPM 5 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_20->SetLogy();
        // c_detint_vs_vpint_20->SetLogx();
        c_detint_vs_vpint_20->cd();
        h_detint_vs_vpint_20->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_20->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_20->Draw();
        c_detint_vs_vpint_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_20.png");
        // c_detint_vs_vpint_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_20.png");
        c_detint_vs_vpint_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_20.png"));
        h_detint_vs_vpint_20->Reset();

        TCanvas* c_detint_vs_vpint_21 = new TCanvas("c_detint_vs_vpint_21", "SiPM 6 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_21->SetLogy();
        // c_detint_vs_vpint_21->SetLogx();
        c_detint_vs_vpint_21->cd();
        h_detint_vs_vpint_21->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_21->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_21->Draw();
        c_detint_vs_vpint_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_21.png");
        // c_detint_vs_vpint_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_21.png");
        c_detint_vs_vpint_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_21.png"));
        h_detint_vs_vpint_21->Reset();

        TCanvas* c_detint_vs_vpint_22 = new TCanvas("c_detint_vs_vpint_22", "SiPM 7 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_22->SetLogy();
        // c_detint_vs_vpint_22->SetLogx();
        c_detint_vs_vpint_22->cd();
        h_detint_vs_vpint_22->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_22->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_22->Draw();
        c_detint_vs_vpint_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_22.png");
        // c_detint_vs_vpint_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_22.png");
        c_detint_vs_vpint_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_22.png"));
        h_detint_vs_vpint_22->Reset();

        TCanvas* c_detint_vs_vpint_23 = new TCanvas("c_detint_vs_vpint_23", "SiPM 8 & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_23->SetLogy();
        // c_detint_vs_vpint_23->SetLogx();
        c_detint_vs_vpint_23->cd();
        h_detint_vs_vpint_23->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_23->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_23->Draw();
        c_detint_vs_vpint_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_23.png");
        // c_detint_vs_vpint_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_23.png");
        c_detint_vs_vpint_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_23.png"));
        h_detint_vs_vpint_23->Reset();

        TCanvas* c_detint_vs_vpint_top = new TCanvas("c_detint_vs_vpint_top", "Top SiPM & PMT Integral Values", 900, 700);
        c_detint_vs_vpint_top->SetLogy();
        // c_detint_vs_vpint_top->SetLogx();
        c_detint_vs_vpint_top->cd();
        h_detint_vs_vpint_top->GetXaxis()->SetTitle("PMT Integral (Ph.e.)");
        h_detint_vs_vpint_top->GetYaxis()->SetTitle("SiPM Integral (ADC)");
        h_detint_vs_vpint_top->Draw();
        c_detint_vs_vpint_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_top.png");
        // c_detint_vs_vpint_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_detint_vs_vpint_top.png");
        c_detint_vs_vpint_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_vs_vpint_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_vs_vpint_top.png"));
        h_detint_vs_vpint_top->Reset();

        TCanvas* c_tB = new TCanvas("c_tB", "TriggerBits Distribution", 1200, 700);
        c_tB->SetLogy();
        c_tB->cd();
        h_tB->GetXaxis()->SetTitle("triggerBits Value");
        h_tB->GetYaxis()->SetTitle("Counts");
        h_tB->Draw();
        c_tB->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_tB.png");
        c_tB->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_tB.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots12636/h_tB.png"));
        h_tB->Reset();

        TCanvas* c_blM_pmt = new TCanvas("c_blM_pmt", "PMT baselineMean Distribution", 1200, 700);
        c_blM_pmt->SetLogy();
        c_blM_pmt->cd();
        h_blM_pmt->GetXaxis()->SetTitle("baselineMean Value");
        h_blM_pmt->GetYaxis()->SetTitle("Counts");
        h_blM_pmt->Draw();
        c_blM_pmt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_blM_pmt.png");
        c_blM_pmt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_blM_pmt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots12636/h_blM_pmt.png"));
        h_blM_pmt->Reset();

        TCanvas* c_blM_vp = new TCanvas("c_blM_vp", "SiPM baselineMean Distribution", 1200, 700);
        c_blM_vp->SetLogy();
        c_blM_vp->cd();
        h_blM_vp->GetXaxis()->SetTitle("baselineMean Value");
        h_blM_vp->GetYaxis()->SetTitle("Counts");
        h_blM_vp->Draw();
        c_blM_vp->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_blM_vp.png");
        c_blM_vp->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_blM_vp.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots12636/h_blM_vp.png"));
        h_blM_vp->Reset();

        TCanvas* c_blM_ev61 = new TCanvas("c_blM_ev61", "Event 61 baselineMean Distribution", 1200, 700);
        c_blM_ev61->SetLogy();
        c_blM_ev61->cd();
        h_blM_ev61->GetXaxis()->SetTitle("baselineMean Value");
        h_blM_ev61->GetYaxis()->SetTitle("Counts");
        h_blM_ev61->Draw();
        c_blM_ev61->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_blM_ev61.png");
        c_blM_ev61->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_blM_ev61.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots12636/h_blM_ev61.png"));
        h_blM_ev61->Reset();

    }

    TCanvas* c_time_bt_ev61 = new TCanvas("c_time_bt_ev61", "Time Between Event 61 Peaks", 1200, 700);
    c_time_bt_ev61->SetLogy();
    c_time_bt_ev61->cd();
    h_time_bt_ev61->GetXaxis()->SetTitle("Delta-T (ns)");
    h_time_bt_ev61->GetYaxis()->SetTitle("Counts");
    h_time_bt_ev61->Draw();
    c_time_bt_ev61->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_time_bt_ev61.png");
    c_time_bt_ev61->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_time_bt_ev61.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_time_bt_ev61.png"));
    h_time_bt_ev61->Reset();

    TCanvas* c_time_bt_minLED = new TCanvas("c_time_bt_minLED", "Time Between Min Bias Peaks", 1200, 700);
    c_time_bt_minLED->SetLogy();
    c_time_bt_minLED->cd();
    h_time_bt_minLED->GetXaxis()->SetTitle("Delta-T (ns)");
    h_time_bt_minLED->GetYaxis()->SetTitle("Counts");
    h_time_bt_minLED->Draw();
    c_time_bt_minLED->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_time_bt_minLED.png");
    c_time_bt_minLED->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_time_bt_minLED.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_time_bt_minLED.png"));
    h_time_bt_minLED->Reset();

    TCanvas* c_time_bt_lowLED = new TCanvas("c_time_bt_lowLED", "Time Between Low Light Peaks", 1200, 700);
    c_time_bt_lowLED->SetLogy();
    c_time_bt_lowLED->cd();
    h_time_bt_lowLED->GetXaxis()->SetTitle("Delta-T (ns)");
    h_time_bt_lowLED->GetYaxis()->SetTitle("Counts");
    h_time_bt_lowLED->Draw();
    c_time_bt_lowLED->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_time_bt_lowLED.png");
    c_time_bt_lowLED->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_time_bt_lowLED.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_time_bt_lowLED.png"));
    h_time_bt_lowLED->Reset();

    TCanvas* c_time_bt_highLED = new TCanvas("c_time_bt_highLED", "Time Between High Light Peaks", 1200, 700);
    c_time_bt_highLED->SetLogy();
    c_time_bt_highLED->cd();
    h_time_bt_highLED->GetXaxis()->SetTitle("Delta-T (ns)");
    h_time_bt_highLED->GetYaxis()->SetTitle("Counts");
    h_time_bt_highLED->Draw();
    c_time_bt_highLED->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_time_bt_highLED.png");
    c_time_bt_highLED->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_time_bt_highLED.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_time_bt_highLED.png"));
    h_time_bt_highLED->Reset();

    TCanvas* c_intall = new TCanvas("c_intall", "All Integral Values", 1200, 700);
    c_intall->SetLogy();
    c_intall->cd();
    h_intall->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intall->GetYaxis()->SetTitle("Counts");
    h_intall->Draw();
    TLatex *lat_intall = new TLatex(0.75, 0.7, Form("RMS = %.2f", rmsValue(intall)));
    lat_intall->SetNDC();
    lat_intall->SetTextColor(1);
    lat_intall->SetTextSize(0.035);
    lat_intall->Draw();
    c_intall->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intall.png");
    c_intall->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_intall.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intall.png"));
    h_intall->Reset();

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
    c_intlowLED->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intlowLED_%i.png", run));
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
    c_inthighLED->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_inthighLED_%i.png", run));
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
    c_intminLED->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intminLED_%i.png", run));
    c_intminLED->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_intminLED.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intminLED.png"));
    h_intminLED->Reset();

    TCanvas* c_MB_HL_dt = new TCanvas("c_MB_HL_dt", "MB-HL dt Values", 1200, 700);
    c_MB_HL_dt->SetLogy();
    c_MB_HL_dt->cd();
    h_MB_HL_dt->GetXaxis()->SetTitle("Delta-T (ns)");
    h_MB_HL_dt->GetYaxis()->SetTitle("Counts");
    h_MB_HL_dt->Draw();
    c_MB_HL_dt->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_MB_HL_dt_%i.png", run));
    c_MB_HL_dt->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_MB_HL_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_MB_HL_dt.png"));
    h_MB_HL_dt->Reset();

    TCanvas* c_intallint = new TCanvas("c_intallint", "All Internally Triggered Integral Values", 1200, 700);
    c_intallint->SetLogy();
    c_intallint->cd();
    h_intallint->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intallint->GetYaxis()->SetTitle("Counts");
    h_intallint->Draw();
    TLatex *lat_intallint = new TLatex(0.75, 0.7, Form("RMS = %.2f", rmsValue(intallint)));
    lat_intallint->SetNDC();
    lat_intallint->SetTextColor(1);
    lat_intallint->SetTextSize(0.035);
    lat_intallint->Draw();
    c_intallint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_intallint.png");
    c_intallint->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_intallint.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_intallint.png"));
    h_intallint->Reset();

    TCanvas* c_intallmultcut = new TCanvas("c_intallmultcut", "All Internally Triggered Integral Values, Multiplicity Cut", 1200, 700);
    c_intallmultcut->SetLogy();
    c_intallmultcut->cd();
    h_intallmultcut->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intallmultcut->GetYaxis()->SetTitle("Counts");
    h_intallmultcut->Draw();
    TLatex *lat_intallmultcut = new TLatex(0.75, 0.7, Form("RMS = %.2f", rmsValue(intallmultcut)));
    lat_intallmultcut->SetNDC();
    lat_intallmultcut->SetTextColor(1);
    lat_intallmultcut->SetTextSize(0.035);
    lat_intallmultcut->Draw();
    c_intallmultcut->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_intallmultcut.png");
    c_intallmultcut->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_intallmultcut.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_intallmultcut.png"));
    h_intallmultcut->Reset();

    TCanvas* c_intallmultmuoncut = new TCanvas("c_intallmultmuoncut", "All Internally Triggered Integral Values, Multiplicity & Muon Cut", 1200, 700);
    c_intallmultmuoncut->SetLogy();
    c_intallmultmuoncut->cd();
    h_intallmultmuoncut->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intallmultmuoncut->GetYaxis()->SetTitle("Counts");
    h_intallmultmuoncut->Draw();
    TLatex *lat_intallmultmuoncut = new TLatex(0.75, 0.7, Form("RMS = %.2f", rmsValue(intallmultmuoncut)));
    lat_intallmultmuoncut->SetNDC();
    lat_intallmultmuoncut->SetTextColor(1);
    lat_intallmultmuoncut->SetTextSize(0.035);
    lat_intallmultmuoncut->Draw();
    c_intallmultmuoncut->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_intallmultmuoncut.png");
    c_intallmultmuoncut->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_intallmultmuoncut.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_intallmultmuoncut.png"));
    h_intallmultmuoncut->Reset();

    auto p_vpeff = new TGraph();
    p_vpeff->SetTitle("Efficiency of Top Veto Panel as a Function of Top Veto Panel Energy Cuts");
    p_vpeff->GetXaxis()->SetTitle("Top Veto Panel Integral Threshold (ADC)");
    p_vpeff->GetYaxis()->SetTitle("Efficiency");
    p_vpeff->AddPoint(400, 0.9997);         // 7654/7656
    p_vpeff->AddPoint(500, 0.9990);         // 7648/7656
    p_vpeff->AddPoint(600, 0.9918);         // 7593/7656
    p_vpeff->AddPoint(700, 0.9653);         // 7390/7656
    p_vpeff->AddPoint(800, 0.9070);         // 6944/7656
    p_vpeff->Draw("AL*");
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/p_vpeff.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/p_vpeff.png"));

    TCanvas* c_amp = new TCanvas("c_amp", "Distribution of All Event Amplitude Values", 1200, 700);
    c_amp->SetLogy();
    c_amp->cd();
    h_amp->GetXaxis()->SetTitle("Amplitude (ADC)");
    h_amp->GetYaxis()->SetTitle("Counts");
    h_amp->Draw();
    c_amp->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_amp.png");
    c_amp->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_amp.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_amp.png"));
    h_amp->Reset();

    TCanvas* c_int = new TCanvas("c_int", "Distribution of All Event Integral Values", 1200, 700);
    c_int->SetLogy();
    c_int->cd();
    h_int->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_int->GetYaxis()->SetTitle("Counts");
    h_int->Draw();
    c_int->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_int.png");
    c_int->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int.png"));
    h_int->Reset();

    h_no_prior_muons_dt->GetXaxis()->SetTitle("Delta-T (ns)");
    h_no_prior_muons_dt->GetYaxis()->SetTitle("Counts");
    h_no_prior_muons_dt->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_no_prior_muons_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_no_prior_muons_dt.png"));
    h_no_prior_muons_dt->Reset();

    TCanvas* c_dt_61_smc = new TCanvas("c_dt_61_smc", "Delta-T", 900, 700);
    // c_dt_61_smc->SetLogy();
    c_dt_61_smc->cd();

    // TF1* func_61 = new TF1("func_61", "[0]+[1]*exp(-x/[2])", 2100, 12100);
    // func_61->SetParNames("Baseline", "Amplitude", "Time Constant");
    // func_61->SetParameters(50, 300, 2100);
    // h_dt_61->Fit("func_61", "R");

    h_dt_61_smc->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_smc->GetYaxis()->SetTitle("Counts");
    // h_dt_61_smc->GetYaxis()->SetRange(0, 200);

    // h_dt_61->Draw();
    // h_dt_61->SetLineColor(kGreen);
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61.png"));

    h_dt_61_smc->Draw("same");
    h_dt_61_smc->SetLineColor(kRed);
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_smc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_smc.png"));

    h_dt_61_spmc->Draw("same");
    h_dt_61_spmc->SetLineColor(kViolet);
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_spmc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_spmc.png"));

    h_dt_61_pmc->Draw("same");
    h_dt_61_pmc->SetLineColor(kBlue);
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_pmc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_pmc.png"));

    // func_61->Draw("same");
    
    TLegend *leg_dt_61 = new TLegend(0.9, 0.65, 0.65, 0.45);
    // leg_dt_61->AddEntry(h_dt_61_nmc, "No Muon Cuts", "l");
    leg_dt_61->AddEntry(h_dt_61_smc, "With Simultaneous Muon Cut", "l");
    leg_dt_61->AddEntry(h_dt_61_pmc, "With Prior Muon Cut", "l");
    leg_dt_61->AddEntry(h_dt_61_spmc, "With Simult. + Prior Muon Cuts", "l");
    leg_dt_61->Draw();
    
    c_dt_61_smc->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_61.png");
    c_dt_61_smc->Close();

    h_dt_61->Reset();
    h_dt_61_vc->Reset();

    TCanvas* c_int_61 = new TCanvas("c_int_61", "Delta-T", 900, 700);
    // c_int_61->SetLogy();
    c_int_61->cd();

    h_int_61->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_int_61->GetYaxis()->SetTitle("Counts");
    
    h_int_61->Draw();
    h_int_61->SetLineColor(kGreen);
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_61.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_61.png"));

    h_int_61_smc->Draw("same");
    h_int_61_smc->SetLineColor(kRed);
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_61_vc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_61_vc.png"));
    
    TLegend *leg_int_61 = new TLegend(0.65, 0.55, 0.85, 0.75);
    leg_int_61->AddEntry(h_int_61, "No Muon Cut", "l");
    leg_int_61->AddEntry(h_int_61_smc, "Simult. Muon Cut", "l");
    leg_int_61->Draw();

    c_int_61->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_int_61.png");
    c_int_61->Close();

    h_int_61->Reset();
    h_int_61_vc->Reset();

    h_dt_61_smc->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_smc->GetYaxis()->SetTitle("Counts");
    h_dt_61_smc->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_61_smc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_smc.png"));
    h_dt_61_smc->Reset();

    h_int_61_smc->GetXaxis()->SetTitle("Integral (ADC)");
    h_int_61_smc->GetYaxis()->SetTitle("Counts");
    h_int_61_smc->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_int_61_smc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_61_smc.png"));
    h_int_61_smc->Reset();

    h_dt_61_spmc->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_spmc->GetYaxis()->SetTitle("Counts");
    h_dt_61_spmc->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_61_spmc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_spmc.png"));
    h_dt_61_spmc->Reset();

    h_int_61_spmc->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_int_61_spmc->GetYaxis()->SetTitle("Counts");
    h_int_61_spmc->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_int_61_spmc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_61_spmc.png"));
    h_int_61_spmc->Reset();

    h_dt_61_pmc->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_pmc->GetYaxis()->SetTitle("Counts");
    h_dt_61_pmc->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_61_pmc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_pmc.png"));
    h_dt_61_pmc->Reset();

    h_int_61_pmc->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_int_61_pmc->GetYaxis()->SetTitle("Counts");
    h_int_61_pmc->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_int_61_pmc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_61_pmc.png"));
    h_int_61_pmc->Reset();

    TCanvas* c_int_61_spmc_outer = new TCanvas("c_int_61_spmc_outer", "Inner vs Outer Detector Integral Distribution, No Muons Before or During", 1200, 700);
    // c_int_61_spmc_outer->SetLogy();
    c_int_61_spmc_outer->cd();
    h_int_61_spmc_outer->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_int_61_spmc_outer->GetYaxis()->SetTitle("Counts");
    h_int_61_spmc_outer->Draw();
    h_int_61_spmc_outer->SetLineColor(kRed);
    h_int_61_spmc_inner->Draw("same");
    h_int_61_spmc_inner->SetLineColor(kGreen);
    TLegend *leg_int_61_spmc_outer = new TLegend(0.65, 0.55, 0.85, 0.75);
    leg_int_61_spmc_outer->AddEntry(h_int_61_spmc_inner, "dt < |2 us|", "l");
    leg_int_61_spmc_outer->AddEntry(h_int_61_spmc_outer, "dt > |2 us|", "l");
    leg_int_61_spmc_outer->Draw();
    c_int_61_spmc_outer->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_int_61_spmc_inner_vs_outer.png");
    c_int_61_spmc_outer->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_61_spmc_inner_vs_outer.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_61_spmc_inner_vs_outer.png"));
    h_int_61_spmc_outer->Reset();

    TCanvas* c_dt_v_int_61 = new TCanvas("c_dt_v_int_61", "Event 61 Detector dt vs Integral Distribution", 1200, 700);
    // c_dt_v_int_61->SetLogy();
    c_dt_v_int_61->cd();
    h_dt_v_int_61->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_v_int_61->GetYaxis()->SetTitle("Integral (Ph.e.)");
    h_dt_v_int_61->SetMarkerStyle(7);
    h_dt_v_int_61->Draw();
    c_dt_v_int_61->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_v_int_61.png");
    c_dt_v_int_61->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_int_61.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_int_61.png"));
    h_dt_v_int_61->Reset();

    TCanvas* c_dt_v_int_61_smc = new TCanvas("c_dt_v_int_61_smc", "Event 61 Veto Panel dt vs Integral Distribution", 1200, 700);
    // c_dt_v_int_61_smc->SetLogy();
    c_dt_v_int_61_smc->cd();
    h_dt_v_int_61_smc->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_v_int_61_smc->GetYaxis()->SetTitle("Integral (Ph.e.)");
    h_dt_v_int_61_smc->SetMarkerStyle(7);
    h_dt_v_int_61_smc->Draw();
    c_dt_v_int_61_smc->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_v_int_61_smc.png");
    c_dt_v_int_61_smc->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_int_61_smc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_int_61_smc.png"));
    h_dt_v_int_61_smc->Reset();

    TCanvas* c_dt_v_int_61_spmc = new TCanvas("c_dt_v_int_61_spmc", "Event 61 Detector dt vs Integral Distribution, No Muons", 1200, 700);
    // c_dt_v_int_61_spmc->SetLogy();
    c_dt_v_int_61_spmc->cd();
    h_dt_v_int_61_spmc->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_v_int_61_spmc->GetYaxis()->SetTitle("Integral (Ph.e.)");
    h_dt_v_int_61_spmc->SetMarkerStyle(7);
    h_dt_v_int_61_spmc->Draw();
    c_dt_v_int_61_spmc->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_v_int_61_spmc.png");
    c_dt_v_int_61_spmc->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_int_61_spmc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_int_61_spmc.png"));
    h_dt_v_int_61_spmc->Reset();

    TCanvas* c_dt_v_int_61_pmc = new TCanvas("c_dt_v_int_61_pmc", "Event 61 Detector dt vs Integral Distribution, No Muons", 1200, 700);
    // c_dt_v_int_61_pmc->SetLogy();
    c_dt_v_int_61_pmc->cd();
    h_dt_v_int_61_pmc->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_v_int_61_pmc->GetYaxis()->SetTitle("Integral (Ph.e.)");
    h_dt_v_int_61_pmc->SetMarkerStyle(7);
    h_dt_v_int_61_pmc->Draw();
    c_dt_v_int_61_pmc->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_v_int_61_pmc.png");
    c_dt_v_int_61_pmc->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_int_61_pmc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_int_61_pmc.png"));
    h_dt_v_int_61_pmc->Reset();

    h_neg_dt_ev61_peak_pos->GetXaxis()->SetTitle("Event 61 Peak Position (ns)");
    h_neg_dt_ev61_peak_pos->GetYaxis()->SetTitle("Counts");
    h_neg_dt_ev61_peak_pos->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_neg_dt_ev61_peak_pos.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_neg_dt_ev61_peak_pos.png"));
    h_neg_dt_ev61_peak_pos->Reset();

    TCanvas* c_tvp_int = new TCanvas("c_tvp_int", "Distribution of Top Veto Panel Integral Values", 1200, 700);
    c_tvp_int->SetLogy();
    c_tvp_int->cd();
    h_tvp_int->GetXaxis()->SetTitle("Integral (ADC)");
    h_tvp_int->GetYaxis()->SetTitle("Counts");
    h_tvp_int->Draw();
    c_tvp_int->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_tvp_int.png");
    c_tvp_int->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_tvp_int.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_tvp_int.png"));
    h_tvp_int->Reset();

    /*TCanvas* c_intneu = new TCanvas("c_intneu", "Distribution of Beam On Event Integral Values", 1200, 700);
    c_intneu->SetLogy();
    c_intneu->cd();
    h_intneu->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intneu->GetYaxis()->SetTitle("Counts");
    h_intneu->Draw();
    c_intneu->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_intneu.png");
    c_intneu->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_intneu.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_intneu.png"));
    h_intneu->Reset();*/

    TCanvas* c_spheLL = new TCanvas("c_spheLL", "Single Channel Low Light Integral Values", 1200, 700);
    c_spheLL->SetLogy();
    c_spheLL->cd();
    h_spheLL->GetXaxis()->SetTitle("Integral (ADC)");
    h_spheLL->GetYaxis()->SetTitle("Counts");
    h_spheLL->Draw();
    c_spheLL->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_spheLL.png");
    c_spheLL->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_spheLL.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_spheLL.png"));
    h_spheLL->Reset();

    TCanvas* c_spheMB = new TCanvas("c_spheMB", "Single Channel Min Bias Integral Values", 1200, 700);
    c_spheMB->SetLogy();
    c_spheMB->cd();
    h_spheMB->GetXaxis()->SetTitle("Integral (ADC)");
    h_spheMB->GetYaxis()->SetTitle("Counts");
    h_spheMB->Draw();
    c_spheMB->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_spheMB.png");
    c_spheMB->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_spheMB.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_spheMB.png"));
    h_spheMB->Reset();

    TCanvas* c_ampOnOff = new TCanvas("c_ampOnOff", "Peak Amplitudes with SNS Beam On vs Off", 1200, 700);
    // c_ampOnOff->SetLogy();
    c_ampOnOff->cd();

    for (int i = 0; i < amp_on_off_val.size(); i++) {

        for (int j = 0; j < amp_on_off_val[i].size(); j++) {

            if (i == 0) {

                h_ampOn->Fill(amp_on_off_val[i][j]);


            }

            else if (i == 1) {

                h_ampOff->Fill(amp_on_off_val[i][j]);

            }

        }

        if (i == 0) {

            h_ampOn->Draw();

            h_ampOn->SetLineColor(kRed);


        }

        else if (i == 1) {

            h_ampOff->Draw("same");

            h_ampOff->SetLineColor(kGreen);

        }

    }

    h_ampOn->GetXaxis()->SetTitle("Amplitude (Ph.e.)");
    h_ampOn->GetYaxis()->SetTitle("Counts");
    c_ampOnOff->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_ampOnOff.png");
    c_ampOnOff->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_ampOnOff.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_ampOnOff.png"));
    h_ampOn->Reset();
    h_ampOff->Reset();

    h_ampP1->GetXaxis()->SetTitle("Amplitude (Ph.e.)");
    h_ampP1->GetYaxis()->SetTitle("Counts");
    h_ampP1->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_ampP1.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_ampP1.png"));
    h_ampP1->Reset();

    h_ampP2->GetXaxis()->SetTitle("Amplitude (Ph.e.)");
    h_ampP2->GetYaxis()->SetTitle("Counts");
    h_ampP2->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_ampP2.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_ampP2.png"));
    h_ampP2->Reset();

    h_var->GetXaxis()->SetTitle("Variance (ns)");
    h_var->GetYaxis()->SetTitle("Counts");
    h_var->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_var.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_var.png"));
    h_var->Reset();

    TCanvas* c_pst = new TCanvas("c_pst", "Peak Start Time Distribution", 1200, 700);
    c_pst->SetLogy();
    c_pst->cd();
    h_pst->GetXaxis()->SetTitle("Pulse Start Time (ns)");
    h_pst->GetYaxis()->SetTitle("Counts");
    h_pst->Draw();
    c_pst->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_pst.png");
    c_pst->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_pst.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_pst.png"));
    h_pst->Reset();

    TCanvas* c_mfpst = new TCanvas("c_mfpst", "Most Frequent Peak Start Time Distribution", 1200, 700);
    c_mfpst->SetLogy();
    c_mfpst->cd();
    h_mfpst->GetXaxis()->SetTitle("'Averaged' Pulse Start Time (ns)");
    h_mfpst->GetYaxis()->SetTitle("Counts");
    h_mfpst->Draw();
    c_mfpst->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_mfpst.png");
    c_mfpst->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_mfpst.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_mfpst.png"));
    h_mfpst->Reset();

    h_pst_v_pint->GetXaxis()->SetTitle("Pulse Start Time (ns)");
    h_pst_v_pint->GetYaxis()->SetTitle("Pulse Integral Value (Ph.e.)");
    h_pst_v_pint->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_pst_v_pint.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_pst_v_pint.png"));
    h_pst_v_pint->Reset();

    h_mfpst_v_mfpint->GetXaxis()->SetTitle("Pulse Start Time (ns)");
    h_mfpst_v_mfpint->GetYaxis()->SetTitle("Pulse Integral Value (Ph.e.)");
    h_mfpst_v_mfpint->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_mfpst_v_mfpint.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_mfpst_v_mfpint.png"));
    h_mfpst_v_mfpint->Reset();

    TCanvas* c_ml = new TCanvas("c_ml", "Multiplicity Distribution", 1200, 700);
    c_ml->SetLogy();
    c_ml->cd();
    h_ml->GetXaxis()->SetTitle("Number of Channels a Peak was Found in");
    h_ml->GetYaxis()->SetTitle("Counts");
    h_ml->Draw();
    c_ml->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_ml.png");
    c_ml->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_ml.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_ml.png"));
    h_ml->Reset();

    TCanvas* c_ml_v_int = new TCanvas("c_ml_v_int", "Multiplicity vs Integral Distribution", 1200, 700);
    // c_ml_v_int->SetLogy();
    c_ml_v_int->cd();
    h_ml_v_int->GetXaxis()->SetTitle("Number of Channels a Peak was Found in");
    h_ml_v_int->GetYaxis()->SetTitle("Integral (Ph.e.)");
    h_ml_v_int->Draw();
    c_ml_v_int->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_ml_v_int.png");
    c_ml_v_int->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_ml_v_int.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots7894/h_ml_v_int.png"));
    h_ml_v_int->Reset();

    h_dtint->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dtint->GetYaxis()->SetTitle("Counts");
    h_dtint->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dtint.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dtint.png"));
    h_dtint->Reset();

    h_dtmin->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dtmin->GetYaxis()->SetTitle("Counts");
    h_dtmin->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dtmin.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dtmin.png"));
    h_dtmin->Reset();

    h_dthlt->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dthlt->GetYaxis()->SetTitle("Counts");
    h_dthlt->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dthlt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dthlt.png"));
    h_dthlt->Reset();

    h_dtllt->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dtllt->GetYaxis()->SetTitle("Counts");
    h_dtllt->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dtllt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dtllt.png"));
    h_dtllt->Reset();

    h_ampdiff->GetXaxis()->SetTitle("Peak Amplitude Difference");
    h_ampdiff->GetYaxis()->SetTitle("Counts");
    h_ampdiff->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_ampdiff.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_ampdiff.png"));
    h_ampdiff->Reset();

    h_ampboth->GetXaxis()->SetTitle("First Peak Amplitude");
    h_ampboth->GetYaxis()->SetTitle("Second Peak Amplitude");
    h_ampboth->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_ampboth.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_ampboth.png"));
    h_ampboth->Reset();

    h_intboth->GetXaxis()->SetTitle("First Peak Integral Values (Ph.e.)");
    h_intboth->GetYaxis()->SetTitle("Second Peak Integral Values (Ph.e.)");
    h_intboth->Draw("col");
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_intboth.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_intboth.png"));
    h_intboth->Reset();

    TCanvas* c_cent_spke_tB_smc = new TCanvas("c_cent_spke_tB_smc", "Central Spike TriggerBits Distribution, Simult. Muon Cuts", 1200, 700);
    c_cent_spke_tB_smc->SetLogy();
    c_cent_spke_tB_smc->cd();
    h_cent_spke_tB_smc->GetXaxis()->SetTitle("triggerBits Value");
    h_cent_spke_tB_smc->GetYaxis()->SetTitle("Counts");
    h_cent_spke_tB_smc->Draw();
    c_cent_spke_tB_smc->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_cent_spke_tB_smc.png");
    c_cent_spke_tB_smc->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_cent_spke_tB_smc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_cent_spke_tB_smc.png"));
    h_cent_spke_tB_smc->Reset();

    TCanvas* c_cent_spke_tB_pmc = new TCanvas("c_cent_spke_tB_pmc", "Central Spike TriggerBits Distribution, Prior Muon Cuts", 1200, 700);
    c_cent_spke_tB_pmc->SetLogy();
    c_cent_spke_tB_pmc->cd();
    h_cent_spke_tB_pmc->GetXaxis()->SetTitle("triggerBits Value");
    h_cent_spke_tB_pmc->GetYaxis()->SetTitle("Counts");
    h_cent_spke_tB_pmc->Draw();
    c_cent_spke_tB_pmc->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_cent_spke_tB_pmc.png");
    c_cent_spke_tB_pmc->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_cent_spke_tB_pmc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_cent_spke_tB_pmc.png"));
    h_cent_spke_tB_pmc->Reset();

    TCanvas* c_cent_spke_tB_spmc = new TCanvas("c_cent_spke_tB_spmc", "Central Spike TriggerBits Distribution, Prior + Simult. Muon Cuts", 1200, 700);
    c_cent_spke_tB_spmc->SetLogy();
    c_cent_spke_tB_spmc->cd();
    h_cent_spke_tB_spmc->GetXaxis()->SetTitle("triggerBits Value");
    h_cent_spke_tB_spmc->GetYaxis()->SetTitle("Counts");
    h_cent_spke_tB_spmc->Draw();
    c_cent_spke_tB_spmc->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_cent_spke_tB_spmc.png");
    c_cent_spke_tB_spmc->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_cent_spke_tB_spmc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_cent_spke_tB_spmc.png"));
    h_cent_spke_tB_spmc->Reset();

    // cout << "\n" << "Size of pulse_time_sep: " << pulse_time_sep.size() << endl;

    // cout << "\n" << "Size of event_num: " << event_num.size() << endl;

    // TH2D* h_dten = new TH2D("h_dten", "Delta-T vs Event Number", 1000, 0, pow(10, 7), 1000, 0, MAX_NUM_ENTRIES);

    // for (int i = 0; i < pulse_time_sep.size(); i++) {
    //     h_dten->Fill(pulse_time_sep[i], event_num[i]);
        // cout << "\n" << "dt value = " << pulse_time_sep[i] << endl;
        // cout << "\n" << "P2 trigger value = " << p2_trig_type[i] << endl;
    // }
    h_dten->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dten->GetYaxis()->SetTitle("Event Number");
    h_dten->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dten.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dten.png"));
    h_dten->Reset();
    
    h_tc->GetXaxis()->SetTitle("Delta-T (ns)");
    h_tc->GetYaxis()->SetTitle("Counts");
    h_tc->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_tc.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_tc.png"));
    h_tc->Reset();

    double veto_det_ratio_1 = num_veto_found_1 * pow(num_large_found_1, -1);

    double veto_det_ratio_2 = num_veto_found_2 * pow(num_large_found_2, -1);

    double veto_det_ratio_3 = num_veto_found_3 * pow(num_large_found_3, -1);

    double veto_det_ratio_4 = num_veto_found_4 * pow(num_large_found_4, -1);

    double veto_det_ratio_5 = num_veto_found_5 * pow(num_large_found_5, -1);

    double veto_det_ratio_6 = num_veto_found_6 * pow(num_large_found_6, -1);

    // double veto_det_ratio_7 = num_veto_found_7 * pow(num_large_found_7, -1);

    // double veto_det_ratio_8 = num_veto_found_8 * pow(num_large_found_8, -1);

    // double veto_det_ratio_9 = num_veto_found_9 * pow(num_large_found_9, -1);

    // double veto_det_ratio_10 = num_veto_found_10 * pow(num_large_found_10, -1);

    cout << "\n" << "1-200 ph.e.:     num_veto_found_1 / num_large_found_1 = " << num_veto_found_1 << " / " << num_large_found_1 << " = " << veto_det_ratio_1 << endl;
    
    cout << "\n" << "201-400 ph.e.:   num_veto_found_2 / num_large_found_2 = " << num_veto_found_2 << " / " << num_large_found_2 << " = " << veto_det_ratio_2 << endl;

    cout << "\n" << "401-600 ph.e.:   num_veto_found_3 / num_large_found_3 = " << num_veto_found_3 << " / " << num_large_found_3 << " = " << veto_det_ratio_3 << endl;

    cout << "\n" << "601-800 ph.e.:   num_veto_found_4 / num_large_found_4 = " << num_veto_found_4 << " / " << num_large_found_4 << " = " << veto_det_ratio_4 << endl;

    cout << "\n" << "801-1000 ph.e.:  num_veto_found_5 / num_large_found_5 = " << num_veto_found_5 << " / " << num_large_found_5 << " = " << veto_det_ratio_5 << endl;

    cout << "\n" << ">1000 ph.e.: num_veto_found_6 / num_large_found_6 = " << num_veto_found_6 << " / " << num_large_found_6 << " = " << veto_det_ratio_6 << endl;

    // cout << "\n" << "1201-1400 ph.e.: num_veto_found_7 / num_large_found_7 = " << num_veto_found_7 << " / " << num_large_found_7 << " = " << veto_det_ratio_7 << endl;

    // cout << "\n" << "1401-1600 ph.e.: num_veto_found_8 / num_large_found_8 = " << num_veto_found_8 << " / " << num_large_found_8 << " = " << veto_det_ratio_8 << endl;

    // cout << "\n" << "1601-1800 ph.e.: num_veto_found_9 / num_large_found_9 = " << num_veto_found_9 << " / " << num_large_found_9 << " = " << veto_det_ratio_9 << endl;

    // cout << "\n" << "1801-2000 ph.e.: num_veto_found_10 / num_large_found_10 = " << num_veto_found_10 << " / " << num_large_found_10 << " = " << veto_det_ratio_10 << endl;

    cout << "\n" << "End of code." << endl;
    
    return 0;

}
