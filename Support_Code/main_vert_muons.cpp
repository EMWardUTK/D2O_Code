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

// std::vector<double> amplitudeToPE{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15.9,16.4,15.65,15.9,16.3,0,10.0,15.3,14.9,0,14.75,14,14.45,15.35 };          // 10/5: run2757

// std::vector<double> amplitudeToPE{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15.4,16.5,15.95,15.3,15.7,11.0,14.95,14.8,16.1,16.0,16.15,16.25 };              // 11/7: run3529

// std::vector<double> amplitudeToPE{ 28,27,30,29,26,11,29,28,30,29,24,29 };             // 12/2: run4144

std::vector<double> amplitudeToPE{ 26,28,30,29,27,12,30,28,30,30,24,30 };                // 12/3: run4176

// std::vector<double> amplitudeToPE{ 25,27,30,29,25,10,27,27,30,30,23,30 };             // 12/4: run4193

// std::vector<double> integralToPE{ 102.068,104.926,103.562,103.482,106.928,69.5778,101.677,100.208,100.963,120.495,102.519,104.107 };                  // 01/23: run2001

std::vector<double> MichelSiPMNoiseAvg1{ 187.887,188.1116667,314.3093333,195.7781667,190.6345,202.5948333,209.708,201.6861667,243.2 };                          // Jul 17, 2023 - Aug 29, 2023

std::vector<double> MichelSiPMNoiseAvg2{ 185.7410476,174.2948095,461.303619,192.1831429,190.1772857,189.2594762,187.5979524,195.2732857,282.2847619 };          // Sep 3, 2023 - Jan 30, 2024

std::vector<double> MichelSiPMNoiseAvg3{ 220.497,202.8251905,545.8977619,220.8942857,220.4852381,222.7695238,220.9086667,231.9756667,300.433381 };              // Feb 9, 2024 - Jul 21, 2024

std::vector<double> MichelSiPMNoiseAvg4{ 222.2114286,213.4123571,231.4738571,222.3042143,220.602,234.8156429,237.2396429,240.2786429,282.2042143 };             // Jul 31, 2024 - Nov 23, 2024

std::vector<double> MichelSiPMNoiseAvg5{ 214.9131176,219.8434706,227.0522353,211.6268824,215.7442353,222.5197059,217.8758824,230.3478235,258.0880588 };         // Jan 12, 2025 - May 24, 2025

std::vector<double> MichelSiPMNoiseAvg6{ 202.2175,214.0718,190.558875,202.427125,209.80825,212.1195,195.1185,225.44875,227.461 };                               // May 31, 2025 - Aug 23, 2025

std::vector<double> MichelSiPMNoiseAvg7{ 201.0280968,219.3897823,216.1231694,207.7725242,209.4667177,215.816371,227.0340161,231.9643629,268.8965081 };          // Aug 11, 2025 - Jun 24, 2026

std::vector<double> MinMuonSiPMCut1{ 400,500,800,1000,300,400,350,250,400 };            // Jul 17, 2023 - Jan 30, 2024

std::vector<double> MinMuonSiPMCut2{ 600,700,900,1100,400,500,450,350,500 };            // Feb 9, 2024 - Dec 15, 2024

std::vector<double> MinMuonSiPMCut3{ 500,600,900,1100,400,400,350,300,500 };            // Jan 12, 2025 - Aug 24, 2025

std::vector<double> MinMuonSiPMCut23{ 850,1050,1300,1475,625,800,800,600,500 };         // Jul 31, 2024 - Dec 15, 2024

std::vector<double> MinMuonSiPMCut4{ 500,600,800,900,400,500,400,400,500 };             // Aug 11, 2025 - Jun 24, 2026

std::vector<double> integralToPEF25S26{ 104.1735026,102.1484265,105.5107949,100.8622726,107.9768889,108.4633248,107.4503248,106.1125983,104.1273932,112.0,105.3396068,104.9230085 };            // Aug 11, 2025 - Jun 24, 2026

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

    TH1D* h_intlowLED = new TH1D("h_intlowLED", "Distribution of Low Light LED Integral Values", 1000, -500, 500);

    TH1D* h_inthighLED = new TH1D("h_inthighLED", "Distribution of High Light LED Integral Values", 80, 0, 800);            //, 11000, -1000, 10000);

    TH1D* h_intminLED = new TH1D("h_intminLED", "Distribution of Minimum Bias LED Integral Values", 1000, -500, 500);

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

    TH1D* h_ml = new TH1D("h_ml", "Distribution of Internally Triggered Event Multiplicity", 17, -2, 15);

    TH2D* h_ml_v_int = new TH2D("h_ml_v_int", "Internally Triggered Event Multiplicity vs Integral Values", 17, -2, 15, 5500, -500, 5000);

    TH1D* h_pst = new TH1D("h_pst", "Distribution of Peak Start Times in Each Channel", 45, 0, 720);

    TH1D* h_mfpst = new TH1D("h_mfpst", "Distribution of Most Frequent Peak Start Times in Each Channel", 45, 0, 720);

    TH2D* h_pst_v_pint = new TH2D("h_pst_v_pint", "Distribution of Peak Start Times in Each Channel vs Peak Integral Values", 45, 0, 720, 100, 0, 1000);

    TH2D* h_mfpst_v_mfpint = new TH2D("h_mfpst_v_mfpint", "Distribution of Most Frequent Peak Start Times vs Most Frequent Peak Integral Values", 45, 0, 720, 100, 0, 1000);

    TH1D* h_var = new TH1D("h_var", "Distribution of Variance in Peak Start Times", 48, 0, 480);

    TH1D* h_spheLL = new TH1D("h_spheLL", "Distribution of Single Channel Low Light Integral Values", 1000, -500, 500);

    TH1D* h_spheMB = new TH1D("h_spheMB", "Distribution of Single Channel Minimum Bias Integral Values", 1000, -500, 500);

    TH1D* h_intmuons = new TH1D("h_intmuons", "Distribution of Muon Integral Values", 100, -500, 5000);

    TH1D* h_intmichels = new TH1D("h_intmichels", "Distribution of Michel Electron Integral Values", 80, 0, 800);

    TH1D* h_intcosmics = new TH1D("h_intcosmics", "Distribution of Cosmic Integral Values", 400, 0, 4000);

    TH1D* h_intcosmics_S25 = new TH1D("h_intcosmics_S25", "Vertical Cosmics Integral Values, with No High Light Mean Correction", 400, 0, 4000);

    TH1D* h_intcosmics_F25 = new TH1D("h_intcosmics_F25", "Vertical Cosmics Integral Values, with No High Light Mean Correction", 400, 0, 4000);

    TH1D* h_intcosmics_S26 = new TH1D("h_intcosmics_S26", "Vertical Cosmics Integral Values, with No High Light Mean Correction", 400, 0, 4000);

    TH1D* h_inttemp = new TH1D("h_inttemp", "Cosmic Integral Values, One Hour Run", 400, 0, 4000);

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

    TH1D* h_vpint_16 = new TH1D("h_vpint_16", "SiPM 1 Vertical Cosmics Integral Values", 150, -400, 2000);

    TH1D* h_vpint_17 = new TH1D("h_vpint_17", "SiPM 2 Vertical Cosmics Integral Values", 150, -400, 2000);

    TH1D* h_vpint_18 = new TH1D("h_vpint_18", "SiPM 3 Vertical Cosmics Integral Values", 150, -400, 2000);

    TH1D* h_vpint_19 = new TH1D("h_vpint_19", "SiPM 4 Vertical Cosmics Integral Values", 150, -400, 2000);

    TH1D* h_vpint_20 = new TH1D("h_vpint_20", "SiPM 5 Vertical Cosmics Integral Values", 150, -400, 2000);

    TH1D* h_vpint_21 = new TH1D("h_vpint_21", "SiPM 6 Vertical Cosmics Integral Values", 150, -400, 2000);

    TH1D* h_vpint_22 = new TH1D("h_vpint_22", "SiPM 7 Vertical Cosmics Integral Values", 150, -400, 2000);

    TH1D* h_vpint_23 = new TH1D("h_vpint_23", "SiPM 8 Vertical Cosmics Integral Values", 150, -400, 2000);

    TH1D* h_vpint_top = new TH1D("h_vpint_top", "SiPM 9 & 10 Vertical Cosmics Integral Values", 500, -500, 7500);

    TH2D* h_detint_vs_vpint_16 = new TH2D("h_detint_vs_vpint_16", "SiPM 1 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_17 = new TH2D("h_detint_vs_vpint_17", "SiPM 2 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_18 = new TH2D("h_detint_vs_vpint_18", "SiPM 3 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_19 = new TH2D("h_detint_vs_vpint_19", "SiPM 4 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_20 = new TH2D("h_detint_vs_vpint_20", "SiPM 5 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_21 = new TH2D("h_detint_vs_vpint_21", "SiPM 6 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_22 = new TH2D("h_detint_vs_vpint_22", "SiPM 7 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_23 = new TH2D("h_detint_vs_vpint_23", "SiPM 8 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH2D* h_detint_vs_vpint_top = new TH2D("h_detint_vs_vpint_top", "SiPM 9 & 10 Integral vs Detector Integral", 1000, 0, 3000, 1000, 0, 6000);

    TH1D* h_dt_low = new TH1D("h_dt_low", "Time Separations Between Peaks, <500phe Muons", 1000, 0, muon_dt_max);

    TH1D* h_intmuons_low = new TH1D("h_intmuons_low", "Muon Integral Values, <500phe Muons", 250, 0, 2500);

    TH1D* h_intmichels_low = new TH1D("h_intmichels_low", "Michel Electron Integral Values, <500phe Muons", 100, 0, 1000);

    TH1D* h_vpintmichels_low = new TH1D("h_vpintmichels_low", "Michel Electron Top Veto Panel Integral Values, <500phe Muons", 100, 0, 1000);

    TH1D* h_dt_med_low = new TH1D("h_dt_med_low", "Time Separations Between Peaks, 500-1000phe Muons", 1000, 0, muon_dt_max);

    TH1D* h_intmuons_med_low = new TH1D("h_intmuons_med_low", "Muon Integral Values, 500-1000phe Muons", 250, 0, 2500);

    TH1D* h_intmichels_med_low = new TH1D("h_intmichels_med_low", "Michel Electron Integral Values, 500-1000phe Muons", 100, 0, 1000);

    TH1D* h_vpintmichels_med_low = new TH1D("h_vpintmichels_med_low", "Michel Electron Top Veto Panel Integral Values, <500phe Muons", 100, 0, 1000);

    TH1D* h_dt_med = new TH1D("h_dt_med", "Time Separations Between Peaks, 1000-1500phe Muons", 1000, 0, muon_dt_max);

    TH1D* h_intmuons_med = new TH1D("h_intmuons_med", "Muon Integral Values, 1000-1500phe Muons", 250, 0, 2500);

    TH1D* h_intmichels_med = new TH1D("h_intmichels_med", "Michel Electron Integral Values, 1000-1500phe Muons", 100, 0, 1000);

    TH1D* h_vpintmichels_med = new TH1D("h_vpintmichels_med", "Michel Electron Top Veto Panel Integral Values, <500phe Muons", 100, 0, 1000);

    TH1D* h_dt_med_high = new TH1D("h_dt_med_high", "Time Separations Between Peaks, 1500-2000phe Muons", 1000, 0, muon_dt_max);

    TH1D* h_intmuons_med_high = new TH1D("h_intmuons_med_high", "Muon Integral Values, 1500-2000phe Muons", 250, 0, 2500);

    TH1D* h_intmichels_med_high = new TH1D("h_intmichels_med_high", "Michel Electron Integral Values, >1500phe Muons", 100, 0, 1000);

    TH1D* h_vpintmichels_med_high = new TH1D("h_vpintmichels_med_high", "Michel Electron Top Veto Panel Integral Values, <500phe Muons", 100, 0, 1000);

    TH1D* h_dt_high = new TH1D("h_dt_high", "Time Separations Between Peaks, >2000phe Muons", 1000, 0, muon_dt_max);

    TH1D* h_intmuons_high = new TH1D("h_intmuons_high", "Muon Integral Values, >2000phe Muons", 250, 0, 2500);

    TH1D* h_intmichels_high = new TH1D("h_intmichels_high", "Michel Electron Integral Values, >2000phe Muons", 100, 0, 1000);

    TH1D* h_vpintmichels_high = new TH1D("h_vpintmichels_high", "Michel Electron Top Veto Panel Integral Values, <500phe Muons", 100, 0, 1000);

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
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/allweeklyruns/allweeklyruns_%i.txt", run));
    // 16151-16174 = 1721442002-1721524801 = 7/19/24-7/20/24 = file ~48 // 16351-16374 = 1722090001-1722172801 = 7/27/24-7/28/24 = file ~49
    // 16647-16671 = 1723155601-1723242001 = 8/08/24-8/09/24 = file 50  // 16855-16878 = 1723870861-1723953661 = 8/17/24-8/18/24 = file 51
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/F24-S25_GR_Runs/F24-S25_GR_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/F25-S26_GR_Runs_4084GR/F25-S26_GR_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/F25-S26_GR_Weekly_4084GR/F25-S26_GR_Weekly_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/S25-S26_GR_Runs_6235GR/S25-S26_GR_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/S25-S26_GR_Weekly_6149GR/S25-S26_GR_Weekly_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/S25F25S26_GR_3_rand_days_%i.txt", run));
    ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/New_F25-S26_GR_%i.txt", run));

    // Use a while loop together with the getline() function to read the file line by line

    int iLine = 0;

    int RunListLines = 0;

    while (getline(ReadRunListFile, intval1)) {

        // Output the text from the file

        runlist.push_back(stod(intval1));

        RunListLines += 1;

        // if (RunListLines == 24) {break;}
        
    }
    
    // Close the file

    ReadRunListFile.close();

    ReadRunListFile.clear();

    // Create a text string, which is used to output the text file

    int HLmeancounter = 0;

    double sumHLmean = 0.0;

    string int_val_1;

    // Read from the text file         

    // To use F24-S25 data, delete "_F25-S26_GR"
    // To use F25-S26 data, add "_F25-S26_GR" before "_1-124"
    // ifstream ReadHLMFile("Vertical_Muon_Peak_Center_1-124_Run_Vals.txt");          // Vertical_Muon_Peak_Center_1-124_Val_Only.txt         // hl_mean_val_all_1.txt          // hl_mean_val_all_temp.txt         // hl_mean_val_all_F24-S25_GR.txt
    ifstream ReadHLMFile("Vertical_Muon_Peak_Center_F25-S26_GR_1-124_Run_Vals.txt");

    // Use a while loop together with the getline() function to read the file line by line

    while (getline(ReadHLMFile, int_val_1)) {

        if (stod(int_val_1) > 0.0) {

            HLmeancounter += 1;

            sumHLmean += stod(int_val_1);

        }
        
    }

    double defaultHLmean = 1813.94;         // = 1613.271109; (F25-S26 Vert Muon Mean Value)           // = 1813.94; (12/1/24 Vert Muon Mean Value)         // = 34.5323; (11/26/24 HL Mean Value)         // sumHLmean / HLmeancounter;           // sumHLmean / HLmeancounter = 1783.874 (Averaged value)

    double currentHLmean = defaultHLmean;

    // Close the file

    ReadHLMFile.close();

    ReadHLMFile.clear();

    // Create and open a text file
    ofstream VMCFile(Form("Vertical_Muon_Count_F25-S26_GR_Run%i.txt", run));
    ofstream VMAFile(Form("Vertical_Muon_Amplitude_F25-S26_GR_Run%i.txt", run));
    ofstream VMPFile(Form("Vertical_Muon_Peak_Center_F25-S26_New_GR_Run%i.txt", run));
    ofstream VMIFile(Form("Vertical_Muon_Peak_Center_Val_Only_F25-S26_New_GR_Run%i.txt", run));

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

        int ev61_min = 1200; int ev61_max = 2000;

        if (run_iterable < 24813) {ev61_min = 1200; ev61_max = 2000;}

        else if (run_iterable >= 24813) {ev61_min = 700; ev61_max = 1000;}

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

        // To use F24-S25 data, delete "_F25-S26_GR"
        // To use F25-S26 data, add "_F25-S26_GR" before "_1-124"
        // ifstream ReadHLMFile("Vertical_Muon_Peak_Center_1-124_Run_Nums.txt");          // Vertical_Muon_Peak_Center_1-124_No_HL_Mean_Correction.txt          // hl_mean_val_look_up_1.txt          // hl_mean_val_look_up_temp.txt           // hl_mean_val_look_up_F24-S25_GR.txt      // hl_mean_val_look_up_51-82_no63-65.txt (Aug 17th, 2024 - May 24th, 2025, excluding beam off period in December 2024, these runs are NOT just golden runs)
        // ifstream ReadHLMFile("Vertical_Muon_Peak_Center_F25-S26_GR_1-124_Run_Nums.txt");
        // Runs 16830 - 32398 (F24-S26)
        ifstream ReadHLMFile("Vertical_Muon_Peak_Center_All_F24-S26_GR_Run_Nums.txt");

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

        // if (run_iterable < 24813) {defaultHLmean = 1813.94;}           // = 1813.94; (12/1/24 Vert Muon Mean Value)

        // else if (run_iterable >= 24813) {defaultHLmean = 1613.271109;}         // = 1613.271109; (F25-S26 Vert Muon Mean Value)

        // Close the file

        ReadHLMFile.close();

        ReadHLMFile.clear();

        // Experimental rotational correction, to correct for Michel integral change over time leftover/created by vertical muon correction

        double rot_up_corr = run_starttime * 0.00000000042694381419 + 0.2600815812;

        double rot_down_corr = run_starttime * -0.00000000060991973456 + 2.0570263125;

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

        else if (run_starttime > 1748775661 && run_starttime < 1754892061) {vp_cut_vals = MichelSiPMNoiseAvg6;}

        else if (run_starttime >= 1754892061) {vp_cut_vals = MichelSiPMNoiseAvg7;}

        // Select appropriate Minimum Muon values

        std::vector<double> min_muon_vals;

        if (run_starttime <= 1706628002) {min_muon_vals = MinMuonSiPMCut1;}

        else if (run_starttime > 1706628002 && run_starttime <= 1734260462) {min_muon_vals = MinMuonSiPMCut2;}

        // else if (run_starttime > 1706628002 && run_starttime <= 1722417601) {min_muon_vals = MinMuonSiPMCut2;}

        // else if (run_starttime > 1722417601 && run_starttime <= 1734260462) {min_muon_vals = MinMuonSiPMCut23;}         // std::vector<double> MinMuonSiPMCut23{ 850,1050,1300,1475,625,800,800,600,500 };          // Jul 31, 2024 - Dec 15, 2024

        else if (run_starttime > 1734260462 && run_starttime < 1754892061) {min_muon_vals = MinMuonSiPMCut3;}

        else if (run_starttime >= 1754892061) {min_muon_vals = MinMuonSiPMCut4;}

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

        double avg_peak_pos_RMS = 0.0;

        double var_peak_pos_RMS = 0.0;

        int ch15_on = 0;

        double ev61_time = 0.0;
        double minLED_time = 0.0;
        double lowLED_time = 0.0;
        double highLED_time = 0.0;

        int counter = 0;
    
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

                if (iChan <= 11) {sphe_int = run_starttime * y_val[iChan] + b_val[iChan];}

                if (run_starttime >= 1754892061 && iChan <= 11) {sphe_int = integralToPEF25S26[iChan];}

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

                    if (Ev61Energy > ev61_min && Ev61Energy < ev61_max) {

                        all_chan_beam = true;

                        ch15_on += 1;

                    }

                    // if (triggerBits == 1 || triggerBits == 3 || triggerBits == 33 || triggerBits == 35) {

                    h_intev61->Fill(Ev61Energy);
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
                                p.energy = pulseEnergy / sphe_int; // sphe_int;

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
            avg_pulse.energy = std::accumulate(all_chan_energy.begin(), all_chan_energy.end(), 0.0);            // / (currentHLmean / defaultHLmean);           // Integrals still sum over secondary peaks discarded in .start and .end values 
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

            if (avg_pulse.beam == true && avg_pulse.number >= 10) {           // && (avg_pulse.trigger >= 32 && avg_pulse.trigger <= 35)

                h_intneu->Fill(avg_pulse.energy);

            }

            if ((avg_pulse.trigger >= 32 && avg_pulse.trigger <= 35) && avg_pulse.number >= 10) {

                amp_val.push_back(avg_pulse.peak);

            }

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

            if (avg_pulse.trigger >= 32 && avg_pulse.trigger <= 35) {

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

                    highLED_time = avg_pulse.start;

                }

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

            }

            else if (avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {

                h_intallint->Fill(avg_pulse.energy);

                intallint.push_back(avg_pulse.energy);

                if (avg_pulse.number >= 10) {

                    h_intallmultcut->Fill(avg_pulse.energy);

                    intallmultcut.push_back(avg_pulse.energy);

                    // if (p_int_16 <= 750 && p_int_17 <= 950 && p_int_18 <= 1200 && p_int_19 <= 1375 && p_int_20 <= 525 && p_int_21 <= 700 && p_int_22 <= 700 && p_int_23 <= 500 && avg_pulse.top_vp_energy <= 450) {

                    if (p_int_16 <= min_muon_vals[0] && p_int_17 <= min_muon_vals[1] && p_int_18 <= min_muon_vals[2] && p_int_19 <= min_muon_vals[3] && p_int_20 <= min_muon_vals[4] && p_int_21 <= min_muon_vals[5] && p_int_22 <= min_muon_vals[6] && p_int_23 <= min_muon_vals[7] && avg_pulse.vp_int_top <= min_muon_vals[8]) {

                        h_intallmultmuoncut->Fill(avg_pulse.energy);

                        intallmultmuoncut.push_back(avg_pulse.energy);

                    }

                }

            }

            // Search for cosmic events

            if (avg_pulse.energy >= 50 && avg_pulse.number >= 10 && avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16 && avg_pulse.start - last_event_time > echo_cut) {

                // if ((p_int_16 < 750 && p_int_17 < 950 && p_int_18 < 1200 && p_int_19 < 1375 && p_int_20 < 525 && p_int_21 < 700 && p_int_22 < 700 && p_int_23 < 500) && avg_pulse.vp_int_top > 450) {

                // if (p_int_16 > vp_cut_vals[0]  p_int_17 > vp_cut_vals[1]  p_int_18 > vp_cut_vals[2]  p_int_19 > vp_cut_vals[3]  p_int_20 > vp_cut_vals[4]  p_int_21 > vp_cut_vals[5]  p_int_22 > vp_cut_vals[6]  p_int_23 > vp_cut_vals[7]  avg_pulse.vp_int_top > vp_cut_vals[8]) {

                if ((p_int_16 < min_muon_vals[0] && p_int_17 < min_muon_vals[1] && p_int_18 < min_muon_vals[2] && p_int_19 < min_muon_vals[3] && p_int_20 < min_muon_vals[4] && p_int_21 < min_muon_vals[5] && p_int_22 < min_muon_vals[6] && p_int_23 < min_muon_vals[7]) && avg_pulse.vp_int_top > min_muon_vals[8]) {

                    cosmic_events += 1;

                    h_inttemp->Fill(avg_pulse.energy);

                    h_intcosmics->Fill(avg_pulse.energy);

                    if (run_iterable >= 22338 && run_iterable <= 22361) {h_intcosmics_S25->Fill(avg_pulse.energy);}

                    if (run_iterable >= 25749 && run_iterable <= 25772) {h_intcosmics_F25->Fill(avg_pulse.energy);}

                    if (run_iterable >= 31291 && run_iterable <= 31314) {h_intcosmics_S26->Fill(avg_pulse.energy);}

                    h_vpint_16->Fill(avg_pulse.vp_int_16);

                    h_vpint_17->Fill(avg_pulse.vp_int_17);

                    h_vpint_18->Fill(avg_pulse.vp_int_18);

                    h_vpint_19->Fill(avg_pulse.vp_int_19);

                    h_vpint_20->Fill(avg_pulse.vp_int_20);

                    h_vpint_21->Fill(avg_pulse.vp_int_21);

                    h_vpint_22->Fill(avg_pulse.vp_int_22);

                    h_vpint_23->Fill(avg_pulse.vp_int_23);

                    h_vpint_top->Fill(avg_pulse.vp_int_top);
                    
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

            if (time_count == 0.0 && (avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) && avg_pulse.number >= 10 && avg_pulse.energy >= 50 && avg_pulse.start - last_event_time > echo_cut) {            // avg_pulse.peak > PULSE_PE_THRESHOLD && 

                // if ((p_int_16 < 750 && p_int_17 < 950 && p_int_18 < 1200 && p_int_19 < 1375 && p_int_20 < 525 && p_int_21 < 700 && p_int_22 < 700 && p_int_23 < 500) && avg_pulse.vp_int_top > 450) {

                if ((p_int_16 < min_muon_vals[0] && p_int_17 < min_muon_vals[1] && p_int_18 < min_muon_vals[2] && p_int_19 < min_muon_vals[3] && p_int_20 < min_muon_vals[4] && p_int_21 < min_muon_vals[5] && p_int_22 < min_muon_vals[6] && p_int_23 < min_muon_vals[7]) && avg_pulse.vp_int_top > min_muon_vals[8]) {
                
                    WithinPulse = true;

                    first_pulse_start_time = avg_pulse.start;

                    MichelVetoActivity = true;

                }

            }

            if (WithinPulse && (avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) && avg_pulse.number >= 10 && avg_pulse.energy >= 50 && avg_pulse.start - last_event_time > echo_cut) {                // Conditions on if statement mean not all counted pulses get stored in mult_pulses!!!
                
                // if ((p_int_16 < 750 && p_int_17 < 950 && p_int_18 < 1200 && p_int_19 < 1375 && p_int_20 < 525 && p_int_21 < 700 && p_int_22 < 700 && p_int_23 < 500) && avg_pulse.vp_int_top > 450) {

                if ((p_int_16 < min_muon_vals[0] && p_int_17 < min_muon_vals[1] && p_int_18 < min_muon_vals[2] && p_int_19 < min_muon_vals[3] && p_int_20 < min_muon_vals[4] && p_int_21 < min_muon_vals[5] && p_int_22 < min_muon_vals[6] && p_int_23 < min_muon_vals[7]) && avg_pulse.vp_int_top > min_muon_vals[8]) {
                
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

                else if (p_int_16 <= min_muon_vals[0] && p_int_17 <= min_muon_vals[1] && p_int_18 <= min_muon_vals[2] && p_int_19 <= min_muon_vals[3] && p_int_20 <= min_muon_vals[4] && p_int_21 <= min_muon_vals[5] && p_int_22 <= min_muon_vals[6] && p_int_23 <= min_muon_vals[7] && avg_pulse.vp_int_top <= min_muon_vals[8]) {

                    mult_pulses.push_back(avg_pulse);

                    time_count = avg_pulse.start - first_pulse_start_time;

                }

                if (time_count == 0.0) {time_count = 0.5;}

                if (time_count > 1 && p_int_16 <= 300 && p_int_17 <= 300 && p_int_18 <= 600 && p_int_19 <= 300 && p_int_20 <= 300 && p_int_21 <= 300 && p_int_22 <= 300 && p_int_23 <= 300 && avg_pulse.vp_int_top <= 300) {MichelVetoActivity = false;}
                
                // else if (run == 7894 && time_count > 1 && p_int_16 <= 100 && p_int_17 <= 100 && p_int_18 <= 200 && p_int_19 <= 100 && p_int_20 <= 100 && p_int_21 <= 100 && p_int_22 <= 100 && p_int_23 <= 100 && avg_pulse.vp_int_top <= 150) {MichelVetoActivity = false;}
                
                // else if (run == 11241 && time_count > 1 && p_int_16 <= 100 && p_int_17 <= 100 && p_int_18 <= 300 && p_int_19 <= 100 && p_int_20 <= 100 && p_int_21 <= 100 && p_int_22 <= 100 && p_int_23 <= 100 && avg_pulse.vp_int_top <= 200) {MichelVetoActivity = false;}
                
                // else if (run == 12636 && time_count > 1 && p_int_16 <= 150 && p_int_17 <= 150 && p_int_18 <= 400 && p_int_19 <= 150 && p_int_20 <= 150 && p_int_21 <= 150 && p_int_22 <= 150 && p_int_23 <= 150 && avg_pulse.vp_int_top <= 250) {MichelVetoActivity = false;}

                // else if (run == 15312 && time_count > 1 && p_int_16 <= 300 && p_int_17 <= 300 && p_int_18 <= 600 && p_int_19 <= 300 && p_int_20 <= 300 && p_int_21 <= 300 && p_int_22 <= 300 && p_int_23 <= 300 && avg_pulse.vp_int_top <= 300) {MichelVetoActivity = false;}

                else {MichelVetoActivity = true;}

            }

            if (time_count > 2000 && time_count <= muon_dt_max) {          // && MichelVetoActivity == false

                // Fill histograms here?

                hist_events += 1;

                tot_michel_energy += mult_pulses.back().energy;

                h_dt->Fill(time_count);
                            
                h_intmuons->Fill(mult_pulses.front().energy);
                            
                h_intmichels->Fill(mult_pulses.back().energy);

                if (mult_pulses.front().energy < 500) {h_dt_low->Fill(time_count); h_intmuons_low->Fill(mult_pulses.front().energy); h_intmichels_low->Fill(mult_pulses.back().energy); h_vpintmichels_low->Fill(mult_pulses.back().vp_int_top);}

                if (mult_pulses.front().energy >= 500 && mult_pulses.front().energy < 1000) {h_dt_med_low->Fill(time_count); h_intmuons_med_low->Fill(mult_pulses.front().energy); h_intmichels_med_low->Fill(mult_pulses.back().energy); h_vpintmichels_med_low->Fill(mult_pulses.back().vp_int_top);}

                if (mult_pulses.front().energy >= 1000 && mult_pulses.front().energy < 1500) {h_dt_med->Fill(time_count); h_intmuons_med->Fill(mult_pulses.front().energy); h_intmichels_med->Fill(mult_pulses.back().energy); h_vpintmichels_med->Fill(mult_pulses.back().vp_int_top);}

                if (mult_pulses.front().energy >= 1500) {h_dt_med_high->Fill(time_count); h_intmuons_med_high->Fill(mult_pulses.front().energy); h_intmichels_med_high->Fill(mult_pulses.back().energy); h_vpintmichels_med_high->Fill(mult_pulses.back().vp_int_top);}

                // if (mult_pulses.front().energy > 2000) {h_dt_high->Fill(time_count); h_intmuons_high->Fill(mult_pulses.front().energy); h_intmichels_high->Fill(mult_pulses.back().energy); h_vpintmichels_high->Fill(mult_pulses.back().vp_int_top);}

                if (mult_pulses.back().energy > 200 && mult_pulses.back().energy < 400) {h_slice_vpintavg->Fill((mult_pulses.back().side_vp_energy + mult_pulses.back().top_vp_energy) / 9);}

                // if (mult_pulses.back().vp_int_18 > 350 && mult_pulses.back().vp_int_18 < 450) {cout << "\n" << "350 - 450 ADC: " << run_iterable << " " << iEnt << endl;}

                // if (mult_pulses.back().vp_int_18 > -450 && mult_pulses.back().vp_int_18 < -350) {cout << "\n" << "-450 - -350 ADC: " << run_iterable << " " << iEnt << endl;}

                // if (mult_pulses.back().vp_int_18 > 1400 && mult_pulses.back().vp_int_18 < 1500) {cout << "\n" << "1200 - 1400 ADC: " << run_iterable << " " << iEnt << endl;}

                h_dt_61_16->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_16 - mult_pulses.front().wf_time); // h_vpint_16->Fill(mult_pulses.front().vp_int_16);

                h_dt_61_17->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_17 - mult_pulses.front().wf_time); // h_vpint_17->Fill(mult_pulses.front().vp_int_17);

                h_dt_61_18->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_18 - mult_pulses.front().wf_time); // h_vpint_18->Fill(mult_pulses.front().vp_int_18);

                h_dt_61_19->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_19 - mult_pulses.front().wf_time); // h_vpint_19->Fill(mult_pulses.front().vp_int_19);

                h_dt_61_20->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_20 - mult_pulses.front().wf_time); // h_vpint_20->Fill(mult_pulses.front().vp_int_20);

                h_dt_61_21->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_21 - mult_pulses.front().wf_time); // h_vpint_21->Fill(mult_pulses.front().vp_int_21);

                h_dt_61_22->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_22 - mult_pulses.front().wf_time); // h_vpint_22->Fill(mult_pulses.front().vp_int_22);

                h_dt_61_23->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_23 - mult_pulses.front().wf_time); // h_vpint_23->Fill(mult_pulses.front().vp_int_23);

                h_dt_61_top->Fill(mult_pulses.back().wf_time + mult_pulses.back().vp_start_top - mult_pulses.front().wf_time); // h_vpint_top->Fill(mult_pulses.front().vp_int_top);

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

                // Write to the file

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

                // if ((p_int_16 < 750 && p_int_17 < 950 && p_int_18 < 1200 && p_int_19 < 1375 && p_int_20 < 525 && p_int_21 < 700 && p_int_22 < 700 && p_int_23 < 500) && avg_pulse.vp_int_top > 450) {

                if ((p_int_16 < min_muon_vals[0] && p_int_17 < min_muon_vals[1] && p_int_18 < min_muon_vals[2] && p_int_19 < min_muon_vals[3] && p_int_20 < min_muon_vals[4] && p_int_21 < min_muon_vals[5] && p_int_22 < min_muon_vals[6] && p_int_23 < min_muon_vals[7]) && avg_pulse.vp_int_top > min_muon_vals[8]) {

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

                if (avg_pulse.energy > 1000 && avg_pulse.energy <= 1200) {num_large_found_6 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_6 += 1;}}
            
                if (avg_pulse.energy > 1200 && avg_pulse.energy <= 1400) {num_large_found_7 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_7 += 1;}}
            
                if (avg_pulse.energy > 1400 && avg_pulse.energy <= 1600) {num_large_found_8 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_8 += 1;}}
            
                if (avg_pulse.energy > 1600 && avg_pulse.energy <= 1800) {num_large_found_9 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_9 += 1;}}
            
                if (avg_pulse.energy > 1800 && avg_pulse.energy <= 2000) {num_large_found_10 += 1; if (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450) {num_veto_found_10 += 1;}}
            
            }

            last_event_time = avg_pulse.start;

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

        TF1* func_inttemp = new TF1("func_inttemp", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))", 1200, 2200);
        double Par_temp[] = {10, 500, 1700, 200};
        func_inttemp->SetParameters(Par_temp);
        func_inttemp->SetParName(0, "Baseline");
        func_inttemp->SetParName(1, "Amplitude");
        func_inttemp->SetParName(2, "PeakCenter");
        func_inttemp->SetParName(3, "StandardDeviation");
        h_inttemp->Fit("func_inttemp", "R");
        double_t am_cosmics = func_inttemp->GetParameter(1);
        double_t pc_cosmics = func_inttemp->GetParameter(2);
        h_inttemp->Reset();

        // Write to the file
        VMCFile << run_iterable << "\t" << cosmic_events << endl;
        VMAFile << run_iterable << "\t" << am_cosmics << endl;
        VMPFile << run_iterable << "\t" << pc_cosmics << endl;
        VMIFile << pc_cosmics << endl;

        cosmic_events = 0;

        ch15_on = 0;

        all_nsTime.clear();

        amp_on_off_val.push_back(amp_val);

        amp_val.clear();

        integralToPE.clear();

        vp_cut_vals.clear();

        u_event_bin_val.clear();

        avg_peak_pos_RMS = getAverage(peak_pos_RMS);

        var_peak_pos_RMS = variance(peak_pos_RMS);

        peak_pos_RMS.clear();

        run_counter += 1;

        if (run_counter == floor(0.5 * RunListLines)) {hour12_st = run_starttime; hour12_run = run_iterable;}

        // else if (run_counter == 24) {run_counter = 0;}

        f->Close();

    } // Run loop

    double avg_mich_en = tot_michel_energy * pow(hist_events, -1);

    hist_events = 0;

    tot_michel_energy = 0.0;

    TCanvas* c_intcosmics = new TCanvas("c_intcosmics", "Cosmic Event Integral Values", 900, 700);
    // c_intcosmics->SetLogy();
    c_intcosmics->cd();

    TF1* func_intcosmics = new TF1("func_intcosmics", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))", 1200, 2200);
    double Par_cosmics[] = {100, 5500, 1700, 200};
    func_intcosmics->SetParameters(Par_cosmics);
    func_intcosmics->SetParName(0, "Baseline");
    func_intcosmics->SetParName(1, "Amplitude");
    func_intcosmics->SetParName(2, "PeakCenter");
    func_intcosmics->SetParName(3, "StandardDeviation");
    h_intcosmics->Fit("func_intcosmics", "R");

    /*double_t am_cosmics = func_intcosmics->GetParameter(1);
    double_t pc_cosmics = func_intcosmics->GetParameter(2);

    VMCFile << hour12_run << "\t" << cosmic_events << endl;
    VMAFile << hour12_run << "\t" << am_cosmics << endl;
    VMPFile << hour12_run << "\t" << pc_cosmics << endl;
    VMIFile << pc_cosmics << endl;*/

    h_intcosmics->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intcosmics->GetYaxis()->SetTitle("Counts");
    h_intcosmics->Draw();
    func_intcosmics->Draw("same");
    c_intcosmics->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_intcosmics_%i.png", run));
    // c_intcosmics->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_intcosmics.png");
    c_intcosmics->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_intcosmics_%i.png", run));
    h_intcosmics->Reset();

    TCanvas* c_intcosmics_F25 = new TCanvas("c_intcosmics_F25", "Vertical Cosmics Integral Values", 1200, 700);
    // c_intcosmics_F25->SetLogy();
    c_intcosmics_F25->cd();
    h_intcosmics_F25->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intcosmics_F25->GetYaxis()->SetTitle("Counts");
    h_intcosmics_S26->Draw("same");
    h_intcosmics_S26->SetLineColor(kGreen);
    h_intcosmics_F25->Draw("same");
    h_intcosmics_F25->SetLineColor(kBlue);
    h_intcosmics_S25->Draw("same");
    h_intcosmics_S25->SetLineColor(kRed);
    TLegend *leg_intcosmics_F25 = new TLegend(0.9, 0.65, 0.65, 0.45);
    leg_intcosmics_F25->AddEntry(h_intcosmics_S25, "Spring 2025", "l");
    leg_intcosmics_F25->AddEntry(h_intcosmics_F25, "Fall 2025", "l");
    leg_intcosmics_F25->AddEntry(h_intcosmics_S26, "Spring 2026", "l");
    leg_intcosmics_F25->Draw();
    c_intcosmics_F25->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_intcosmics_3days_%i.png", run));
    // c_intcosmics_F25->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_intcosmics_3days.png");
    c_intcosmics_F25->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_intcosmics_3days.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_intcosmics_3days.png"));
    h_intcosmics_F25->Reset();

    TCanvas* c_vpint_16 = new TCanvas("c_vpint_16", "Event 61 and SiPM dt Integral Values", 900, 700);
    c_vpint_16->SetLogy();
    c_vpint_16->cd();
    h_vpint_16->GetXaxis()->SetTitle("Integral (ADC)");
    h_vpint_16->GetYaxis()->SetTitle("Counts");
    h_vpint_16->Draw();
    c_vpint_16->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_vpint_16_%i.png", run));
    // c_vpint_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_vpint_16.png");
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
    c_vpint_17->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_vpint_17_%i.png", run));
    // c_vpint_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_vpint_17.png");
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
    c_vpint_18->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_vpint_18_%i.png", run));
    // c_vpint_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_vpint_18.png");
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
    c_vpint_19->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_vpint_19_%i.png", run));
    // c_vpint_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_vpint_19.png");
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
    c_vpint_20->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_vpint_20_%i.png", run));
    // c_vpint_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_vpint_20.png");
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
    c_vpint_21->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_vpint_21_%i.png", run));
    // c_vpint_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_vpint_21.png");
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
    c_vpint_22->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_vpint_22_%i.png", run));
    // c_vpint_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_vpint_22.png");
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
    c_vpint_23->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_vpint_23_%i.png", run));
    // c_vpint_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_vpint_23.png");
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
    c_vpint_top->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_vpint_top_%i.png", run));
    // c_vpint_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_vpint_top.png");
    c_vpint_top->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_top.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_top.png"));
    h_vpint_top->Reset();

    /*h_intmichels_low->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intmichels_low->GetYaxis()->SetTitle("Counts");
    h_intmichels_low->Scale( 1./h_intmichels_low->Integral());
    h_intmichels_low->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/jskweres/OutputFiles/plots/h_intmichels_low.png"));
    h_intmichels_low->Reset();

    h_intmichels_med_low->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intmichels_med_low->GetYaxis()->SetTitle("Counts");
    h_intmichels_med_low->Scale( 1./h_intmichels_med_low->Integral());
    h_intmichels_med_low->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/jskweres/OutputFiles/plots/h_intmichels_med_low.png"));
    h_intmichels_med_low->Reset();

    h_intmichels_med->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intmichels_med->GetYaxis()->SetTitle("Counts");
    h_intmichels_med->Scale( 1./h_intmichels_med->Integral());
    h_intmichels_med->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/jskweres/OutputFiles/plots/h_intmichels_med.png"));
    h_intmichels_med->Reset();

    h_intmichels_med_high->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intmichels_med_high->GetYaxis()->SetTitle("Counts");
    h_intmichels_med_high->Scale( 1./h_intmichels_med_high->Integral());
    h_intmichels_med_high->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/jskweres/OutputFiles/plots/h_intmichels_med_high.png"));
    h_intmichels_med_high->Reset();

    h_intmichels_high->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_intmichels_high->GetYaxis()->SetTitle("Counts");
    h_intmichels_high->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/jskweres/OutputFiles/plots/h_intmichels_high.png"));
    h_intmichels_high->Reset();

    h_vpintmichels_low->GetXaxis()->SetTitle("Integral (ADC)");
    h_vpintmichels_low->GetYaxis()->SetTitle("Counts");
    h_vpintmichels_low->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/jskweres/OutputFiles/plots/h_vpintmichels_low.png"));
    h_vpintmichels_low->Reset();

    h_vpintmichels_med_low->GetXaxis()->SetTitle("Integral (ADC)");
    h_vpintmichels_med_low->GetYaxis()->SetTitle("Counts");
    h_vpintmichels_med_low->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/jskweres/OutputFiles/plots/h_vpintmichels_med_low.png"));
    h_vpintmichels_med_low->Reset();

    h_vpintmichels_med->GetXaxis()->SetTitle("Integral (ADC)");
    h_vpintmichels_med->GetYaxis()->SetTitle("Counts");
    h_vpintmichels_med->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/jskweres/OutputFiles/plots/h_vpintmichels_med.png"));
    h_vpintmichels_med->Reset();

    h_vpintmichels_med_high->GetXaxis()->SetTitle("Integral (ADC)");
    h_vpintmichels_med_high->GetYaxis()->SetTitle("Counts");
    h_vpintmichels_med_high->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/jskweres/OutputFiles/plots/h_vpintmichels_med_high.png"));
    h_vpintmichels_med_high->Reset();

    h_vpintmichels_high->GetXaxis()->SetTitle("Integral (ADC)");
    h_vpintmichels_high->GetYaxis()->SetTitle("Counts");
    h_vpintmichels_high->Draw();
    gPad->Print(Form("/data7/coherent/data/d2o/jskweres/OutputFiles/plots/h_vpintmichels_high.png"));
    h_vpintmichels_high->Reset();*/

    // Close the file
    VMCFile.close();

    VMCFile.clear();
    
    return 0;

}