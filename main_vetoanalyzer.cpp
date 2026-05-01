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
struct pulse {
  int file_num;             /* Entry # from inputFile TTree T */
  int entry;                /* Entry # from inputFile TTree T */
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
  long long last_muon_time; /* Time value of most recent detected muon event */
  long long last_lead_time; /* Time value of most recent detected muon stopping in lead event */
  long long last_det_time;  /* Time value of most recent detected detector event */
  bool issue;               /* Flag to keep track of unusual eventTree entries */
  long long ev61_start_time;/* Start time of Event 61 pulse */
  bool ev61_significant;    /* Flag to keep track of whether Event 61 event is also potentially significant */
  long long wf_time;         /* nsTime of event */
  long long vp_start_16;    /* Veto panel 16 start time */
  long long vp_start_17;    /* Veto panel 17 start time */
  long long vp_start_18;    /* Veto panel 18 start time */
  long long vp_start_19;    /* Veto panel 19 start time */
  long long vp_start_20;    /* Veto panel 20 start time */
  long long vp_start_21;    /* Veto panel 21 start time */
  long long vp_start_22;    /* Veto panel 22 start time */
  long long vp_start_23;    /* Veto panel 23 start time */
  long long vp_start_24;    /* Veto panel 24 start time */
  long long vp_start_25;    /* Veto panel 25 start time */
  long long vp_start_top;   /* Top veto panel start time */
  long long vp_int_16;      /* Veto panel 16 integral value */
  long long vp_int_17;      /* Veto panel 17 integral value */
  long long vp_int_18;      /* Veto panel 18 integral value */
  long long vp_int_19;      /* Veto panel 19 integral value */
  long long vp_int_20;      /* Veto panel 20 integral value */
  long long vp_int_21;      /* Veto panel 21 integral value */
  long long vp_int_22;      /* Veto panel 22 integral value */
  long long vp_int_23;      /* Veto panel 23 integral value */
  long long vp_int_24;      /* Veto panel 24 integral value */
  long long vp_int_25;      /* Veto panel 25 integral value */
  long long vp_int_top;     /* Top veto panel integral value */
  double vp_ttp_amp_ratio;  /* Waveform pulse tail-to-peak amplitude ratio */
};

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

    /* Initialize vectors and variables */

    std::vector<struct pulse> pulses;

    long long sample_range = 3 * pow(10, 4);            // pow(10, 4);

    bool bt_HL_n_LL = false;

    int num_events_low_dt = 0;

    int num_events_low_dt_cent_spike = 0;

    int num_events_low_dt_no_spike_neg = 0;

    int num_events_low_dt_no_spike_pos = 0;

    int take_sample_count = 0;

    int hist_events = 0;

    int BRN_count = 0;

    int background_count = 0;

    int counter = 0;

    /* Initialize histograms */                 // NUMBER OF BINS SHOULD BE A MULTIPLE OF THE SMALLEST AXIS UNIT!!! (OR LENGTH OF AXIS DIVIDED BY SMALLEST AXIS UNIT?)

    TH1D* h_dt_61 = new TH1D("h_dt_61", "Delta-t Between Event 61 and Detector Peaks", 1250, - sample_range, -10000);                  // 1 bin = 16 ns           // 550, 1200, sample_range);           // 1250, - sample_range, sample_range);

    TH1D* h_dt_61_16 = new TH1D("h_dt_61_16", "Delta-t Between Event 61 and SiPM 1 Peaks", 1250, - sample_range, -10000);              // 1 bin = 16 ns           // 562, 1000, 9992);

    TH1D* h_dt_61_17 = new TH1D("h_dt_61_17", "Delta-t Between Event 61 and SiPM 2 Peaks", 1250, - sample_range, -10000);              // 1 bin = 16 ns

    TH1D* h_dt_61_18 = new TH1D("h_dt_61_18", "Delta-t Between Event 61 and SiPM 3 Peaks", 1250, - sample_range, -10000);              // 1 bin = 16 ns

    TH1D* h_dt_61_19 = new TH1D("h_dt_61_19", "Delta-t Between Event 61 and SiPM 4 Peaks", 1250, - sample_range, -10000);              // 1 bin = 16 ns

    TH1D* h_dt_61_20 = new TH1D("h_dt_61_20", "Delta-t Between Event 61 and SiPM 5 Peaks", 1250, - sample_range, -10000);              // 1 bin = 16 ns

    TH1D* h_dt_61_21 = new TH1D("h_dt_61_21", "Delta-t Between Event 61 and SiPM 6 Peaks", 1250, - sample_range, -10000);              // 1 bin = 16 ns

    TH1D* h_dt_61_22 = new TH1D("h_dt_61_22", "Delta-t Between Event 61 and SiPM 7 Peaks", 1250, - sample_range, -10000);              // 1 bin = 16 ns

    TH1D* h_dt_61_23 = new TH1D("h_dt_61_23", "Delta-t Between Event 61 and SiPM 8 Peaks", 1250, - sample_range, -10000);              // 1 bin = 16 ns

    TH1D* h_dt_61_24 = new TH1D("h_dt_61_24", "Delta-t Between Event 61 and SiPM 9 Peaks", 1250, - sample_range, -10000);              // 1 bin = 16 ns

    TH1D* h_dt_61_25 = new TH1D("h_dt_61_25", "Delta-t Between Event 61 and SiPM 10 Peaks", 1250, - sample_range, -10000);             // 1 bin = 16 ns

    TH1D* h_dt_61_top = new TH1D("h_dt_61_top", "Delta-t Between Event 61 and SiPM 9 & 10 Peaks", 1250, - sample_range, -10000);       // 1 bin = 16 ns

    TH1D* h_vpint_16 = new TH1D("h_vpint_16", "Low dt SiPM 1 Integral Values", 275, -400, 4000);                                             // 1 bin = 16 ns

    TH1D* h_vpint_17 = new TH1D("h_vpint_17", "Low dt SiPM 2 Integral Values", 275, -400, 4000);                                             // 1 bin = 16 ns

    TH1D* h_vpint_18 = new TH1D("h_vpint_18", "Low dt SiPM 3 Integral Values", 275, -400, 4000);                                             // 1 bin = 16 ns

    TH1D* h_vpint_19 = new TH1D("h_vpint_19", "Low dt SiPM 4 Integral Values", 275, -400, 4000);                                             // 1 bin = 16 ns

    TH1D* h_vpint_20 = new TH1D("h_vpint_20", "Low dt SiPM 5 Integral Values", 275, -400, 4000);                                             // 1 bin = 16 ns

    TH1D* h_vpint_21 = new TH1D("h_vpint_21", "Low dt SiPM 6 Integral Values", 275, -400, 4000);                                             // 1 bin = 16 ns

    TH1D* h_vpint_22 = new TH1D("h_vpint_22", "Low dt SiPM 7 Integral Values", 275, -400, 4000);                                             // 1 bin = 16 ns

    TH1D* h_vpint_23 = new TH1D("h_vpint_23", "Low dt SiPM 8 Integral Values", 275, -400, 4000);                                             // 1 bin = 16 ns

    TH1D* h_vpint_24 = new TH1D("h_vpint_24", "Low dt SiPM 9 Integral Values", 275, -400, 4000);                                             // 1 bin = 16 ns

    TH1D* h_vpint_25 = new TH1D("h_vpint_25", "Low dt SiPM 10 Integral Values", 275, -400, 4000);                                            // 1 bin = 16 ns

    TH1D* h_vpint_top = new TH1D("h_vpint_top", "Low dt SiPM 9 & 10 Integral Values", 275, -400, 4000);                                      // 1 bin = 16 ns

    TH2D* h_dt_v_detint = new TH2D("h_dt_v_detint", "Delta-T vs Detector Integral Value", 1250, - sample_range, -10000, 250, 0, 4000);                                            // 1 bin = 16 ns

    TH2D* h_dt_v_tmuon = new TH2D("h_dt_v_tmuon", "Delta-T vs Most Recent Muon Time Difference", 1250, - sample_range, -10000, 1000, 0, 2 * pow(10, 7));                          // 1 bin = 16 ns

    TH2D* h_dt_v_tlead = new TH2D("h_dt_v_tlead", "Delta-T vs Most Recent Muon Stopping in Lead Time Difference", 1250, - sample_range, -10000, 1000, 0, 2 * pow(10, 7));         // 1 bin = 32 ns

    TH1D* h_int = new TH1D("h_int", "Distribution of All Event Integral Values Between High Light and Low Light LED Events", 200, 0, 2000);

    TH1D* h_LLint = new TH1D("h_LLint", "Distribution of All Low Light LED Integral Values", 200, 0, 2000);

    TH1D* h_HLint = new TH1D("h_HLint", "Distribution of All High Light LED Integral Values", 200, 0, 2000);

    TH1D* h_eventint = new TH1D("h_eventint", "Distribution of All Internally Triggered Detector Integral Values", 200, 0, 2000);

    TH1D* h_eventvpint = new TH1D("h_eventvpint", "Distribution of All Internally Triggered Veto Panel Integral Values", 200, 0, 2000);

    TH1D* h_ev61int = new TH1D("h_ev61int", "Distribution of All Event 61 Detector Integral Values", 200, 0, 2000);

    TH1D* h_ev61vpint = new TH1D("h_ev61vpint", "Distribution of All Event 61 Veto Panel Integral Values", 200, 0, 2000);

    TH1D* h_34detint = new TH1D("h_34detint", "Distribution of All triggerBits = 34 Detector Integral Values", 200, 0, 2000);

    TH1D* h_34vpint = new TH1D("h_34vpint", "Distribution of All triggerBits = 34 Veto Panel Integral Values", 200, 0, 2000);

    TH1D* h_detint = new TH1D("h_detint", "Distribution of All Internally Triggered Detector Integral Values", 200, 0, 2000);

    TH1D* h_detint_low_dt = new TH1D("h_detint_low_dt", "Distribution of Low dt Detector Integral Values", 80, 0, 800);

    TH1D* h_int_cent_spike = new TH1D("h_int_cent_spike", "Distribution of All Central Spike Detector Integral Values", 200, 0, 2000);

    TH1D* h_int_pos_low_dt = new TH1D("h_int_pos_low_dt", "Distribution of Low dt Detector Integral Values", 200, 0, 2000);

    TH1D* h_int_neg_low_dt = new TH1D("h_int_neg_low_dt", "Distribution of Low dt Detector Integral Values", 200, 0, 2000);

    TH1D* h_avpint_low_dt = new TH1D("h_avpint_low_dt", "Distribution of Low dt Veto Panel Integral Values", 240, -400, 2000);

    TH1D* h_avpint_cent_spike = new TH1D("h_avpint_cent_spike", "Distribution of All Central Spike Veto Panel Integral Values", 240, -400, 2000);

    TH1D* h_avpint_pos_low_dt = new TH1D("h_avpint_pos_low_dt", "Distribution of Low dt Veto Panel Integral Values", 240, -400, 2000);

    TH1D* h_avpint_neg_low_dt = new TH1D("h_avpint_neg_low_dt", "Distribution of Low dt Veto Panel Integral Values", 240, -400, 2000);

    TH1D* h_svpint_cent_spike = new TH1D("h_svpint_cent_spike", "Distribution of All Central Spike Side Veto Panel Integral Values", 275, -400, 4000);

    TH1D* h_svpint_pos_low_dt = new TH1D("h_svpint_pos_low_dt", "Distribution of Low dt Side Veto Panel Integral Values", 275, -400, 4000);

    TH1D* h_svpint_neg_low_dt = new TH1D("h_svpint_neg_low_dt", "Distribution of Low dt Side Veto Panel Integral Values", 275, -400, 4000);

    TH1D* h_tvpint_cent_spike = new TH1D("h_tvpint_cent_spike", "Distribution of All Central Spike Top Veto Panel Integral Values", 275, -400, 4000);

    TH1D* h_tvpint_pos_low_dt = new TH1D("h_tvpint_pos_low_dt", "Distribution of Low dt Top Veto Panel Integral Values", 275, -400, 4000);

    TH1D* h_tvpint_neg_low_dt = new TH1D("h_tvpint_neg_low_dt", "Distribution of Low dt Top Veto Panel Integral Values", 275, -400, 4000);

    TH1D* h_vec_size = new TH1D("h_vec_size", "Distribution of Number of Internally Triggered Events Between High Light and Low Light LED Events", 15, 0, 15);

    TH1D* h_vpint_16_BRNpeak = new TH1D("h_vpint_16_BRNpeak", "BRN Peak SiPM 1 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_17_BRNpeak = new TH1D("h_vpint_17_BRNpeak", "BRN Peak SiPM 2 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_18_BRNpeak = new TH1D("h_vpint_18_BRNpeak", "BRN Peak SiPM 3 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_19_BRNpeak = new TH1D("h_vpint_19_BRNpeak", "BRN Peak SiPM 4 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_20_BRNpeak = new TH1D("h_vpint_20_BRNpeak", "BRN Peak SiPM 5 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_21_BRNpeak = new TH1D("h_vpint_21_BRNpeak", "BRN Peak SiPM 6 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_22_BRNpeak = new TH1D("h_vpint_22_BRNpeak", "BRN Peak SiPM 7 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_23_BRNpeak = new TH1D("h_vpint_23_BRNpeak", "BRN Peak SiPM 8 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_top_BRNpeak = new TH1D("h_vpint_top_BRNpeak", "BRN Peak SiPM 9 & 10 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_16_background = new TH1D("h_vpint_16_background", "Background SiPM 1 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_17_background = new TH1D("h_vpint_17_background", "Background SiPM 2 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_18_background = new TH1D("h_vpint_18_background", "Background SiPM 3 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_19_background = new TH1D("h_vpint_19_background", "Background SiPM 4 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_20_background = new TH1D("h_vpint_20_background", "Background SiPM 5 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_21_background = new TH1D("h_vpint_21_background", "Background SiPM 6 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_22_background = new TH1D("h_vpint_22_background", "Background SiPM 7 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_23_background = new TH1D("h_vpint_23_background", "Background SiPM 8 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_top_background = new TH1D("h_vpint_top_background", "Background SiPM 9 & 10 Integral Values", 275, -400, 4000);

    TH1D* h_BRNint_16 = new TH1D("h_BRNint_16", "SiPM 1 ONLY BRN Integral Values", 275, -400, 4000);

    TH1D* h_BRNint_17 = new TH1D("h_BRNint_17", "SiPM 2 ONLY BRN Integral Values", 275, -400, 4000);

    TH1D* h_BRNint_18 = new TH1D("h_BRNint_18", "SiPM 3 ONLY BRN Integral Values", 275, -400, 4000);

    TH1D* h_BRNint_19 = new TH1D("h_BRNint_19", "SiPM 4 ONLY BRN Integral Values", 275, -400, 4000);

    TH1D* h_BRNint_20 = new TH1D("h_BRNint_20", "SiPM 5 ONLY BRN Integral Values", 275, -400, 4000);

    TH1D* h_BRNint_21 = new TH1D("h_BRNint_21", "SiPM 6 ONLY BRN Integral Values", 275, -400, 4000);

    TH1D* h_BRNint_22 = new TH1D("h_BRNint_22", "SiPM 7 ONLY BRN Integral Values", 275, -400, 4000);

    TH1D* h_BRNint_23 = new TH1D("h_BRNint_23", "SiPM 8 ONLY BRN Integral Values", 275, -400, 4000);

    TH1D* h_BRNint_top = new TH1D("h_BRNint_top", "SiPM 9 & 10 ONLY BRN Integral Values", 275, -400, 4000);

    TH1D* h_svpint_not_16 = new TH1D("h_svpint_not_16", "Side Veto Panel Integral Values, Excluding SiPM1", 275, -400, 4000);

    TH1D* h_svpint_not_17 = new TH1D("h_svpint_not_17", "Side Veto Panel Integral Values, Excluding SiPM2", 275, -400, 4000);

    TH1D* h_svpint_not_18 = new TH1D("h_svpint_not_18", "Side Veto Panel Integral Values, Excluding SiPM3", 275, -400, 4000);

    TH1D* h_svpint_not_19 = new TH1D("h_svpint_not_19", "Side Veto Panel Integral Values, Excluding SiPM4", 275, -400, 4000);

    TH1D* h_svpint_not_20 = new TH1D("h_svpint_not_20", "Side Veto Panel Integral Values, Excluding SiPM5", 275, -400, 4000);

    TH1D* h_svpint_not_21 = new TH1D("h_svpint_not_21", "Side Veto Panel Integral Values, Excluding SiPM6", 275, -400, 4000);

    TH1D* h_svpint_not_22 = new TH1D("h_svpint_not_22", "Side Veto Panel Integral Values, Excluding SiPM7", 275, -400, 4000);

    TH1D* h_svpint_not_23 = new TH1D("h_svpint_not_23", "Side Veto Panel Integral Values, Excluding SiPM8", 275, -400, 4000);

    TH1D* h_tvpint_not_top = new TH1D("h_tvpint_not_top", "Top Veto Panel Integral Values, Excluding SiPM9&10", 275, -400, 4000);

    TH1D* h_vpint_16_all = new TH1D("h_vpint_16_all", "Low dt SiPM 1 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_17_all = new TH1D("h_vpint_17_all", "Low dt SiPM 2 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_18_all = new TH1D("h_vpint_18_all", "Low dt SiPM 3 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_19_all = new TH1D("h_vpint_19_all", "Low dt SiPM 4 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_20_all = new TH1D("h_vpint_20_all", "Low dt SiPM 5 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_21_all = new TH1D("h_vpint_21_all", "Low dt SiPM 6 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_22_all = new TH1D("h_vpint_22_all", "Low dt SiPM 7 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_23_all = new TH1D("h_vpint_23_all", "Low dt SiPM 8 Integral Values", 275, -400, 4000);

    TH1D* h_vpint_top_all = new TH1D("h_vpint_top_all", "Low dt SiPM 9 & 10 Integral Values", 275, -400, 4000);

    TH1D* h_BRN_tail_amp = new TH1D("h_BRN_tail_amp", "BRN Pulse Tail-to-Peak Amplitude Ratio", 1000, 0, 1);

    TH1D* h_background_tail_amp = new TH1D("h_background_tail_amp", "Background Pulse Tail-to-Peak Amplitude Ratio", 1000, 0, 1);

    TFile *f; 
    //your root file location here
        if(gSystem->AccessPathName(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/int_trig_events/Veto_Batch_Analysis_Run%i.root", run))){
        //    if (gSystem->AccessPathName(Form("/data9/coherent/data/d2o/processedData/Veto_Analysis_24hrs_Run%i.root", run))) {
        //    if(gSystem->AccessPathName(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/Veto_Analysis_24hrs_Run%i.root", run))){
        cout << "Could not open file" << endl;
        return -1; 
    } else{
        f = new TFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/int_trig_events/Veto_Batch_Analysis_Run%i.root", run));
        //        f = new TFile(Form("/data9/coherent/data/d2o/processedData/Veto_Analysis_24hrs_Run%i.root", run));
        //        f = new TFile(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/Veto_Analysis_24hrs_Run%i.root", run));
    }

    TTree* t = (TTree*)f->Get("eventTree");

    // Declaration of leaf types
    Int_t           fnum;
    Int_t           entry;
    Long64_t        st;
    Long64_t        ed;
    Double_t        pk;
    Double_t        ey;
    Int_t           nr;
    Bool_t          se;
    Bool_t          bm;
    Int_t           tr;
    Int_t           lh;
    Double_t        svpe;
    Double_t        tvpe;
    Double_t        avpe;
    Long64_t        lmt;
    Long64_t        llt;
    Long64_t        ldt;
    Bool_t          issue;
    Long64_t        e61st;
    Bool_t          e61sig;
    Long64_t        wftime;
    Long64_t        vptime16;
    Long64_t        vptime17;
    Long64_t        vptime18;
    Long64_t        vptime19;
    Long64_t        vptime20;
    Long64_t        vptime21;
    Long64_t        vptime22;
    Long64_t        vptime23;
    Long64_t        vptime24;
    Long64_t        vptime25;
    Long64_t        vptimetop;
    Long64_t        vpint16;
    Long64_t        vpint17;
    Long64_t        vpint18;
    Long64_t        vpint19;
    Long64_t        vpint20;
    Long64_t        vpint21;
    Long64_t        vpint22;
    Long64_t        vpint23;
    Long64_t        vpint24;
    Long64_t        vpint25;
    Long64_t        vpinttop;
    Double_t        vpar;

    // List of branches
    TBranch        *b_fnum;
    TBranch        *b_entry;
    TBranch        *b_st;
    TBranch        *b_ed;
    TBranch        *b_pk;
    TBranch        *b_ey;
    TBranch        *b_nr;
    TBranch        *b_se;
    TBranch        *b_bm;
    TBranch        *b_tr;
    TBranch        *b_lh;
    TBranch        *b_svpe;
    TBranch        *b_tvpe;
    TBranch        *b_avpe;
    TBranch        *b_lmt;
    TBranch        *b_llt;
    TBranch        *b_ldt;
    TBranch        *b_issue;
    TBranch        *b_e61st;
    TBranch        *b_e61sig;
    TBranch        *b_wftime;
    TBranch        *b_vptime16;
    TBranch        *b_vptime17;
    TBranch        *b_vptime18;
    TBranch        *b_vptime19;
    TBranch        *b_vptime20;
    TBranch        *b_vptime21;
    TBranch        *b_vptime22;
    TBranch        *b_vptime23;
    TBranch        *b_vptime24;
    TBranch        *b_vptime25;
    TBranch        *b_vptimetop;
    TBranch        *b_vpint16;
    TBranch        *b_vpint17;
    TBranch        *b_vpint18;
    TBranch        *b_vpint19;
    TBranch        *b_vpint20;
    TBranch        *b_vpint21;
    TBranch        *b_vpint22;
    TBranch        *b_vpint23;
    TBranch        *b_vpint24;
    TBranch        *b_vpint25;
    TBranch        *b_vpinttop;
    TBranch        *b_vpar;

    t->SetBranchAddress("fnum", &fnum, &b_fnum);
    t->SetBranchAddress("entry", &entry, &b_entry);
    t->SetBranchAddress("st", &st, &b_st);
    t->SetBranchAddress("ed", &ed, &b_ed);
    t->SetBranchAddress("pk", &pk, &b_pk);
    t->SetBranchAddress("ey", &ey, &b_ey);
    t->SetBranchAddress("nr", &nr, &b_nr);
    t->SetBranchAddress("se", &se, &b_se);
    t->SetBranchAddress("bm", &bm, &b_bm);
    t->SetBranchAddress("tr", &tr, &b_tr);
    t->SetBranchAddress("lh", &lh, &b_lh);
    t->SetBranchAddress("svpe", &svpe, &b_svpe);
    t->SetBranchAddress("tvpe", &tvpe, &b_tvpe);
    t->SetBranchAddress("avpe", &avpe, &b_avpe);
    t->SetBranchAddress("lmt", &lmt, &b_lmt);
    t->SetBranchAddress("llt", &llt, &b_llt);
    t->SetBranchAddress("ldt", &ldt, &b_ldt);
    t->SetBranchAddress("issue", &issue, &b_issue);
    t->SetBranchAddress("e61st", &e61st, &b_e61st);
    t->SetBranchAddress("e61sig", &e61sig, &b_e61sig);
    t->SetBranchAddress("wftime", &wftime, &b_wftime);
    t->SetBranchAddress("vptime16", &vptime16, &b_vptime16);
    t->SetBranchAddress("vptime17", &vptime17, &b_vptime17);
    t->SetBranchAddress("vptime18", &vptime18, &b_vptime18);
    t->SetBranchAddress("vptime19", &vptime19, &b_vptime19);
    t->SetBranchAddress("vptime20", &vptime20, &b_vptime20);
    t->SetBranchAddress("vptime21", &vptime21, &b_vptime21);
    t->SetBranchAddress("vptime22", &vptime22, &b_vptime22);
    t->SetBranchAddress("vptime23", &vptime23, &b_vptime23);
    t->SetBranchAddress("vptime24", &vptime24, &b_vptime24);
    t->SetBranchAddress("vptime25", &vptime25, &b_vptime25);
    t->SetBranchAddress("vptimetop", &vptimetop, &b_vptimetop);
    t->SetBranchAddress("vpint16", &vpint16, &b_vpint16);
    t->SetBranchAddress("vpint17", &vpint17, &b_vpint17);
    t->SetBranchAddress("vpint18", &vpint18, &b_vpint18);
    t->SetBranchAddress("vpint19", &vpint19, &b_vpint19);
    t->SetBranchAddress("vpint20", &vpint20, &b_vpint20);
    t->SetBranchAddress("vpint21", &vpint21, &b_vpint21);
    t->SetBranchAddress("vpint22", &vpint22, &b_vpint22);
    t->SetBranchAddress("vpint23", &vpint23, &b_vpint23);
    t->SetBranchAddress("vpint24", &vpint24, &b_vpint24);
    t->SetBranchAddress("vpint25", &vpint25, &b_vpint25);
    t->SetBranchAddress("vpinttop", &vpinttop, &b_vpinttop);
    t->SetBranchAddress("vpar", &vpar, &b_vpar);

    for(int iEvent = 0; iEvent < t->GetEntries(); iEvent++){

        Long64_t tentry = t->LoadTree(iEvent);

        b_fnum->GetEntry(tentry);
        b_entry->GetEntry(tentry);
        b_st->GetEntry(tentry);
        b_ed->GetEntry(tentry);
        b_pk->GetEntry(tentry);
        b_ey->GetEntry(tentry);
        b_nr->GetEntry(tentry);
        b_se->GetEntry(tentry);
        b_bm->GetEntry(tentry);
        b_tr->GetEntry(tentry);
        b_lh->GetEntry(tentry);
        b_svpe->GetEntry(tentry);
        b_tvpe->GetEntry(tentry);
        b_avpe->GetEntry(tentry);
        b_lmt->GetEntry(tentry);
        b_llt->GetEntry(tentry);
        b_ldt->GetEntry(tentry);
        b_issue->GetEntry(tentry);
        b_e61st->GetEntry(tentry);
        b_e61sig->GetEntry(tentry);
        b_wftime->GetEntry(tentry);
        b_vptime16->GetEntry(tentry);
        b_vptime17->GetEntry(tentry);
        b_vptime18->GetEntry(tentry);
        b_vptime19->GetEntry(tentry);
        b_vptime20->GetEntry(tentry);
        b_vptime21->GetEntry(tentry);
        b_vptime22->GetEntry(tentry);
        b_vptime23->GetEntry(tentry);
        b_vptime24->GetEntry(tentry);
        b_vptime25->GetEntry(tentry);
        b_vptimetop->GetEntry(tentry);
        b_vpint16->GetEntry(tentry);
        b_vpint17->GetEntry(tentry);
        b_vpint18->GetEntry(tentry);
        b_vpint19->GetEntry(tentry);
        b_vpint20->GetEntry(tentry);
        b_vpint21->GetEntry(tentry);
        b_vpint22->GetEntry(tentry);
        b_vpint23->GetEntry(tentry);
        b_vpint24->GetEntry(tentry);
        b_vpint25->GetEntry(tentry);
        b_vpinttop->GetEntry(tentry);
        b_vpar->GetEntry(tentry);

        struct pulse avg_pulse;

        avg_pulse.file_num = fnum;
        avg_pulse.entry = entry;
        avg_pulse.start = st;
        avg_pulse.end = ed;
        avg_pulse.peak = pk;
        avg_pulse.energy = ey;
        avg_pulse.number = nr;
        avg_pulse.single = se;
        avg_pulse.beam = bm;
        avg_pulse.trigger = tr;
        avg_pulse.length = lh;
        avg_pulse.side_vp_energy = svpe;
        avg_pulse.top_vp_energy = tvpe;
        avg_pulse.all_vp_energy = avpe;
        avg_pulse.last_muon_time = lmt;
        avg_pulse.last_lead_time = llt;
        avg_pulse.last_det_time = ldt;
        avg_pulse.issue = issue;
        avg_pulse.ev61_start_time = e61st;
        avg_pulse.ev61_significant = e61sig;
        avg_pulse.wf_time = wftime;
        avg_pulse.vp_start_16 = vptime16;
        avg_pulse.vp_start_17 = vptime17;
        avg_pulse.vp_start_18 = vptime18;
        avg_pulse.vp_start_19 = vptime19;
        avg_pulse.vp_start_20 = vptime20;
        avg_pulse.vp_start_21 = vptime21;
        avg_pulse.vp_start_22 = vptime22;
        avg_pulse.vp_start_23 = vptime23;
        avg_pulse.vp_start_24 = vptime24;
        avg_pulse.vp_start_25 = vptime25;
        avg_pulse.vp_start_top = vptimetop;
        avg_pulse.vp_int_16 = vpint16;
        avg_pulse.vp_int_17 = vpint17;
        avg_pulse.vp_int_18 = vpint18;
        avg_pulse.vp_int_19 = vpint19;
        avg_pulse.vp_int_20 = vpint20;
        avg_pulse.vp_int_21 = vpint21;
        avg_pulse.vp_int_22 = vpint22;
        avg_pulse.vp_int_23 = vpint23;
        avg_pulse.vp_int_24 = vpint24;
        avg_pulse.vp_int_25 = vpint25;
        avg_pulse.vp_int_top = vpinttop;
        avg_pulse.vp_ttp_amp_ratio = vpar;

        /*
        
        Only looking for events between high light LED events and low light LED events
        
        */

        // After every high light LED event, "turn on" internally triggered event storing

        if (avg_pulse.trigger == 8) {bt_HL_n_LL = true; h_HLint->Fill(avg_pulse.energy);}

        // Store internally triggered events in vector

        if (bt_HL_n_LL && avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {pulses.push_back(avg_pulse); h_eventint->Fill(avg_pulse.energy); h_eventvpint->Fill(avg_pulse.all_vp_energy);}
        
        // After every low light LED event, "turn off" internally triggered event storing and search through compiled vector of events

        if (avg_pulse.trigger == 16) {

            // Iterate over all events in compiled vector

            for (size_t iVec = 0; iVec < pulses.size(); iVec++) {

                if (pulses[iVec].beam) {

                    take_sample_count += 1;

                    if (iVec > 0 && (pulses[iVec - 1].trigger != 1 && pulses[iVec - 1].trigger != 33 && pulses[iVec - 1].trigger != 35)) {

                        if ((pulses[iVec].ev61_start_time - pulses[iVec - 1].wf_time) <= sample_range && (pulses[iVec].ev61_start_time - pulses[iVec - 1].wf_time) > 0) {

                            // Plot negative dt values

                            if ((pulses[iVec - 1].start - pulses[iVec].ev61_start_time) >= -10000 && (pulses[iVec - 1].start - pulses[iVec].ev61_start_time) <= 10000) {num_events_low_dt += 1; h_detint_low_dt->Fill(pulses[iVec - 1].energy); h_avpint_low_dt->Fill(pulses[iVec - 1].all_vp_energy);}
                
                            if ((pulses[iVec - 1].start - pulses[iVec].ev61_start_time) > -500 && (pulses[iVec - 1].start - pulses[iVec].ev61_start_time) < 500) {num_events_low_dt_cent_spike += 1; h_int_cent_spike->Fill(pulses[iVec - 1].energy); h_avpint_cent_spike->Fill(pulses[iVec - 1].all_vp_energy); h_svpint_cent_spike->Fill(pulses[iVec - 1].side_vp_energy); h_tvpint_cent_spike->Fill(pulses[iVec - 1].top_vp_energy);}          // cout << "\n" << "Run, Event, triggerBits = " << iRun << " " << iEnt + 1 << " " << avg_pulse.trigger << endl;}

                            else if ((pulses[iVec - 1].start - pulses[iVec].ev61_start_time) <= -500 && (pulses[iVec - 1].start - pulses[iVec].ev61_start_time) >= -10000) {num_events_low_dt_no_spike_neg += 1; h_int_neg_low_dt->Fill(pulses[iVec - 1].energy); h_avpint_neg_low_dt->Fill(pulses[iVec - 1].all_vp_energy); h_svpint_neg_low_dt->Fill(pulses[iVec - 1].side_vp_energy); h_tvpint_neg_low_dt->Fill(pulses[iVec - 1].top_vp_energy);}

                            else if ((pulses[iVec - 1].start - pulses[iVec].ev61_start_time) >= 500 && (pulses[iVec - 1].start - pulses[iVec].ev61_start_time) <= 10000) {num_events_low_dt_no_spike_pos += 1; h_int_pos_low_dt->Fill(pulses[iVec - 1].energy); h_avpint_pos_low_dt->Fill(pulses[iVec - 1].all_vp_energy); h_svpint_pos_low_dt->Fill(pulses[iVec - 1].side_vp_energy); h_tvpint_pos_low_dt->Fill(pulses[iVec - 1].top_vp_energy);}

                            hist_events += 1;
                            
                            // h_dt_61->Fill(pulses[iVec - 1].start - pulses[iVec].ev61_start_time);

                            if (pulses[iVec - 1].vp_start_16 != 0) {h_dt_61_16->Fill(pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_16 - pulses[iVec].ev61_start_time); h_vpint_16->Fill(pulses[iVec - 1].vp_int_16); if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_16 - pulses[iVec].ev61_start_time > -19500 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_16 - pulses[iVec].ev61_start_time < -18500) {h_vpint_16_BRNpeak->Fill(pulses[iVec - 1].vp_int_16);} else if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_16 - pulses[iVec].ev61_start_time > -15000 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_16 - pulses[iVec].ev61_start_time < -5000) {h_vpint_16_background->Fill(pulses[iVec - 1].vp_int_16);}}

                            if (pulses[iVec - 1].vp_start_17 != 0) {h_dt_61_17->Fill(pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_17 - pulses[iVec].ev61_start_time); h_vpint_17->Fill(pulses[iVec - 1].vp_int_17); if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_17 - pulses[iVec].ev61_start_time > -19500 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_17 - pulses[iVec].ev61_start_time < -18500) {h_vpint_17_BRNpeak->Fill(pulses[iVec - 1].vp_int_17);} else if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_17 - pulses[iVec].ev61_start_time > -15000 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_17 - pulses[iVec].ev61_start_time < -5000) {h_vpint_17_background->Fill(pulses[iVec - 1].vp_int_17);}}

                            if (pulses[iVec - 1].vp_start_18 != 0) {h_dt_61_18->Fill(pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_18 - pulses[iVec].ev61_start_time); h_vpint_18->Fill(pulses[iVec - 1].vp_int_18); if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_18 - pulses[iVec].ev61_start_time > -19500 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_18 - pulses[iVec].ev61_start_time < -18500) {h_vpint_18_BRNpeak->Fill(pulses[iVec - 1].vp_int_18);} else if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_18 - pulses[iVec].ev61_start_time > -15000 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_18 - pulses[iVec].ev61_start_time < -5000) {h_vpint_18_background->Fill(pulses[iVec - 1].vp_int_18);}}

                            if (pulses[iVec - 1].vp_start_19 != 0) {h_dt_61_19->Fill(pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_19 - pulses[iVec].ev61_start_time); h_vpint_19->Fill(pulses[iVec - 1].vp_int_19); if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_19 - pulses[iVec].ev61_start_time > -19500 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_19 - pulses[iVec].ev61_start_time < -18500) {h_vpint_19_BRNpeak->Fill(pulses[iVec - 1].vp_int_19);} else if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_19 - pulses[iVec].ev61_start_time > -15000 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_19 - pulses[iVec].ev61_start_time < -5000) {h_vpint_19_background->Fill(pulses[iVec - 1].vp_int_19);}}

                            if (pulses[iVec - 1].vp_start_20 != 0) {h_dt_61_20->Fill(pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_20 - pulses[iVec].ev61_start_time); h_vpint_20->Fill(pulses[iVec - 1].vp_int_20); if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_20 - pulses[iVec].ev61_start_time > -19500 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_20 - pulses[iVec].ev61_start_time < -18500) {h_vpint_20_BRNpeak->Fill(pulses[iVec - 1].vp_int_20);} else if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_20 - pulses[iVec].ev61_start_time > -15000 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_20 - pulses[iVec].ev61_start_time < -5000) {h_vpint_20_background->Fill(pulses[iVec - 1].vp_int_20);}}

                            if (pulses[iVec - 1].vp_start_21 != 0) {h_dt_61_21->Fill(pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_21 - pulses[iVec].ev61_start_time); h_vpint_21->Fill(pulses[iVec - 1].vp_int_21); if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_21 - pulses[iVec].ev61_start_time > -19500 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_21 - pulses[iVec].ev61_start_time < -18500) {h_vpint_21_BRNpeak->Fill(pulses[iVec - 1].vp_int_21);} else if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_21 - pulses[iVec].ev61_start_time > -15000 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_21 - pulses[iVec].ev61_start_time < -5000) {h_vpint_21_background->Fill(pulses[iVec - 1].vp_int_21);}}

                            if (pulses[iVec - 1].vp_start_22 != 0) {h_dt_61_22->Fill(pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_22 - pulses[iVec].ev61_start_time); h_vpint_22->Fill(pulses[iVec - 1].vp_int_22); if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_22 - pulses[iVec].ev61_start_time > -19500 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_22 - pulses[iVec].ev61_start_time < -18500) {h_vpint_22_BRNpeak->Fill(pulses[iVec - 1].vp_int_22);} else if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_22 - pulses[iVec].ev61_start_time > -15000 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_22 - pulses[iVec].ev61_start_time < -5000) {h_vpint_22_background->Fill(pulses[iVec - 1].vp_int_22);}}

                            if (pulses[iVec - 1].vp_start_23 != 0) {h_dt_61_23->Fill(pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_23 - pulses[iVec].ev61_start_time); h_vpint_23->Fill(pulses[iVec - 1].vp_int_23); if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_23 - pulses[iVec].ev61_start_time > -19500 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_23 - pulses[iVec].ev61_start_time < -18500) {h_vpint_23_BRNpeak->Fill(pulses[iVec - 1].vp_int_23);} else if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_23 - pulses[iVec].ev61_start_time > -15000 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_23 - pulses[iVec].ev61_start_time < -5000) {h_vpint_23_background->Fill(pulses[iVec - 1].vp_int_23);}}

                            if (pulses[iVec - 1].vp_start_24 != 0) {h_dt_61_24->Fill(pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_24 - pulses[iVec].ev61_start_time); h_vpint_24->Fill(pulses[iVec - 1].vp_int_24);}

                            if (pulses[iVec - 1].vp_start_25 != 0) {h_dt_61_25->Fill(pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_25 - pulses[iVec].ev61_start_time); h_vpint_25->Fill(pulses[iVec - 1].vp_int_25);}

                            if (pulses[iVec - 1].vp_start_top != 0) {h_dt_61_top->Fill(pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_top - pulses[iVec].ev61_start_time); h_vpint_top->Fill(pulses[iVec - 1].vp_int_top); if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_top - pulses[iVec].ev61_start_time > -19500 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_top - pulses[iVec].ev61_start_time < -18500) {h_vpint_top_BRNpeak->Fill(pulses[iVec - 1].vp_int_top);} else if (pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_top - pulses[iVec].ev61_start_time > -15000 && pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_top - pulses[iVec].ev61_start_time < -5000) {h_vpint_top_background->Fill(pulses[iVec - 1].vp_int_top);}}

                            if ((pulses[iVec - 1].vp_start_17 != 0 || pulses[iVec - 1].vp_start_18 != 0 || pulses[iVec - 1].vp_start_19 != 0 || pulses[iVec - 1].vp_start_20 != 0 || pulses[iVec - 1].vp_start_21 != 0 || pulses[iVec - 1].vp_start_22 != 0 || pulses[iVec - 1].vp_start_23 != 0 || pulses[iVec - 1].vp_start_top != 0) && pulses[iVec - 1].vp_start_16 == 0) {h_svpint_not_16->Fill(pulses[iVec - 1].vp_int_16);}         // .side_vp_energy);}

                            if ((pulses[iVec - 1].vp_start_16 != 0 || pulses[iVec - 1].vp_start_18 != 0 || pulses[iVec - 1].vp_start_19 != 0 || pulses[iVec - 1].vp_start_20 != 0 || pulses[iVec - 1].vp_start_21 != 0 || pulses[iVec - 1].vp_start_22 != 0 || pulses[iVec - 1].vp_start_23 != 0 || pulses[iVec - 1].vp_start_top != 0) && pulses[iVec - 1].vp_start_17 == 0) {h_svpint_not_17->Fill(pulses[iVec - 1].vp_int_17);}

                            if ((pulses[iVec - 1].vp_start_16 != 0 || pulses[iVec - 1].vp_start_17 != 0 || pulses[iVec - 1].vp_start_19 != 0 || pulses[iVec - 1].vp_start_20 != 0 || pulses[iVec - 1].vp_start_21 != 0 || pulses[iVec - 1].vp_start_22 != 0 || pulses[iVec - 1].vp_start_23 != 0 || pulses[iVec - 1].vp_start_top != 0) && pulses[iVec - 1].vp_start_18 == 0) {h_svpint_not_18->Fill(pulses[iVec - 1].vp_int_18);}

                            if ((pulses[iVec - 1].vp_start_16 != 0 || pulses[iVec - 1].vp_start_17 != 0 || pulses[iVec - 1].vp_start_18 != 0 || pulses[iVec - 1].vp_start_20 != 0 || pulses[iVec - 1].vp_start_21 != 0 || pulses[iVec - 1].vp_start_22 != 0 || pulses[iVec - 1].vp_start_23 != 0 || pulses[iVec - 1].vp_start_top != 0) && pulses[iVec - 1].vp_start_19 == 0) {h_svpint_not_19->Fill(pulses[iVec - 1].vp_int_19);}

                            if ((pulses[iVec - 1].vp_start_16 != 0 || pulses[iVec - 1].vp_start_17 != 0 || pulses[iVec - 1].vp_start_18 != 0 || pulses[iVec - 1].vp_start_19 != 0 || pulses[iVec - 1].vp_start_21 != 0 || pulses[iVec - 1].vp_start_22 != 0 || pulses[iVec - 1].vp_start_23 != 0 || pulses[iVec - 1].vp_start_top != 0) && pulses[iVec - 1].vp_start_20 == 0) {h_svpint_not_20->Fill(pulses[iVec - 1].vp_int_20);}
                            
                            if ((pulses[iVec - 1].vp_start_16 != 0 || pulses[iVec - 1].vp_start_17 != 0 || pulses[iVec - 1].vp_start_18 != 0 || pulses[iVec - 1].vp_start_19 != 0 || pulses[iVec - 1].vp_start_20 != 0 || pulses[iVec - 1].vp_start_22 != 0 || pulses[iVec - 1].vp_start_23 != 0 || pulses[iVec - 1].vp_start_top != 0) && pulses[iVec - 1].vp_start_21 == 0) {h_svpint_not_21->Fill(pulses[iVec - 1].vp_int_21);}

                            if ((pulses[iVec - 1].vp_start_16 != 0 || pulses[iVec - 1].vp_start_17 != 0 || pulses[iVec - 1].vp_start_18 != 0 || pulses[iVec - 1].vp_start_19 != 0 || pulses[iVec - 1].vp_start_20 != 0 || pulses[iVec - 1].vp_start_21 != 0 || pulses[iVec - 1].vp_start_23 != 0 || pulses[iVec - 1].vp_start_top != 0) && pulses[iVec - 1].vp_start_22 == 0) {h_svpint_not_22->Fill(pulses[iVec - 1].vp_int_22);}

                            if ((pulses[iVec - 1].vp_start_16 != 0 || pulses[iVec - 1].vp_start_17 != 0 || pulses[iVec - 1].vp_start_18 != 0 || pulses[iVec - 1].vp_start_19 != 0 || pulses[iVec - 1].vp_start_20 != 0 || pulses[iVec - 1].vp_start_21 != 0 || pulses[iVec - 1].vp_start_22 != 0 || pulses[iVec - 1].vp_start_top != 0) && pulses[iVec - 1].vp_start_23 == 0) {h_svpint_not_23->Fill(pulses[iVec - 1].vp_int_23);}

                            if ((pulses[iVec - 1].vp_start_16 != 0 || pulses[iVec - 1].vp_start_17 != 0 || pulses[iVec - 1].vp_start_18 != 0 || pulses[iVec - 1].vp_start_19 != 0 || pulses[iVec - 1].vp_start_20 != 0 || pulses[iVec - 1].vp_start_21 != 0 || pulses[iVec - 1].vp_start_22 != 0 || pulses[iVec - 1].vp_start_23 != 0) && pulses[iVec - 1].vp_start_top == 0) {h_tvpint_not_top->Fill(pulses[iVec - 1].vp_int_top);}

                            h_vpint_16_all->Fill(pulses[iVec - 1].vp_int_16);

                            h_vpint_17_all->Fill(pulses[iVec - 1].vp_int_17);

                            h_vpint_18_all->Fill(pulses[iVec - 1].vp_int_18); 

                            h_vpint_19_all->Fill(pulses[iVec - 1].vp_int_19); 

                            h_vpint_20_all->Fill(pulses[iVec - 1].vp_int_20);

                            h_vpint_21_all->Fill(pulses[iVec - 1].vp_int_21);

                            h_vpint_22_all->Fill(pulses[iVec - 1].vp_int_22);

                            h_vpint_23_all->Fill(pulses[iVec - 1].vp_int_23);

                            h_vpint_top_all->Fill(pulses[iVec - 1].vp_int_top);                             

                            h_dt_v_tmuon->Fill(pulses[iVec - 1].start - pulses[iVec].ev61_start_time, pulses[iVec - 1].start - pulses[iVec - 1].last_muon_time);

                            h_dt_v_tlead->Fill(pulses[iVec - 1].start - pulses[iVec].ev61_start_time, pulses[iVec - 1].start - pulses[iVec - 1].last_lead_time);

                            h_detint->Fill(pulses[iVec - 1].energy);
                            
                            h_dt_v_detint->Fill(pulses[iVec - 1].start - pulses[iVec].ev61_start_time, pulses[iVec - 1].energy);

                            if (pulses[iVec - 1].vp_int_21 > 1200 && counter < 1000) {
                                
                                counter += 1;

                                cout << "\n" << "NEGATIVE dt Value Event:" << endl;

                                cout << "\n" << "Ev61 Event: Run, Event Number, triggerBits, SiPM6 Integral Value = " << pulses[iVec].file_num << " " << pulses[iVec].entry << " " << pulses[iVec].trigger << " " << pulses[iVec].vp_int_21 << endl;

                                cout << "\n" << "VP Event: Run, Event Number, triggerBits, SiPM6 Integral Value, dt Value = " << pulses[iVec - 1].file_num << " " << pulses[iVec - 1].entry << " " << pulses[iVec - 1].trigger << " " << pulses[iVec - 1].vp_int_21 << " " << pulses[iVec - 1].wf_time + pulses[iVec - 1].vp_start_21 - pulses[iVec].ev61_start_time << endl;
                            
                            }

                        }

                    }

                    if (pulses[iVec].ev61_significant) {

                        // Plot low dt values

                        if ((pulses[iVec].start - pulses[iVec].ev61_start_time) >= -10000 && (pulses[iVec].start - pulses[iVec].ev61_start_time) <= 10000) {num_events_low_dt += 1; h_detint_low_dt->Fill(pulses[iVec].energy); h_avpint_low_dt->Fill(pulses[iVec].all_vp_energy);}
            
                        if ((pulses[iVec].start - pulses[iVec].ev61_start_time) > -500 && (pulses[iVec].start - pulses[iVec].ev61_start_time) < 500) {num_events_low_dt_cent_spike += 1; h_int_cent_spike->Fill(pulses[iVec].energy); h_avpint_cent_spike->Fill(pulses[iVec].all_vp_energy); h_svpint_cent_spike->Fill(pulses[iVec].side_vp_energy); h_tvpint_cent_spike->Fill(pulses[iVec].top_vp_energy);}          // cout << "\n" << "Run, Event, triggerBits = " << iRun << " " << iEnt + 1 << " " << avg_pulse.trigger << endl;}

                        else if ((pulses[iVec].start - pulses[iVec].ev61_start_time) <= -500 && (pulses[iVec].start - pulses[iVec].ev61_start_time) >= -10000) {num_events_low_dt_no_spike_neg += 1; h_int_neg_low_dt->Fill(pulses[iVec].energy); h_avpint_neg_low_dt->Fill(pulses[iVec].all_vp_energy); h_svpint_neg_low_dt->Fill(pulses[iVec].side_vp_energy); h_tvpint_neg_low_dt->Fill(pulses[iVec].top_vp_energy);}

                        else if ((pulses[iVec].start - pulses[iVec].ev61_start_time) >= 500 && (pulses[iVec].start - pulses[iVec].ev61_start_time) <= 10000) {num_events_low_dt_no_spike_pos += 1; h_int_pos_low_dt->Fill(pulses[iVec].energy); h_avpint_pos_low_dt->Fill(pulses[iVec].all_vp_energy); h_svpint_pos_low_dt->Fill(pulses[iVec].side_vp_energy); h_tvpint_pos_low_dt->Fill(pulses[iVec].top_vp_energy);}

                        hist_events += 1;
                        
                        // h_dt_61->Fill(pulses[iVec].start - pulses[iVec].ev61_start_time);
                        /*
                        if (pulses[iVec].vp_start_16 != 0) {h_dt_61_16->Fill(pulses[iVec].wf_time + pulses[iVec].vp_start_16 - pulses[iVec].ev61_start_time); h_vpint_16->Fill(pulses[iVec].vp_int_16);}

                        if (pulses[iVec].vp_start_17 != 0) {h_dt_61_17->Fill(pulses[iVec].wf_time + pulses[iVec].vp_start_17 - pulses[iVec].ev61_start_time); h_vpint_17->Fill(pulses[iVec].vp_int_17);}

                        if (pulses[iVec].vp_start_18 != 0) {h_dt_61_18->Fill(pulses[iVec].wf_time + pulses[iVec].vp_start_18 - pulses[iVec].ev61_start_time); h_vpint_18->Fill(pulses[iVec].vp_int_18);}

                        if (pulses[iVec].vp_start_19 != 0) {h_dt_61_19->Fill(pulses[iVec].wf_time + pulses[iVec].vp_start_19 - pulses[iVec].ev61_start_time); h_vpint_19->Fill(pulses[iVec].vp_int_19);}

                        if (pulses[iVec].vp_start_20 != 0) {h_dt_61_20->Fill(pulses[iVec].wf_time + pulses[iVec].vp_start_20 - pulses[iVec].ev61_start_time); h_vpint_20->Fill(pulses[iVec].vp_int_20);}

                        if (pulses[iVec].vp_start_21 != 0) {h_dt_61_21->Fill(pulses[iVec].wf_time + pulses[iVec].vp_start_21 - pulses[iVec].ev61_start_time); h_vpint_21->Fill(pulses[iVec].vp_int_21);}

                        if (pulses[iVec].vp_start_22 != 0) {h_dt_61_22->Fill(pulses[iVec].wf_time + pulses[iVec].vp_start_22 - pulses[iVec].ev61_start_time); h_vpint_22->Fill(pulses[iVec].vp_int_22);}

                        if (pulses[iVec].vp_start_23 != 0) {h_dt_61_23->Fill(pulses[iVec].wf_time + pulses[iVec].vp_start_23 - pulses[iVec].ev61_start_time); h_vpint_23->Fill(pulses[iVec].vp_int_23);}

                        if (pulses[iVec].vp_start_top != 0) {h_dt_61_top->Fill(pulses[iVec].wf_time + pulses[iVec].vp_start_top - pulses[iVec].ev61_start_time); h_vpint_top->Fill(pulses[iVec].vp_int_top);}
                        */
                        h_dt_v_tmuon->Fill(pulses[iVec].start - pulses[iVec].ev61_start_time, pulses[iVec].start - pulses[iVec].last_muon_time);

                        h_dt_v_tlead->Fill(pulses[iVec].start - pulses[iVec].ev61_start_time, pulses[iVec].start - pulses[iVec].last_lead_time);

                        h_detint->Fill(pulses[iVec].energy);

                        h_dt_v_detint->Fill(pulses[iVec].start - pulses[iVec].ev61_start_time, pulses[iVec].energy);

                    }

                    if (iVec < pulses.size() - 1 && (pulses[iVec + 1].trigger != 1 && pulses[iVec + 1].trigger != 33 && pulses[iVec + 1].trigger != 35)) {

                        if ((pulses[iVec + 1].wf_time - pulses[iVec].ev61_start_time) <= sample_range && (pulses[iVec + 1].wf_time - pulses[iVec].ev61_start_time) > 0) {

                            // Plot positive dt values

                            if ((pulses[iVec + 1].start - pulses[iVec].ev61_start_time) >= -10000 && (pulses[iVec + 1].start - pulses[iVec].ev61_start_time) <= 10000) {num_events_low_dt += 1; h_detint_low_dt->Fill(pulses[iVec + 1].energy); h_avpint_low_dt->Fill(pulses[iVec + 1].all_vp_energy);}
                
                            if ((pulses[iVec + 1].start - pulses[iVec].ev61_start_time) > -500 && (pulses[iVec + 1].start - pulses[iVec].ev61_start_time) < 500) {num_events_low_dt_cent_spike += 1; h_int_cent_spike->Fill(pulses[iVec + 1].energy); h_avpint_cent_spike->Fill(pulses[iVec + 1].all_vp_energy); h_svpint_cent_spike->Fill(pulses[iVec + 1].side_vp_energy); h_tvpint_cent_spike->Fill(pulses[iVec + 1].top_vp_energy);}          // cout << "\n" << "Run, Event, triggerBits = " << iRun << " " << iEnt + 1 << " " << avg_pulse.trigger << endl;}

                            else if ((pulses[iVec + 1].start - pulses[iVec].ev61_start_time) <= -500 && (pulses[iVec + 1].start - pulses[iVec].ev61_start_time) >= -10000) {num_events_low_dt_no_spike_neg += 1; h_int_neg_low_dt->Fill(pulses[iVec + 1].energy); h_avpint_neg_low_dt->Fill(pulses[iVec + 1].all_vp_energy); h_svpint_neg_low_dt->Fill(pulses[iVec + 1].side_vp_energy); h_tvpint_neg_low_dt->Fill(pulses[iVec + 1].top_vp_energy);}

                            else if ((pulses[iVec + 1].start - pulses[iVec].ev61_start_time) >= 500 && (pulses[iVec + 1].start - pulses[iVec].ev61_start_time) <= 10000) {num_events_low_dt_no_spike_pos += 1; h_int_pos_low_dt->Fill(pulses[iVec + 1].energy); h_avpint_pos_low_dt->Fill(pulses[iVec + 1].all_vp_energy); h_svpint_pos_low_dt->Fill(pulses[iVec + 1].side_vp_energy); h_tvpint_pos_low_dt->Fill(pulses[iVec + 1].top_vp_energy);}

                            hist_events += 1;
                            
                            // h_dt_61->Fill(pulses[iVec + 1].start - pulses[iVec].ev61_start_time);

                            if (pulses[iVec + 1].vp_start_16 != 0) {h_dt_61_16->Fill(pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_16 - pulses[iVec].ev61_start_time); h_vpint_16->Fill(pulses[iVec + 1].vp_int_16);} // if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_16 - pulses[iVec].ev61_start_time > 2000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_16 - pulses[iVec].ev61_start_time < 2500) {h_vpint_16_BRNpeak->Fill(pulses[iVec + 1].vp_int_16);} else if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_16 - pulses[iVec].ev61_start_time > 3000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_16 - pulses[iVec].ev61_start_time < 8000) {h_vpint_16_background->Fill(pulses[iVec + 1].vp_int_16);}}

                            if (pulses[iVec + 1].vp_start_17 != 0) {h_dt_61_17->Fill(pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_17 - pulses[iVec].ev61_start_time); h_vpint_17->Fill(pulses[iVec + 1].vp_int_17);} // if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_17 - pulses[iVec].ev61_start_time > 2000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_17 - pulses[iVec].ev61_start_time < 2500) {h_vpint_17_BRNpeak->Fill(pulses[iVec + 1].vp_int_17);} else if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_17 - pulses[iVec].ev61_start_time > 3000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_17 - pulses[iVec].ev61_start_time < 8000) {h_vpint_17_background->Fill(pulses[iVec + 1].vp_int_17);}}            // else if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_17 - pulses[iVec].ev61_start_time > 9000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_17 - pulses[iVec].ev61_start_time < 9200) {cout << "\n" << "SiPM 2 Unexplained dt Peak Event, Run & Event Number = " << pulses[iVec + 1].file_num << " " << pulses[iVec + 1].entry << endl;}}

                            if (pulses[iVec + 1].vp_start_18 != 0) {h_dt_61_18->Fill(pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_18 - pulses[iVec].ev61_start_time); h_vpint_18->Fill(pulses[iVec + 1].vp_int_18);} // if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_18 - pulses[iVec].ev61_start_time > 2000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_18 - pulses[iVec].ev61_start_time < 2500) {h_vpint_18_BRNpeak->Fill(pulses[iVec + 1].vp_int_18);} else if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_18 - pulses[iVec].ev61_start_time > 3000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_18 - pulses[iVec].ev61_start_time < 8000) {h_vpint_18_background->Fill(pulses[iVec + 1].vp_int_18);}}

                            if (pulses[iVec + 1].vp_start_19 != 0) {h_dt_61_19->Fill(pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_19 - pulses[iVec].ev61_start_time); h_vpint_19->Fill(pulses[iVec + 1].vp_int_19);} // if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_19 - pulses[iVec].ev61_start_time > 2000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_19 - pulses[iVec].ev61_start_time < 2500) {h_vpint_19_BRNpeak->Fill(pulses[iVec + 1].vp_int_19);} else if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_19 - pulses[iVec].ev61_start_time > 3000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_19 - pulses[iVec].ev61_start_time < 8000) {h_vpint_19_background->Fill(pulses[iVec + 1].vp_int_19);}}

                            if (pulses[iVec + 1].vp_start_20 != 0) {h_dt_61_20->Fill(pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_20 - pulses[iVec].ev61_start_time); h_vpint_20->Fill(pulses[iVec + 1].vp_int_20);} // if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_20 - pulses[iVec].ev61_start_time > 2000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_20 - pulses[iVec].ev61_start_time < 2500) {h_vpint_20_BRNpeak->Fill(pulses[iVec + 1].vp_int_20);} else if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_20 - pulses[iVec].ev61_start_time > 3000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_20 - pulses[iVec].ev61_start_time < 8000) {h_vpint_20_background->Fill(pulses[iVec + 1].vp_int_20);}}

                            if (pulses[iVec + 1].vp_start_21 != 0) {h_dt_61_21->Fill(pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_21 - pulses[iVec].ev61_start_time); h_vpint_21->Fill(pulses[iVec + 1].vp_int_21);} // if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_21 - pulses[iVec].ev61_start_time > 2000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_21 - pulses[iVec].ev61_start_time < 2500) {h_vpint_21_BRNpeak->Fill(pulses[iVec + 1].vp_int_21); BRN_count += 1; h_BRN_tail_amp->Fill(pulses[iVec + 1].vp_ttp_amp_ratio);} else if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_21 - pulses[iVec].ev61_start_time > 3000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_21 - pulses[iVec].ev61_start_time < 8000) {h_vpint_21_background->Fill(pulses[iVec + 1].vp_int_21); background_count += 1; h_background_tail_amp->Fill(pulses[iVec + 1].vp_ttp_amp_ratio);}} // if (pulses[iVec + 1].vp_int_21 > 500 && pulses[iVec + 1].vp_int_21 < 550) {cout << "\n" << "SiPM 6 Integral Peak Event, Run & Event Number = " << pulses[iVec + 1].file_num << " " << pulses[iVec + 1].entry << endl;}}

                            if (pulses[iVec + 1].vp_start_22 != 0) {h_dt_61_22->Fill(pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_22 - pulses[iVec].ev61_start_time); h_vpint_22->Fill(pulses[iVec + 1].vp_int_22);} // if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_22 - pulses[iVec].ev61_start_time > 2000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_22 - pulses[iVec].ev61_start_time < 2500) {h_vpint_22_BRNpeak->Fill(pulses[iVec + 1].vp_int_22);} else if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_22 - pulses[iVec].ev61_start_time > 3000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_22 - pulses[iVec].ev61_start_time < 8000) {h_vpint_22_background->Fill(pulses[iVec + 1].vp_int_22);}}

                            if (pulses[iVec + 1].vp_start_23 != 0) {h_dt_61_23->Fill(pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_23 - pulses[iVec].ev61_start_time); h_vpint_23->Fill(pulses[iVec + 1].vp_int_23);} // if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_23 - pulses[iVec].ev61_start_time > 2000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_23 - pulses[iVec].ev61_start_time < 2500) {h_vpint_23_BRNpeak->Fill(pulses[iVec + 1].vp_int_23);} else if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_23 - pulses[iVec].ev61_start_time > 3000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_23 - pulses[iVec].ev61_start_time < 8000) {h_vpint_23_background->Fill(pulses[iVec + 1].vp_int_23);}}

                            if (pulses[iVec + 1].vp_start_24 != 0) {h_dt_61_24->Fill(pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_24 - pulses[iVec].ev61_start_time); h_vpint_24->Fill(pulses[iVec + 1].vp_int_24);}

                            if (pulses[iVec + 1].vp_start_25 != 0) {h_dt_61_25->Fill(pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_25 - pulses[iVec].ev61_start_time); h_vpint_25->Fill(pulses[iVec + 1].vp_int_25);}
                            
                            if (pulses[iVec + 1].vp_start_top != 0) {h_dt_61_top->Fill(pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_top - pulses[iVec].ev61_start_time); h_vpint_top->Fill(pulses[iVec + 1].vp_int_top);} // if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_top - pulses[iVec].ev61_start_time > 2000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_top - pulses[iVec].ev61_start_time < 2500) {h_vpint_top_BRNpeak->Fill(pulses[iVec + 1].vp_int_top);} else if (pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_top - pulses[iVec].ev61_start_time > 3000 && pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_top - pulses[iVec].ev61_start_time < 8000) {h_vpint_top_background->Fill(pulses[iVec + 1].vp_int_top);}} // if (pulses[iVec + 1].vp_int_top > 300 && pulses[iVec + 1].vp_int_top < 350) {cout << "\n" << "SiPM 9 & 10 Integral Peak Event, Run & Event Number = " << pulses[iVec + 1].file_num << " " << pulses[iVec + 1].entry << endl;}}

                            if ((pulses[iVec + 1].vp_start_17 != 0 || pulses[iVec + 1].vp_start_18 != 0 || pulses[iVec + 1].vp_start_19 != 0 || pulses[iVec + 1].vp_start_20 != 0 || pulses[iVec + 1].vp_start_21 != 0 || pulses[iVec + 1].vp_start_22 != 0 || pulses[iVec + 1].vp_start_23 != 0 || pulses[iVec + 1].vp_start_top != 0) && pulses[iVec + 1].vp_start_16 == 0) {h_svpint_not_16->Fill(pulses[iVec + 1].vp_int_16);}         // .side_vp_energy);}

                            if ((pulses[iVec + 1].vp_start_16 != 0 || pulses[iVec + 1].vp_start_18 != 0 || pulses[iVec + 1].vp_start_19 != 0 || pulses[iVec + 1].vp_start_20 != 0 || pulses[iVec + 1].vp_start_21 != 0 || pulses[iVec + 1].vp_start_22 != 0 || pulses[iVec + 1].vp_start_23 != 0 || pulses[iVec + 1].vp_start_top != 0) && pulses[iVec + 1].vp_start_17 == 0) {h_svpint_not_17->Fill(pulses[iVec + 1].vp_int_17);}

                            if ((pulses[iVec + 1].vp_start_16 != 0 || pulses[iVec + 1].vp_start_17 != 0 || pulses[iVec + 1].vp_start_19 != 0 || pulses[iVec + 1].vp_start_20 != 0 || pulses[iVec + 1].vp_start_21 != 0 || pulses[iVec + 1].vp_start_22 != 0 || pulses[iVec + 1].vp_start_23 != 0 || pulses[iVec + 1].vp_start_top != 0) && pulses[iVec + 1].vp_start_18 == 0) {h_svpint_not_18->Fill(pulses[iVec + 1].vp_int_18);}

                            if ((pulses[iVec + 1].vp_start_16 != 0 || pulses[iVec + 1].vp_start_17 != 0 || pulses[iVec + 1].vp_start_18 != 0 || pulses[iVec + 1].vp_start_20 != 0 || pulses[iVec + 1].vp_start_21 != 0 || pulses[iVec + 1].vp_start_22 != 0 || pulses[iVec + 1].vp_start_23 != 0 || pulses[iVec + 1].vp_start_top != 0) && pulses[iVec + 1].vp_start_19 == 0) {h_svpint_not_19->Fill(pulses[iVec + 1].vp_int_19);}

                            if ((pulses[iVec + 1].vp_start_16 != 0 || pulses[iVec + 1].vp_start_17 != 0 || pulses[iVec + 1].vp_start_18 != 0 || pulses[iVec + 1].vp_start_19 != 0 || pulses[iVec + 1].vp_start_21 != 0 || pulses[iVec + 1].vp_start_22 != 0 || pulses[iVec + 1].vp_start_23 != 0 || pulses[iVec + 1].vp_start_top != 0) && pulses[iVec + 1].vp_start_20 == 0) {h_svpint_not_20->Fill(pulses[iVec + 1].vp_int_20);}
                            
                            if ((pulses[iVec + 1].vp_start_16 != 0 || pulses[iVec + 1].vp_start_17 != 0 || pulses[iVec + 1].vp_start_18 != 0 || pulses[iVec + 1].vp_start_19 != 0 || pulses[iVec + 1].vp_start_20 != 0 || pulses[iVec + 1].vp_start_22 != 0 || pulses[iVec + 1].vp_start_23 != 0 || pulses[iVec + 1].vp_start_top != 0) && pulses[iVec + 1].vp_start_21 == 0) {h_svpint_not_21->Fill(pulses[iVec + 1].vp_int_21);}

                            if ((pulses[iVec + 1].vp_start_16 != 0 || pulses[iVec + 1].vp_start_17 != 0 || pulses[iVec + 1].vp_start_18 != 0 || pulses[iVec + 1].vp_start_19 != 0 || pulses[iVec + 1].vp_start_20 != 0 || pulses[iVec + 1].vp_start_21 != 0 || pulses[iVec + 1].vp_start_23 != 0 || pulses[iVec + 1].vp_start_top != 0) && pulses[iVec + 1].vp_start_22 == 0) {h_svpint_not_22->Fill(pulses[iVec + 1].vp_int_22);}

                            if ((pulses[iVec + 1].vp_start_16 != 0 || pulses[iVec + 1].vp_start_17 != 0 || pulses[iVec + 1].vp_start_18 != 0 || pulses[iVec + 1].vp_start_19 != 0 || pulses[iVec + 1].vp_start_20 != 0 || pulses[iVec + 1].vp_start_21 != 0 || pulses[iVec + 1].vp_start_22 != 0 || pulses[iVec + 1].vp_start_top != 0) && pulses[iVec + 1].vp_start_23 == 0) {h_svpint_not_23->Fill(pulses[iVec + 1].vp_int_23);}

                            if ((pulses[iVec + 1].vp_start_16 != 0 || pulses[iVec + 1].vp_start_17 != 0 || pulses[iVec + 1].vp_start_18 != 0 || pulses[iVec + 1].vp_start_19 != 0 || pulses[iVec + 1].vp_start_20 != 0 || pulses[iVec + 1].vp_start_21 != 0 || pulses[iVec + 1].vp_start_22 != 0 || pulses[iVec + 1].vp_start_23 != 0) && pulses[iVec + 1].vp_start_top == 0) {h_tvpint_not_top->Fill(pulses[iVec + 1].vp_int_top);}            // if (pulses[iVec + 1].vp_int_top > 800) {cout << "\n" << "SiPM9&10 Untagged Muon Event, Run & Event Number = " << pulses[iVec + 1].file_num << " " << pulses[iVec + 1].entry << endl;}}
                            
                            h_vpint_16_all->Fill(pulses[iVec + 1].vp_int_16);

                            h_vpint_17_all->Fill(pulses[iVec + 1].vp_int_17);

                            h_vpint_18_all->Fill(pulses[iVec + 1].vp_int_18); 

                            h_vpint_19_all->Fill(pulses[iVec + 1].vp_int_19); 

                            h_vpint_20_all->Fill(pulses[iVec + 1].vp_int_20);

                            h_vpint_21_all->Fill(pulses[iVec + 1].vp_int_21);

                            h_vpint_22_all->Fill(pulses[iVec + 1].vp_int_22);

                            h_vpint_23_all->Fill(pulses[iVec + 1].vp_int_23);

                            h_vpint_top_all->Fill(pulses[iVec + 1].vp_int_top);
                            
                            h_dt_v_tmuon->Fill(pulses[iVec + 1].start - pulses[iVec].ev61_start_time, pulses[iVec + 1].start - pulses[iVec + 1].last_muon_time);

                            h_dt_v_tlead->Fill(pulses[iVec + 1].start - pulses[iVec].ev61_start_time, pulses[iVec + 1].start - pulses[iVec + 1].last_lead_time);

                            h_detint->Fill(pulses[iVec + 1].energy);

                            h_dt_v_detint->Fill(pulses[iVec + 1].start - pulses[iVec].ev61_start_time, pulses[iVec + 1].energy);

                            if (pulses[iVec - 1].vp_int_21 > 1200 && counter < 1000) {
                                
                                counter += 1;

                                cout << "\n" << "POSITIVE dt Value Event:" << endl;

                                cout << "\n" << "Ev61 Event: Run, Event Number, triggerBits, SiPM6 Integral Value = " << pulses[iVec].file_num << " " << pulses[iVec].entry << " " << pulses[iVec].trigger << " " << pulses[iVec].vp_int_21 << endl;

                                cout << "\n" << "VP Event: Run, Event Number, triggerBits, SiPM6 Integral Value, dt Value = " << pulses[iVec + 1].file_num << " " << pulses[iVec + 1].entry << " " << pulses[iVec + 1].trigger << " " << pulses[iVec + 1].vp_int_21 << " " << pulses[iVec + 1].wf_time + pulses[iVec + 1].vp_start_21 - pulses[iVec].ev61_start_time << endl;
                            
                            }

                        }

                    }

                }

            }

            // Reset search parameter after finding low light LED event
            
            bt_HL_n_LL = false; h_LLint->Fill(avg_pulse.energy); h_vec_size->Fill(pulses.size()); pulses.clear();
            
        }
        
        if (avg_pulse.beam) {h_ev61int->Fill(avg_pulse.energy); h_ev61vpint->Fill(avg_pulse.all_vp_energy);}

        if (avg_pulse.trigger == 34) {h_34detint->Fill(avg_pulse.energy); h_34vpint->Fill(avg_pulse.all_vp_energy);}

        h_int->Fill(avg_pulse.energy);

    } // Event loop

    h_vpint_16_BRNpeak->Sumw2(); h_vpint_16_background->Sumw2(); h_vpint_16_background->Scale(0.1);

    h_vpint_17_BRNpeak->Sumw2(); h_vpint_17_background->Sumw2(); h_vpint_17_background->Scale(0.1);

    h_vpint_18_BRNpeak->Sumw2(); h_vpint_18_background->Sumw2(); h_vpint_18_background->Scale(0.1);

    h_vpint_19_BRNpeak->Sumw2(); h_vpint_19_background->Sumw2(); h_vpint_19_background->Scale(0.1);

    h_vpint_20_BRNpeak->Sumw2(); h_vpint_20_background->Sumw2(); h_vpint_20_background->Scale(0.1);

    h_vpint_21_BRNpeak->Sumw2(); h_vpint_21_background->Sumw2(); h_vpint_21_background->Scale(0.1);

    h_vpint_22_BRNpeak->Sumw2(); h_vpint_22_background->Sumw2(); h_vpint_22_background->Scale(0.1);

    h_vpint_23_BRNpeak->Sumw2(); h_vpint_23_background->Sumw2(); h_vpint_23_background->Scale(0.1);

    h_vpint_top_BRNpeak->Sumw2(); h_vpint_top_background->Sumw2(); h_vpint_top_background->Scale(0.1);

    h_BRNint_16->Add(h_vpint_16_BRNpeak, h_vpint_16_background, 1, -1);

    h_BRNint_17->Add(h_vpint_17_BRNpeak, h_vpint_17_background, 1, -1);

    h_BRNint_18->Add(h_vpint_18_BRNpeak, h_vpint_18_background, 1, -1);

    h_BRNint_19->Add(h_vpint_19_BRNpeak, h_vpint_19_background, 1, -1);

    h_BRNint_20->Add(h_vpint_20_BRNpeak, h_vpint_20_background, 1, -1);

    h_BRNint_21->Add(h_vpint_21_BRNpeak, h_vpint_21_background, 1, -1);         // h_BRNint_21->Sumw2();

    h_BRNint_22->Add(h_vpint_22_BRNpeak, h_vpint_22_background, 1, -1);

    h_BRNint_23->Add(h_vpint_23_BRNpeak, h_vpint_23_background, 1, -1);

    h_BRNint_top->Add(h_vpint_top_BRNpeak, h_vpint_top_background, 1, -1);

    Double_t A_error = h_vpint_21_BRNpeak->GetBinError(62);
    Double_t B_error = h_vpint_21_background->GetBinError(62);
    Double_t C_error = h_BRNint_21->GetBinError(62);
    cout << "\n" << "Error of bin 62 for plots A, B and C = " << A_error << ", " << B_error << ", " << C_error << endl;
    
    if (run != 11241 && run != 12636) {

        cout << "\n" << "Number of low dt events, Runs 7894-7918 = " << num_events_low_dt << endl;

        cout << "\n" << "Number of central spike events, Runs 7894-7918 = " << num_events_low_dt_cent_spike << endl;

        cout << "\n" << "Number of negative low dt events, excluding central spike, Runs 7894-7918 = " << num_events_low_dt_no_spike_neg << endl;

        cout << "\n" << "Number of positive low dt events, excluding central spike, Runs 7894-7918 = " << num_events_low_dt_no_spike_pos << endl;

        cout << "\n" << "Number of Event 61 events found, Runs 7894-7918 = " << take_sample_count << endl;

        cout << "\n" << "Number of histogram events within 10 us of central sample, Runs 7894-7918 = " << hist_events << endl;

        cout << "\n" << "Number of SiPM6 BRN events = " << BRN_count - background_count << "\n" << endl;

    }

    if (run == 11608) {

        cout << "\n" << "Number of low dt events, Runs 11608-11632 = " << num_events_low_dt << endl;

        cout << "\n" << "Number of central spike events, Runs 11608-11632 = " << num_events_low_dt_cent_spike << endl;

        cout << "\n" << "Number of negative low dt events, excluding central spike, Runs 11608-11632 = " << num_events_low_dt_no_spike_neg << endl;

        cout << "\n" << "Number of positive low dt events, excluding central spike, Runs 11608-11632 = " << num_events_low_dt_no_spike_pos << endl;

        cout << "\n" << "Number of Event 61 events found, Runs 11608-11632 = " << take_sample_count << endl;

        cout << "\n" << "Number of histogram events within 10 us of central sample, Runs 11608-11632 = " << hist_events << "\n" << endl;

    }

    if (run == 12636) {

        cout << "\n" << "Number of low dt events, Runs 12636-12660 = " << num_events_low_dt << endl;

        cout << "\n" << "Number of central spike events, Runs 12636-12660 = " << num_events_low_dt_cent_spike << endl;

        cout << "\n" << "Number of negative low dt events, excluding central spike, Runs 12636-12660 = " << num_events_low_dt_no_spike_neg << endl;

        cout << "\n" << "Number of positive low dt events, excluding central spike, Runs 12636-12660 = " << num_events_low_dt_no_spike_pos << endl;

        cout << "\n" << "Number of Event 61 events found, Runs 12636-12660 = " << take_sample_count << endl;

        cout << "\n" << "Number of histogram events within 10 us of central sample, Runs 12636-12660 = " << hist_events << "\n" << endl;

    }

    /* Fill & plot histograms */

    if (run != 11241 && run != 12636) {

        TCanvas* c_dt_61 = new TCanvas("c_dt_61", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61->SetLogy();
        c_dt_61->cd();
        h_dt_61->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61->GetYaxis()->SetTitle("Counts");
        h_dt_61->Draw();
        c_dt_61->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61.png");
        // c_dt_61->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61.png");
        c_dt_61->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61.png"));
        h_dt_61->Reset();

        TCanvas* c_dt_61_16 = new TCanvas("c_dt_61_16", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_16->SetLogy();
        c_dt_61_16->cd();
        TF1* func16 = new TF1("func16", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -22000, -16000);           // -22000, -16000); // 1000, 10000);
        func16->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func16->SetParameters(5, 20, -19000, 100);          // 5, 20, -19000, 100);          // 5, 10, 2200, 200);
        h_dt_61_16->Fit("func16", "R");
        h_dt_61_16->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_16->GetYaxis()->SetTitle("Counts");
        h_dt_61_16->Draw("E1");
        func16->Draw("same");
        func16->SetLineColor(2);

        Double_t param16_1 = func16->GetParameter(0);
        Double_t param16_2 = func16->GetParameter(1);
        Double_t param16_3 = func16->GetParameter(2);
        Double_t param16_4 = func16->GetParameter(3);
        Double_t error16_1 = func16->GetParError(0);
        Double_t error16_2 = func16->GetParError(1);
        Double_t error16_3 = func16->GetParError(2);
        Double_t error16_4 = func16->GetParError(3);
        Double_t chi2_16 = func16->GetChisquare();

        TLatex *lat_dt_16_1 = new TLatex(0.3, 0.88, Form("Baseline = %.2f +/- %.2f", param16_1, error16_1));
        lat_dt_16_1->SetNDC();
        lat_dt_16_1->SetTextColor(1);
        lat_dt_16_1->SetTextSize(0.03);
        lat_dt_16_1->Draw();

        TLatex *lat_dt_16_2 = new TLatex(0.3, 0.84, Form("Amplitude = %.2f +/- %.2f", param16_2, error16_2));
        lat_dt_16_2->SetNDC();
        lat_dt_16_2->SetTextColor(1);
        lat_dt_16_2->SetTextSize(0.03);
        lat_dt_16_2->Draw();

        TLatex *lat_dt_16_3 = new TLatex(0.3, 0.80, Form("Peak Center = %.2f +/- %.2f", param16_3, error16_3));
        lat_dt_16_3->SetNDC();
        lat_dt_16_3->SetTextColor(1);
        lat_dt_16_3->SetTextSize(0.03);
        lat_dt_16_3->Draw();

        TLatex *lat_dt_16_4 = new TLatex(0.3, 0.76, Form("Standard Deviation = %.2f +/- %.2f", param16_4, error16_4));
        lat_dt_16_4->SetNDC();
        lat_dt_16_4->SetTextColor(1);
        lat_dt_16_4->SetTextSize(0.03);
        lat_dt_16_4->Draw();

        TLatex *lat_dt_16_5 = new TLatex(0.3, 0.72, Form("Chi-Squared = %.2f", chi2_16));
        lat_dt_16_5->SetNDC();
        lat_dt_16_5->SetTextColor(1);
        lat_dt_16_5->SetTextSize(0.03);
        lat_dt_16_5->Draw();

        c_dt_61_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_16.png");
        // c_dt_61_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_16.png");
        c_dt_61_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_16.png"));
        h_dt_61_16->Reset();

        TCanvas* c_dt_61_17 = new TCanvas("c_dt_61_17", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_17->SetLogy();
        c_dt_61_17->cd();
        TF1* func17 = new TF1("func17", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -22000, -16000);
        func17->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func17->SetParameters(15, 40, -19000, 100);          // 15, 40, -19000, 100);          // 15, 25, 2200, 200);
        h_dt_61_17->Fit("func17", "R");
        h_dt_61_17->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_17->GetYaxis()->SetTitle("Counts");
        h_dt_61_17->Draw("E1");
        func17->Draw("same");
        func17->SetLineColor(2);

        Double_t param17_1 = func17->GetParameter(0);
        Double_t param17_2 = func17->GetParameter(1);
        Double_t param17_3 = func17->GetParameter(2);
        Double_t param17_4 = func17->GetParameter(3);
        Double_t error17_1 = func17->GetParError(0);
        Double_t error17_2 = func17->GetParError(1);
        Double_t error17_3 = func17->GetParError(2);
        Double_t error17_4 = func17->GetParError(3);
        Double_t chi2_17 = func17->GetChisquare();

        TLatex *lat_dt_17_1 = new TLatex(0.3, 0.88, Form("Baseline = %.2f +/- %.2f", param17_1, error17_1));
        lat_dt_17_1->SetNDC();
        lat_dt_17_1->SetTextColor(1);
        lat_dt_17_1->SetTextSize(0.03);
        lat_dt_17_1->Draw();

        TLatex *lat_dt_17_2 = new TLatex(0.3, 0.84, Form("Amplitude = %.2f +/- %.2f", param17_2, error17_2));
        lat_dt_17_2->SetNDC();
        lat_dt_17_2->SetTextColor(1);
        lat_dt_17_2->SetTextSize(0.03);
        lat_dt_17_2->Draw();

        TLatex *lat_dt_17_3 = new TLatex(0.3, 0.80, Form("Peak Center = %.2f +/- %.2f", param17_3, error17_3));
        lat_dt_17_3->SetNDC();
        lat_dt_17_3->SetTextColor(1);
        lat_dt_17_3->SetTextSize(0.03);
        lat_dt_17_3->Draw();

        TLatex *lat_dt_17_4 = new TLatex(0.3, 0.76, Form("Standard Deviation = %.2f +/- %.2f", param17_4, error17_4));
        lat_dt_17_4->SetNDC();
        lat_dt_17_4->SetTextColor(1);
        lat_dt_17_4->SetTextSize(0.03);
        lat_dt_17_4->Draw();

        TLatex *lat_dt_17_5 = new TLatex(0.3, 0.72, Form("Chi-Squared = %.2f", chi2_17));
        lat_dt_17_5->SetNDC();
        lat_dt_17_5->SetTextColor(1);
        lat_dt_17_5->SetTextSize(0.03);
        lat_dt_17_5->Draw();

        c_dt_61_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_17.png");
        // c_dt_61_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_17.png");
        c_dt_61_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_17.png"));
        h_dt_61_17->Reset();

        TCanvas* c_dt_61_18 = new TCanvas("c_dt_61_18", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_18->SetLogy();
        c_dt_61_18->cd();
        TF1* func18 = new TF1("func18", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -22000, -16000);
        func18->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func18->SetParameters(50, 5, -19000, 100);          // 50, 5, -19000, 100);          // 35, 5, 2200, 200);
        h_dt_61_18->Fit("func18", "R");
        h_dt_61_18->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_18->GetYaxis()->SetTitle("Counts");
        // h_dt_61_18->SetMaximum(100);
        h_dt_61_18->Draw("E1");
        func18->Draw("same");
        func18->SetLineColor(2);
        
        Double_t param18_1 = func18->GetParameter(0);
        Double_t param18_2 = func18->GetParameter(1);
        Double_t param18_3 = func18->GetParameter(2);
        Double_t param18_4 = func18->GetParameter(3);
        Double_t error18_1 = func18->GetParError(0);
        Double_t error18_2 = func18->GetParError(1);
        Double_t error18_3 = func18->GetParError(2);
        Double_t error18_4 = func18->GetParError(3);
        Double_t chi2_18 = func18->GetChisquare();

        TLatex *lat_dt_18_1 = new TLatex(0.3, 0.88, Form("Baseline = %.2f +/- %.2f", param18_1, error18_1));
        lat_dt_18_1->SetNDC();
        lat_dt_18_1->SetTextColor(1);
        lat_dt_18_1->SetTextSize(0.03);
        lat_dt_18_1->Draw();

        TLatex *lat_dt_18_2 = new TLatex(0.3, 0.84, Form("Amplitude = %.2f +/- %.2f", param18_2, error18_2));
        lat_dt_18_2->SetNDC();
        lat_dt_18_2->SetTextColor(1);
        lat_dt_18_2->SetTextSize(0.03);
        lat_dt_18_2->Draw();

        TLatex *lat_dt_18_3 = new TLatex(0.3, 0.80, Form("Peak Center = %.2f +/- %.2f", param18_3, error18_3));
        lat_dt_18_3->SetNDC();
        lat_dt_18_3->SetTextColor(1);
        lat_dt_18_3->SetTextSize(0.03);
        lat_dt_18_3->Draw();

        TLatex *lat_dt_18_4 = new TLatex(0.3, 0.76, Form("Standard Deviation = %.2f +/- %.2f", param18_4, error18_4));
        lat_dt_18_4->SetNDC();
        lat_dt_18_4->SetTextColor(1);
        lat_dt_18_4->SetTextSize(0.03);
        lat_dt_18_4->Draw();

        TLatex *lat_dt_18_5 = new TLatex(0.3, 0.72, Form("Chi-Squared = %.2f", chi2_18));
        lat_dt_18_5->SetNDC();
        lat_dt_18_5->SetTextColor(1);
        lat_dt_18_5->SetTextSize(0.03);
        lat_dt_18_5->Draw();

        c_dt_61_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_18.png");
        // c_dt_61_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_18.png");
        c_dt_61_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_18.png"));
        h_dt_61_18->Reset();

        TCanvas* c_dt_61_19 = new TCanvas("c_dt_61_19", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_19->SetLogy();
        c_dt_61_19->cd();
        TF1* func19 = new TF1("func19", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -22000, -16000);
        func19->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func19->SetParameters(30, 5, -19000, 100);          // 30, 5, -19000, 100);          // 25, 5, 2200, 200);
        h_dt_61_19->Fit("func19", "R");
        h_dt_61_19->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_19->GetYaxis()->SetTitle("Counts");
        // h_dt_61_19->SetMaximum(85);
        h_dt_61_19->Draw("E1");
        func19->Draw("same");
        func19->SetLineColor(2);
        
        Double_t param19_1 = func19->GetParameter(0);
        Double_t param19_2 = func19->GetParameter(1);
        Double_t param19_3 = func19->GetParameter(2);
        Double_t param19_4 = func19->GetParameter(3);
        Double_t error19_1 = func19->GetParError(0);
        Double_t error19_2 = func19->GetParError(1);
        Double_t error19_3 = func19->GetParError(2);
        Double_t error19_4 = func19->GetParError(3);
        Double_t chi2_19 = func19->GetChisquare();

        TLatex *lat_dt_19_1 = new TLatex(0.3, 0.88, Form("Baseline = %.2f +/- %.2f", param19_1, error19_1));
        lat_dt_19_1->SetNDC();
        lat_dt_19_1->SetTextColor(1);
        lat_dt_19_1->SetTextSize(0.03);
        lat_dt_19_1->Draw();

        TLatex *lat_dt_19_2 = new TLatex(0.3, 0.84, Form("Amplitude = %.2f +/- %.2f", param19_2, error19_2));
        lat_dt_19_2->SetNDC();
        lat_dt_19_2->SetTextColor(1);
        lat_dt_19_2->SetTextSize(0.03);
        lat_dt_19_2->Draw();

        TLatex *lat_dt_19_3 = new TLatex(0.3, 0.80, Form("Peak Center = %.2f +/- %.2f", param19_3, error19_3));
        lat_dt_19_3->SetNDC();
        lat_dt_19_3->SetTextColor(1);
        lat_dt_19_3->SetTextSize(0.03);
        lat_dt_19_3->Draw();

        TLatex *lat_dt_19_4 = new TLatex(0.3, 0.76, Form("Standard Deviation = %.2f +/- %.2f", param19_4, error19_4));
        lat_dt_19_4->SetNDC();
        lat_dt_19_4->SetTextColor(1);
        lat_dt_19_4->SetTextSize(0.03);
        lat_dt_19_4->Draw();

        TLatex *lat_dt_19_5 = new TLatex(0.3, 0.72, Form("Chi-Squared = %.2f", chi2_19));
        lat_dt_19_5->SetNDC();
        lat_dt_19_5->SetTextColor(1);
        lat_dt_19_5->SetTextSize(0.03);
        lat_dt_19_5->Draw();

        c_dt_61_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_19.png");
        // c_dt_61_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_19.png");
        c_dt_61_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_19.png"));
        h_dt_61_19->Reset();

        TCanvas* c_dt_61_20 = new TCanvas("c_dt_61_20", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_20->SetLogy();
        c_dt_61_20->cd();
        TF1* func20 = new TF1("func20", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -22000, -16000);
        func20->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func20->SetParameters(10, 15, -19000, 100);          // 10, 15, -19000, 100);          // 5, 10, 2200, 200);
        h_dt_61_20->Fit("func20", "R");
        h_dt_61_20->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_20->GetYaxis()->SetTitle("Counts");
        h_dt_61_20->Draw("E1");
        func20->Draw("same");
        func20->SetLineColor(2);
        
        Double_t param20_1 = func20->GetParameter(0);
        Double_t param20_2 = func20->GetParameter(1);
        Double_t param20_3 = func20->GetParameter(2);
        Double_t param20_4 = func20->GetParameter(3);
        Double_t error20_1 = func20->GetParError(0);
        Double_t error20_2 = func20->GetParError(1);
        Double_t error20_3 = func20->GetParError(2);
        Double_t error20_4 = func20->GetParError(3);
        Double_t chi2_20 = func20->GetChisquare();

        TLatex *lat_dt_20_1 = new TLatex(0.3, 0.88, Form("Baseline = %.2f +/- %.2f", param20_1, error20_1));
        lat_dt_20_1->SetNDC();
        lat_dt_20_1->SetTextColor(1);
        lat_dt_20_1->SetTextSize(0.03);
        lat_dt_20_1->Draw();

        TLatex *lat_dt_20_2 = new TLatex(0.3, 0.84, Form("Amplitude = %.2f +/- %.2f", param20_2, error20_2));
        lat_dt_20_2->SetNDC();
        lat_dt_20_2->SetTextColor(1);
        lat_dt_20_2->SetTextSize(0.03);
        lat_dt_20_2->Draw();

        TLatex *lat_dt_20_3 = new TLatex(0.3, 0.80, Form("Peak Center = %.2f +/- %.2f", param20_3, error20_3));
        lat_dt_20_3->SetNDC();
        lat_dt_20_3->SetTextColor(1);
        lat_dt_20_3->SetTextSize(0.03);
        lat_dt_20_3->Draw();

        TLatex *lat_dt_20_4 = new TLatex(0.3, 0.76, Form("Standard Deviation = %.2f +/- %.2f", param20_4, error20_4));
        lat_dt_20_4->SetNDC();
        lat_dt_20_4->SetTextColor(1);
        lat_dt_20_4->SetTextSize(0.03);
        lat_dt_20_4->Draw();

        TLatex *lat_dt_20_5 = new TLatex(0.3, 0.72, Form("Chi-Squared = %.2f", chi2_20));
        lat_dt_20_5->SetNDC();
        lat_dt_20_5->SetTextColor(1);
        lat_dt_20_5->SetTextSize(0.03);
        lat_dt_20_5->Draw();

        c_dt_61_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_20.png");
        // c_dt_61_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_20.png");
        c_dt_61_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_20.png"));
        h_dt_61_20->Reset();

        TCanvas* c_dt_61_21 = new TCanvas("c_dt_61_21", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_21->SetLogy();
        c_dt_61_21->cd();
        TF1* func21 = new TF1("func21", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -22000, -16000);
        func21->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func21->SetParameters(20, 55, -19000, 100);          // 20, 55, -19000, 100);          // 15, 30, 2200, 200);
        h_dt_61_21->Fit("func21", "R");
        h_dt_61_21->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_21->GetYaxis()->SetTitle("Counts");
        h_dt_61_21->Draw("E1");
        func21->Draw("same");
        func21->SetLineColor(2);
                
        Double_t param21_1 = func21->GetParameter(0);
        Double_t param21_2 = func21->GetParameter(1);
        Double_t param21_3 = func21->GetParameter(2);
        Double_t param21_4 = func21->GetParameter(3);
        Double_t error21_1 = func21->GetParError(0);
        Double_t error21_2 = func21->GetParError(1);
        Double_t error21_3 = func21->GetParError(2);
        Double_t error21_4 = func21->GetParError(3);
        Double_t chi2_21 = func21->GetChisquare();

        TLatex *lat_dt_21_1 = new TLatex(0.3, 0.88, Form("Baseline = %.2f +/- %.2f", param21_1, error21_1));
        lat_dt_21_1->SetNDC();
        lat_dt_21_1->SetTextColor(1);
        lat_dt_21_1->SetTextSize(0.03);
        lat_dt_21_1->Draw();

        TLatex *lat_dt_21_2 = new TLatex(0.3, 0.84, Form("Amplitude = %.2f +/- %.2f", param21_2, error21_2));
        lat_dt_21_2->SetNDC();
        lat_dt_21_2->SetTextColor(1);
        lat_dt_21_2->SetTextSize(0.03);
        lat_dt_21_2->Draw();

        TLatex *lat_dt_21_3 = new TLatex(0.3, 0.80, Form("Peak Center = %.2f +/- %.2f", param21_3, error21_3));
        lat_dt_21_3->SetNDC();
        lat_dt_21_3->SetTextColor(1);
        lat_dt_21_3->SetTextSize(0.03);
        lat_dt_21_3->Draw();

        TLatex *lat_dt_21_4 = new TLatex(0.3, 0.76, Form("Standard Deviation = %.2f +/- %.2f", param21_4, error21_4));
        lat_dt_21_4->SetNDC();
        lat_dt_21_4->SetTextColor(1);
        lat_dt_21_4->SetTextSize(0.03);
        lat_dt_21_4->Draw();

        TLatex *lat_dt_21_5 = new TLatex(0.3, 0.72, Form("Chi-Squared = %.2f", chi2_21));
        lat_dt_21_5->SetNDC();
        lat_dt_21_5->SetTextColor(1);
        lat_dt_21_5->SetTextSize(0.03);
        lat_dt_21_5->Draw();

        c_dt_61_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_21.png");
        // c_dt_61_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_21.png");
        c_dt_61_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_21.png"));
        h_dt_61_21->Reset();

        TCanvas* c_dt_61_22 = new TCanvas("c_dt_61_22", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_22->SetLogy();
        c_dt_61_22->cd();
        TF1* func22 = new TF1("func22", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -22000, -16000);
        func22->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func22->SetParameters(20, 40, -19000, 100);          // 20, 40, -19000, 100);          // 15, 25, 2200, 200);
        h_dt_61_22->Fit("func22", "R");
        h_dt_61_22->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_22->GetYaxis()->SetTitle("Counts");
        h_dt_61_22->Draw("E1");
        func22->Draw("same");
        func22->SetLineColor(2);
                
        Double_t param22_1 = func22->GetParameter(0);
        Double_t param22_2 = func22->GetParameter(1);
        Double_t param22_3 = func22->GetParameter(2);
        Double_t param22_4 = func22->GetParameter(3);
        Double_t error22_1 = func22->GetParError(0);
        Double_t error22_2 = func22->GetParError(1);
        Double_t error22_3 = func22->GetParError(2);
        Double_t error22_4 = func22->GetParError(3);
        Double_t chi2_22 = func22->GetChisquare();

        TLatex *lat_dt_22_1 = new TLatex(0.3, 0.88, Form("Baseline = %.2f +/- %.2f", param22_1, error22_1));
        lat_dt_22_1->SetNDC();
        lat_dt_22_1->SetTextColor(1);
        lat_dt_22_1->SetTextSize(0.03);
        lat_dt_22_1->Draw();

        TLatex *lat_dt_22_2 = new TLatex(0.3, 0.84, Form("Amplitude = %.2f +/- %.2f", param22_2, error22_2));
        lat_dt_22_2->SetNDC();
        lat_dt_22_2->SetTextColor(1);
        lat_dt_22_2->SetTextSize(0.03);
        lat_dt_22_2->Draw();

        TLatex *lat_dt_22_3 = new TLatex(0.3, 0.80, Form("Peak Center = %.2f +/- %.2f", param22_3, error22_3));
        lat_dt_22_3->SetNDC();
        lat_dt_22_3->SetTextColor(1);
        lat_dt_22_3->SetTextSize(0.03);
        lat_dt_22_3->Draw();

        TLatex *lat_dt_22_4 = new TLatex(0.3, 0.76, Form("Standard Deviation = %.2f +/- %.2f", param22_4, error22_4));
        lat_dt_22_4->SetNDC();
        lat_dt_22_4->SetTextColor(1);
        lat_dt_22_4->SetTextSize(0.03);
        lat_dt_22_4->Draw();

        TLatex *lat_dt_22_5 = new TLatex(0.3, 0.72, Form("Chi-Squared = %.2f", chi2_22));
        lat_dt_22_5->SetNDC();
        lat_dt_22_5->SetTextColor(1);
        lat_dt_22_5->SetTextSize(0.03);
        lat_dt_22_5->Draw();

        c_dt_61_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_22.png");
        // c_dt_61_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_22.png");
        c_dt_61_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_22.png"));
        h_dt_61_22->Reset();

        TCanvas* c_dt_61_23 = new TCanvas("c_dt_61_23", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_23->SetLogy();
        c_dt_61_23->cd();
        TF1* func23 = new TF1("func23", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -22000, -16000);
        func23->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func23->SetParameters(15, 15, -19000, 100);          // 15, 15, -19000, 100);          // 10, 5, 2200, 200);
        h_dt_61_23->Fit("func23", "R");
        h_dt_61_23->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_23->GetYaxis()->SetTitle("Counts");
        // h_dt_61_23->SetMaximum(70);
        h_dt_61_23->Draw("E1");
        func23->Draw("same");
        func23->SetLineColor(2);
                
        Double_t param23_1 = func23->GetParameter(0);
        Double_t param23_2 = func23->GetParameter(1);
        Double_t param23_3 = func23->GetParameter(2);
        Double_t param23_4 = func23->GetParameter(3);
        Double_t error23_1 = func23->GetParError(0);
        Double_t error23_2 = func23->GetParError(1);
        Double_t error23_3 = func23->GetParError(2);
        Double_t error23_4 = func23->GetParError(3);
        Double_t chi2_23 = func23->GetChisquare();

        TLatex *lat_dt_23_1 = new TLatex(0.3, 0.88, Form("Baseline = %.2f +/- %.2f", param23_1, error23_1));
        lat_dt_23_1->SetNDC();
        lat_dt_23_1->SetTextColor(1);
        lat_dt_23_1->SetTextSize(0.03);
        lat_dt_23_1->Draw();

        TLatex *lat_dt_23_2 = new TLatex(0.3, 0.84, Form("Amplitude = %.2f +/- %.2f", param23_2, error23_2));
        lat_dt_23_2->SetNDC();
        lat_dt_23_2->SetTextColor(1);
        lat_dt_23_2->SetTextSize(0.03);
        lat_dt_23_2->Draw();

        TLatex *lat_dt_23_3 = new TLatex(0.3, 0.80, Form("Peak Center = %.2f +/- %.2f", param23_3, error23_3));
        lat_dt_23_3->SetNDC();
        lat_dt_23_3->SetTextColor(1);
        lat_dt_23_3->SetTextSize(0.03);
        lat_dt_23_3->Draw();

        TLatex *lat_dt_23_4 = new TLatex(0.3, 0.76, Form("Standard Deviation = %.2f +/- %.2f", param23_4, error23_4));
        lat_dt_23_4->SetNDC();
        lat_dt_23_4->SetTextColor(1);
        lat_dt_23_4->SetTextSize(0.03);
        lat_dt_23_4->Draw();

        TLatex *lat_dt_23_5 = new TLatex(0.3, 0.72, Form("Chi-Squared = %.2f", chi2_23));
        lat_dt_23_5->SetNDC();
        lat_dt_23_5->SetTextColor(1);
        lat_dt_23_5->SetTextSize(0.03);
        lat_dt_23_5->Draw();

        c_dt_61_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_23.png");
        // c_dt_61_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_23.png");
        c_dt_61_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_23.png"));
        h_dt_61_23->Reset();

        TCanvas* c_dt_61_24 = new TCanvas("c_dt_61_24", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_24->SetLogy();
        c_dt_61_24->cd();
        TF1* func24 = new TF1("func24", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -22000, -16000);
        func24->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func24->SetParameters(5, 20, -19000, 100);          // 5, 20, -19000, 100);          // 160, 20, 2200, 200);
        h_dt_61_24->Fit("func24", "R");
        h_dt_61_24->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_24->GetYaxis()->SetTitle("Counts");
        // h_dt_61_24->SetMaximum(280);
        h_dt_61_24->Draw("E1");
        func24->Draw("same");
        func24->SetLineColor(2);
                
        Double_t param24_1 = func24->GetParameter(0);
        Double_t param24_2 = func24->GetParameter(1);
        Double_t param24_3 = func24->GetParameter(2);
        Double_t param24_4 = func24->GetParameter(3);
        Double_t error24_1 = func24->GetParError(0);
        Double_t error24_2 = func24->GetParError(1);
        Double_t error24_3 = func24->GetParError(2);
        Double_t error24_4 = func24->GetParError(3);
        Double_t chi2_24 = func24->GetChisquare();

        TLatex *lat_dt_24_1 = new TLatex(0.3, 0.88, Form("Baseline = %.2f +/- %.2f", param24_1, error24_1));
        lat_dt_24_1->SetNDC();
        lat_dt_24_1->SetTextColor(1);
        lat_dt_24_1->SetTextSize(0.03);
        lat_dt_24_1->Draw();

        TLatex *lat_dt_24_2 = new TLatex(0.3, 0.84, Form("Amplitude = %.2f +/- %.2f", param24_2, error24_2));
        lat_dt_24_2->SetNDC();
        lat_dt_24_2->SetTextColor(1);
        lat_dt_24_2->SetTextSize(0.03);
        lat_dt_24_2->Draw();

        TLatex *lat_dt_24_3 = new TLatex(0.3, 0.80, Form("Peak Center = %.2f +/- %.2f", param24_3, error24_3));
        lat_dt_24_3->SetNDC();
        lat_dt_24_3->SetTextColor(1);
        lat_dt_24_3->SetTextSize(0.03);
        lat_dt_24_3->Draw();

        TLatex *lat_dt_24_4 = new TLatex(0.3, 0.76, Form("Standard Deviation = %.2f +/- %.2f", param24_4, error24_4));
        lat_dt_24_4->SetNDC();
        lat_dt_24_4->SetTextColor(1);
        lat_dt_24_4->SetTextSize(0.03);
        lat_dt_24_4->Draw();

        TLatex *lat_dt_24_5 = new TLatex(0.3, 0.72, Form("Chi-Squared = %.2f", chi2_24));
        lat_dt_24_5->SetNDC();
        lat_dt_24_5->SetTextColor(1);
        lat_dt_24_5->SetTextSize(0.03);
        lat_dt_24_5->Draw();

        c_dt_61_24->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_24.png");
        // c_dt_61_24->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_24.png");
        c_dt_61_24->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_24.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_24.png"));
        h_dt_61_24->Reset();

        TCanvas* c_dt_61_25 = new TCanvas("c_dt_61_25", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_25->SetLogy();
        c_dt_61_25->cd();
        TF1* func25 = new TF1("func25", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -22000, -16000);
        func25->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func25->SetParameters(5, 20, -19000, 100);          // 5, 20, -19000, 100);          // 160, 20, 2200, 200);
        h_dt_61_25->Fit("func25", "R");
        h_dt_61_25->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_25->GetYaxis()->SetTitle("Counts");
        // h_dt_61_25->SetMaximum(280);
        h_dt_61_25->Draw("E1");
        func25->Draw("same");
        func25->SetLineColor(2);
                
        Double_t param25_1 = func25->GetParameter(0);
        Double_t param25_2 = func25->GetParameter(1);
        Double_t param25_3 = func25->GetParameter(2);
        Double_t param25_4 = func25->GetParameter(3);
        Double_t error25_1 = func25->GetParError(0);
        Double_t error25_2 = func25->GetParError(1);
        Double_t error25_3 = func25->GetParError(2);
        Double_t error25_4 = func25->GetParError(3);
        Double_t chi2_25 = func25->GetChisquare();

        TLatex *lat_dt_25_1 = new TLatex(0.3, 0.88, Form("Baseline = %.2f +/- %.2f", param25_1, error25_1));
        lat_dt_25_1->SetNDC();
        lat_dt_25_1->SetTextColor(1);
        lat_dt_25_1->SetTextSize(0.03);
        lat_dt_25_1->Draw();

        TLatex *lat_dt_25_2 = new TLatex(0.3, 0.84, Form("Amplitude = %.2f +/- %.2f", param25_2, error25_2));
        lat_dt_25_2->SetNDC();
        lat_dt_25_2->SetTextColor(1);
        lat_dt_25_2->SetTextSize(0.03);
        lat_dt_25_2->Draw();

        TLatex *lat_dt_25_3 = new TLatex(0.3, 0.80, Form("Peak Center = %.2f +/- %.2f", param25_3, error25_3));
        lat_dt_25_3->SetNDC();
        lat_dt_25_3->SetTextColor(1);
        lat_dt_25_3->SetTextSize(0.03);
        lat_dt_25_3->Draw();

        TLatex *lat_dt_25_4 = new TLatex(0.3, 0.76, Form("Standard Deviation = %.2f +/- %.2f", param25_4, error25_4));
        lat_dt_25_4->SetNDC();
        lat_dt_25_4->SetTextColor(1);
        lat_dt_25_4->SetTextSize(0.03);
        lat_dt_25_4->Draw();

        TLatex *lat_dt_25_5 = new TLatex(0.3, 0.72, Form("Chi-Squared = %.2f", chi2_25));
        lat_dt_25_5->SetNDC();
        lat_dt_25_5->SetTextColor(1);
        lat_dt_25_5->SetTextSize(0.03);
        lat_dt_25_5->Draw();

        c_dt_61_25->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_25.png");
        // c_dt_61_25->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_25.png");
        c_dt_61_25->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_25.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_25.png"));
        h_dt_61_25->Reset();

        TCanvas* c_dt_61_top = new TCanvas("c_dt_61_top", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_top->SetLogy();
        c_dt_61_top->cd();
        TF1* functop = new TF1("functop", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -22000, -16000);
        functop->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        functop->SetParameters(150, 50, -19000, 100);         // 150, 50, -19000, 100);          // 130, 30, 2100, 200);
        h_dt_61_top->Fit("functop", "R");
        h_dt_61_top->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_top->GetYaxis()->SetTitle("Counts");
        // h_dt_61_top->SetMaximum(280);
        h_dt_61_top->Draw("E1");
        functop->Draw("same");
        functop->SetLineColor(2);
                
        Double_t paramtop_1 = functop->GetParameter(0);
        Double_t paramtop_2 = functop->GetParameter(1);
        Double_t paramtop_3 = functop->GetParameter(2);
        Double_t paramtop_4 = functop->GetParameter(3);
        Double_t errortop_1 = functop->GetParError(0);
        Double_t errortop_2 = functop->GetParError(1);
        Double_t errortop_3 = functop->GetParError(2);
        Double_t errortop_4 = functop->GetParError(3);
        Double_t chi2_top = functop->GetChisquare();

        TLatex *lat_dt_top_1 = new TLatex(0.3, 0.88, Form("Baseline = %.2f +/- %.2f", paramtop_1, errortop_1));
        lat_dt_top_1->SetNDC();
        lat_dt_top_1->SetTextColor(1);
        lat_dt_top_1->SetTextSize(0.03);
        lat_dt_top_1->Draw();

        TLatex *lat_dt_top_2 = new TLatex(0.3, 0.84, Form("Amplitude = %.2f +/- %.2f", paramtop_2, errortop_2));
        lat_dt_top_2->SetNDC();
        lat_dt_top_2->SetTextColor(1);
        lat_dt_top_2->SetTextSize(0.03);
        lat_dt_top_2->Draw();

        TLatex *lat_dt_top_3 = new TLatex(0.3, 0.80, Form("Peak Center = %.2f +/- %.2f", paramtop_3, errortop_3));
        lat_dt_top_3->SetNDC();
        lat_dt_top_3->SetTextColor(1);
        lat_dt_top_3->SetTextSize(0.03);
        lat_dt_top_3->Draw();

        TLatex *lat_dt_top_4 = new TLatex(0.3, 0.76, Form("Standard Deviation = %.2f +/- %.2f", paramtop_4, errortop_4));
        lat_dt_top_4->SetNDC();
        lat_dt_top_4->SetTextColor(1);
        lat_dt_top_4->SetTextSize(0.03);
        lat_dt_top_4->Draw();

        TLatex *lat_dt_top_5 = new TLatex(0.3, 0.72, Form("Chi-Squared = %.2f", chi2_top));
        lat_dt_top_5->SetNDC();
        lat_dt_top_5->SetTextColor(1);
        lat_dt_top_5->SetTextSize(0.03);
        lat_dt_top_5->Draw();

        c_dt_61_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_top.png");
        // c_dt_61_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_61_top.png");
        c_dt_61_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_top.png"));
        h_dt_61_top->Reset();

        TCanvas* c_vpint_16 = new TCanvas("c_vpint_16", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_16->SetLogy();
        c_vpint_16->cd();
        h_vpint_16_all->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_16_all->GetYaxis()->SetTitle("Counts");
        h_vpint_16_all->Draw("same");
        h_vpint_16_all->SetLineColor(kBlue);
        h_svpint_not_16->Draw("same");
        h_svpint_not_16->SetLineColor(kRed);
        h_vpint_16->Draw("same");
        h_vpint_16->SetLineColor(kGreen);
        TLegend *leg_vpint_16 = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_vpint_16->AddEntry(h_vpint_16, "ONLY SiPM1 Events", "l");
        leg_vpint_16->AddEntry(h_svpint_not_16, "NOT SiPM1 Events", "l");
        leg_vpint_16->AddEntry(h_vpint_16_all, "ALL SiPM1 Events", "l");
        leg_vpint_16->Draw();
        c_vpint_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_16.png");
        // c_vpint_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_16.png");
        c_vpint_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_16.png"));
        h_vpint_16->Reset();

        TCanvas* c_vpint_17 = new TCanvas("c_vpint_17", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_17->SetLogy();
        c_vpint_17->cd();
        h_vpint_17_all->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_17_all->GetYaxis()->SetTitle("Counts");
        h_vpint_17_all->Draw("same");
        h_vpint_17_all->SetLineColor(kBlue);
        h_svpint_not_17->Draw("same");
        h_svpint_not_17->SetLineColor(kRed);
        h_vpint_17->Draw("same");
        h_vpint_17->SetLineColor(kGreen);
        TLegend *leg_vpint_17 = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_vpint_17->AddEntry(h_vpint_17, "ONLY SiPM2 Events", "l");
        leg_vpint_17->AddEntry(h_svpint_not_17, "NOT SiPM2 Events", "l");
        leg_vpint_17->AddEntry(h_vpint_17_all, "ALL SiPM2 Events", "l");
        leg_vpint_17->Draw();
        c_vpint_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_17.png");
        // c_vpint_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_17.png");
        c_vpint_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_17.png"));
        h_vpint_17->Reset();

        TCanvas* c_vpint_18 = new TCanvas("c_vpint_18", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_18->SetLogy();
        c_vpint_18->cd();
        h_vpint_18_all->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_18_all->GetYaxis()->SetTitle("Counts");
        h_vpint_18_all->Draw("same");
        h_vpint_18_all->SetLineColor(kBlue);
        h_svpint_not_18->Draw("same");
        h_svpint_not_18->SetLineColor(kRed);
        h_vpint_18->Draw("same");
        h_vpint_18->SetLineColor(kGreen);
        TLegend *leg_vpint_18 = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_vpint_18->AddEntry(h_vpint_18, "ONLY SiPM3 Events", "l");
        leg_vpint_18->AddEntry(h_svpint_not_18, "NOT SiPM3 Events", "l");
        leg_vpint_18->AddEntry(h_vpint_18_all, "ALL SiPM3 Events", "l");
        leg_vpint_18->Draw();
        c_vpint_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_18.png");
        // c_vpint_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_18.png");
        c_vpint_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_18.png"));
        h_vpint_18->Reset();

        TCanvas* c_vpint_19 = new TCanvas("c_vpint_19", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_19->SetLogy();
        c_vpint_19->cd();
        h_vpint_19_all->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_19_all->GetYaxis()->SetTitle("Counts");
        h_vpint_19_all->Draw("same");
        h_vpint_19_all->SetLineColor(kBlue);
        h_svpint_not_19->Draw("same");
        h_svpint_not_19->SetLineColor(kRed);
        h_vpint_19->Draw("same");
        h_vpint_19->SetLineColor(kGreen);
        TLegend *leg_vpint_19 = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_vpint_19->AddEntry(h_vpint_19, "ONLY SiPM4 Events", "l");
        leg_vpint_19->AddEntry(h_svpint_not_19, "NOT SiPM4 Events", "l");
        leg_vpint_19->AddEntry(h_vpint_19_all, "ALL SiPM4 Events", "l");
        leg_vpint_19->Draw();
        c_vpint_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_19.png");
        // c_vpint_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_19.png");
        c_vpint_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_19.png"));
        h_vpint_19->Reset();

        TCanvas* c_vpint_20 = new TCanvas("c_vpint_20", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_20->SetLogy();
        c_vpint_20->cd();
        h_vpint_20_all->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_20_all->GetYaxis()->SetTitle("Counts");
        h_vpint_20_all->Draw("same");
        h_vpint_20_all->SetLineColor(kBlue);
        h_svpint_not_20->Draw("same");
        h_svpint_not_20->SetLineColor(kRed);
        h_vpint_20->Draw("same");
        h_vpint_20->SetLineColor(kGreen);
        TLegend *leg_vpint_20 = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_vpint_20->AddEntry(h_vpint_20, "ONLY SiPM5 Events", "l");
        leg_vpint_20->AddEntry(h_svpint_not_20, "NOT SiPM5 Events", "l");
        leg_vpint_20->AddEntry(h_vpint_20_all, "ALL SiPM5 Events", "l");
        leg_vpint_20->Draw();
        c_vpint_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_20.png");
        // c_vpint_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_20.png");
        c_vpint_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_20.png"));
        h_vpint_20->Reset();
        h_svpint_not_20->Reset();

        TCanvas* c_vpint_21 = new TCanvas("c_vpint_21", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_21->SetLogy();
        c_vpint_21->cd();
        h_vpint_21_all->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_21_all->GetYaxis()->SetTitle("Counts");
        h_vpint_21_all->Draw("same");
        h_vpint_21_all->SetLineColor(kBlue);
        h_svpint_not_21->Draw("same");
        h_svpint_not_21->SetLineColor(kRed);
        h_vpint_21->Draw("same");
        h_vpint_21->SetLineColor(kGreen);
        TLegend *leg_vpint_21 = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_vpint_21->AddEntry(h_vpint_21, "ONLY SiPM6 Events", "l");
        leg_vpint_21->AddEntry(h_svpint_not_21, "NOT SiPM6 Events", "l");
        leg_vpint_21->AddEntry(h_vpint_21_all, "ALL SiPM6 Events", "l");
        leg_vpint_21->Draw();
        c_vpint_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_21.png");
        // c_vpint_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_21.png");
        c_vpint_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_21.png"));
        h_vpint_21->Reset();
        h_svpint_not_21->Reset();

        TCanvas* c_vpint_22 = new TCanvas("c_vpint_22", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_22->SetLogy();
        c_vpint_22->cd();
        h_vpint_22_all->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_22_all->GetYaxis()->SetTitle("Counts");
        h_vpint_22_all->Draw("same");
        h_vpint_22_all->SetLineColor(kBlue);
        h_svpint_not_22->Draw("same");
        h_svpint_not_22->SetLineColor(kRed);
        h_vpint_22->Draw("same");
        h_vpint_22->SetLineColor(kGreen);
        TLegend *leg_vpint_22 = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_vpint_22->AddEntry(h_vpint_22, "ONLY SiPM7 Events", "l");
        leg_vpint_22->AddEntry(h_svpint_not_22, "NOT SiPM7 Events", "l");
        leg_vpint_22->AddEntry(h_vpint_22_all, "ALL SiPM7 Events", "l");
        leg_vpint_22->Draw();
        c_vpint_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_22.png");
        // c_vpint_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_22.png");
        c_vpint_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_22.png"));
        h_vpint_22->Reset();
        h_svpint_not_22->Reset();

        TCanvas* c_vpint_23 = new TCanvas("c_vpint_23", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_23->SetLogy();
        c_vpint_23->cd();
        h_vpint_23_all->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_23_all->GetYaxis()->SetTitle("Counts");
        h_vpint_23_all->Draw("same");
        h_vpint_23_all->SetLineColor(kBlue);
        h_svpint_not_23->Draw("same");
        h_svpint_not_23->SetLineColor(kRed);
        h_vpint_23->Draw("same");
        h_vpint_23->SetLineColor(kGreen);
        TLegend *leg_vpint_23 = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_vpint_23->AddEntry(h_vpint_23, "ONLY SiPM8 Events", "l");
        leg_vpint_23->AddEntry(h_svpint_not_23, "NOT SiPM8 Events", "l");
        leg_vpint_23->AddEntry(h_vpint_23_all, "ALL SiPM8 Events", "l");
        leg_vpint_23->Draw();
        c_vpint_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_23.png");
        // c_vpint_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_23.png");
        c_vpint_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_23.png"));
        h_vpint_23->Reset();
        h_svpint_not_23->Reset();

        TCanvas* c_vpint_24 = new TCanvas("c_vpint_24", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_24->SetLogy();
        c_vpint_24->cd();
        h_vpint_24->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_24->GetYaxis()->SetTitle("Counts");
        h_vpint_24->Draw();
        c_vpint_24->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_24.png");
        // c_vpint_24->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_24.png");
        c_vpint_24->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_24.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_24.png"));
        h_vpint_24->Reset();

        TCanvas* c_vpint_25 = new TCanvas("c_vpint_25", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_25->SetLogy();
        c_vpint_25->cd();
        h_vpint_25->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_25->GetYaxis()->SetTitle("Counts");
        h_vpint_25->Draw();
        c_vpint_25->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_25.png");
        // c_vpint_25->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_25.png");
        c_vpint_25->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_25.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_25.png"));
        h_vpint_25->Reset();

        TCanvas* c_vpint_top = new TCanvas("c_vpint_top", "Event 61 and SiPM dt Integral Values", 900, 700);
        c_vpint_top->SetLogy();
        c_vpint_top->cd();
        h_vpint_top_all->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_top_all->GetYaxis()->SetTitle("Counts");
        h_vpint_top_all->Draw("same");
        h_vpint_top_all->SetLineColor(kBlue);
        h_tvpint_not_top->Draw("same");
        h_tvpint_not_top->SetLineColor(kRed);
        h_vpint_top->Draw("same");
        h_vpint_top->SetLineColor(kGreen);
        TLegend *leg_vpint_top = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_vpint_top->AddEntry(h_vpint_top, "ONLY SiPM9&10 Events", "l");
        leg_vpint_top->AddEntry(h_tvpint_not_top, "NOT SiPM9&10 Events", "l");
        leg_vpint_top->AddEntry(h_vpint_top_all, "ALL SiPM9&10 Events", "l");
        leg_vpint_top->Draw();
        c_vpint_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_top.png");
        // c_vpint_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_top.png");
        c_vpint_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_top.png"));
        h_vpint_top->Reset();
        h_tvpint_not_top->Reset();

        TCanvas* c_BRN_tail_amp = new TCanvas("c_BRN_tail_amp", "BRN Pulse Tail Amplitude", 900, 700);
        c_BRN_tail_amp->SetLogy();
        c_BRN_tail_amp->cd();
        h_BRN_tail_amp->GetXaxis()->SetTitle("Amplitude (ADC)");
        h_BRN_tail_amp->GetYaxis()->SetTitle("Counts");
        h_BRN_tail_amp->Draw();
        c_BRN_tail_amp->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_BRN_tail_amp.png");
        // c_BRN_tail_amp->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_BRN_tail_amp.png");
        c_BRN_tail_amp->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_BRN_tail_amp.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_BRN_tail_amp.png"));
        h_BRN_tail_amp->Reset();

        TCanvas* c_background_tail_amp = new TCanvas("c_background_tail_amp", "Background Pulse Tail Amplitude", 900, 700);
        c_background_tail_amp->SetLogy();
        c_background_tail_amp->cd();
        h_background_tail_amp->GetXaxis()->SetTitle("Amplitude (ADC)");
        h_background_tail_amp->GetYaxis()->SetTitle("Counts");
        h_background_tail_amp->Draw();
        c_background_tail_amp->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_background_tail_amp.png");
        // c_background_tail_amp->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_background_tail_amp.png");
        c_background_tail_amp->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_background_tail_amp.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_background_tail_amp.png"));
        h_background_tail_amp->Reset();

        TCanvas* c_vpint_16_BRNpeak = new TCanvas("c_vpint_16_BRNpeak", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_16_BRNpeak->SetLogy();
        c_vpint_16_BRNpeak->cd();
        h_vpint_16_BRNpeak->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_16_BRNpeak->GetYaxis()->SetTitle("Counts");
        h_vpint_16_BRNpeak->Draw();
        c_vpint_16_BRNpeak->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_16_BRNpeak.png");
        // c_vpint_16_BRNpeak->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_16_BRNpeak.png");
        c_vpint_16_BRNpeak->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_16_BRNpeak.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_16_BRNpeak.png"));
        h_vpint_16_BRNpeak->Reset();

        TCanvas* c_vpint_17_BRNpeak = new TCanvas("c_vpint_17_BRNpeak", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_17_BRNpeak->SetLogy();
        c_vpint_17_BRNpeak->cd();
        h_vpint_17_BRNpeak->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_17_BRNpeak->GetYaxis()->SetTitle("Counts");
        h_vpint_17_BRNpeak->Draw();
        c_vpint_17_BRNpeak->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_17_BRNpeak.png");
        // c_vpint_17_BRNpeak->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_17_BRNpeak.png");
        c_vpint_17_BRNpeak->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_17_BRNpeak.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_17_BRNpeak.png"));
        h_vpint_17_BRNpeak->Reset();

        TCanvas* c_vpint_18_BRNpeak = new TCanvas("c_vpint_18_BRNpeak", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_18_BRNpeak->SetLogy();
        c_vpint_18_BRNpeak->cd();
        h_vpint_18_BRNpeak->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_18_BRNpeak->GetYaxis()->SetTitle("Counts");
        h_vpint_18_BRNpeak->Draw();
        c_vpint_18_BRNpeak->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_18_BRNpeak.png");
        // c_vpint_18_BRNpeak->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_18_BRNpeak.png");
        c_vpint_18_BRNpeak->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_18_BRNpeak.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_18_BRNpeak.png"));
        h_vpint_18_BRNpeak->Reset();

        TCanvas* c_vpint_19_BRNpeak = new TCanvas("c_vpint_19_BRNpeak", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_19_BRNpeak->SetLogy();
        c_vpint_19_BRNpeak->cd();
        h_vpint_19_BRNpeak->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_19_BRNpeak->GetYaxis()->SetTitle("Counts");
        h_vpint_19_BRNpeak->Draw();
        c_vpint_19_BRNpeak->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_19_BRNpeak.png");
        // c_vpint_19_BRNpeak->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_19_BRNpeak.png");
        c_vpint_19_BRNpeak->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_19_BRNpeak.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_19_BRNpeak.png"));
        h_vpint_19_BRNpeak->Reset();

        TCanvas* c_vpint_20_BRNpeak = new TCanvas("c_vpint_20_BRNpeak", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_20_BRNpeak->SetLogy();
        c_vpint_20_BRNpeak->cd();
        h_vpint_20_BRNpeak->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_20_BRNpeak->GetYaxis()->SetTitle("Counts");
        h_vpint_20_BRNpeak->Draw();
        c_vpint_20_BRNpeak->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_20_BRNpeak.png");
        // c_vpint_20_BRNpeak->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_20_BRNpeak.png");
        c_vpint_20_BRNpeak->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_20_BRNpeak.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_20_BRNpeak.png"));
        h_vpint_20_BRNpeak->Reset();

        TCanvas* c_vpint_21_BRNpeak = new TCanvas("c_vpint_21_BRNpeak", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_21_BRNpeak->SetLogy();
        c_vpint_21_BRNpeak->cd();
        h_vpint_21_BRNpeak->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_21_BRNpeak->GetYaxis()->SetTitle("Counts");
        h_vpint_21_BRNpeak->Draw("same");
        h_vpint_21_BRNpeak->SetLineColor(kGreen);
        h_vpint_21_background->Draw("same");
        h_vpint_21_background->SetLineColor(kRed);
        TLegend *leg_vpint_21_BRNs = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_vpint_21_BRNs->AddEntry(h_vpint_21_BRNpeak, "BRN Peak Sample", "l");
        leg_vpint_21_BRNs->AddEntry(h_vpint_21_background, "Background Sample", "l");
        leg_vpint_21_BRNs->Draw();
        c_vpint_21_BRNpeak->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_21_BRNpeak.png");
        // c_vpint_21_BRNpeak->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_21_BRNpeak.png");
        c_vpint_21_BRNpeak->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_21_BRNpeak.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_21_BRNpeak.png"));
        h_vpint_21_BRNpeak->Reset();

        TCanvas* c_vpint_22_BRNpeak = new TCanvas("c_vpint_22_BRNpeak", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_22_BRNpeak->SetLogy();
        c_vpint_22_BRNpeak->cd();
        h_vpint_22_BRNpeak->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_22_BRNpeak->GetYaxis()->SetTitle("Counts");
        h_vpint_22_BRNpeak->Draw();
        c_vpint_22_BRNpeak->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_22_BRNpeak.png");
        // c_vpint_22_BRNpeak->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_22_BRNpeak.png");
        c_vpint_22_BRNpeak->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_22_BRNpeak.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_22_BRNpeak.png"));
        h_vpint_22_BRNpeak->Reset();

        TCanvas* c_vpint_23_BRNpeak = new TCanvas("c_vpint_23_BRNpeak", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_23_BRNpeak->SetLogy();
        c_vpint_23_BRNpeak->cd();
        h_vpint_23_BRNpeak->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_23_BRNpeak->GetYaxis()->SetTitle("Counts");
        h_vpint_23_BRNpeak->Draw();
        c_vpint_23_BRNpeak->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_23_BRNpeak.png");
        // c_vpint_23_BRNpeak->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_23_BRNpeak.png");
        c_vpint_23_BRNpeak->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_23_BRNpeak.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_23_BRNpeak.png"));
        h_vpint_23_BRNpeak->Reset();

        TCanvas* c_vpint_top_BRNpeak = new TCanvas("c_vpint_top_BRNpeak", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_top_BRNpeak->SetLogy();
        c_vpint_top_BRNpeak->cd();
        h_vpint_top_BRNpeak->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_top_BRNpeak->GetYaxis()->SetTitle("Counts");
        h_vpint_top_BRNpeak->Draw();
        c_vpint_top_BRNpeak->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_top_BRNpeak.png");
        // c_vpint_top_BRNpeak->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_top_BRNpeak.png");
        c_vpint_top_BRNpeak->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_top_BRNpeak.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_top_BRNpeak.png"));
        h_vpint_top_BRNpeak->Reset();

        TCanvas* c_vpint_16_background = new TCanvas("c_vpint_16_background", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_16_background->SetLogy();
        c_vpint_16_background->cd();
        h_vpint_16_background->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_16_background->GetYaxis()->SetTitle("Counts");
        h_vpint_16_background->Draw();
        c_vpint_16_background->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_16_background.png");
        // c_vpint_16_background->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_16_background.png");
        c_vpint_16_background->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_16_background.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_16_background.png"));
        h_vpint_16_background->Reset();

        TCanvas* c_vpint_17_background = new TCanvas("c_vpint_17_background", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_17_background->SetLogy();
        c_vpint_17_background->cd();
        h_vpint_17_background->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_17_background->GetYaxis()->SetTitle("Counts");
        h_vpint_17_background->Draw();
        c_vpint_17_background->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_17_background.png");
        // c_vpint_17_background->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_17_background.png");
        c_vpint_17_background->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_17_background.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_17_background.png"));
        h_vpint_17_background->Reset();

        TCanvas* c_vpint_18_background = new TCanvas("c_vpint_18_background", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_18_background->SetLogy();
        c_vpint_18_background->cd();
        h_vpint_18_background->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_18_background->GetYaxis()->SetTitle("Counts");
        h_vpint_18_background->Draw();
        c_vpint_18_background->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_18_background.png");
        // c_vpint_18_background->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_18_background.png");
        c_vpint_18_background->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_18_background.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_18_background.png"));
        h_vpint_18_background->Reset();

        TCanvas* c_vpint_19_background = new TCanvas("c_vpint_19_background", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_19_background->SetLogy();
        c_vpint_19_background->cd();
        h_vpint_19_background->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_19_background->GetYaxis()->SetTitle("Counts");
        h_vpint_19_background->Draw();
        c_vpint_19_background->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_19_background.png");
        // c_vpint_19_background->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_19_background.png");
        c_vpint_19_background->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_19_background.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_19_background.png"));
        h_vpint_19_background->Reset();

        TCanvas* c_vpint_20_background = new TCanvas("c_vpint_20_background", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_20_background->SetLogy();
        c_vpint_20_background->cd();
        h_vpint_20_background->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_20_background->GetYaxis()->SetTitle("Counts");
        h_vpint_20_background->Draw();
        c_vpint_20_background->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_20_background.png");
        // c_vpint_20_background->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_20_background.png");
        c_vpint_20_background->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_20_background.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_20_background.png"));
        h_vpint_20_background->Reset();

        TCanvas* c_vpint_21_background = new TCanvas("c_vpint_21_background", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_21_background->SetLogy();
        c_vpint_21_background->cd();
        h_vpint_21_background->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_21_background->GetYaxis()->SetTitle("Counts");
        h_vpint_21_background->Draw();
        c_vpint_21_background->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_21_background.png");
        // c_vpint_21_background->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_21_background.png");
        c_vpint_21_background->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_21_background.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_21_background.png"));
        h_vpint_21_background->Reset();

        TCanvas* c_vpint_22_background = new TCanvas("c_vpint_22_background", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_22_background->SetLogy();
        c_vpint_22_background->cd();
        h_vpint_22_background->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_22_background->GetYaxis()->SetTitle("Counts");
        h_vpint_22_background->Draw();
        c_vpint_22_background->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_22_background.png");
        // c_vpint_22_background->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_22_background.png");
        c_vpint_22_background->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_22_background.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_22_background.png"));
        h_vpint_22_background->Reset();

        TCanvas* c_vpint_23_background = new TCanvas("c_vpint_23_background", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_23_background->SetLogy();
        c_vpint_23_background->cd();
        h_vpint_23_background->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_23_background->GetYaxis()->SetTitle("Counts");
        h_vpint_23_background->Draw();
        c_vpint_23_background->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_23_background.png");
        // c_vpint_23_background->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_23_background.png");
        c_vpint_23_background->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_23_background.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_23_background.png"));
        h_vpint_23_background->Reset();

        TCanvas* c_vpint_top_background = new TCanvas("c_vpint_top_background", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_top_background->SetLogy();
        c_vpint_top_background->cd();
        h_vpint_top_background->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_top_background->GetYaxis()->SetTitle("Counts");
        h_vpint_top_background->Draw();
        c_vpint_top_background->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vpint_top_background.png");
        // c_vpint_top_background->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vpint_top_background.png");
        c_vpint_top_background->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_top_background.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_top_background.png"));
        h_vpint_top_background->Reset();

        TCanvas* c_BRNint_16 = new TCanvas("c_BRNint_16", "BRN Integral Values", 900, 700);
        // c_BRNint_16->SetLogy();
        c_BRNint_16->cd();
        h_BRNint_16->GetXaxis()->SetTitle("Integral (ADC)");
        h_BRNint_16->GetYaxis()->SetTitle("Counts");
        h_BRNint_16->Draw("E1");
        c_BRNint_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_BRNint_16.png");
        // c_BRNint_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_BRNint_16.png");
        c_BRNint_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_BRNint_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_BRNint_16.png"));
        h_BRNint_16->Reset();

        TCanvas* c_BRNint_17 = new TCanvas("c_BRNint_17", "BRN Integral Values", 900, 700);
        // c_BRNint_17->SetLogy();
        c_BRNint_17->cd();
        h_BRNint_17->GetXaxis()->SetTitle("Integral (ADC)");
        h_BRNint_17->GetYaxis()->SetTitle("Counts");
        h_BRNint_17->Draw("E1");
        c_BRNint_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_BRNint_17.png");
        // c_BRNint_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_BRNint_17.png");
        c_BRNint_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_BRNint_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_BRNint_17.png"));
        h_BRNint_17->Reset();

        TCanvas* c_BRNint_18 = new TCanvas("c_BRNint_18", "BRN Integral Values", 900, 700);
        // c_BRNint_18->SetLogy();
        c_BRNint_18->cd();
        h_BRNint_18->GetXaxis()->SetTitle("Integral (ADC)");
        h_BRNint_18->GetYaxis()->SetTitle("Counts");
        h_BRNint_18->Draw("E1");
        c_BRNint_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_BRNint_18.png");
        // c_BRNint_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_BRNint_18.png");
        c_BRNint_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_BRNint_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_BRNint_18.png"));
        h_BRNint_18->Reset();

        TCanvas* c_BRNint_19 = new TCanvas("c_BRNint_19", "BRN Integral Values", 900, 700);
        // c_BRNint_19->SetLogy();
        c_BRNint_19->cd();
        h_BRNint_19->GetXaxis()->SetTitle("Integral (ADC)");
        h_BRNint_19->GetYaxis()->SetTitle("Counts");
        h_BRNint_19->Draw("E1");
        c_BRNint_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_BRNint_19.png");
        // c_BRNint_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_BRNint_19.png");
        c_BRNint_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_BRNint_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_BRNint_19.png"));
        h_BRNint_19->Reset();

        TCanvas* c_BRNint_20 = new TCanvas("c_BRNint_20", "BRN Integral Values", 900, 700);
        // c_BRNint_20->SetLogy();
        c_BRNint_20->cd();
        h_BRNint_20->GetXaxis()->SetTitle("Integral (ADC)");
        h_BRNint_20->GetYaxis()->SetTitle("Counts");
        h_BRNint_20->Draw("E1");
        c_BRNint_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_BRNint_20.png");
        // c_BRNint_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_BRNint_20.png");
        c_BRNint_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_BRNint_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_BRNint_20.png"));
        h_BRNint_20->Reset();

        TCanvas* c_BRNint_21 = new TCanvas("c_BRNint_21", "BRN Integral Values", 900, 700);
        // c_BRNint_21->SetLogy();
        c_BRNint_21->cd();
        h_BRNint_21->GetXaxis()->SetTitle("Integral (ADC)");
        h_BRNint_21->GetYaxis()->SetTitle("Counts");
        h_BRNint_21->Draw("E1");
        c_BRNint_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_BRNint_21.png");
        // c_BRNint_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_BRNint_21.png");
        c_BRNint_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_BRNint_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_BRNint_21.png"));
        h_BRNint_21->Reset();

        TCanvas* c_BRNint_22 = new TCanvas("c_BRNint_22", "BRN Integral Values", 900, 700);
        // c_BRNint_22->SetLogy();
        c_BRNint_22->cd();
        h_BRNint_22->GetXaxis()->SetTitle("Integral (ADC)");
        h_BRNint_22->GetYaxis()->SetTitle("Counts");
        h_BRNint_22->Draw("E1");
        c_BRNint_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_BRNint_22.png");
        // c_BRNint_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_BRNint_22.png");
        c_BRNint_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_BRNint_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_BRNint_22.png"));
        h_BRNint_22->Reset();

        TCanvas* c_BRNint_23 = new TCanvas("c_BRNint_23", "BRN Integral Values", 900, 700);
        // c_BRNint_23->SetLogy();
        c_BRNint_23->cd();
        h_BRNint_23->GetXaxis()->SetTitle("Integral (ADC)");
        h_BRNint_23->GetYaxis()->SetTitle("Counts");
        h_BRNint_23->Draw("E1");
        c_BRNint_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_BRNint_23.png");
        // c_BRNint_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_BRNint_23.png");
        c_BRNint_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_BRNint_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_BRNint_23.png"));
        h_BRNint_23->Reset();

        TCanvas* c_BRNint_top = new TCanvas("c_BRNint_top", "BRN Integral Values", 900, 700);
        // c_BRNint_top->SetLogy();
        c_BRNint_top->cd();
        h_BRNint_top->GetXaxis()->SetTitle("Integral (ADC)");
        h_BRNint_top->GetYaxis()->SetTitle("Counts");
        h_BRNint_top->Draw("E1");
        c_BRNint_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_BRNint_top.png");
        // c_BRNint_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_BRNint_top.png");
        c_BRNint_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_BRNint_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_BRNint_top.png"));
        h_BRNint_top->Reset();

        TCanvas* c_detint = new TCanvas("c_detint", "Distribution of All Ev61Det Integral Values", 1200, 700);
        c_detint->SetLogy();
        c_detint->cd();
        h_detint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_detint->GetYaxis()->SetTitle("Counts");
        h_detint->Draw();
        c_detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint.png");
        // c_detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_detint.png");
        c_detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint.png"));
        h_detint->Reset();

        TCanvas* c_detint_low_dt = new TCanvas("c_detint_low_dt", "Distribution of Low dt Ev61Det Integral Values", 1200, 700);
        c_detint_low_dt->SetLogy();
        c_detint_low_dt->cd();
        h_detint_low_dt->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_detint_low_dt->GetYaxis()->SetTitle("Counts");
        h_detint_low_dt->Draw();
        c_detint_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint_low_dt.png");
        // c_detint_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_detint_low_dt.png");
        c_detint_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_low_dt.png"));
        h_detint_low_dt->Reset();

        TCanvas* c_dt_v_detint = new TCanvas("c_dt_v_detint", "Event 61 Detector dt vs Integral Value", 1200, 700);
        // c_dt_v_detint->SetLogy();
        c_dt_v_detint->cd();
        h_dt_v_detint->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_detint->GetYaxis()->SetTitle("Integral (Ph.e.)");
        h_dt_v_detint->SetMarkerStyle(7);
        h_dt_v_detint->Draw();
        c_dt_v_detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_v_detint.png");
        // c_dt_v_detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_v_detint.png");
        c_dt_v_detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_detint.png"));
        h_dt_v_detint->Reset();

        TCanvas* c_int = new TCanvas("c_int", "Distribution of All Event Integral Values", 1200, 700);
        c_int->SetLogy();
        c_int->cd();
        h_int->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_int->GetYaxis()->SetTitle("Counts");
        h_int->Draw();
        c_int->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_int.png");
        // c_int->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_int.png");
        c_int->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int.png"));
        h_int->Reset();

        TCanvas* c_eventint = new TCanvas("c_eventint", "Distribution of All Int. Trig. Det. Integral Values", 1200, 700);
        c_eventint->SetLogy();
        c_eventint->cd();
        h_eventint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_eventint->GetYaxis()->SetTitle("Counts");
        h_eventint->Draw();
        c_eventint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_eventint.png");
        // c_eventint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_eventint.png");
        c_eventint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_eventint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_eventint.png"));
        h_eventint->Reset();

        TCanvas* c_eventvpint = new TCanvas("c_eventvpint", "Distribution of All Int. Trig. VP Integral Values", 1200, 700);
        c_eventvpint->SetLogy();
        c_eventvpint->cd();
        h_eventvpint->GetXaxis()->SetTitle("Integral (ADC)");
        h_eventvpint->GetYaxis()->SetTitle("Counts");
        h_eventvpint->Draw();
        c_eventvpint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_eventvpint.png");
        // c_eventvpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_eventvpint.png");
        c_eventvpint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_eventvpint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_eventvpint.png"));
        h_eventvpint->Reset();

        TCanvas* c_LLint = new TCanvas("c_LLint", "Distribution of All Low Light Integral Values", 1200, 700);
        c_LLint->SetLogy();
        c_LLint->cd();
        h_LLint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_LLint->GetYaxis()->SetTitle("Counts");
        h_LLint->Draw();
        c_LLint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_LLint.png");
        // c_LLint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_LLint.png");
        c_LLint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_LLint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_LLint.png"));
        h_LLint->Reset();

        TCanvas* c_HLint = new TCanvas("c_HLint", "Distribution of All High Light Integral Values", 1200, 700);
        c_HLint->SetLogy();
        c_HLint->cd();
        h_HLint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_HLint->GetYaxis()->SetTitle("Counts");
        h_HLint->Draw();
        c_HLint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_HLint.png");
        // c_HLint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_HLint.png");
        c_HLint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_HLint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_HLint.png"));
        h_HLint->Reset();

        TCanvas* c_ev61int = new TCanvas("c_ev61int", "Distribution of All Ev61 Det. Integral Values", 1200, 700);
        c_ev61int->SetLogy();
        c_ev61int->cd();
        h_ev61int->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_ev61int->GetYaxis()->SetTitle("Counts");
        h_ev61int->Draw();
        c_ev61int->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_ev61int.png");
        // c_ev61int->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_ev61int.png");
        c_ev61int->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_ev61int.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_ev61int.png"));
        h_ev61int->Reset();

        TCanvas* c_ev61vpint = new TCanvas("c_ev61vpint", "Distribution of All Ev61 VP Integral Values", 1200, 700);
        c_ev61vpint->SetLogy();
        c_ev61vpint->cd();
        h_ev61vpint->GetXaxis()->SetTitle("Integral (ADC)");
        h_ev61vpint->GetYaxis()->SetTitle("Counts");
        h_ev61vpint->Draw();
        c_ev61vpint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_ev61vpint.png");
        // c_ev61vpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_ev61vpint.png");
        c_ev61vpint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_ev61vpint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_ev61vpint.png"));
        h_ev61vpint->Reset();

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

        TCanvas* c_int_cent_spike = new TCanvas("c_int_cent_spike", "Distribution of All Central Spike Integral Values", 1200, 700);
        c_int_cent_spike->SetLogy();
        c_int_cent_spike->cd();
        h_int_cent_spike->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_int_cent_spike->GetYaxis()->SetTitle("Counts");
        h_int_cent_spike->Draw();
        c_int_cent_spike->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_int_cent_spike.png");
        // c_int_cent_spike->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_int_cent_spike.png");
        c_int_cent_spike->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_cent_spike.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_cent_spike.png"));
        h_int_cent_spike->Reset();

        TCanvas* c_int_neg_low_dt = new TCanvas("c_int_neg_low_dt", "Distribution of All Low dt Detector Integral Values", 1200, 700);
        c_int_neg_low_dt->SetLogy();
        c_int_neg_low_dt->cd();
        h_int_neg_low_dt->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_int_neg_low_dt->GetYaxis()->SetTitle("Counts");
        h_int_neg_low_dt->Draw("same");
        h_int_neg_low_dt->SetLineColor(kRed);
        h_int_pos_low_dt->Draw("same");
        h_int_pos_low_dt->SetLineColor(kGreen);
        TLegend *leg_int_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_int_neg_low_dt->AddEntry(h_int_pos_low_dt, "Positive Low dt Events", "l");
        leg_int_neg_low_dt->AddEntry(h_int_neg_low_dt, "Negative Low dt Events", "l");
        leg_int_neg_low_dt->Draw();
        c_int_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_int_low_dt.png");
        // c_int_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_int_low_dt.png");
        c_int_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_neg_low_dt.png"));
        h_int_pos_low_dt->Reset();
        h_int_neg_low_dt->Reset();

        TCanvas* c_svpint_neg_low_dt = new TCanvas("c_svpint_neg_low_dt", "Distribution of All Low dt Side Veto Panel Integral Values", 1200, 700);
        c_svpint_neg_low_dt->SetLogy();
        c_svpint_neg_low_dt->cd();
        h_svpint_cent_spike->GetXaxis()->SetTitle("Integral (ADC)");
        h_svpint_cent_spike->GetYaxis()->SetTitle("Counts");
        h_svpint_cent_spike->Draw("same");
        h_svpint_cent_spike->SetLineColor(kBlue);
        h_svpint_neg_low_dt->Draw("same");
        h_svpint_neg_low_dt->SetLineColor(kRed);
        h_svpint_pos_low_dt->Draw("same");
        h_svpint_pos_low_dt->SetLineColor(kGreen);
        TLegend *leg_svpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_svpint_neg_low_dt->AddEntry(h_svpint_cent_spike, "Central Spike Events", "l");
        leg_svpint_neg_low_dt->AddEntry(h_svpint_pos_low_dt, "Positive Low dt Events", "l");
        leg_svpint_neg_low_dt->AddEntry(h_svpint_neg_low_dt, "Negative Low dt Events", "l");
        leg_svpint_neg_low_dt->Draw();
        c_svpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_svpint_low_dt.png");
        // c_svpint_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_svpint_low_dt.png");
        c_svpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_svpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_svpint_neg_low_dt.png"));
        h_svpint_pos_low_dt->Reset();
        h_svpint_neg_low_dt->Reset();

        TCanvas* c_tvpint_neg_low_dt = new TCanvas("c_tvpint_neg_low_dt", "Distribution of All Low dt Top Veto Panel Integral Values", 1200, 700);
        c_tvpint_neg_low_dt->SetLogy();
        c_tvpint_neg_low_dt->cd();
        h_tvpint_cent_spike->GetXaxis()->SetTitle("Integral (ADC)");
        h_tvpint_cent_spike->GetYaxis()->SetTitle("Counts");
        h_tvpint_cent_spike->Draw("same");
        h_tvpint_cent_spike->SetLineColor(kBlue);
        h_tvpint_neg_low_dt->Draw("same");
        h_tvpint_neg_low_dt->SetLineColor(kRed);
        h_tvpint_pos_low_dt->Draw("same");
        h_tvpint_pos_low_dt->SetLineColor(kGreen);
        TLegend *leg_tvpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_tvpint_neg_low_dt->AddEntry(h_tvpint_cent_spike, "Central Spike Events", "l");
        leg_tvpint_neg_low_dt->AddEntry(h_tvpint_pos_low_dt, "Positive Low dt Events", "l");
        leg_tvpint_neg_low_dt->AddEntry(h_tvpint_neg_low_dt, "Negative Low dt Events", "l");
        leg_tvpint_neg_low_dt->Draw();
        c_tvpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_tvpint_low_dt.png");
        // c_tvpint_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_tvpint_low_dt.png");
        c_tvpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_tvpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_tvpint_neg_low_dt.png"));
        h_tvpint_pos_low_dt->Reset();
        h_tvpint_neg_low_dt->Reset();

        TCanvas* c_avpint_neg_low_dt = new TCanvas("c_avpint_neg_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
        c_avpint_neg_low_dt->SetLogy();
        c_avpint_neg_low_dt->cd();
        // TF1* func1 = new TF1("func1", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", 200, 800);
        // func1->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        // func1->SetParameters(1, 10, 0, 50);
        // h_avpint_pos_low_dt->Fit("func1", "R");
        h_avpint_cent_spike->GetXaxis()->SetTitle("Integral (ADC)");
        h_avpint_cent_spike->GetYaxis()->SetTitle("Counts");
        // h_avpint_cent_spike->SetMaximum(1200);
        h_avpint_cent_spike->Draw("same");
        h_avpint_cent_spike->SetLineColor(kBlue);
        h_avpint_neg_low_dt->Draw("same");
        h_avpint_neg_low_dt->SetLineColor(kRed);
        h_avpint_pos_low_dt->Draw("same");
        h_avpint_pos_low_dt->SetLineColor(kGreen);
        // func1->Draw("same");
        // func1->SetLineColor(1);
        TLegend *leg_avpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_avpint_neg_low_dt->AddEntry(h_avpint_cent_spike, "Central Spike Events", "l");
        leg_avpint_neg_low_dt->AddEntry(h_avpint_pos_low_dt, "Positive Low dt Events", "l");
        leg_avpint_neg_low_dt->AddEntry(h_avpint_neg_low_dt, "Negative Low dt Events", "l");
        leg_avpint_neg_low_dt->Draw();
        c_avpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_avpint_low_dt.png");
        // c_avpint_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_avpint_low_dt.png");
        c_avpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_neg_low_dt.png"));
        h_avpint_pos_low_dt->Reset();
        h_avpint_neg_low_dt->Reset();

        TCanvas* c_avpint_low_dt = new TCanvas("c_avpint_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
        // c_avpint_low_dt->SetLogy();
        c_avpint_low_dt->cd();
        TF1* func2 = new TF1("func2", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", 0, 800);
        func2->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func2->SetParameters(100, 10000, 300, 300);
        h_avpint_low_dt->Fit("func2", "R");
        h_avpint_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
        h_avpint_low_dt->GetYaxis()->SetTitle("Counts");
        // h_avpint_low_dt->GetYaxis()->SetRange(0, 50);
        h_avpint_low_dt->Draw();
        func2->Draw("same");
        func2->SetLineColor(2);
        c_avpint_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_avpint_low_dt_1c.png");
        // c_avpint_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_avpint_low_dt_1c.png");
        c_avpint_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_low_dt.png"));
        h_avpint_low_dt->Reset();

        TCanvas* c_dt_v_tmuon = new TCanvas("c_dt_v_tmuon", "Event 61 Detector dt vs Muon dt", 1200, 700);
        // c_dt_v_tmuon->SetLogy();
        c_dt_v_tmuon->cd();
        h_dt_v_tmuon->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_tmuon->GetYaxis()->SetTitle("Time Since Last Muon (ns)");
        h_dt_v_tmuon->SetMarkerStyle(7);
        h_dt_v_tmuon->Draw();
        c_dt_v_tmuon->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_v_tmuon.png");
        // c_dt_v_tmuon->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_dt_v_tmuon.png");
        c_dt_v_tmuon->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_tmuon.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_tmuon.png"));
        h_dt_v_tmuon->Reset();

        TCanvas* c_dt_v_tlead = new TCanvas("c_dt_v_tlead", "Event 61 Detector dt vs Lead dt", 1200, 700);
        // c_dt_v_tlead->SetLogy();
        c_dt_v_tlead->cd();
        h_dt_v_tlead->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_tlead->GetYaxis()->SetTitle("Time Since Last Lead Muon (ns)");
        h_dt_v_tlead->SetMarkerStyle(7);
        h_dt_v_tlead->Draw();
        c_dt_v_tlead->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_v_tlead.png");
        // c_dt_v_tlead->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_dt_v_tlead.png");
        c_dt_v_tlead->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_tlead.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_tlead.png"));
        h_dt_v_tlead->Reset();

        TCanvas* c_vec_size = new TCanvas("c_vec_size", "Distribution of pulses.size() Values", 1200, 700);
        c_vec_size->SetLogy();
        c_vec_size->cd();
        h_vec_size->GetXaxis()->SetTitle("Number of Events");
        h_vec_size->GetYaxis()->SetTitle("Counts");
        h_vec_size->Draw();
        c_vec_size->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_vec_size.png");
        // c_vec_size->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots7894/h_vec_size.png");
        c_vec_size->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vec_size.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vec_size.png"));
        h_vec_size->Reset();

    }

    if (run == 11608) {

        TCanvas* c_dt_61 = new TCanvas("c_dt_61", "Delta-T Between Event 61 and Detector Events", 900, 700);
        c_dt_61->SetLogy();
        c_dt_61->cd();
        h_dt_61->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61->GetYaxis()->SetTitle("Counts");
        h_dt_61->Draw();
        c_dt_61->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61.png");
        // c_dt_61->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_dt_61.png");
        c_dt_61->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61.png"));
        h_dt_61->Reset();

        TCanvas* c_dt_61_16 = new TCanvas("c_dt_61_16", "Delta-T Between Event 61 and SiPM Events", 900, 700);
        // c_dt_61_16->SetLogy();
        c_dt_61_16->cd();
        h_dt_61_16->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_16->GetYaxis()->SetTitle("Counts");
        h_dt_61_16->Draw();
        c_dt_61_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_16.png");
        // c_dt_61_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_dt_61_16.png");
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
        c_dt_61_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_17.png");
        // c_dt_61_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_dt_61_17.png");
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
        c_dt_61_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_18.png");
        // c_dt_61_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_dt_61_18.png");
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
        c_dt_61_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_19.png");
        // c_dt_61_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_dt_61_19.png");
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
        c_dt_61_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_20.png");
        // c_dt_61_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_dt_61_20.png");
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
        c_dt_61_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_21.png");
        // c_dt_61_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_dt_61_21.png");
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
        c_dt_61_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_22.png");
        // c_dt_61_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_dt_61_22.png");
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
        c_dt_61_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_23.png");
        // c_dt_61_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_dt_61_23.png");
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
        c_dt_61_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_top.png");
        // c_dt_61_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_dt_61_top.png");
        c_dt_61_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_top.png"));
        h_dt_61_top->Reset();
        
        TCanvas* c_vpint_16 = new TCanvas("c_vpint_16", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_16->SetLogy();
        c_vpint_16->cd();
        h_vpint_16->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_16->GetYaxis()->SetTitle("Counts");
        h_vpint_16->Draw();
        c_vpint_16->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_vpint_16.png");
        // c_vpint_16->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_vpint_16.png");
        c_vpint_16->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_16.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_16.png"));
        h_vpint_16->Reset();

        TCanvas* c_vpint_17 = new TCanvas("c_vpint_17", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_17->SetLogy();
        c_vpint_17->cd();
        h_vpint_17->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_17->GetYaxis()->SetTitle("Counts");
        h_vpint_17->Draw();
        c_vpint_17->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_vpint_17.png");
        // c_vpint_17->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_vpint_17.png");
        c_vpint_17->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_17.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_17.png"));
        h_vpint_17->Reset();

        TCanvas* c_vpint_18 = new TCanvas("c_vpint_18", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_18->SetLogy();
        c_vpint_18->cd();
        h_vpint_18->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_18->GetYaxis()->SetTitle("Counts");
        h_vpint_18->Draw();
        c_vpint_18->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_vpint_18.png");
        // c_vpint_18->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_vpint_18.png");
        c_vpint_18->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_18.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_18.png"));
        h_vpint_18->Reset();

        TCanvas* c_vpint_19 = new TCanvas("c_vpint_19", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_19->SetLogy();
        c_vpint_19->cd();
        h_vpint_19->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_19->GetYaxis()->SetTitle("Counts");
        h_vpint_19->Draw();
        c_vpint_19->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_vpint_19.png");
        // c_vpint_19->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_vpint_19.png");
        c_vpint_19->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_19.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_19.png"));
        h_vpint_19->Reset();

        TCanvas* c_vpint_20 = new TCanvas("c_vpint_20", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_20->SetLogy();
        c_vpint_20->cd();
        h_vpint_20->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_20->GetYaxis()->SetTitle("Counts");
        h_vpint_20->Draw();
        c_vpint_20->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_vpint_20.png");
        // c_vpint_20->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_vpint_20.png");
        c_vpint_20->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_20.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_20.png"));
        h_vpint_20->Reset();

        TCanvas* c_vpint_21 = new TCanvas("c_vpint_21", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_21->SetLogy();
        c_vpint_21->cd();
        h_vpint_21->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_21->GetYaxis()->SetTitle("Counts");
        h_vpint_21->Draw();
        c_vpint_21->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_vpint_21.png");
        // c_vpint_21->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_vpint_21.png");
        c_vpint_21->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_21.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_21.png"));
        h_vpint_21->Reset();

        TCanvas* c_vpint_22 = new TCanvas("c_vpint_22", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_22->SetLogy();
        c_vpint_22->cd();
        h_vpint_22->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_22->GetYaxis()->SetTitle("Counts");
        h_vpint_22->Draw();
        c_vpint_22->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_vpint_22.png");
        // c_vpint_22->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_vpint_22.png");
        c_vpint_22->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_22.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_22.png"));
        h_vpint_22->Reset();

        TCanvas* c_vpint_23 = new TCanvas("c_vpint_23", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_23->SetLogy();
        c_vpint_23->cd();
        h_vpint_23->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_23->GetYaxis()->SetTitle("Counts");
        h_vpint_23->Draw();
        c_vpint_23->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_vpint_23.png");
        // c_vpint_23->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_vpint_23.png");
        c_vpint_23->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_23.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_23.png"));
        h_vpint_23->Reset();

        TCanvas* c_vpint_top = new TCanvas("c_vpint_top", "Event 61 and SiPM dt Integral Values", 900, 700);
        // c_vpint_top->SetLogy();
        c_vpint_top->cd();
        h_vpint_top->GetXaxis()->SetTitle("Integral (ADC)");
        h_vpint_top->GetYaxis()->SetTitle("Counts");
        h_vpint_top->Draw();
        c_vpint_top->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_vpint_top.png");
        // c_vpint_top->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_vpint_top.png");
        c_vpint_top->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vpint_top.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vpint_top.png"));
        h_vpint_top->Reset();

        TCanvas* c_detint = new TCanvas("c_detint", "Distribution of All Ev61Det Integral Values", 1200, 700);
        c_detint->SetLogy();
        c_detint->cd();
        h_detint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_detint->GetYaxis()->SetTitle("Counts");
        h_detint->Draw();
        c_detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_detint.png");
        // c_detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_detint.png");
        c_detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint.png"));
        h_detint->Reset();

        TCanvas* c_detint_low_dt = new TCanvas("c_detint_low_dt", "Distribution of Low dt Ev61Det Integral Values", 1200, 700);
        c_detint_low_dt->SetLogy();
        c_detint_low_dt->cd();
        h_detint_low_dt->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_detint_low_dt->GetYaxis()->SetTitle("Counts");
        h_detint_low_dt->Draw();
        c_detint_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_detint_low_dt.png");
        // c_detint_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_detint_low_dt.png");
        c_detint_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_low_dt.png"));
        h_detint_low_dt->Reset();

        TCanvas* c_dt_v_detint = new TCanvas("c_dt_v_detint", "Event 61 Detector dt vs Integral Value", 1200, 700);
        // c_dt_v_detint->SetLogy();
        c_dt_v_detint->cd();
        h_dt_v_detint->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_detint->GetYaxis()->SetTitle("Integral (Ph.e.)");
        h_dt_v_detint->SetMarkerStyle(7);
        h_dt_v_detint->Draw();
        c_dt_v_detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_v_detint.png");
        // c_dt_v_detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_dt_v_detint.png");
        c_dt_v_detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_detint.png"));
        h_dt_v_detint->Reset();

        TCanvas* c_int = new TCanvas("c_int", "Distribution of All Event Integral Values", 1200, 700);
        c_int->SetLogy();
        c_int->cd();
        h_int->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_int->GetYaxis()->SetTitle("Counts");
        h_int->Draw();
        c_int->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_int.png");
        // c_int->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_int.png");
        c_int->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int.png"));
        h_int->Reset();

        TCanvas* c_eventint = new TCanvas("c_eventint", "Distribution of All Int. Trig. Det. Integral Values", 1200, 700);
        c_eventint->SetLogy();
        c_eventint->cd();
        h_eventint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_eventint->GetYaxis()->SetTitle("Counts");
        h_eventint->Draw();
        c_eventint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_eventint.png");
        // c_eventint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_eventint.png");
        c_eventint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_eventint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_eventint.png"));
        h_eventint->Reset();

        TCanvas* c_eventvpint = new TCanvas("c_eventvpint", "Distribution of All Int. Trig. VP Integral Values", 1200, 700);
        c_eventvpint->SetLogy();
        c_eventvpint->cd();
        h_eventvpint->GetXaxis()->SetTitle("Integral (ADC)");
        h_eventvpint->GetYaxis()->SetTitle("Counts");
        h_eventvpint->Draw();
        c_eventvpint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_eventvpint.png");
        // c_eventvpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_eventvpint.png");
        c_eventvpint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_eventvpint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_eventvpint.png"));
        h_eventvpint->Reset();

        TCanvas* c_LLint = new TCanvas("c_LLint", "Distribution of All Low Light Integral Values", 1200, 700);
        c_LLint->SetLogy();
        c_LLint->cd();
        h_LLint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_LLint->GetYaxis()->SetTitle("Counts");
        h_LLint->Draw();
        c_LLint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_LLint.png");
        // c_LLint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_LLint.png");
        c_LLint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_LLint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_LLint.png"));
        h_LLint->Reset();

        TCanvas* c_HLint = new TCanvas("c_HLint", "Distribution of All High Light Integral Values", 1200, 700);
        c_HLint->SetLogy();
        c_HLint->cd();
        h_HLint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_HLint->GetYaxis()->SetTitle("Counts");
        h_HLint->Draw();
        c_HLint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_HLint.png");
        // c_HLint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_HLint.png");
        c_HLint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_HLint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_HLint.png"));
        h_HLint->Reset();

        TCanvas* c_ev61int = new TCanvas("c_ev61int", "Distribution of All Ev61 Det. Integral Values", 1200, 700);
        c_ev61int->SetLogy();
        c_ev61int->cd();
        h_ev61int->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_ev61int->GetYaxis()->SetTitle("Counts");
        h_ev61int->Draw();
        c_ev61int->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_ev61int.png");
        // c_ev61int->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_ev61int.png");
        c_ev61int->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_ev61int.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_ev61int.png"));
        h_ev61int->Reset();

        TCanvas* c_ev61vpint = new TCanvas("c_ev61vpint", "Distribution of All Ev61 VP Integral Values", 1200, 700);
        c_ev61vpint->SetLogy();
        c_ev61vpint->cd();
        h_ev61vpint->GetXaxis()->SetTitle("Integral (ADC)");
        h_ev61vpint->GetYaxis()->SetTitle("Counts");
        h_ev61vpint->Draw();
        c_ev61vpint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_ev61vpint.png");
        // c_ev61vpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_ev61vpint.png");
        c_ev61vpint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_ev61vpint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_ev61vpint.png"));
        h_ev61vpint->Reset();

        TCanvas* c_34detint = new TCanvas("c_34detint", "Distribution of All tB34 Det. Integral Values", 1200, 700);
        c_34detint->SetLogy();
        c_34detint->cd();
        h_34detint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_34detint->GetYaxis()->SetTitle("Counts");
        h_34detint->Draw();
        c_34detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_34detint.png");
        // c_34detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_34detint.png");
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
        c_34vpint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_34vpint.png");
        // c_34vpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_34vpint.png");
        c_34vpint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_34vpint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_34vpint.png"));
        h_34vpint->Reset();

        TCanvas* c_int_cent_spike = new TCanvas("c_int_cent_spike", "Distribution of All Central Spike Integral Values", 1200, 700);
        c_int_cent_spike->SetLogy();
        c_int_cent_spike->cd();
        h_int_cent_spike->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_int_cent_spike->GetYaxis()->SetTitle("Counts");
        h_int_cent_spike->Draw();
        c_int_cent_spike->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_int_cent_spike.png");
        // c_int_cent_spike->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_int_cent_spike.png");
        c_int_cent_spike->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_cent_spike.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_cent_spike.png"));
        h_int_cent_spike->Reset();

        TCanvas* c_int_neg_low_dt = new TCanvas("c_int_neg_low_dt", "Distribution of All Low dt Detector Integral Values", 1200, 700);
        c_int_neg_low_dt->SetLogy();
        c_int_neg_low_dt->cd();
        h_int_neg_low_dt->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_int_neg_low_dt->GetYaxis()->SetTitle("Counts");
        h_int_neg_low_dt->Draw("same");
        h_int_neg_low_dt->SetLineColor(kRed);
        h_int_pos_low_dt->Draw("same");
        h_int_pos_low_dt->SetLineColor(kGreen);
        TLegend *leg_int_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_int_neg_low_dt->AddEntry(h_int_pos_low_dt, "Positive Low dt Events", "l");
        leg_int_neg_low_dt->AddEntry(h_int_neg_low_dt, "Negative Low dt Events", "l");
        leg_int_neg_low_dt->Draw();
        c_int_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_int_low_dt.png");
        // c_int_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_int_low_dt.png");
        c_int_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_neg_low_dt.png"));
        h_int_pos_low_dt->Reset();
        h_int_neg_low_dt->Reset();

        TCanvas* c_svpint_neg_low_dt = new TCanvas("c_svpint_neg_low_dt", "Distribution of All Low dt Side Veto Panel Integral Values", 1200, 700);
        c_svpint_neg_low_dt->SetLogy();
        c_svpint_neg_low_dt->cd();
        h_svpint_cent_spike->GetXaxis()->SetTitle("Integral (ADC)");
        h_svpint_cent_spike->GetYaxis()->SetTitle("Counts");
        h_svpint_cent_spike->Draw("same");
        h_svpint_cent_spike->SetLineColor(kBlue);
        h_svpint_neg_low_dt->Draw("same");
        h_svpint_neg_low_dt->SetLineColor(kRed);
        h_svpint_pos_low_dt->Draw("same");
        h_svpint_pos_low_dt->SetLineColor(kGreen);
        TLegend *leg_svpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_svpint_neg_low_dt->AddEntry(h_svpint_cent_spike, "Central Spike Events", "l");
        leg_svpint_neg_low_dt->AddEntry(h_svpint_pos_low_dt, "Positive Low dt Events", "l");
        leg_svpint_neg_low_dt->AddEntry(h_svpint_neg_low_dt, "Negative Low dt Events", "l");
        leg_svpint_neg_low_dt->Draw();
        c_svpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_svpint_low_dt.png");
        // c_svpint_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_svpint_low_dt.png");
        c_svpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_svpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_svpint_neg_low_dt.png"));
        h_svpint_pos_low_dt->Reset();
        h_svpint_neg_low_dt->Reset();

        TCanvas* c_tvpint_neg_low_dt = new TCanvas("c_tvpint_neg_low_dt", "Distribution of All Low dt Top Veto Panel Integral Values", 1200, 700);
        c_tvpint_neg_low_dt->SetLogy();
        c_tvpint_neg_low_dt->cd();
        h_tvpint_cent_spike->GetXaxis()->SetTitle("Integral (ADC)");
        h_tvpint_cent_spike->GetYaxis()->SetTitle("Counts");
        h_tvpint_cent_spike->Draw("same");
        h_tvpint_cent_spike->SetLineColor(kBlue);
        h_tvpint_neg_low_dt->Draw("same");
        h_tvpint_neg_low_dt->SetLineColor(kRed);
        h_tvpint_pos_low_dt->Draw("same");
        h_tvpint_pos_low_dt->SetLineColor(kGreen);
        TLegend *leg_tvpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_tvpint_neg_low_dt->AddEntry(h_tvpint_cent_spike, "Central Spike Events", "l");
        leg_tvpint_neg_low_dt->AddEntry(h_tvpint_pos_low_dt, "Positive Low dt Events", "l");
        leg_tvpint_neg_low_dt->AddEntry(h_tvpint_neg_low_dt, "Negative Low dt Events", "l");
        leg_tvpint_neg_low_dt->Draw();
        c_tvpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_tvpint_low_dt.png");
        // c_tvpint_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_tvpint_low_dt.png");
        c_tvpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_tvpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_tvpint_neg_low_dt.png"));
        h_tvpint_pos_low_dt->Reset();
        h_tvpint_neg_low_dt->Reset();

        TCanvas* c_avpint_neg_low_dt = new TCanvas("c_avpint_neg_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
        c_avpint_neg_low_dt->SetLogy();
        c_avpint_neg_low_dt->cd();
        // TF1* func1 = new TF1("func1", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", 200, 800);
        // func1->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        // func1->SetParameters(1, 10, 0, 50);
        // h_avpint_pos_low_dt->Fit("func1", "R");
        h_avpint_cent_spike->GetXaxis()->SetTitle("Integral (ADC)");
        h_avpint_cent_spike->GetYaxis()->SetTitle("Counts");
        // h_avpint_cent_spike->SetMaximum(1200);
        h_avpint_cent_spike->Draw("same");
        h_avpint_cent_spike->SetLineColor(kBlue);
        h_avpint_neg_low_dt->Draw("same");
        h_avpint_neg_low_dt->SetLineColor(kRed);
        h_avpint_pos_low_dt->Draw("same");
        h_avpint_pos_low_dt->SetLineColor(kGreen);
        // func1->Draw("same");
        // func1->SetLineColor(1);
        TLegend *leg_avpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_avpint_neg_low_dt->AddEntry(h_avpint_cent_spike, "Central Spike Events", "l");
        leg_avpint_neg_low_dt->AddEntry(h_avpint_pos_low_dt, "Positive Low dt Events", "l");
        leg_avpint_neg_low_dt->AddEntry(h_avpint_neg_low_dt, "Negative Low dt Events", "l");
        leg_avpint_neg_low_dt->Draw();
        c_avpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_avpint_low_dt.png");
        // c_avpint_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_avpint_low_dt.png");
        c_avpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_neg_low_dt.png"));
        h_avpint_pos_low_dt->Reset();
        h_avpint_neg_low_dt->Reset();

        TCanvas* c_avpint_low_dt = new TCanvas("c_avpint_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
        // c_avpint_low_dt->SetLogy();
        c_avpint_low_dt->cd();
        TF1* func2 = new TF1("func2", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -100, 350);
        func2->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func2->SetParameters(1, 10, 0, 50);
        h_avpint_low_dt->Fit("func2", "R");
        h_avpint_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
        h_avpint_low_dt->GetYaxis()->SetTitle("Counts");
        // h_avpint_low_dt->GetYaxis()->SetRange(0, 50);
        h_avpint_low_dt->Draw();
        func2->Draw("same");
        func2->SetLineColor(2);
        c_avpint_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_avpint_low_dt_1c.png");
        // c_avpint_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_avpint_low_dt_1c.png");
        c_avpint_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_low_dt.png"));
        h_avpint_low_dt->Reset();

        TCanvas* c_dt_v_tmuon = new TCanvas("c_dt_v_tmuon", "Event 61 Detector dt vs Muon dt", 1200, 700);
        // c_dt_v_tmuon->SetLogy();
        c_dt_v_tmuon->cd();
        h_dt_v_tmuon->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_tmuon->GetYaxis()->SetTitle("Time Since Last Muon (ns)");
        h_dt_v_tmuon->SetMarkerStyle(7);
        h_dt_v_tmuon->Draw();
        c_dt_v_tmuon->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_v_tmuon.png");
        // c_dt_v_tmuon->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_dt_v_tmuon.png");
        c_dt_v_tmuon->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_tmuon.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_tmuon.png"));
        h_dt_v_tmuon->Reset();

        TCanvas* c_vec_size = new TCanvas("c_vec_size", "Distribution of pulses.size() Values", 1200, 700);
        c_vec_size->SetLogy();
        c_vec_size->cd();
        h_vec_size->GetXaxis()->SetTitle("Number of Events");
        h_vec_size->GetYaxis()->SetTitle("Counts");
        h_vec_size->Draw();
        c_vec_size->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_vec_size.png");
        // c_vec_size->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots11608/h_vec_size.png");
        c_vec_size->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vec_size.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vec_size.png"));
        h_vec_size->Reset();

    }

    if (run == 12636) {

        TCanvas* c_dt_61 = new TCanvas("c_dt_61", "Delta-T Between Event 61 and Detector Events", 900, 700);
        c_dt_61->SetLogy();
        c_dt_61->cd();
        h_dt_61->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61->GetYaxis()->SetTitle("Counts");
        h_dt_61->Draw();
        c_dt_61->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61.png");
        // c_dt_61->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_dt_61.png");
        c_dt_61->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61.png"));
        h_dt_61->Reset();

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
        // c_vpint_16->SetLogy();
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
        // c_vpint_17->SetLogy();
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
        // c_vpint_18->SetLogy();
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
        // c_vpint_19->SetLogy();
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
        // c_vpint_20->SetLogy();
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
        // c_vpint_21->SetLogy();
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
        // c_vpint_22->SetLogy();
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
        // c_vpint_23->SetLogy();
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
        // c_vpint_top->SetLogy();
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

        TCanvas* c_detint = new TCanvas("c_detint", "Distribution of All Ev61Det Integral Values", 1200, 700);
        c_detint->SetLogy();
        c_detint->cd();
        h_detint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_detint->GetYaxis()->SetTitle("Counts");
        h_detint->Draw();
        c_detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint.png");
        // c_detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_detint.png");
        c_detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint.png"));
        h_detint->Reset();

        TCanvas* c_detint_low_dt = new TCanvas("c_detint_low_dt", "Distribution of Low dt Ev61Det Integral Values", 1200, 700);
        c_detint_low_dt->SetLogy();
        c_detint_low_dt->cd();
        h_detint_low_dt->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_detint_low_dt->GetYaxis()->SetTitle("Counts");
        h_detint_low_dt->Draw();
        c_detint_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint_low_dt.png");
        // c_detint_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_detint_low_dt.png");
        c_detint_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_low_dt.png"));
        h_detint_low_dt->Reset();

        TCanvas* c_dt_v_detint = new TCanvas("c_dt_v_detint", "Event 61 Detector dt vs Integral Value", 1200, 700);
        // c_dt_v_detint->SetLogy();
        c_dt_v_detint->cd();
        h_dt_v_detint->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_detint->GetYaxis()->SetTitle("Integral (Ph.e.)");
        h_dt_v_detint->SetMarkerStyle(7);
        h_dt_v_detint->Draw();
        c_dt_v_detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_v_detint.png");
        // c_dt_v_detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_dt_v_detint.png");
        c_dt_v_detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_detint.png"));
        h_dt_v_detint->Reset();

        TCanvas* c_int = new TCanvas("c_int", "Distribution of All Event Integral Values", 1200, 700);
        c_int->SetLogy();
        c_int->cd();
        h_int->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_int->GetYaxis()->SetTitle("Counts");
        h_int->Draw();
        c_int->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_int.png");
        // c_int->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_int.png");
        c_int->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int.png"));
        h_int->Reset();

        TCanvas* c_eventint = new TCanvas("c_eventint", "Distribution of All Int. Trig. Det. Integral Values", 1200, 700);
        c_eventint->SetLogy();
        c_eventint->cd();
        h_eventint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_eventint->GetYaxis()->SetTitle("Counts");
        h_eventint->Draw();
        c_eventint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_eventint.png");
        // c_eventint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_eventint.png");
        c_eventint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_eventint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_eventint.png"));
        h_eventint->Reset();

        TCanvas* c_eventvpint = new TCanvas("c_eventvpint", "Distribution of All Int. Trig. VP Integral Values", 1200, 700);
        c_eventvpint->SetLogy();
        c_eventvpint->cd();
        h_eventvpint->GetXaxis()->SetTitle("Integral (ADC)");
        h_eventvpint->GetYaxis()->SetTitle("Counts");
        h_eventvpint->Draw();
        c_eventvpint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_eventvpint.png");
        // c_eventvpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_eventvpint.png");
        c_eventvpint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_eventvpint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_eventvpint.png"));
        h_eventvpint->Reset();

        TCanvas* c_LLint = new TCanvas("c_LLint", "Distribution of All Low Light Integral Values", 1200, 700);
        c_LLint->SetLogy();
        c_LLint->cd();
        h_LLint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_LLint->GetYaxis()->SetTitle("Counts");
        h_LLint->Draw();
        c_LLint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_LLint.png");
        // c_LLint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_LLint.png");
        c_LLint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_LLint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_LLint.png"));
        h_LLint->Reset();

        TCanvas* c_HLint = new TCanvas("c_HLint", "Distribution of All High Light Integral Values", 1200, 700);
        c_HLint->SetLogy();
        c_HLint->cd();
        h_HLint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_HLint->GetYaxis()->SetTitle("Counts");
        h_HLint->Draw();
        c_HLint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_HLint.png");
        // c_HLint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_HLint.png");
        c_HLint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_HLint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_HLint.png"));
        h_HLint->Reset();

        TCanvas* c_ev61int = new TCanvas("c_ev61int", "Distribution of All Ev61 Det. Integral Values", 1200, 700);
        c_ev61int->SetLogy();
        c_ev61int->cd();
        h_ev61int->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_ev61int->GetYaxis()->SetTitle("Counts");
        h_ev61int->Draw();
        c_ev61int->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_ev61int.png");
        // c_ev61int->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_ev61int.png");
        c_ev61int->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_ev61int.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_ev61int.png"));
        h_ev61int->Reset();

        TCanvas* c_ev61vpint = new TCanvas("c_ev61vpint", "Distribution of All Ev61 VP Integral Values", 1200, 700);
        c_ev61vpint->SetLogy();
        c_ev61vpint->cd();
        h_ev61vpint->GetXaxis()->SetTitle("Integral (ADC)");
        h_ev61vpint->GetYaxis()->SetTitle("Counts");
        h_ev61vpint->Draw();
        c_ev61vpint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_ev61vpint.png");
        // c_ev61vpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_ev61vpint.png");
        c_ev61vpint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_ev61vpint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_ev61vpint.png"));
        h_ev61vpint->Reset();

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

        TCanvas* c_int_cent_spike = new TCanvas("c_int_cent_spike", "Distribution of All Central Spike Integral Values", 1200, 700);
        c_int_cent_spike->SetLogy();
        c_int_cent_spike->cd();
        h_int_cent_spike->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_int_cent_spike->GetYaxis()->SetTitle("Counts");
        h_int_cent_spike->Draw();
        c_int_cent_spike->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_int_cent_spike.png");
        // c_int_cent_spike->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_int_cent_spike.png");
        c_int_cent_spike->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_cent_spike.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_cent_spike.png"));
        h_int_cent_spike->Reset();

        TCanvas* c_int_neg_low_dt = new TCanvas("c_int_neg_low_dt", "Distribution of All Low dt Detector Integral Values", 1200, 700);
        c_int_neg_low_dt->SetLogy();
        c_int_neg_low_dt->cd();
        h_int_neg_low_dt->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_int_neg_low_dt->GetYaxis()->SetTitle("Counts");
        h_int_neg_low_dt->Draw("same");
        h_int_neg_low_dt->SetLineColor(kRed);
        h_int_pos_low_dt->Draw("same");
        h_int_pos_low_dt->SetLineColor(kGreen);
        TLegend *leg_int_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_int_neg_low_dt->AddEntry(h_int_pos_low_dt, "Positive Low dt Events", "l");
        leg_int_neg_low_dt->AddEntry(h_int_neg_low_dt, "Negative Low dt Events", "l");
        leg_int_neg_low_dt->Draw();
        c_int_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_int_low_dt.png");
        // c_int_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_int_low_dt.png");
        c_int_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_neg_low_dt.png"));
        h_int_pos_low_dt->Reset();
        h_int_neg_low_dt->Reset();

        TCanvas* c_svpint_neg_low_dt = new TCanvas("c_svpint_neg_low_dt", "Distribution of All Low dt Side Veto Panel Integral Values", 1200, 700);
        c_svpint_neg_low_dt->SetLogy();
        c_svpint_neg_low_dt->cd();
        h_svpint_cent_spike->GetXaxis()->SetTitle("Integral (ADC)");
        h_svpint_cent_spike->GetYaxis()->SetTitle("Counts");
        h_svpint_cent_spike->Draw("same");
        h_svpint_cent_spike->SetLineColor(kBlue);
        h_svpint_neg_low_dt->Draw("same");
        h_svpint_neg_low_dt->SetLineColor(kRed);
        h_svpint_pos_low_dt->Draw("same");
        h_svpint_pos_low_dt->SetLineColor(kGreen);
        TLegend *leg_svpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_svpint_neg_low_dt->AddEntry(h_svpint_cent_spike, "Central Spike Events", "l");
        leg_svpint_neg_low_dt->AddEntry(h_svpint_pos_low_dt, "Positive Low dt Events", "l");
        leg_svpint_neg_low_dt->AddEntry(h_svpint_neg_low_dt, "Negative Low dt Events", "l");
        leg_svpint_neg_low_dt->Draw();
        c_svpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_svpint_low_dt.png");
        // c_svpint_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_svpint_low_dt.png");
        c_svpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_svpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_svpint_neg_low_dt.png"));
        h_svpint_pos_low_dt->Reset();
        h_svpint_neg_low_dt->Reset();

        TCanvas* c_tvpint_neg_low_dt = new TCanvas("c_tvpint_neg_low_dt", "Distribution of All Low dt Top Veto Panel Integral Values", 1200, 700);
        c_tvpint_neg_low_dt->SetLogy();
        c_tvpint_neg_low_dt->cd();
        h_tvpint_cent_spike->GetXaxis()->SetTitle("Integral (ADC)");
        h_tvpint_cent_spike->GetYaxis()->SetTitle("Counts");
        h_tvpint_cent_spike->Draw("same");
        h_tvpint_cent_spike->SetLineColor(kBlue);
        h_tvpint_neg_low_dt->Draw("same");
        h_tvpint_neg_low_dt->SetLineColor(kRed);
        h_tvpint_pos_low_dt->Draw("same");
        h_tvpint_pos_low_dt->SetLineColor(kGreen);
        TLegend *leg_tvpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_tvpint_neg_low_dt->AddEntry(h_tvpint_cent_spike, "Central Spike Events", "l");
        leg_tvpint_neg_low_dt->AddEntry(h_tvpint_pos_low_dt, "Positive Low dt Events", "l");
        leg_tvpint_neg_low_dt->AddEntry(h_tvpint_neg_low_dt, "Negative Low dt Events", "l");
        leg_tvpint_neg_low_dt->Draw();
        c_tvpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_tvpint_low_dt.png");
        // c_tvpint_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_tvpint_low_dt.png");
        c_tvpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_tvpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_tvpint_neg_low_dt.png"));
        h_tvpint_pos_low_dt->Reset();
        h_tvpint_neg_low_dt->Reset();

        TCanvas* c_avpint_neg_low_dt = new TCanvas("c_avpint_neg_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
        c_avpint_neg_low_dt->SetLogy();
        c_avpint_neg_low_dt->cd();
        // TF1* func1 = new TF1("func1", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", 200, 800);
        // func1->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        // func1->SetParameters(1, 10, 0, 50);
        // h_avpint_pos_low_dt->Fit("func1", "R");
        h_avpint_cent_spike->GetXaxis()->SetTitle("Integral (ADC)");
        h_avpint_cent_spike->GetYaxis()->SetTitle("Counts");
        // h_avpint_cent_spike->SetMaximum(1200);
        h_avpint_cent_spike->Draw("same");
        h_avpint_cent_spike->SetLineColor(kBlue);
        h_avpint_neg_low_dt->Draw("same");
        h_avpint_neg_low_dt->SetLineColor(kRed);
        h_avpint_pos_low_dt->Draw("same");
        h_avpint_pos_low_dt->SetLineColor(kGreen);
        // func1->Draw("same");
        // func1->SetLineColor(1);
        TLegend *leg_avpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_avpint_neg_low_dt->AddEntry(h_avpint_cent_spike, "Central Spike Events", "l");
        leg_avpint_neg_low_dt->AddEntry(h_avpint_pos_low_dt, "Positive Low dt Events", "l");
        leg_avpint_neg_low_dt->AddEntry(h_avpint_neg_low_dt, "Negative Low dt Events", "l");
        leg_avpint_neg_low_dt->Draw();
        c_avpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_avpint_low_dt.png");
        // c_avpint_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_avpint_low_dt.png");
        c_avpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_neg_low_dt.png"));
        h_avpint_pos_low_dt->Reset();
        h_avpint_neg_low_dt->Reset();

        TCanvas* c_avpint_low_dt = new TCanvas("c_avpint_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
        // c_avpint_low_dt->SetLogy();
        c_avpint_low_dt->cd();
        TF1* func2 = new TF1("func2", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -100, 350);
        func2->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func2->SetParameters(1, 10, 0, 50);
        h_avpint_low_dt->Fit("func2", "R");
        h_avpint_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
        h_avpint_low_dt->GetYaxis()->SetTitle("Counts");
        // h_avpint_low_dt->GetYaxis()->SetRange(0, 50);
        h_avpint_low_dt->Draw();
        func2->Draw("same");
        func2->SetLineColor(2);
        c_avpint_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_avpint_low_dt_1c.png");
        // c_avpint_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_avpint_low_dt_1c.png");
        c_avpint_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_low_dt.png"));
        h_avpint_low_dt->Reset();

        TCanvas* c_dt_v_tmuon = new TCanvas("c_dt_v_tmuon", "Event 61 Detector dt vs Muon dt", 1200, 700);
        // c_dt_v_tmuon->SetLogy();
        c_dt_v_tmuon->cd();
        h_dt_v_tmuon->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_tmuon->GetYaxis()->SetTitle("Time Since Last Muon (ns)");
        h_dt_v_tmuon->SetMarkerStyle(7);
        h_dt_v_tmuon->Draw();
        c_dt_v_tmuon->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_v_tmuon.png");
        // c_dt_v_tmuon->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_dt_v_tmuon.png");
        c_dt_v_tmuon->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_tmuon.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_tmuon.png"));
        h_dt_v_tmuon->Reset();

        TCanvas* c_vec_size = new TCanvas("c_vec_size", "Distribution of pulses.size() Values", 1200, 700);
        c_vec_size->SetLogy();
        c_vec_size->cd();
        h_vec_size->GetXaxis()->SetTitle("Number of Events");
        h_vec_size->GetYaxis()->SetTitle("Counts");
        h_vec_size->Draw();
        c_vec_size->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_vec_size.png");
        // c_vec_size->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots12636/h_vec_size.png");
        c_vec_size->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vec_size.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vec_size.png"));
        h_vec_size->Reset();

    }

    cout << "\n" << "End of code." << "\n" << endl;

    return 0;

}