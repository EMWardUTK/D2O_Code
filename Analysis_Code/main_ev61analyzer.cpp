#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPad.h>
#include <TTreeIndex.h>
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
#include "TMarker.h"
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
  bool issue;               /* Flag to keep track of unusual eventTree entries */
  long long ev61_start_time;/* Start time of Event 61 pulse */
  bool ev61_significant;    /* Flag to keep track of whether Event 61 event is also potentially significant */
  int selection_cut;        /* Tacks which selection cut(s) a given event satisfies */
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
    cout << "\n" << "Run: " << run << " Number of selection cuts: " << last_run << endl;

    /* Initialize vectors and variables */

    std::vector<struct pulse> pulses;

    int selection_cut_level = last_run;

    long long sample_range = 2 * pow(10, 4);

    double plot_min = 1000;          // -19960;

    double plot_max = 9960;          // -11000;

    double num_bins = 56;

    double sample_bins = 250;

    double ev61_offset = 2000;

    bool bt_HL_n_LL = false;

    double last_file_num = 0.0;

    double cumulative_MWhrs = 0.0;

    int num_events_low_dt = 0;

    int num_events_low_dt_cent_spike = 0;

    int num_events_low_dt_no_spike_neg = 0;

    int num_events_low_dt_no_spike_pos = 0;

    int take_sample_count = 0;

    int hist_events = 0;

    int run_events = 0;

    double tot_detint = 0;

    bool fall_is_true_spring_is_false = true;

    if (run < 66) {fall_is_true_spring_is_false = true;}

    else if (run >= 66) {fall_is_true_spring_is_false = false;}

    if (fall_is_true_spring_is_false) {

        sample_range = 4 * pow(10, 4);

        sample_bins = 500;

        plot_min = 1000;

        plot_max = 9960;

        ev61_offset = 2000;

    }

    else if (!fall_is_true_spring_is_false) {

        sample_range = 4 * pow(10, 4);

        sample_bins = 500;

        plot_min = -19960;

        plot_max = -11000;

        ev61_offset = -19000;

    }

    // sample_range = pow(10, 5);

    // sample_bins = 1250;

    // cout << "\n" << "Number of selection cuts: " << selection_cut_level << "\n" << endl;

    /* Initialize histograms */                 // NUMBER OF BINS SHOULD BE A MULTIPLE OF THE SMALLEST AXIS UNIT!!! (OR LENGTH OF AXIS DIVIDED BY SMALLEST AXIS UNIT?)

    TH1D* h_dt_61 = new TH1D("h_dt_61", "Delta-T Between Event 61 and Detector Peaks", num_bins, 0, 10000);                                                                 // 1 bin = 160 ns

    TH1D* h_dt_61_1 = new TH1D("h_dt_61_1", "Delta-T Between Event 61 and Detector Peaks", 625, - 1 * pow(10, 4), 1 * pow(10, 4));                                          // 1 bin = 32 ns

    TH1D* h_dt_61_2 = new TH1D("h_dt_61_2", "Delta-T Between Event 61 and Detector Peaks", 65, 0, 1040);                                         // 1 bin = 32 ns

    TH1D* h_dt_61_3 = new TH1D("h_dt_61_3", "Delta-T Between Event 61 and Detector Peaks", 3125, - 5 * pow(10, 4), 5 * pow(10, 4));                                         // 1 bin = 32 ns

    TH1D* h_dt_61_4 = new TH1D("h_dt_61_4", "Delta-T Between Event 61 and Detector Peaks", 1250, - 10 * pow(10, 4), 10 * pow(10, 4));                                       // 1 bin = 160 ns

    TH2D* h_dt_v_detint = new TH2D("h_dt_v_detint", "Delta-T vs Detector Integral Value", sample_bins, - sample_range, sample_range, 50, 0, 500);                           // 1 bin = 160 ns

    TH2D* h_dt_v_tmuon = new TH2D("h_dt_v_tmuon", "Delta-T vs Most Recent Muon Time Difference", 1250, - 2 * sample_range, 2 * sample_range, 1000, 0, 2 * pow(10, 7));      // 1 bin = 32 ns

    TH1D* h_int = new TH1D("h_int", "Distribution of All Event Integral Values Between High Light and Low Light LED Events", 200, 0, 2000);

    TH1D* h_LLint = new TH1D("h_LLint", "Distribution of All Low Light LED Integral Values", 200, 0, 2000);

    TH1D* h_HLint = new TH1D("h_HLint", "Distribution of All High Light LED Integral Values", 200, 0, 2000);

    TH1D* h_eventint = new TH1D("h_eventint", "Distribution of All Internally Triggered Detector Integral Values", 200, 0, 2000);

    TH1D* h_eventvpint = new TH1D("h_eventvpint", "Distribution of All Internally Triggered Veto Panel Integral Values", 200, 0, 2000);

    TH1D* h_ev61int = new TH1D("h_ev61int", "Distribution of All Event 61 Detector Integral Values", 200, 0, 2000);

    TH1D* h_ev61vpint = new TH1D("h_ev61vpint", "Distribution of All Event 61 Veto Panel Integral Values", 200, 0, 2000);

    TH1D* h_34detint = new TH1D("h_34detint", "Distribution of All triggerBits = 34 Detector Integral Values", 200, 0, 2000);

    TH1D* h_34vpint = new TH1D("h_34vpint", "Distribution of All triggerBits = 34 Veto Panel Integral Values", 200, 0, 2000);

    TH1D* h_runnum = new TH1D("h_runnum", "Beam Spill Run Number Distribution (0 - 10 us)", 100, 14000, 24000);

    TH2D* h_runnum_v_dt = new TH2D("h_runnum_v_dt", "Beam Spill Run Number vs Delta-T Between Event 61 and Detector Peaks", 100, 14000, 24000, sample_bins, - sample_range, sample_range);

    TH2D* h_runnum_v_detint = new TH2D("h_runnum_v_detint", "Beam Spill Run Number vs Detector Integral Value", 100, 14000, 24000, 50, 0, 500);

    TH1D* h_detint = new TH1D("h_detint", "Beam Spill Detector Integral Values", 50, 0, 500);

    TH1D* h_detint_1 = new TH1D("h_detint_1", "Beam Spill PMT Integral Values, All Selection Cuts", 80, 0, 800);

    TH1D* h_detint_2 = new TH1D("h_detint_2", "Beam Spill PMT Integral Values, All Selection Cuts", 80, 0, 800);

    TH1D* h_detint_3 = new TH1D("h_detint_3", "Beam Spill PMT Integral Values, All Selection Cuts", 80, 0, 800);

    TH1D* h_detint_4 = new TH1D("h_detint_4", "Beam Spill PMT Integral Values, All Selection Cuts", 80, 0, 800);

    TH1D* h_detint_5 = new TH1D("h_detint_5", "Beam Spill PMT Integral Values, All Selection Cuts", 80, 0, 800);

    TH1D* h_detint_6 = new TH1D("h_detint_6", "Beam Spill PMT Integral Values, All Selection Cuts", 80, 0, 800);

    TH1D* h_detint_R1 = new TH1D("h_detint_R1", "Far Negative Beam Spill PMT Integral Values (-20us to -5us)", 50, 0, 500);

    TH1D* h_detint_R2 = new TH1D("h_detint_R2", "Close Negative Beam Spill PMT Integral Values (-0.5us to 0us)", 50, 0, 500);

    TH1D* h_detint_R3 = new TH1D("h_detint_R3", "BRN Spike Beam Spill PMT Integral Values (0.4us to 0.9us)", 50, 0, 500);

    TH1D* h_detint_R4 = new TH1D("h_detint_R4", "Close Positive Beam Spill PMT Integral Values (1.5us to 4.5us)", 50, 0, 500);

    TH1D* h_detint_R5 = new TH1D("h_detint_R5", "Far Positive Beam Spill PMT Integral Values (5us to 20us)", 50, 0, 500);

    TH1D* h_detint_R6 = new TH1D("h_detint_R6", "1us Before S25 Ev61 PMT Integral Values (18us to 19us)", 50, 0, 500);

    TH1D* h_detint_R7 = new TH1D("h_detint_R7", "PMT Integral Values 0.5us Before and After S25 Ev61 Waveform", 50, 0, 500);            // 1us After S25 Ev61 PMT Integral Values (22us to 23us)

    TH1D* h_detint_low_dt = new TH1D("h_detint_low_dt", "Distribution of Low dt Detector Integral Values", 80, 0, 800);

    TH1D* h_int_cent_spike = new TH1D("h_int_cent_spike", "Distribution of All Central Spike Detector Integral Values", 200, 0, 2000);

    TH1D* h_int_pos_low_dt = new TH1D("h_int_pos_low_dt", "Distribution of Low dt Detector Integral Values", 200, 0, 2000);

    TH1D* h_int_neg_low_dt = new TH1D("h_int_neg_low_dt", "Distribution of Low dt Detector Integral Values", 200, 0, 2000);

    TH1D* h_avpint_low_dt = new TH1D("h_avpint_low_dt", "Distribution of Low dt Veto Panel Integral Values", 60, -400, 2000);

    TH1D* h_avpint_cent_spike = new TH1D("h_avpint_cent_spike", "Distribution of All Central Spike Veto Panel Integral Values", 60, -400, 2000);

    TH1D* h_avpint_pos_low_dt = new TH1D("h_avpint_pos_low_dt", "Distribution of Low dt Veto Panel Integral Values", 60, -400, 2000);

    TH1D* h_avpint_neg_low_dt = new TH1D("h_avpint_neg_low_dt", "Distribution of Low dt Veto Panel Integral Values", 60, -400, 2000);

    TH1D* h_svpint_cent_spike = new TH1D("h_svpint_cent_spike", "Distribution of All Central Spike Side Veto Panel Integral Values", 240, -400, 2000);

    TH1D* h_svpint_pos_low_dt = new TH1D("h_svpint_pos_low_dt", "Distribution of Low dt Side Veto Panel Integral Values", 240, -400, 2000);

    TH1D* h_svpint_neg_low_dt = new TH1D("h_svpint_neg_low_dt", "Distribution of Low dt Side Veto Panel Integral Values", 240, -400, 2000);

    TH1D* h_tvpint_cent_spike = new TH1D("h_tvpint_cent_spike", "Distribution of All Central Spike Top Veto Panel Integral Values", 60, -200, 400);

    TH1D* h_tvpint_pos_low_dt = new TH1D("h_tvpint_pos_low_dt", "Distribution of Low dt Top Veto Panel Integral Values", 60, -200, 400);

    TH1D* h_tvpint_neg_low_dt = new TH1D("h_tvpint_neg_low_dt", "Distribution of Low dt Top Veto Panel Integral Values", 60, -200, 400);

    TH1D* h_vec_size = new TH1D("h_vec_size", "Distribution of Number of Internally Triggered Events Between High Light and Low Light LED Events", 5, 0, 5);

    TH1D* h_vpint_16 = new TH1D("h_vpint_16", "Beam Spill SiPM 1 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_17 = new TH1D("h_vpint_17", "Beam Spill SiPM 2 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_18 = new TH1D("h_vpint_18", "Beam Spill SiPM 3 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_19 = new TH1D("h_vpint_19", "Beam Spill SiPM 4 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_20 = new TH1D("h_vpint_20", "Beam Spill SiPM 5 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_21 = new TH1D("h_vpint_21", "Beam Spill SiPM 6 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_22 = new TH1D("h_vpint_22", "Beam Spill SiPM 7 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_23 = new TH1D("h_vpint_23", "Beam Spill SiPM 8 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_top = new TH1D("h_vpint_top", "Beam Spill SiPM 9 & 10 Integral Values", 350, -500, 3000);

    auto p_ks_test = new TGraph();

    TFile *f; 
    //your root file location here
        if(gSystem->AccessPathName(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/int_trig_events/Unblinded_Beam_Spill_Data_S25-S26_GR_%i.root", run))){            // Add/remove "EC1us_" to change to/from 1us echo cut for Spring 2025
        //    if (gSystem->AccessPathName(Form("/data9/coherent/data/d2o/processedData/Event61_Analysis_24hrs_Run%i.root", run))) {
        //    if(gSystem->AccessPathName(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/Event61_Analysis_24hrs_Run%i.root", run))){
        cout << "Could not open file" << endl;
        return -1; 
    } else{
        f = new TFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/int_trig_events/Unblinded_Beam_Spill_Data_S25-S26_GR_%i.root", run));
        //        f = new TFile(Form("/data9/coherent/data/d2o/processedData/Event61_Analysis_24hrs_Run%i.root", run));
        //        f = new TFile(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/Event61_Analysis_24hrs_Run%i.root", run));
    }

    // Create a text string, which is used to output the text file

    std::vector<double> MWhrPerRun_Runs;

    std::vector<double> MWhrPerRun_MWhr;
    
    std::vector<double> num_vec;

    string num_str;

    string intval;

    // Read from the text file

    // ifstream ReadMWHFile("F24-S25_MWhrPerRun.txt");
    ifstream ReadMWHFile("F24-S26_MWhrPerRun.txt");

    // Use a while loop together with the getline() function to read the file line by line

    while (getline(ReadMWHFile, intval)) {

        num_vec.clear();

        num_str.clear();

        intval.push_back('\t');

        for (char ch : intval) {
            
            if (ch != '\t') {num_str.push_back(ch);}
            
            else if (ch == '\t') {num_vec.push_back(stod(num_str)); num_str.clear();}
        
        }

        MWhrPerRun_Runs.push_back(num_vec[0]);

        MWhrPerRun_MWhr.push_back(num_vec[1]);
        
    }

    // Close the file

    ReadMWHFile.close();

    ReadMWHFile.clear();

    TTree* t = (TTree*)f->Get("eventTree");

    // Declaration of leaf types
    Int_t           fnum;
    Int_t           entry;
    Long64_t        utime;
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
    Long64_t        lct;
    Long64_t        ldt;
    Bool_t          issue;
    Long64_t        e61st;
    Bool_t          e61sig;
    Int_t           selcut;
    Long64_t        wftime;
    Long64_t        vptime16;
    Long64_t        vptime17;
    Long64_t        vptime18;
    Long64_t        vptime19;
    Long64_t        vptime20;
    Long64_t        vptime21;
    Long64_t        vptime22;
    Long64_t        vptime23;
    Long64_t        vptimetop;
    Long64_t        vpint16;
    Long64_t        vpint17;
    Long64_t        vpint18;
    Long64_t        vpint19;
    Long64_t        vpint20;
    Long64_t        vpint21;
    Long64_t        vpint22;
    Long64_t        vpint23;
    Long64_t        vpinttop;

    // List of branches
    TBranch        *b_fnum;
    TBranch        *b_entry;
    TBranch        *b_utime;
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
    TBranch        *b_lct;
    TBranch        *b_ldt;
    TBranch        *b_issue;
    TBranch        *b_e61st;
    TBranch        *b_e61sig;
    TBranch        *b_selcut;
    TBranch        *b_wftime;
    TBranch        *b_vptime16;
    TBranch        *b_vptime17;
    TBranch        *b_vptime18;
    TBranch        *b_vptime19;
    TBranch        *b_vptime20;
    TBranch        *b_vptime21;
    TBranch        *b_vptime22;
    TBranch        *b_vptime23;
    TBranch        *b_vptimetop;
    TBranch        *b_vpint16;
    TBranch        *b_vpint17;
    TBranch        *b_vpint18;
    TBranch        *b_vpint19;
    TBranch        *b_vpint20;
    TBranch        *b_vpint21;
    TBranch        *b_vpint22;
    TBranch        *b_vpint23;
    TBranch        *b_vpinttop;

    t->SetBranchAddress("fnum", &fnum, &b_fnum);
    t->SetBranchAddress("entry", &entry, &b_entry);
    t->SetBranchAddress("utime", &utime, &b_utime);
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
    t->SetBranchAddress("lct", &lct, &b_lct);
    t->SetBranchAddress("ldt", &ldt, &b_ldt);
    t->SetBranchAddress("issue", &issue, &b_issue);
    t->SetBranchAddress("e61st", &e61st, &b_e61st);
    t->SetBranchAddress("e61sig", &e61sig, &b_e61sig);
    t->SetBranchAddress("selcut", &selcut, &b_selcut);
    t->SetBranchAddress("wftime", &wftime, &b_wftime);
    t->SetBranchAddress("vptime16", &vptime16, &b_vptime16);
    t->SetBranchAddress("vptime17", &vptime17, &b_vptime17);
    t->SetBranchAddress("vptime18", &vptime18, &b_vptime18);
    t->SetBranchAddress("vptime19", &vptime19, &b_vptime19);
    t->SetBranchAddress("vptime20", &vptime20, &b_vptime20);
    t->SetBranchAddress("vptime21", &vptime21, &b_vptime21);
    t->SetBranchAddress("vptime22", &vptime22, &b_vptime22);
    t->SetBranchAddress("vptime23", &vptime23, &b_vptime23);
    t->SetBranchAddress("vptimetop", &vptimetop, &b_vptimetop);
    t->SetBranchAddress("vpint16", &vpint16, &b_vpint16);
    t->SetBranchAddress("vpint17", &vpint17, &b_vpint17);
    t->SetBranchAddress("vpint18", &vpint18, &b_vpint18);
    t->SetBranchAddress("vpint19", &vpint19, &b_vpint19);
    t->SetBranchAddress("vpint20", &vpint20, &b_vpint20);
    t->SetBranchAddress("vpint21", &vpint21, &b_vpint21);
    t->SetBranchAddress("vpint22", &vpint22, &b_vpint22);
    t->SetBranchAddress("vpint23", &vpint23, &b_vpint23);
    t->SetBranchAddress("vpinttop", &vpinttop, &b_vpinttop);

    // A virtual index so the code always iterates over the TTree in the correct order
    t->BuildIndex("fnum", "entry");
    TTreeIndex *evIndex = (TTreeIndex*)t->GetTreeIndex();
    Long64_t *sortedEntries = evIndex->GetIndex();

    // Create and open a text file
    ofstream UBFile(Form("Unblinded_BeamSpill_dt_Int_%i.txt", run));
    ofstream KSFile(Form("Unblinded_BeamSpill_KS_Test_%i.txt", run));

    for(int iEvent = 0; iEvent < t->GetEntries(); iEvent++){

        Long64_t tentry = t->LoadTree(sortedEntries[iEvent]);
        
        b_fnum->GetEntry(tentry);
        b_entry->GetEntry(tentry);
        b_utime->GetEntry(tentry);
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
        b_lct->GetEntry(tentry);
        b_ldt->GetEntry(tentry);
        b_issue->GetEntry(tentry);
        b_e61st->GetEntry(tentry);
        b_e61sig->GetEntry(tentry);
        b_selcut->GetEntry(tentry);
        b_wftime->GetEntry(tentry);
        b_vptime16->GetEntry(tentry);
        b_vptime17->GetEntry(tentry);
        b_vptime18->GetEntry(tentry);
        b_vptime19->GetEntry(tentry);
        b_vptime20->GetEntry(tentry);
        b_vptime21->GetEntry(tentry);
        b_vptime22->GetEntry(tentry);
        b_vptime23->GetEntry(tentry);
        b_vptimetop->GetEntry(tentry);
        b_vpint16->GetEntry(tentry);
        b_vpint17->GetEntry(tentry);
        b_vpint18->GetEntry(tentry);
        b_vpint19->GetEntry(tentry);
        b_vpint20->GetEntry(tentry);
        b_vpint21->GetEntry(tentry);
        b_vpint22->GetEntry(tentry);
        b_vpint23->GetEntry(tentry);
        b_vpinttop->GetEntry(tentry);

        struct pulse avg_pulse;

        avg_pulse.file_num = fnum;
        avg_pulse.entry = entry;
        avg_pulse.unix_time = utime;
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
        avg_pulse.last_cosmics_time = lct;
        avg_pulse.last_det_time = ldt;
        avg_pulse.issue = issue;
        avg_pulse.ev61_start_time = e61st;
        avg_pulse.ev61_significant = e61sig;
        avg_pulse.selection_cut = selcut;

        avg_pulse.wf_time = wftime;
        avg_pulse.vp_start_16 = vptime16;
        avg_pulse.vp_start_17 = vptime17;
        avg_pulse.vp_start_18 = vptime18;
        avg_pulse.vp_start_19 = vptime19;
        avg_pulse.vp_start_20 = vptime20;
        avg_pulse.vp_start_21 = vptime21;
        avg_pulse.vp_start_22 = vptime22;
        avg_pulse.vp_start_23 = vptime23;
        avg_pulse.vp_start_top = vptimetop;
        avg_pulse.vp_int_16 = vpint16;
        avg_pulse.vp_int_17 = vpint17;
        avg_pulse.vp_int_18 = vpint18;
        avg_pulse.vp_int_19 = vpint19;
        avg_pulse.vp_int_20 = vpint20;
        avg_pulse.vp_int_21 = vpint21;
        avg_pulse.vp_int_22 = vpint22;
        avg_pulse.vp_int_23 = vpint23;
        avg_pulse.vp_int_top = vpinttop;

        int low_bound_1 = -20000;
        int high_bound_1 = -5000;
        int low_bound_2 = -500;
        int high_bound_2 = 0;
        int low_bound_3 = 400;
        int high_bound_3 = 900;
        int low_bound_4 = 1500;
        int high_bound_4 = 4500;
        int low_bound_5 = 5000;
        int high_bound_5 = 20000;
        int low_bound_6 = 18000;
        int high_bound_6 = 19000;
        int low_bound_7 = 22000;
        int high_bound_7 = 23000;
        int low_bound_8 = -2700;
        int high_bound_8 = -2200;
        int low_bound_9 = -500;
        int high_bound_9 = 0;

        int neg_peak_low_bound = -500;
        int neg_peak_high_bound = 0;
        int before_ev61_low_bound = -2700;
        int before_ev61_high_bound = -2200;
        int after_ev61_low_bound = -500;
        int after_ev61_high_bound = 0;

        int time_shift = 0;

        // if (avg_pulse.file_num > 28000) {continue;}          // Dividing line between F25 & S26

        if (avg_pulse.file_num < 20000) {fall_is_true_spring_is_false = true; ev61_offset = 2000 - 625 - 200 + time_shift;}

        else if (avg_pulse.file_num >= 20000 && avg_pulse.file_num < 24000) {fall_is_true_spring_is_false = false; ev61_offset = -19000 - 625 - 662 + time_shift;}

        else if (avg_pulse.file_num >= 24000 && avg_pulse.file_num < 28000) {fall_is_true_spring_is_false = false; ev61_offset = -19000 - 625 - 662 + 67 + time_shift;}

        else if (avg_pulse.file_num >= 28000) {fall_is_true_spring_is_false = false; ev61_offset = -19000 - 625 - 662 + 18 + time_shift;}

        if (fall_is_true_spring_is_false) {

            plot_min = 1000;

            plot_max = 9960;

            neg_peak_low_bound = -350;

            neg_peak_high_bound = -50;

            before_ev61_low_bound = -3000; // low_bound_8;
            
            before_ev61_high_bound = -2500; // high_bound_8;
            
            after_ev61_low_bound = -500; // low_bound_9;
            
            after_ev61_high_bound = 0; // high_bound_9;

        }

        else if (!fall_is_true_spring_is_false) {

            plot_min = -19960;

            plot_max = -11000;

            neg_peak_low_bound = -800;

            neg_peak_high_bound = -500;

            before_ev61_low_bound = 18500; // low_bound_6;
            
            before_ev61_high_bound = 19000; // high_bound_6;
            
            after_ev61_low_bound = 21000; // low_bound_7;
            
            after_ev61_high_bound = 21500; // high_bound_7;

        }

        /*
        
        Only looking for events between high light LED events and low light LED events
        
        */

        // After every high light LED event, "turn on" internally triggered event storing

        if (avg_pulse.trigger == 8) {bt_HL_n_LL = true; h_HLint->Fill(avg_pulse.energy);}

        // Store internally triggered events in vector

        if (bt_HL_n_LL && (avg_pulse.selection_cut >= selection_cut_level || avg_pulse.beam) && avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {pulses.push_back(avg_pulse); h_eventint->Fill(avg_pulse.energy); h_eventvpint->Fill(avg_pulse.all_vp_energy);}
        
        // After every low light LED event, "turn off" internally triggered event storing and search through compiled vector of events

        if (avg_pulse.trigger == 16) {

            // Iterate over all events in compiled vector

            for (size_t iVec = 0; iVec < pulses.size(); iVec++) {

                if (pulses[iVec].beam) {

                    take_sample_count += 1;

                    if (iVec > 0 && (pulses[iVec - 1].trigger != 1 && pulses[iVec - 1].trigger != 3 && pulses[iVec - 1].trigger != 33 && pulses[iVec - 1].trigger != 35)) {

                        long long eventstarttime1 = pulses[iVec - 1].start - ev61_offset;

                        if ((eventstarttime1 - pulses[iVec].ev61_start_time) >= -1 * sample_range && (eventstarttime1 - pulses[iVec].ev61_start_time) < 0 - ev61_offset) {

                            // Plot negative dt values

                            if ((eventstarttime1 - pulses[iVec].ev61_start_time) >= -10000 && (eventstarttime1 - pulses[iVec].ev61_start_time) <= 10000) {num_events_low_dt += 1; h_detint_low_dt->Fill(pulses[iVec - 1].energy); h_avpint_low_dt->Fill(pulses[iVec - 1].all_vp_energy);}
                
                            if ((eventstarttime1 - pulses[iVec].ev61_start_time) > -500 && (eventstarttime1 - pulses[iVec].ev61_start_time) < 500) {num_events_low_dt_cent_spike += 1; h_int_cent_spike->Fill(pulses[iVec - 1].energy); h_avpint_cent_spike->Fill(pulses[iVec - 1].all_vp_energy); h_svpint_cent_spike->Fill(pulses[iVec - 1].side_vp_energy); h_tvpint_cent_spike->Fill(pulses[iVec - 1].top_vp_energy);}          // cout << "\n" << "Run, Event, triggerBits = " << iRun << " " << iEnt + 1 << " " << avg_pulse.trigger << endl;}

                            else if ((eventstarttime1 - pulses[iVec].ev61_start_time) <= -500 && (eventstarttime1 - pulses[iVec].ev61_start_time) >= -10000) {num_events_low_dt_no_spike_neg += 1; h_int_neg_low_dt->Fill(pulses[iVec - 1].energy); h_avpint_neg_low_dt->Fill(pulses[iVec - 1].all_vp_energy); h_svpint_neg_low_dt->Fill(pulses[iVec - 1].side_vp_energy); h_tvpint_neg_low_dt->Fill(pulses[iVec - 1].top_vp_energy);}

                            else if ((eventstarttime1 - pulses[iVec].ev61_start_time) >= 500 && (eventstarttime1 - pulses[iVec].ev61_start_time) <= 10000) {num_events_low_dt_no_spike_pos += 1; h_int_pos_low_dt->Fill(pulses[iVec - 1].energy); h_avpint_pos_low_dt->Fill(pulses[iVec - 1].all_vp_energy); h_svpint_pos_low_dt->Fill(pulses[iVec - 1].side_vp_energy); h_tvpint_pos_low_dt->Fill(pulses[iVec - 1].top_vp_energy);}

                            if (eventstarttime1 - pulses[iVec].ev61_start_time >= 1000 && eventstarttime1 - pulses[iVec].ev61_start_time <= 5000) {

                                hist_events += 1;

                                run_events += 1;

                                tot_detint += pulses[iVec - 1].energy;

                                h_runnum->Fill(pulses[iVec - 1].file_num);

                                h_dt_61->Fill(eventstarttime1 - pulses[iVec].ev61_start_time);

                                h_detint->Fill(pulses[iVec - 1].energy);

                                if (pulses[iVec - 1].selection_cut >= 1) {h_detint_1->Fill(pulses[iVec - 1].energy);}

                                if (pulses[iVec - 1].selection_cut >= 2) {h_detint_2->Fill(pulses[iVec - 1].energy);}

                                if (pulses[iVec - 1].selection_cut >= 3) {h_detint_3->Fill(pulses[iVec - 1].energy);}

                                if (pulses[iVec - 1].selection_cut >= 4) {h_detint_4->Fill(pulses[iVec - 1].energy);}

                                if (pulses[iVec - 1].selection_cut >= 5) {h_detint_5->Fill(pulses[iVec - 1].energy);}

                                if (pulses[iVec - 1].selection_cut >= 6) {h_detint_6->Fill(pulses[iVec - 1].energy);}

                            }

                            h_runnum_v_dt->Fill(pulses[iVec - 1].file_num, eventstarttime1 - pulses[iVec].ev61_start_time);

                            h_runnum_v_detint->Fill(pulses[iVec - 1].file_num, pulses[iVec - 1].energy);
                            
                            h_dt_61_1->Fill(eventstarttime1 - pulses[iVec].ev61_start_time);

                            h_dt_61_2->Fill(eventstarttime1 - pulses[iVec].ev61_start_time);

                            h_dt_61_3->Fill(eventstarttime1 - pulses[iVec].ev61_start_time);

                            h_dt_61_4->Fill(eventstarttime1 - pulses[iVec].ev61_start_time);

                            h_dt_v_tmuon->Fill(eventstarttime1 - pulses[iVec].ev61_start_time, eventstarttime1 - pulses[iVec - 1].last_cosmics_time);

                            h_dt_v_detint->Fill(eventstarttime1 - pulses[iVec].ev61_start_time, pulses[iVec - 1].energy);

                            UBFile << eventstarttime1 - pulses[iVec].ev61_start_time << "\t" << pulses[iVec - 1].energy << "\n";

                            h_vpint_16->Fill(pulses[iVec - 1].vp_int_16); // if (pulses[iVec].vp_int_16 > 150) {h_vpint_16->Fill(pulses[iVec].vp_int_16);}          // pulses[iVec].vp_start_16 != 0

                            h_vpint_17->Fill(pulses[iVec - 1].vp_int_17); // if (pulses[iVec].vp_int_17 > 150) {h_vpint_17->Fill(pulses[iVec].vp_int_17);}

                            h_vpint_18->Fill(pulses[iVec - 1].vp_int_18); // if (pulses[iVec].vp_int_18 > 400) {h_vpint_18->Fill(pulses[iVec].vp_int_18);}

                            h_vpint_19->Fill(pulses[iVec - 1].vp_int_19); // if (pulses[iVec].vp_int_19 > 150) {h_vpint_19->Fill(pulses[iVec].vp_int_19);}

                            h_vpint_20->Fill(pulses[iVec - 1].vp_int_20); // if (pulses[iVec].vp_int_20 > 150) {h_vpint_20->Fill(pulses[iVec].vp_int_20);}

                            h_vpint_21->Fill(pulses[iVec - 1].vp_int_21); // if (pulses[iVec].vp_int_21 > 150) {h_vpint_21->Fill(pulses[iVec].vp_int_21);}

                            h_vpint_22->Fill(pulses[iVec - 1].vp_int_22); // if (pulses[iVec].vp_int_22 > 150) {h_vpint_22->Fill(pulses[iVec].vp_int_22);}

                            h_vpint_23->Fill(pulses[iVec - 1].vp_int_23); // if (pulses[iVec].vp_int_23 > 150) {h_vpint_23->Fill(pulses[iVec].vp_int_23);}

                            h_vpint_top->Fill(pulses[iVec - 1].vp_int_top); // if (pulses[iVec].vp_int_top > 250) {h_vpint_top->Fill(pulses[iVec].vp_int_top);}

                            if (eventstarttime1 - pulses[iVec].ev61_start_time > low_bound_1 && eventstarttime1 - pulses[iVec].ev61_start_time < high_bound_1) {h_detint_R1->Fill(pulses[iVec - 1].energy);}

                            if (eventstarttime1 - pulses[iVec].ev61_start_time > low_bound_2 && eventstarttime1 - pulses[iVec].ev61_start_time < high_bound_2) {h_detint_R2->Fill(pulses[iVec - 1].energy); }

                            if (eventstarttime1 - pulses[iVec].ev61_start_time > low_bound_3 && eventstarttime1 - pulses[iVec].ev61_start_time < high_bound_3) {h_detint_R3->Fill(pulses[iVec - 1].energy);}

                            if (eventstarttime1 - pulses[iVec].ev61_start_time > low_bound_4 && eventstarttime1 - pulses[iVec].ev61_start_time < high_bound_4) {h_detint_R4->Fill(pulses[iVec - 1].energy);}

                            if (eventstarttime1 - pulses[iVec].ev61_start_time > low_bound_5 && eventstarttime1 - pulses[iVec].ev61_start_time < high_bound_5) {h_detint_R5->Fill(pulses[iVec - 1].energy);}

                            if (eventstarttime1 - pulses[iVec].ev61_start_time > before_ev61_low_bound && eventstarttime1 - pulses[iVec].ev61_start_time < before_ev61_high_bound) {h_detint_R6->Fill(pulses[iVec - 1].energy);}

                            if (eventstarttime1 - pulses[iVec].ev61_start_time > after_ev61_low_bound && eventstarttime1 - pulses[iVec].ev61_start_time < after_ev61_high_bound) {h_detint_R7->Fill(pulses[iVec - 1].energy);}

                            if (eventstarttime1 - pulses[iVec].ev61_start_time > before_ev61_low_bound && eventstarttime1 - pulses[iVec].ev61_start_time < before_ev61_high_bound && pulses[iVec - 1].energy > 150 && pulses[iVec - 1].energy < 250) {cout << "\n" << "Before Ev61: Run Number, Event Number, triggerBits, Peak Time, Event 61 Time, dt, Energy = " << pulses[iVec - 1].file_num << "\t" << pulses[iVec - 1].entry << "\t" << pulses[iVec - 1].trigger << "\t" << eventstarttime1 << "\t" << pulses[iVec].ev61_start_time << "\t" << eventstarttime1 - pulses[iVec].ev61_start_time << "\t" << pulses[iVec - 1].energy << endl;}

                            if (eventstarttime1 - pulses[iVec].ev61_start_time > after_ev61_low_bound && eventstarttime1 - pulses[iVec].ev61_start_time < after_ev61_high_bound && pulses[iVec - 1].energy > 150 && pulses[iVec - 1].energy < 250) {cout << "\n" << "After Ev61: Run Number, Event Number, triggerBits, Peak Time, Event 61 Time, dt, Energy = " << pulses[iVec - 1].file_num << "\t" << pulses[iVec - 1].entry << "\t" << pulses[iVec - 1].trigger << "\t" << eventstarttime1 << "\t" << pulses[iVec].ev61_start_time << "\t" << eventstarttime1 - pulses[iVec].ev61_start_time << "\t" << pulses[iVec - 1].energy << endl;}

                            // if (pulses[iVec - 1].energy > 250 && pulses[iVec - 1].energy < 500) {cout << "\n" << "Run Number = " << pulses[iVec - 1].file_num << ", Event Number = " << pulses[iVec - 1].entry << ", triggerBits = " << pulses[iVec - 1].trigger << ", Peak Time = " << eventstarttime1 << ", Event 61 Time = " << pulses[iVec].ev61_start_time << ", dt = " << eventstarttime1 - pulses[iVec].ev61_start_time << ", detint = " << pulses[iVec - 1].energy << endl;}

                        }

                    }

                    if (pulses[iVec].ev61_significant && pulses[iVec].selection_cut >= selection_cut_level) {

                        long long eventstarttime2 = pulses[iVec].start - ev61_offset;

                        // Plot low dt values

                        /*if ((eventstarttime2 - pulses[iVec].ev61_start_time) >= -10000 && (eventstarttime2 - pulses[iVec].ev61_start_time) <= 10000) {num_events_low_dt += 1; h_detint_low_dt->Fill(pulses[iVec].energy); h_avpint_low_dt->Fill(pulses[iVec].all_vp_energy);}
            
                        if ((eventstarttime2 - pulses[iVec].ev61_start_time) > -500 && (eventstarttime2 - pulses[iVec].ev61_start_time) < 500) {num_events_low_dt_cent_spike += 1; h_int_cent_spike->Fill(pulses[iVec].energy); h_avpint_cent_spike->Fill(pulses[iVec].all_vp_energy); h_svpint_cent_spike->Fill(pulses[iVec].side_vp_energy); h_tvpint_cent_spike->Fill(pulses[iVec].top_vp_energy);}          // cout << "\n" << "Run, Event, triggerBits = " << iRun << " " << iEnt + 1 << " " << avg_pulse.trigger << endl;}

                        else if ((eventstarttime2 - pulses[iVec].ev61_start_time) <= -500 && (eventstarttime2 - pulses[iVec].ev61_start_time) >= -10000) {num_events_low_dt_no_spike_neg += 1; h_int_neg_low_dt->Fill(pulses[iVec].energy); h_avpint_neg_low_dt->Fill(pulses[iVec].all_vp_energy); h_svpint_neg_low_dt->Fill(pulses[iVec].side_vp_energy); h_tvpint_neg_low_dt->Fill(pulses[iVec].top_vp_energy);}

                        else if ((eventstarttime2 - pulses[iVec].ev61_start_time) >= 500 && (eventstarttime2 - pulses[iVec].ev61_start_time) <= 10000) {num_events_low_dt_no_spike_pos += 1; h_int_pos_low_dt->Fill(pulses[iVec].energy); h_avpint_pos_low_dt->Fill(pulses[iVec].all_vp_energy); h_svpint_pos_low_dt->Fill(pulses[iVec].side_vp_energy); h_tvpint_pos_low_dt->Fill(pulses[iVec].top_vp_energy);}

                        hist_events += 1;

                        run_events += 1;

                        tot_detint += pulses[iVec].energy;

                        h_runnum_v_dt->Fill(pulses[iVec].file_num, eventstarttime2 - pulses[iVec].ev61_start_time);

                        h_runnum_v_detint->Fill(pulses[iVec].file_num, pulses[iVec].energy);
                        
                        h_dt_61->Fill(eventstarttime2 - pulses[iVec].ev61_start_time);
                        
                        h_dt_61_1->Fill(eventstarttime2 - pulses[iVec].ev61_start_time);

                        h_dt_61_2->Fill(eventstarttime2 - pulses[iVec].ev61_start_time);

                        h_dt_61_3->Fill(eventstarttime2 - pulses[iVec].ev61_start_time);

                        h_dt_61_4->Fill(eventstarttime2 - pulses[iVec].ev61_start_time);

                        h_dt_v_tmuon->Fill(eventstarttime2 - pulses[iVec].ev61_start_time, eventstarttime2 - pulses[iVec].last_cosmics_time);

                        h_detint->Fill(pulses[iVec].energy);

                        h_dt_v_detint->Fill(eventstarttime2 - pulses[iVec].ev61_start_time, pulses[iVec].energy);

                        UBFile << eventstarttime2 - pulses[iVec].ev61_start_time << "\t" << pulses[iVec].energy << "\n";

                        h_vpint_16->Fill(pulses[iVec].vp_int_16); // if (pulses[iVec].vp_int_16 > 150) {h_vpint_16->Fill(pulses[iVec].vp_int_16);}          // pulses[iVec].vp_start_16 != 0

                        h_vpint_17->Fill(pulses[iVec].vp_int_17); // if (pulses[iVec].vp_int_17 > 150) {h_vpint_17->Fill(pulses[iVec].vp_int_17);}

                        h_vpint_18->Fill(pulses[iVec].vp_int_18); // if (pulses[iVec].vp_int_18 > 400) {h_vpint_18->Fill(pulses[iVec].vp_int_18);}

                        h_vpint_19->Fill(pulses[iVec].vp_int_19); // if (pulses[iVec].vp_int_19 > 150) {h_vpint_19->Fill(pulses[iVec].vp_int_19);}

                        h_vpint_20->Fill(pulses[iVec].vp_int_20); // if (pulses[iVec].vp_int_20 > 150) {h_vpint_20->Fill(pulses[iVec].vp_int_20);}

                        h_vpint_21->Fill(pulses[iVec].vp_int_21); // if (pulses[iVec].vp_int_21 > 150) {h_vpint_21->Fill(pulses[iVec].vp_int_21);}

                        h_vpint_22->Fill(pulses[iVec].vp_int_22); // if (pulses[iVec].vp_int_22 > 150) {h_vpint_22->Fill(pulses[iVec].vp_int_22);}

                        h_vpint_23->Fill(pulses[iVec].vp_int_23); // if (pulses[iVec].vp_int_23 > 150) {h_vpint_23->Fill(pulses[iVec].vp_int_23);}

                        h_vpint_top->Fill(pulses[iVec].vp_int_top); // if (pulses[iVec].vp_int_top > 250) {h_vpint_top->Fill(pulses[iVec].vp_int_top);}*/

                    }

                    if (iVec < pulses.size() - 1 && (pulses[iVec + 1].trigger != 1 && pulses[iVec + 1].trigger != 3 && pulses[iVec + 1].trigger != 33 && pulses[iVec + 1].trigger != 35)) {

                        long long eventstarttime3 = pulses[iVec + 1].start - ev61_offset;

                        if ((eventstarttime3 - pulses[iVec].ev61_start_time) <= sample_range && (eventstarttime3 - pulses[iVec].ev61_start_time) > 0 - ev61_offset) {

                            // Plot positive dt values

                            if ((eventstarttime3 - pulses[iVec].ev61_start_time) >= -10000 && (eventstarttime3 - pulses[iVec].ev61_start_time) <= 10000) {num_events_low_dt += 1; h_detint_low_dt->Fill(pulses[iVec + 1].energy); h_avpint_low_dt->Fill(pulses[iVec + 1].all_vp_energy);}
                
                            if ((eventstarttime3 - pulses[iVec].ev61_start_time) > -500 && (eventstarttime3 - pulses[iVec].ev61_start_time) < 500) {num_events_low_dt_cent_spike += 1; h_int_cent_spike->Fill(pulses[iVec + 1].energy); h_avpint_cent_spike->Fill(pulses[iVec + 1].all_vp_energy); h_svpint_cent_spike->Fill(pulses[iVec + 1].side_vp_energy); h_tvpint_cent_spike->Fill(pulses[iVec + 1].top_vp_energy);}          // cout << "\n" << "Run, Event, triggerBits = " << iRun << " " << iEnt + 1 << " " << avg_pulse.trigger << endl;}

                            else if ((eventstarttime3 - pulses[iVec].ev61_start_time) <= -500 && (eventstarttime3 - pulses[iVec].ev61_start_time) >= -10000) {num_events_low_dt_no_spike_neg += 1; h_int_neg_low_dt->Fill(pulses[iVec + 1].energy); h_avpint_neg_low_dt->Fill(pulses[iVec + 1].all_vp_energy); h_svpint_neg_low_dt->Fill(pulses[iVec + 1].side_vp_energy); h_tvpint_neg_low_dt->Fill(pulses[iVec + 1].top_vp_energy);}

                            else if ((eventstarttime3 - pulses[iVec].ev61_start_time) >= 500 && (eventstarttime3 - pulses[iVec].ev61_start_time) <= 10000) {num_events_low_dt_no_spike_pos += 1; h_int_pos_low_dt->Fill(pulses[iVec + 1].energy); h_avpint_pos_low_dt->Fill(pulses[iVec + 1].all_vp_energy); h_svpint_pos_low_dt->Fill(pulses[iVec + 1].side_vp_energy); h_tvpint_pos_low_dt->Fill(pulses[iVec + 1].top_vp_energy);}

                            if (eventstarttime3 - pulses[iVec].ev61_start_time >= 1000 && eventstarttime3 - pulses[iVec].ev61_start_time <= 5000) {

                                hist_events += 1;

                                run_events += 1;

                                tot_detint += pulses[iVec + 1].energy;

                                h_runnum->Fill(pulses[iVec + 1].file_num);

                                h_dt_61->Fill(eventstarttime3 - pulses[iVec].ev61_start_time);

                                h_detint->Fill(pulses[iVec + 1].energy);

                                if (pulses[iVec + 1].selection_cut >= 1) {h_detint_1->Fill(pulses[iVec + 1].energy);}

                                if (pulses[iVec + 1].selection_cut >= 2) {h_detint_2->Fill(pulses[iVec + 1].energy);}

                                if (pulses[iVec + 1].selection_cut >= 3) {h_detint_3->Fill(pulses[iVec + 1].energy);}

                                if (pulses[iVec + 1].selection_cut >= 4) {h_detint_4->Fill(pulses[iVec + 1].energy);}

                                if (pulses[iVec + 1].selection_cut >= 5) {h_detint_5->Fill(pulses[iVec + 1].energy);}

                                if (pulses[iVec + 1].selection_cut >= 6) {h_detint_6->Fill(pulses[iVec + 1].energy);}

                            }

                            h_runnum_v_dt->Fill(pulses[iVec + 1].file_num, eventstarttime3 - pulses[iVec].ev61_start_time);

                            h_runnum_v_detint->Fill(pulses[iVec + 1].file_num, pulses[iVec + 1].energy);
                            
                            h_dt_61_1->Fill(eventstarttime3 - pulses[iVec].ev61_start_time);

                            h_dt_61_2->Fill(eventstarttime3 - pulses[iVec].ev61_start_time);

                            h_dt_61_3->Fill(eventstarttime3 - pulses[iVec].ev61_start_time);

                            h_dt_61_4->Fill(eventstarttime3 - pulses[iVec].ev61_start_time);

                            h_dt_v_tmuon->Fill(eventstarttime3 - pulses[iVec].ev61_start_time, eventstarttime3 - pulses[iVec + 1].last_cosmics_time);

                            h_dt_v_detint->Fill(eventstarttime3 - pulses[iVec].ev61_start_time, pulses[iVec + 1].energy);

                            UBFile << eventstarttime3 - pulses[iVec].ev61_start_time << "\t" << pulses[iVec + 1].energy << "\n";

                            h_vpint_16->Fill(pulses[iVec + 1].vp_int_16); // if (pulses[iVec].vp_int_16 > 150) {h_vpint_16->Fill(pulses[iVec].vp_int_16);}          // pulses[iVec].vp_start_16 != 0

                            h_vpint_17->Fill(pulses[iVec + 1].vp_int_17); // if (pulses[iVec].vp_int_17 > 150) {h_vpint_17->Fill(pulses[iVec].vp_int_17);}

                            h_vpint_18->Fill(pulses[iVec + 1].vp_int_18); // if (pulses[iVec].vp_int_18 > 400) {h_vpint_18->Fill(pulses[iVec].vp_int_18);}

                            h_vpint_19->Fill(pulses[iVec + 1].vp_int_19); // if (pulses[iVec].vp_int_19 > 150) {h_vpint_19->Fill(pulses[iVec].vp_int_19);}

                            h_vpint_20->Fill(pulses[iVec + 1].vp_int_20); // if (pulses[iVec].vp_int_20 > 150) {h_vpint_20->Fill(pulses[iVec].vp_int_20);}

                            h_vpint_21->Fill(pulses[iVec + 1].vp_int_21); // if (pulses[iVec].vp_int_21 > 150) {h_vpint_21->Fill(pulses[iVec].vp_int_21);}

                            h_vpint_22->Fill(pulses[iVec + 1].vp_int_22); // if (pulses[iVec].vp_int_22 > 150) {h_vpint_22->Fill(pulses[iVec].vp_int_22);}

                            h_vpint_23->Fill(pulses[iVec + 1].vp_int_23); // if (pulses[iVec].vp_int_23 > 150) {h_vpint_23->Fill(pulses[iVec].vp_int_23);}

                            h_vpint_top->Fill(pulses[iVec + 1].vp_int_top); // if (pulses[iVec].vp_int_top > 250) {h_vpint_top->Fill(pulses[iVec].vp_int_top);}

                            if (eventstarttime3 - pulses[iVec].ev61_start_time > low_bound_1 && eventstarttime3 - pulses[iVec].ev61_start_time < high_bound_1) {h_detint_R1->Fill(pulses[iVec + 1].energy);}

                            if (eventstarttime3 - pulses[iVec].ev61_start_time > low_bound_2 && eventstarttime3 - pulses[iVec].ev61_start_time < high_bound_2) {h_detint_R2->Fill(pulses[iVec + 1].energy);}

                            if (eventstarttime3 - pulses[iVec].ev61_start_time > low_bound_3 && eventstarttime3 - pulses[iVec].ev61_start_time < high_bound_3) {h_detint_R3->Fill(pulses[iVec + 1].energy);}

                            if (eventstarttime3 - pulses[iVec].ev61_start_time > low_bound_4 && eventstarttime3 - pulses[iVec].ev61_start_time < high_bound_4) {h_detint_R4->Fill(pulses[iVec + 1].energy);}

                            if (eventstarttime3 - pulses[iVec].ev61_start_time > low_bound_5 && eventstarttime3 - pulses[iVec].ev61_start_time < high_bound_5) {h_detint_R5->Fill(pulses[iVec + 1].energy);}

                            if (eventstarttime3 - pulses[iVec].ev61_start_time > before_ev61_low_bound && eventstarttime3 - pulses[iVec].ev61_start_time < before_ev61_high_bound) {h_detint_R6->Fill(pulses[iVec + 1].energy);}

                            if (eventstarttime3 - pulses[iVec].ev61_start_time > after_ev61_low_bound && eventstarttime3 - pulses[iVec].ev61_start_time < after_ev61_high_bound) {h_detint_R7->Fill(pulses[iVec + 1].energy);}

                            if (eventstarttime3 - pulses[iVec].ev61_start_time > before_ev61_low_bound && eventstarttime3 - pulses[iVec].ev61_start_time < before_ev61_high_bound && pulses[iVec + 1].energy > 150 && pulses[iVec + 1].energy < 250) {cout << "\n" << "Before Ev61: Run Number, Event Number, triggerBits, Peak Time, Event 61 Time, dt, Energy = " << pulses[iVec + 1].file_num << "\t" << pulses[iVec + 1].entry << "\t" << pulses[iVec + 1].trigger << "\t" << eventstarttime3 << "\t" << pulses[iVec].ev61_start_time << "\t" << eventstarttime3 - pulses[iVec].ev61_start_time << "\t" << pulses[iVec + 1].energy << endl;}

                            if (eventstarttime3 - pulses[iVec].ev61_start_time > after_ev61_low_bound && eventstarttime3 - pulses[iVec].ev61_start_time < after_ev61_high_bound && pulses[iVec + 1].energy > 150 && pulses[iVec + 1].energy < 250) {cout << "\n" << "After Ev61: Run Number, Event Number, triggerBits, Peak Time, Event 61 Time, dt, Energy = " << pulses[iVec + 1].file_num << "\t" << pulses[iVec + 1].entry << "\t" << pulses[iVec + 1].trigger << "\t" << eventstarttime3 << "\t" << pulses[iVec].ev61_start_time << "\t" << eventstarttime3 - pulses[iVec].ev61_start_time << "\t" << pulses[iVec + 1].energy << endl;}

                            // if (pulses[iVec + 1].energy > 250 && pulses[iVec + 1].energy < 500) {cout << "\n" << "Run Number = " << pulses[iVec + 1].file_num << ", Event Number = " << pulses[iVec + 1].entry << ", triggerBits = " << pulses[iVec + 1].trigger << ", Peak Time = " << eventstarttime3 << ", Event 61 Time = " << pulses[iVec].ev61_start_time << ", dt = " << eventstarttime3 - pulses[iVec].ev61_start_time << ", detint = " << pulses[iVec + 1].energy << endl;}

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

        if (avg_pulse.file_num != last_file_num) {

            auto run_num_it = find(MWhrPerRun_Runs.begin(), MWhrPerRun_Runs.end(), avg_pulse.file_num);

            int run_num_index = distance(MWhrPerRun_Runs.begin(), run_num_it);

            cumulative_MWhrs += MWhrPerRun_MWhr[run_num_index];
            
            p_ks_test->AddPoint(cumulative_MWhrs, run_events);

            KSFile << avg_pulse.file_num << "\t" << MWhrPerRun_MWhr[run_num_index] << "\t" << run_events << "\n";

            // if (avg_pulse.file_num > 19260 && avg_pulse.file_num <= 19265) {cout << "\n" << "Run number, cumulative MWhrs = " << avg_pulse.file_num << " " << cumulative_MWhrs << endl;}

            // if (avg_pulse.file_num > 20435 && avg_pulse.file_num <= 20440) {cout << "\n" << "Run number, cumulative MWhrs = " << avg_pulse.file_num << " " << cumulative_MWhrs << endl;}

            run_events = 0;
        
        }

        last_file_num = avg_pulse.file_num;

    } // Event loop

    UBFile.close();

    UBFile.clear();

    double avg_detint = tot_detint * pow(hist_events, -1);

    // cout << "\n" << "Number of low dt events, Runs 7894-7918 = " << num_events_low_dt << endl;

    // cout << "\n" << "Number of central spike events, Runs 7894-7918 = " << num_events_low_dt_cent_spike << endl;

    // cout << "\n" << "Number of negative low dt events, excluding central spike, Runs 7894-7918 = " << num_events_low_dt_no_spike_neg << endl;

    // cout << "\n" << "Number of positive low dt events, excluding central spike, Runs 7894-7918 = " << num_events_low_dt_no_spike_pos << endl;

    // cout << "\n" << "Number of Event 61 events found, Runs 7894-7918 = " << take_sample_count << endl;

    // cout << "\n" << "Number of histogram events within 10 us of central sample, Runs 7894-7918 = " << hist_events << "\n" << endl;

    /* Fill & plot histograms */

    p_ks_test->SetTitle("KS Test: Cumulative MWhr per Beam Spill Event Rate");
    p_ks_test->GetXaxis()->SetTitle("Cumulative MWhr");
    p_ks_test->GetYaxis()->SetTitle("Cumulative Beam Spill Event Rate");
    p_ks_test->Draw("A*");
    // p_ks_test->GetXaxis()->SetLimits(2600, 3600);
    // p_ks_test->SetMinimum(50);
    // p_ks_test->SetMaximum(75);
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/p_ks_test_%i.png", run));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/p_ks_test_%i.png", run));

    TCanvas* c_dt_61 = new TCanvas("c_dt_61", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61->SetLogy();
    c_dt_61->cd();
    h_dt_61->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61->GetYaxis()->SetTitle("Counts");
    h_dt_61->Draw();
    c_dt_61->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_61_%i.png", run));
    // c_dt_61->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_dt_61.png");
    c_dt_61->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61.png"));
    h_dt_61->Reset();

    TCanvas* c_dt_61_1 = new TCanvas("c_dt_61_1", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_1->SetLogy();
    c_dt_61_1->cd();
    h_dt_61_1->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_1->GetYaxis()->SetTitle("Counts");
    h_dt_61_1->Draw();
    c_dt_61_1->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_61_1_%i.png", run));
    // c_dt_61_1->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_dt_61_1.png");
    c_dt_61_1->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_1.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_1.png"));
    h_dt_61_1->Reset();

    TCanvas* c_dt_61_2 = new TCanvas("c_dt_61_2", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_2->SetLogy();
    c_dt_61_2->cd();

    TF1* func3 = new TF1("func3", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", 200, 1000);

    double Par3[] = {10, 100, 625, 150};
    func3->SetParameters(Par3);
    func3->SetParName(0, "Baseline");
    func3->SetParName(1, "Amplitude");
    func3->SetParName(2, "Peak Center");
    func3->SetParName(3, "Standard Deviation");
    h_dt_61_2->Fit("func3", "R");

    h_dt_61_2->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_2->GetYaxis()->SetTitle("Counts");
    h_dt_61_2->Draw();

    Double_t param30 = func3->GetParameter(0);
    Double_t param31 = func3->GetParameter(1);
    Double_t param32 = func3->GetParameter(2);
    Double_t param33 = func3->GetParameter(3);

    Double_t error30 = func3->GetParError(0);
    Double_t error31 = func3->GetParError(1);
    Double_t error32 = func3->GetParError(2);
    Double_t error33 = func3->GetParError(3);

    TLatex *lat_dt_0 = new TLatex(0.2, 0.85, Form("Baseline = %.2f #pm %.2f", param30, error30));
    lat_dt_0->SetNDC();
    lat_dt_0->SetTextColor(1);
    lat_dt_0->SetTextSize(0.02);
    lat_dt_0->Draw("same");

    TLatex *lat_dt_1 = new TLatex(0.2, 0.8, Form("Amplitude = %.2f #pm %.2f", param31, error31));
    lat_dt_1->SetNDC();
    lat_dt_1->SetTextColor(1);
    lat_dt_1->SetTextSize(0.02);
    lat_dt_1->Draw("same");

    TLatex *lat_dt_2 = new TLatex(0.2, 0.75, Form("Peak Center = %.2f #pm %.2f", param32, error32));
    lat_dt_2->SetNDC();
    lat_dt_2->SetTextColor(1);
    lat_dt_2->SetTextSize(0.02);
    lat_dt_2->Draw("same");

    TLatex *lat_dt_3 = new TLatex(0.2, 0.7, Form("Standard Deviation = %.2f #pm %.2f", param33, error33));
    lat_dt_3->SetNDC();
    lat_dt_3->SetTextColor(1);
    lat_dt_3->SetTextSize(0.02);
    lat_dt_3->Draw("same");

    func3->Draw("same");
    func3->SetLineColor(2);
    c_dt_61_2->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_61_2_%i.png", run));
    // c_dt_61_2->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_dt_61_2.png");
    c_dt_61_2->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_2.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_2.png"));
    h_dt_61_2->Reset();

    TCanvas* c_dt_61_3 = new TCanvas("c_dt_61_3", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_3->SetLogy();
    c_dt_61_3->cd();
    h_dt_61_3->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_3->GetYaxis()->SetTitle("Counts");
    h_dt_61_3->Draw();
    c_dt_61_3->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_61_3_%i.png", run));
    // c_dt_61_3->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_dt_61_3.png");
    c_dt_61_3->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_3.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_3.png"));
    h_dt_61_3->Reset();

    TCanvas* c_dt_61_4 = new TCanvas("c_dt_61_4", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_4->SetLogy();
    c_dt_61_4->cd();
    h_dt_61_4->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_4->GetYaxis()->SetTitle("Counts");
    h_dt_61_4->Draw();
    c_dt_61_4->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_61_4_%i.png", run));
    // c_dt_61_4->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_dt_61_4.png");
    c_dt_61_4->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_4.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_4.png"));
    h_dt_61_4->Reset();

    TCanvas* c_detint = new TCanvas("c_detint", "Distribution of All Ev61Det Integral Values", 1200, 700);
    c_detint->SetLogy();
    c_detint->cd();
    h_detint->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_detint->GetYaxis()->SetTitle("Counts");
    // h_detint->SetMinimum(0);
    // h_detint->SetMaximum(5000);
    // h_detint->SetLineColor(kGreen);
    h_detint->Draw();
    TLatex *lat_detint = new TLatex(0.6, 0.7, Form("Average beam spill energy = %.2f", avg_detint));
    lat_detint->SetNDC();
    lat_detint->SetTextColor(1);
    lat_detint->SetTextSize(0.04);
    lat_detint->Draw();
    c_detint->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_detint_%i.png", run));
    // c_detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_detint.png");
    c_detint->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint.png"));
    h_detint->Reset();

    TCanvas* c_detint_1 = new TCanvas("c_detint_1", "Distribution of All Ev61Det Integral Values", 1200, 700);
    c_detint_1->SetLogy();
    c_detint_1->cd();
    h_detint_1->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_detint_1->GetYaxis()->SetTitle("Counts");
    h_detint_1->SetMinimum(1);
    h_detint_1->SetMaximum(50);
    h_detint_1->Draw("same");
    h_detint_1->SetLineColor(1);
    h_detint_2->Draw("same");
    h_detint_2->SetLineColor(2);
    h_detint_3->Draw("same");
    h_detint_3->SetLineColor(3);
    h_detint_4->Draw("same");
    h_detint_4->SetLineColor(4);
    h_detint_5->Draw("same");
    h_detint_5->SetLineColor(5);
    h_detint_6->Draw("same");
    h_detint_6->SetLineColor(6);
    TLegend *leg_detint_1 = new TLegend(0.9, 0.65, 0.65, 0.45);
    leg_detint_1->AddEntry(h_detint_1, "Selection Cut One", "l");
    leg_detint_1->AddEntry(h_detint_2, "Selection Cut Two", "l");
    leg_detint_1->AddEntry(h_detint_3, "Selection Cut Three", "l");
    leg_detint_1->AddEntry(h_detint_4, "Selection Cut Four", "l");
    leg_detint_1->AddEntry(h_detint_5, "Selection Cut Five", "l");
    leg_detint_1->AddEntry(h_detint_6, "Selection Cut Six", "l");
    leg_detint_1->Draw();
    c_detint_1->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_detint_SC_%i.png", run));
    // c_detint_1->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_detint_SC.png");
    c_detint_1->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_SC.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_SC.png"));
    h_detint_1->Reset();

    TCanvas* c_detint_R1 = new TCanvas("c_detint_R1", "Beam Spill Integral Values", 1200, 700);
    // c_detint_R1->SetLogy();
    c_detint_R1->cd();
    h_detint_R1->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_detint_R1->GetYaxis()->SetTitle("Counts");
    // h_detint_R1->SetMinimum(1);
    // h_detint_R1->SetMaximum(50);
    h_detint_R1->Draw();
    h_detint_R1->SetLineColor(1);
    c_detint_R1->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_detint_R1_%i.png", run));
    // c_detint_R1->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_detint_R1.png");
    c_detint_R1->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_R1.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_R1.png"));
    h_detint_R1->Reset();

    TCanvas* c_detint_R2 = new TCanvas("c_detint_R2", "Beam Spill Integral Values", 1200, 700);
    // c_detint_R2->SetLogy();
    c_detint_R2->cd();
    h_detint_R2->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_detint_R2->GetYaxis()->SetTitle("Counts");
    // h_detint_R2->SetMinimum(1);
    // h_detint_R2->SetMaximum(50);
    h_detint_R2->Draw();
    h_detint_R2->SetLineColor(1);
    c_detint_R2->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_detint_R2_%i.png", run));
    // c_detint_R2->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_detint_R2.png");
    c_detint_R2->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_R2.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_R2.png"));
    h_detint_R2->Reset();

    TCanvas* c_detint_R3 = new TCanvas("c_detint_R3", "Beam Spill Integral Values", 1200, 700);
    // c_detint_R3->SetLogy();
    c_detint_R3->cd();
    h_detint_R3->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_detint_R3->GetYaxis()->SetTitle("Counts");
    // h_detint_R3->SetMinimum(1);
    // h_detint_R3->SetMaximum(50);
    h_detint_R3->Draw();
    h_detint_R3->SetLineColor(1);
    c_detint_R3->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_detint_R3_%i.png", run));
    // c_detint_R3->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_detint_R3.png");
    c_detint_R3->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_R3.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_R3.png"));
    h_detint_R3->Reset();

    TCanvas* c_detint_R4 = new TCanvas("c_detint_R4", "Beam Spill Integral Values", 1200, 700);
    // c_detint_R4->SetLogy();
    c_detint_R4->cd();
    h_detint_R4->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_detint_R4->GetYaxis()->SetTitle("Counts");
    // h_detint_R4->SetMinimum(1);
    // h_detint_R4->SetMaximum(50);
    h_detint_R4->Draw();
    h_detint_R4->SetLineColor(1);
    c_detint_R4->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_detint_R4_%i.png", run));
    // c_detint_R4->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_detint_R4.png");
    c_detint_R4->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_R4.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_R4.png"));
    h_detint_R4->Reset();

    TCanvas* c_detint_R5 = new TCanvas("c_detint_R5", "Beam Spill Integral Values", 1200, 700);
    // c_detint_R5->SetLogy();
    c_detint_R5->cd();
    h_detint_R5->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_detint_R5->GetYaxis()->SetTitle("Counts");
    // h_detint_R5->SetMinimum(1);
    // h_detint_R5->SetMaximum(50);
    h_detint_R5->Draw();
    h_detint_R5->SetLineColor(1);
    c_detint_R5->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_detint_R5_%i.png", run));
    // c_detint_R5->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_detint_R5.png");
    c_detint_R5->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_R5.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_R5.png"));
    h_detint_R5->Reset();

    TCanvas* c_detint_R7 = new TCanvas("c_detint_R7", "Beam Spill Integral Values", 1200, 700);
    // c_detint_R7->SetLogy();
    c_detint_R7->cd();
    h_detint_R7->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_detint_R7->GetYaxis()->SetTitle("Counts");
    // h_detint_R7->SetMinimum(0);
    // h_detint_R7->SetMaximum(100);
    h_detint_R7->Draw("same");
    h_detint_R7->SetLineColor(2);
    h_detint_R6->Draw("same");
    h_detint_R6->SetLineColor(4);
    TLegend *leg_detint_R7 = new TLegend(0.9, 0.65, 0.55, 0.45);
    leg_detint_R7->AddEntry(h_detint_R6, "0.5us Before Ev61 (18.5us to 19us)", "l");
    leg_detint_R7->AddEntry(h_detint_R7, "0.5us After Ev61 (21us to 21.5us)", "l");
    leg_detint_R7->Draw();
    TLatex *lat_detint_R6 = new TLatex(0.3, 0.85, Form("Number of events before Ev61 = %.0f", h_detint_R6->Integral()));
    lat_detint_R6->SetNDC();
    lat_detint_R6->SetTextColor(1);
    lat_detint_R6->SetTextSize(0.04);
    lat_detint_R6->Draw("same");
    TLatex *lat_detint_R7 = new TLatex(0.3, 0.8, Form("Number of events after Ev61 = %.0f", h_detint_R7->Integral()));
    lat_detint_R7->SetNDC();
    lat_detint_R7->SetTextColor(1);
    lat_detint_R7->SetTextSize(0.04);
    lat_detint_R7->Draw("same");
    TLatex *lat_detint_R8 = new TLatex(0.3, 0.75, "All events within 0.5us time window.");
    lat_detint_R8->SetNDC();
    lat_detint_R8->SetTextColor(1);
    lat_detint_R8->SetTextSize(0.04);
    lat_detint_R8->Draw("same");
    c_detint_R7->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_detint_R67_%i.png", run));
    // c_detint_R7->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_detint_R7.png");
    c_detint_R7->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_R7.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_R7.png"));
    h_detint_R7->Reset();

    TCanvas* c_detint_low_dt = new TCanvas("c_detint_low_dt", "Distribution of Low dt Ev61Det Integral Values", 1200, 700);
    //c_detint_low_dt->SetLogy();
    c_detint_low_dt->cd();
    h_detint_low_dt->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_detint_low_dt->GetYaxis()->SetTitle("Counts");
    // h_detint_low_dt->SetMaximum(420);
    // h_detint_low_dt->SetLineColor(kRed);
    h_detint_low_dt->Draw();
    c_detint_low_dt->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_detint_low_dt_%i.png", run));
    // c_detint_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_detint_low_dt.png");
    c_detint_low_dt->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_low_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_low_dt.png"));
    h_detint_low_dt->Reset();

    TCanvas* c_dt_v_detint = new TCanvas("c_dt_v_detint", "Event 61 Detector dt vs Integral Value", 1200, 700);
    // c_dt_v_detint->SetLogy();
    c_dt_v_detint->cd();
    h_dt_v_detint->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_v_detint->GetYaxis()->SetTitle("Integral (Ph.e.)");
    // h_dt_v_detint->SetMarkerSize(2);
    // h_dt_v_detint->SetMarkerStyle(20);
    // gStyle->SetMarkerSize(3.0);
    h_dt_v_detint->Draw("colz");
    c_dt_v_detint->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_v_detint_%i.png", run));
    // c_dt_v_detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_dt_v_detint.png");
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
    c_int->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_int_%i.png", run));
    // c_int->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_int.png");
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
    c_eventint->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_eventint_%i.png", run));
    // c_eventint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_eventint.png");
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
    c_eventvpint->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_eventvpint_%i.png", run));
    // c_eventvpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_eventvpint.png");
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
    c_LLint->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_LLint_%i.png", run));
    // c_LLint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_LLint.png");
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
    c_HLint->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_HLint_%i.png", run));
    // c_HLint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_HLint.png");
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
    c_ev61int->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_ev61int_%i.png", run));
    // c_ev61int->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_ev61int.png");
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
    c_ev61vpint->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_ev61vpint_%i.png", run));
    // c_ev61vpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_ev61vpint.png");
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
    c_34detint->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_34detint_%i.png", run));
    // c_34detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_34detint.png");
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
    c_34vpint->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_34vpint_%i.png", run));
    // c_34vpint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_34vpint.png");
    c_34vpint->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_34vpint.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_34vpint.png"));
    h_34vpint->Reset();

    TCanvas* c_int_cent_spike = new TCanvas("c_int_cent_spike", "Distribution of All Central Spike Integral Values", 1200, 700);
    //c_int_cent_spike->SetLogy();
    c_int_cent_spike->cd();
    h_int_cent_spike->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_int_cent_spike->GetYaxis()->SetTitle("Counts");
    h_int_cent_spike->Draw();
    c_int_cent_spike->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_int_cent_spike_%i.png", run));
    // c_int_cent_spike->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_int_cent_spike.png");
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
    c_int_neg_low_dt->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_int_low_dt_%i.png", run));
    // c_int_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_int_low_dt.png");
    c_int_neg_low_dt->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_neg_low_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_neg_low_dt.png"));
    h_int_pos_low_dt->Reset();
    h_int_neg_low_dt->Reset();

    TCanvas* c_svpint_neg_low_dt = new TCanvas("c_svpint_neg_low_dt", "Distribution of All Low dt Side Veto Panel Integral Values", 1200, 700);
    c_svpint_neg_low_dt->SetLogy();
    c_svpint_neg_low_dt->cd();
    h_svpint_neg_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
    h_svpint_neg_low_dt->GetYaxis()->SetTitle("Counts");
    h_svpint_neg_low_dt->Draw("same");
    h_svpint_neg_low_dt->SetLineColor(kRed);
    h_svpint_pos_low_dt->Draw("same");
    h_svpint_pos_low_dt->SetLineColor(kGreen);
    h_svpint_cent_spike->Draw("same");
    h_svpint_cent_spike->SetLineColor(kBlue);
    TLegend *leg_svpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
    leg_svpint_neg_low_dt->AddEntry(h_svpint_cent_spike, "Central Spike Events", "l");
    leg_svpint_neg_low_dt->AddEntry(h_svpint_pos_low_dt, "Positive Low dt Events", "l");
    leg_svpint_neg_low_dt->AddEntry(h_svpint_neg_low_dt, "Negative Low dt Events", "l");
    leg_svpint_neg_low_dt->Draw();
    c_svpint_neg_low_dt->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_svpint_low_dt_%i.png", run));
    // c_svpint_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_svpint_low_dt.png");
    c_svpint_neg_low_dt->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_svpint_neg_low_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_svpint_neg_low_dt.png"));
    h_svpint_pos_low_dt->Reset();
    h_svpint_neg_low_dt->Reset();

    TCanvas* c_tvpint_neg_low_dt = new TCanvas("c_tvpint_neg_low_dt", "Distribution of All Low dt Top Veto Panel Integral Values", 1200, 700);
    c_tvpint_neg_low_dt->SetLogy();
    c_tvpint_neg_low_dt->cd();
    h_tvpint_neg_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
    h_tvpint_neg_low_dt->GetYaxis()->SetTitle("Counts");
    h_tvpint_neg_low_dt->Draw("same");
    h_tvpint_neg_low_dt->SetLineColor(kRed);
    h_tvpint_pos_low_dt->Draw("same");
    h_tvpint_pos_low_dt->SetLineColor(kGreen);
    h_tvpint_cent_spike->Draw("same");
    h_tvpint_cent_spike->SetLineColor(kBlue);
    TLegend *leg_tvpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
    leg_tvpint_neg_low_dt->AddEntry(h_tvpint_cent_spike, "Central Spike Events", "l");
    leg_tvpint_neg_low_dt->AddEntry(h_tvpint_pos_low_dt, "Positive Low dt Events", "l");
    leg_tvpint_neg_low_dt->AddEntry(h_tvpint_neg_low_dt, "Negative Low dt Events", "l");
    leg_tvpint_neg_low_dt->Draw();
    c_tvpint_neg_low_dt->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_tvpint_low_dt_%i.png", run));
    // c_tvpint_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_tvpint_low_dt.png");
    c_tvpint_neg_low_dt->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_tvpint_neg_low_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_tvpint_neg_low_dt.png"));
    h_tvpint_pos_low_dt->Reset();
    h_tvpint_neg_low_dt->Reset();

    TCanvas* c_avpint_neg_low_dt = new TCanvas("c_avpint_neg_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
    // c_avpint_neg_low_dt->SetLogy();
    c_avpint_neg_low_dt->cd();
    TF1* func1 = new TF1("func1", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -200, 200);
    func1->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
    func1->SetParameters(1, 10, 0, 50);
    h_avpint_pos_low_dt->Fit("func1", "R");
    h_avpint_neg_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
    h_avpint_neg_low_dt->GetYaxis()->SetTitle("Counts");
    // h_avpint_neg_low_dt->GetYaxis()->SetRange(0, 10);
    h_avpint_neg_low_dt->Draw();
    h_avpint_neg_low_dt->SetLineColor(kRed);
    h_avpint_pos_low_dt->Draw("same");
    h_avpint_pos_low_dt->SetLineColor(kGreen);
    h_avpint_cent_spike->Draw("same");
    h_avpint_cent_spike->SetLineColor(kBlue);
    func1->Draw("same");
    func1->SetLineColor(1);
    TLegend *leg_avpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
    leg_avpint_neg_low_dt->AddEntry(h_avpint_cent_spike, "Central Spike Events", "l");
    leg_avpint_neg_low_dt->AddEntry(h_avpint_pos_low_dt, "Positive Low dt Events", "l");
    leg_avpint_neg_low_dt->AddEntry(h_avpint_neg_low_dt, "Negative Low dt Events", "l");
    leg_avpint_neg_low_dt->Draw();
    c_avpint_neg_low_dt->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_avpint_low_dt_%i.png", run));
    // c_avpint_neg_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_avpint_low_dt.png");
    c_avpint_neg_low_dt->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_neg_low_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_neg_low_dt.png"));
    h_avpint_pos_low_dt->Reset();
    h_avpint_neg_low_dt->Reset();

    TCanvas* c_avpint_low_dt = new TCanvas("c_avpint_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
    // c_avpint_low_dt->SetLogy();
    c_avpint_low_dt->cd();
    TF1* func2 = new TF1("func2", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -300, 300);
    func2->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
    func2->SetParameters(1, 10, 0, 50);
    h_avpint_low_dt->Fit("func2", "R");
    h_avpint_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
    h_avpint_low_dt->GetYaxis()->SetTitle("Counts");
    // h_avpint_low_dt->GetYaxis()->SetRange(0, 50);
    h_avpint_low_dt->Draw();
    func2->Draw("same");
    func2->SetLineColor(2);
    c_avpint_low_dt->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_avpint_low_dt_1c_%i.png", run));
    // c_avpint_low_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_avpint_low_dt_1c.png");
    c_avpint_low_dt->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_low_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_low_dt.png"));
    h_avpint_low_dt->Reset();

    TCanvas* c_runnum = new TCanvas("c_runnum", "Beam Spill Run Numbers", 1200, 700);
    // c_runnum->SetLogy();
    c_runnum->cd();
    h_runnum->GetXaxis()->SetTitle("Run Number");
    h_runnum->GetYaxis()->SetTitle("Counts");
    // h_runnum->SetMinimum(1);
    // h_runnum->SetMaximum(50);
    h_runnum->Draw();
    h_runnum->SetLineColor(1);
    c_runnum->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_runnum_%i.png", run));
    // c_runnum->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_runnum.png");
    c_runnum->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_runnum.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_runnum.png"));
    h_runnum->Reset();

    TCanvas* c_runnum_v_dt = new TCanvas("c_runnum_v_dt", "Beam Spill Run Number vs dt", 1200, 700);
    // c_runnum_v_dt->SetLogy();
    c_runnum_v_dt->cd();
    h_runnum_v_dt->GetXaxis()->SetTitle("Run Number");
    h_runnum_v_dt->GetYaxis()->SetTitle("Delta-T (ns)");
    // h_runnum_v_dt->SetMinimum(1);
    // h_runnum_v_dt->SetMaximum(50);
    h_runnum_v_dt->Draw();
    // h_runnum_v_dt->SetLineColor(1);
    c_runnum_v_dt->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_runnum_v_dt_%i.png", run));
    // c_runnum_v_dt->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_runnum_v_dt.png");
    c_runnum_v_dt->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_runnum_v_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_runnum_v_dt.png"));
    h_runnum_v_dt->Reset();

    TCanvas* c_runnum_v_detint = new TCanvas("c_runnum_v_detint", "Beam Spill Run Number vs Integral", 1200, 700);
    // c_runnum_v_detint->SetLogy();
    c_runnum_v_detint->cd();
    h_runnum_v_detint->GetXaxis()->SetTitle("Run Number");
    h_runnum_v_detint->GetYaxis()->SetTitle("Integral (Ph.e.)");
    // h_runnum_v_detint->SetMinimum(1);
    // h_runnum_v_detint->SetMaximum(50);
    h_runnum_v_detint->Draw();
    // h_runnum_v_detint->SetLineColor(1);
    c_runnum_v_detint->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_runnum_v_detint_%i.png", run));
    // c_runnum_v_detint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_runnum_v_detint.png");
    c_runnum_v_detint->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_runnum_v_detint.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_runnum_v_detint.png"));
    h_runnum_v_detint->Reset();

    TCanvas* c_dt_v_tmuon = new TCanvas("c_dt_v_tmuon", "Event 61 Detector dt vs Muon dt", 1200, 700);
    // c_dt_v_tmuon->SetLogy();
    c_dt_v_tmuon->cd();
    h_dt_v_tmuon->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_v_tmuon->GetYaxis()->SetTitle("Time Since Last Muon (ns)");
    h_dt_v_tmuon->SetMarkerStyle(7);
    h_dt_v_tmuon->Draw();
    c_dt_v_tmuon->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_v_tmuon_%i.png", run));
    // c_dt_v_tmuon->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_dt_v_tmuon.png");
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
    c_vec_size->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_vec_size_%i.png", run));
    // c_vec_size->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_vec_size.png");
    c_vec_size->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_vec_size.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_vec_size.png"));
    h_vec_size->Reset();

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

    cout << "\n" << "End of code." << "\n" << endl;

    return 0;

}