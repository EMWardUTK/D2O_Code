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
#include "TGraphErrors.h"
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
  long long unix_time;      /* Unix starttime value from inputFile TTree T */
  long long start;          /* Universal start time of pulse (10% peak) in waveform (ns) */
  long long end;            /* Universal end time of pulse (reach baseline) in waveform (ns) */
  double peak;              /* Max amplitude of pulse (photo-electrons) */
  double energy;            /* Energy (integral) of pulse (photo-electrons) */
  int number;               /* Number of channels in which we see pulse (photo-electrons) */
  bool single;              /* Is pulse (photo-electrons) timing consistent across all channels */
  bool beam;                /* Tracks whether beam is on or off */
  bool recent_ev61;         /* Tacks whether most recent Event 61 event had a large waveform peak or not */
  int trigger;              /* Tracks whether trigger is external (2) or internal (16) */
  int length;               /* Length of waveform in number of bins */
  double side_vp_energy;    /* Energy (integral) of pulse (photo-electrons) in SIDE veto panels */
  double top_vp_energy;     /* Energy (integral) of pulse (photo-electrons) in TOP veto panel */
  double all_vp_energy;     /* Energy (integral) of pulse (photo-electrons) in ALL veto panels */
  long long last_cosmics_time; /* Time value of most recent detected cosmic particle event */
  long long last_det_time;  /* Time value of most recent detected detector event */
  int selection_cut;        /* Tacks which selection cut(s) a given event satisfies */
  bool issue;               /* Flag to keep track of unusual eventTree entries */
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

/** @brief These values should change for each PMT */
const char* outputFilePath = "";

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

    std::vector<int> bg_samples(100);

    std::vector<struct pulse> pulses;

    int selection_cut_level = last_run;

    long long sample_range = 2 * pow(10, 4);

    long long sample_rate = pow(10, 5);

    long long LL_time = 0.0;

    long long HL_time = 0.0;

    long long iTime = 0.0;

    long long hour12_st = 0.0;

    long long last_utime = 0.0;

    double plot_min = 1000;          // -19960;

    double plot_max = 9960;          // -11000;

    double num_bins = 56;

    double sample_bins = 250;

    bool bt_LL_n_HL = false;

    int num_events_low_dt = 0;

    int num_events_low_dt_cent_spike = 0;

    int num_events_low_dt_no_spike_neg = 0;

    int num_events_low_dt_no_spike_pos = 0;

    int run_counter = 0;

    int take_sample_count = 0;

    int sample_window_count = 0;

    int sample_window_events = 0;

    int hist_events = 0;

    double tot_detint = 0;

    bool fall_is_true_spring_is_false = false;

    // if (run < 66) {fall_is_true_spring_is_false = true;}

    // else if (run >= 66) {fall_is_true_spring_is_false = false;}

    if (fall_is_true_spring_is_false) {

        sample_range = 2 * pow(10, 4);

        sample_bins = 250;

        plot_min = 1000;

        plot_max = 9960;

    }

    else if (!fall_is_true_spring_is_false) {

        sample_range = 3 * pow(10, 4);

        sample_bins = 375;

        plot_min = -19960;

        plot_max = -11000;

    }

    // cout << "\n" << "Number of selection cuts: " << selection_cut_level << "\n" << endl;

    /* Initialize histograms */                 // NUMBER OF BINS SHOULD BE A MULTIPLE OF THE SMALLEST AXIS UNIT!!! (OR LENGTH OF AXIS DIVIDED BY SMALLEST AXIS UNIT?)

    TH1D* h_dt_61 = new TH1D("h_dt_61", "Delta-T Between Sample Center and Background Peaks", num_bins, plot_min, plot_max);                                             // 1 bin = 160 ns

    TH1D* h_dt_61_1 = new TH1D("h_dt_61_1", "Delta-T Between Sample Center and Background Peaks", 625, - 1 * pow(10, 4), 1 * pow(10, 4));                                // 1 bin = 32 ns

    TH1D* h_dt_61_2 = new TH1D("h_dt_61_2", "Delta-T Between Sample Center and Background Peaks", 1250, - 2 * pow(10, 4), 2 * pow(10, 4));                               // 1 bin = 32 ns

    TH1D* h_dt_61_3 = new TH1D("h_dt_61_3", "Delta-T Between Sample Center and Background Peaks", 3125, - 5 * pow(10, 4), 5 * pow(10, 4));                               // 1 bin = 32 ns

    TH1D* h_dt_61_4 = new TH1D("h_dt_61_4", "Delta-T Between Sample Center and Background Peaks", 1250, - 10 * pow(10, 4), 10 * pow(10, 4));                             // 1 bin = 160 ns

    TH2D* h_dt_v_detint = new TH2D("h_dt_v_detint", "Delta-T vs Detector Integral Value", sample_bins, - sample_range, sample_range, 100, 0, 2000);                      // 1 bin = 160 ns

    TH2D* h_dt_v_tmuon = new TH2D("h_dt_v_tmuon", "Delta-T vs Most Recent Muon Time Difference", 1250, - 2 * sample_range, 2 * sample_range, 1000, 0, 2 * pow(10, 7));   // 1 bin = 32 ns

    TH2D* h_dt_v_tdet = new TH2D("h_dt_v_tdet", "Delta-T vs Most Recent Michel Time Difference", 1250, - 2 * sample_range, 2 * sample_range, 1000, 0, 2 * pow(10, 7));   // 1 bin = 32 ns

    TH1D* h_int = new TH1D("h_int", "Distribution of All Event Integral Values Between High Light and Low Light LED Events", 200, 0, 2000);

    TH1D* h_LLint = new TH1D("h_LLint", "Distribution of All Low Light LED Integral Values", 200, 0, 2000);

    TH1D* h_HLint = new TH1D("h_HLint", "Distribution of All High Light LED Integral Values", 200, 0, 2000);

    TH1D* h_eventint = new TH1D("h_eventint", "Distribution of All Internally Triggered Detector Integral Values", 200, 0, 2000);

    TH1D* h_eventvpint = new TH1D("h_eventvpint", "Distribution of All Internally Triggered Veto Panel Integral Values", 200, 0, 2000);

    TH1D* h_34detint = new TH1D("h_34detint", "Distribution of All triggerBits = 34 Detector Integral Values", 200, 0, 2000);

    TH1D* h_34vpint = new TH1D("h_34vpint", "Distribution of All triggerBits = 34 Veto Panel Integral Values", 200, 0, 2000);

    TH1D* h_detint = new TH1D("h_detint", "Background PMT Integral Values", 50, 0, 500);

    TH1D* h_detint_1 = new TH1D("h_detint_1", "Background PMT Integral Values, All Selection Cuts", 80, 0, 800);

    TH1D* h_detint_2 = new TH1D("h_detint_2", "Background PMT Integral Values, All Selection Cuts", 80, 0, 800);

    TH1D* h_detint_3 = new TH1D("h_detint_3", "Background PMT Integral Values, All Selection Cuts", 80, 0, 800);

    TH1D* h_detint_4 = new TH1D("h_detint_4", "Background PMT Integral Values, All Selection Cuts", 80, 0, 800);

    TH1D* h_detint_5 = new TH1D("h_detint_5", "Background PMT Integral Values, All Selection Cuts", 80, 0, 800);

    TH1D* h_detint_6 = new TH1D("h_detint_6", "Background PMT Integral Values, All Selection Cuts", 80, 0, 800);

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

    TH1D* h_vec_size = new TH1D("h_vec_size", "Distribution of Number of Internally Triggered Events Between Low Light and High Light LED Events", 5, 0, 5);

    // TH1D* h_oldbkgdint = new TH1D("h_oldbkgdint", "Distribution of All Background Detector Integral Values", 400, 0, 4000);

    TH1D* h_slice_vpintavg = new TH1D("h_slice_vpintavg", "Distribution of Averaged Veto Panel Integral Values from Background Between 200-400 Ph.e.", 100, 0, 100);

    TH1D* h_vpint_16 = new TH1D("h_vpint_16", "Background SiPM 1 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_17 = new TH1D("h_vpint_17", "Background SiPM 2 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_18 = new TH1D("h_vpint_18", "Background SiPM 3 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_19 = new TH1D("h_vpint_19", "Background SiPM 4 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_20 = new TH1D("h_vpint_20", "Background SiPM 5 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_21 = new TH1D("h_vpint_21", "Background SiPM 6 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_22 = new TH1D("h_vpint_22", "Background SiPM 7 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_23 = new TH1D("h_vpint_23", "Background SiPM 8 Integral Values", 350, -500, 3000);

    TH1D* h_vpint_top = new TH1D("h_vpint_top", "Background SiPM 9 & 10 Integral Values", 350, -500, 3000);

    // TH2D* h_sample_events = new TH2D("h_sample_events", "Background Statistics as a Function of Distance from Event 61", 100, 0, 100, 5, 0, 5);

    /*std::vector<double> BackgroundIntVal;

    // Create a text string, which is used to output the text file
    string int_val;

    // Read from the text file
    // ifstream ReadBIFile(Form("Background_Integrals_Run%i.txt", run));
    ifstream ReadBIFile("Background_Integrals_Run7894.txt");

    // Use a while loop together with the getline() function to read the file line by line
    while (getline(ReadBIFile, int_val)) {

        // Output the text from the file
        BackgroundIntVal.push_back(stod(int_val));

    }

    // Close the file
    ReadBIFile.close();

    for (size_t iVec = 0; iVec < BackgroundIntVal.size(); iVec++) {h_oldbkgdint->Fill(BackgroundIntVal[iVec]);}

    BackgroundIntVal.clear();*/
    
    TFile *f; 
    //your root file location here
        if(gSystem->AccessPathName(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/int_trig_events/Event_Analysis_AllGR_%i.root", run))){           // Event_Analysis_24hrs_Run         // rootfiles_vprm/EventAnalysis_VP16_Run
        //    if (gSystem->AccessPathName(Form("/data9/coherent/data/d2o/processedData/Event_Analysis_24hrs_Run%i.root", run))) {
        //    if(gSystem->AccessPathName(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/Event_Analysis_24hrs_Run%i.root", run))){
        cout << "Could not open file" << endl;
        return -1; 
    } else{
        f = new TFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/int_trig_events/Event_Analysis_AllGR_%i.root", run));         // Event_Analysis_24hrs_Run         // rootfiles_vprm/EventAnalysis_VP16_Run
        //        f = new TFile(Form("/data9/coherent/data/d2o/processedData/Event_Analysis_24hrs_Run%i.root", run));
        //        f = new TFile(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/Event_Analysis_24hrs_Run%i.root", run));
    }

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
    Bool_t          re;
    Int_t           tr;
    Int_t           lh;
    Double_t        svpe;
    Double_t        tvpe;
    Double_t        avpe;
    Long64_t        lct;
    Long64_t        ldt;
    Int_t           selcut;
    Bool_t          issue;
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
    TBranch        *b_re;
    TBranch        *b_tr;
    TBranch        *b_lh;
    TBranch        *b_svpe;
    TBranch        *b_tvpe;
    TBranch        *b_avpe;
    TBranch        *b_lct;
    TBranch        *b_ldt;
    TBranch        *b_selcut;
    TBranch        *b_issue;
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
    t->SetBranchAddress("re", &re, &b_re);
    t->SetBranchAddress("tr", &tr, &b_tr);
    t->SetBranchAddress("lh", &lh, &b_lh);
    t->SetBranchAddress("svpe", &svpe, &b_svpe);
    t->SetBranchAddress("tvpe", &tvpe, &b_tvpe);
    t->SetBranchAddress("avpe", &avpe, &b_avpe);
    t->SetBranchAddress("lct", &lct, &b_lct);
    t->SetBranchAddress("ldt", &ldt, &b_ldt);
    t->SetBranchAddress("selcut", &selcut, &b_selcut);
    t->SetBranchAddress("issue", &issue, &b_issue);
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

    // Create output root file of processed waveform data
    TString outputFile;
    if (strcmp(outputFilePath, "") == 0) {
        outputFile.Form("%s.root", Form("BackgroundPlots%iGoldenRuns", run));
    } else {
        outputFile.Form("%s/%s.root", outputFilePath, Form("BackgroundPlots%iGoldenRuns", run));
    }
    TFile *fileOut = new TFile(outputFile, "RECREATE");
    TTree *backgroundTree = new TTree("backgroundTree", "backgroundTree");
    backgroundTree->SetDirectory(fileOut);

    // Data to keep track of for each backgroundTree entry
    Long64_t *br_bg;
    backgroundTree->Branch("bg", &br_bg, "bg/L");

    // Create and open a text file
    ofstream BIFile(Form("Background_Integrals_Run%i.txt", run));
    ofstream BERFile(Form("Background_Event_Rate_%i.txt", run));
    ofstream BITFile(Form("Background_Int_dt_%i.txt", run));

    for(int iEvent = 0; iEvent < t->GetEntries(); iEvent++){

        Long64_t tentry = t->LoadTree(iEvent);

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
        b_re->GetEntry(tentry);
        b_tr->GetEntry(tentry);
        b_lh->GetEntry(tentry);
        b_svpe->GetEntry(tentry);
        b_tvpe->GetEntry(tentry);
        b_avpe->GetEntry(tentry);
        b_lct->GetEntry(tentry);
        b_ldt->GetEntry(tentry);
        b_selcut->GetEntry(tentry);
        b_issue->GetEntry(tentry);
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
        avg_pulse.recent_ev61 = re;
        avg_pulse.trigger = tr;
        avg_pulse.length = lh;
        avg_pulse.side_vp_energy = svpe;
        avg_pulse.top_vp_energy = tvpe;
        avg_pulse.all_vp_energy = avpe;
        avg_pulse.last_cosmics_time = lct;
        avg_pulse.last_det_time = ldt;
        avg_pulse.selection_cut = selcut;
        avg_pulse.issue = issue;
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

        /*
        
        Only looking for events between low light LED events and high light LED events
        
        */

        // After every low light LED event, "turn on" internally triggered event storing

        if (avg_pulse.trigger == 16 && avg_pulse.recent_ev61) {bt_LL_n_HL = true; LL_time = ceil(avg_pulse.start); h_LLint->Fill(avg_pulse.energy);}

        // Store internally triggered events in vector

        if (bt_LL_n_HL && avg_pulse.selection_cut >= selection_cut_level && avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {pulses.push_back(avg_pulse); h_eventint->Fill(avg_pulse.energy); h_eventvpint->Fill(avg_pulse.all_vp_energy);}

        // After every high light LED event, "turn off" internally triggered event storing and search through compiled vector of events

        if (avg_pulse.trigger == 8) {

            HL_time = floor(avg_pulse.start);

            sample_rate = round((HL_time - LL_time) / 101);

            // Scan through time in 100 us intervals, from time of last low light event to time of current high light event

            iTime = LL_time;

            sample_window_count = 0;

            while (sample_window_count < 100) {
                
                iTime += sample_rate;

                take_sample_count += 1;

                sample_window_events = 0;

                // Iterate over all events in compiled vector

                for (size_t iVec = 0; iVec < pulses.size(); iVec++) {

                    // Check if each event is within time range of current time counter value

                    if (pulses[iVec].start >= iTime - sample_range && pulses[iVec].start <= iTime + sample_range) {

                        sample_window_events += 1;

                        // Plot histograms

                        if ((pulses[iVec].start - iTime) >= -10000 && (pulses[iVec].start - iTime) <= 10000) {num_events_low_dt += 1; h_detint_low_dt->Fill(pulses[iVec].energy); h_avpint_low_dt->Fill(pulses[iVec].all_vp_energy);}
            
                        if ((pulses[iVec].start - iTime) > -500 && (pulses[iVec].start - iTime) < 500) {num_events_low_dt_cent_spike += 1; h_int_cent_spike->Fill(pulses[iVec].energy); h_avpint_cent_spike->Fill(pulses[iVec].all_vp_energy); h_svpint_cent_spike->Fill(pulses[iVec].side_vp_energy); h_tvpint_cent_spike->Fill(pulses[iVec].top_vp_energy);}

                        else if ((pulses[iVec].start - iTime) <= -500 && (pulses[iVec].start - iTime) >= -10000) {num_events_low_dt_no_spike_neg += 1; h_int_neg_low_dt->Fill(pulses[iVec].energy); h_avpint_neg_low_dt->Fill(pulses[iVec].all_vp_energy); h_svpint_neg_low_dt->Fill(pulses[iVec].side_vp_energy); h_tvpint_neg_low_dt->Fill(pulses[iVec].top_vp_energy);}

                        else if ((pulses[iVec].start - iTime) >= 500 && (pulses[iVec].start - iTime) <= 10000) {num_events_low_dt_no_spike_pos += 1; h_int_pos_low_dt->Fill(pulses[iVec].energy); h_avpint_pos_low_dt->Fill(pulses[iVec].all_vp_energy); h_svpint_pos_low_dt->Fill(pulses[iVec].side_vp_energy); h_tvpint_pos_low_dt->Fill(pulses[iVec].top_vp_energy);}

                        if (pulses[iVec].energy > 200 && pulses[iVec].energy < 400) {h_slice_vpintavg->Fill(avg_pulse.all_vp_energy / 9);}

                        if (pulses[iVec].start - iTime >= plot_min && pulses[iVec].start - iTime <= plot_max) {

                            hist_events += 1;

                            tot_detint += pulses[iVec].energy;

                            h_dt_61->Fill(pulses[iVec].start - iTime);

                            // h_detint->Fill(pulses[iVec].energy);

                            if (pulses[iVec].selection_cut >= 1) {h_detint_1->Fill(pulses[iVec].energy);}

                            if (pulses[iVec].selection_cut >= 2) {h_detint_2->Fill(pulses[iVec].energy);}

                            if (pulses[iVec].selection_cut >= 3) {h_detint_3->Fill(pulses[iVec].energy);}

                            if (pulses[iVec].selection_cut >= 4) {h_detint_4->Fill(pulses[iVec].energy);}

                            if (pulses[iVec].selection_cut >= 5) {h_detint_5->Fill(pulses[iVec].energy);}

                            if (pulses[iVec].selection_cut >= 6) {h_detint_6->Fill(pulses[iVec].energy);}

                            Long64_t bg = pulses[iVec].energy;

                            backgroundTree->SetBranchAddress("bg", &bg);

                            backgroundTree->Fill();
                            
                        }
                        
                        h_dt_61_1->Fill(pulses[iVec].start - iTime);

                        h_dt_61_2->Fill(pulses[iVec].start - iTime);

                        h_dt_61_3->Fill(pulses[iVec].start - iTime);

                        h_dt_61_4->Fill(pulses[iVec].start - iTime);

                        h_dt_v_tmuon->Fill(pulses[iVec].start - iTime, pulses[iVec].start - pulses[iVec].last_cosmics_time);

                        h_dt_v_tdet->Fill(pulses[iVec].start - iTime, pulses[iVec].start - pulses[iVec].last_det_time);

                        h_detint->Fill(pulses[iVec].energy);

                        h_dt_v_detint->Fill(pulses[iVec].start - iTime, pulses[iVec].energy);

                        BITFile << pulses[iVec].start - iTime << "\t" << pulses[iVec].energy << "\n";

                        h_vpint_16->Fill(pulses[iVec].vp_int_16); // if (pulses[iVec].vp_int_16 > 150) {h_vpint_16->Fill(pulses[iVec].vp_int_16);}          // pulses[iVec].vp_start_16 != 0

                        h_vpint_17->Fill(pulses[iVec].vp_int_17); // if (pulses[iVec].vp_int_17 > 150) {h_vpint_17->Fill(pulses[iVec].vp_int_17);}

                        h_vpint_18->Fill(pulses[iVec].vp_int_18); // if (pulses[iVec].vp_int_18 > 400) {h_vpint_18->Fill(pulses[iVec].vp_int_18);}

                        h_vpint_19->Fill(pulses[iVec].vp_int_19); // if (pulses[iVec].vp_int_19 > 150) {h_vpint_19->Fill(pulses[iVec].vp_int_19);}

                        h_vpint_20->Fill(pulses[iVec].vp_int_20); // if (pulses[iVec].vp_int_20 > 150) {h_vpint_20->Fill(pulses[iVec].vp_int_20);}

                        h_vpint_21->Fill(pulses[iVec].vp_int_21); // if (pulses[iVec].vp_int_21 > 150) {h_vpint_21->Fill(pulses[iVec].vp_int_21);}

                        h_vpint_22->Fill(pulses[iVec].vp_int_22); // if (pulses[iVec].vp_int_22 > 150) {h_vpint_22->Fill(pulses[iVec].vp_int_22);}

                        h_vpint_23->Fill(pulses[iVec].vp_int_23); // if (pulses[iVec].vp_int_23 > 150) {h_vpint_23->Fill(pulses[iVec].vp_int_23);}

                        h_vpint_top->Fill(pulses[iVec].vp_int_top); // if (pulses[iVec].vp_int_top > 250) {h_vpint_top->Fill(pulses[iVec].vp_int_top);}

                        // cout << "\n" << "Run Number = " << run << ", Event Number = " << pulses[iVec].entry + 1 << ", triggerBits = " << pulses[iVec].trigger << ", Integral = " << pulses[iVec].energy << endl;

                    }

                }

                // cout << "\n" << "sample_window_count = " << sample_window_count << ", sample_window_events = " << sample_window_events << endl;

                // h_sample_events->Fill(sample_window_count, sample_window_events);

                bg_samples[sample_window_count] += sample_window_events;

                sample_window_count += 1;

            }

            // Reset search parameter after finding high light LED event
            
            bt_LL_n_HL = false; h_HLint->Fill(avg_pulse.energy); h_vec_size->Fill(pulses.size()); pulses.clear();
            
        }

        if (avg_pulse.trigger == 34) {h_34detint->Fill(avg_pulse.energy); h_34vpint->Fill(avg_pulse.all_vp_energy);}

        h_int->Fill(avg_pulse.energy);

        if (avg_pulse.unix_time != last_utime) {run_counter += 1;}

        if (run_counter == 12) {hour12_st = avg_pulse.unix_time;}

        last_utime = avg_pulse.unix_time;

    } // Event loop

    double avg_detint = tot_detint * pow(hist_events, -1);

    // Write to the file

    BIFile << hour12_st << "\t" << avg_detint << "\n";

    BERFile << hour12_st << "\t" << hist_events << "\n";

    // Close the file

    BIFile.close();

    BIFile.clear();

    BERFile.close();

    BERFile.clear();

    BITFile.close();

    BITFile.clear();

    // h_detint->Add(h_oldbkgdint, -1);

    // h_detint->Scale(0.01);

    backgroundTree->Write();

    cout << "\n" << "Number of low dt events, Runs " << run << " = " << num_events_low_dt << endl;

    cout << "\n" << "Number of central spike events, Runs " << run << " = " << num_events_low_dt_cent_spike << endl;

    cout << "\n" << "Number of negative low dt events, excluding central spike, Runs " << run << " = " << num_events_low_dt_no_spike_neg << endl;

    cout << "\n" << "Number of positive low dt events, excluding central spike, Runs " << run << " = " << num_events_low_dt_no_spike_pos << endl;

    cout << "\n" << "Number of central samples taken, Runs " << run << " = " << take_sample_count << endl;

    cout << "\n" << "Number of histogram events within 10 us of central sample, Runs " << run << " = " << hist_events << endl;

    cout << "\n" << "Number of runs in ROOT file, Runs " << run << " = " << run_counter << "\n" << endl;

    /* Fill & plot histograms */

    TCanvas* c_dt_61 = new TCanvas("c_dt_61", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61->SetLogy();
    c_dt_61->cd();
    h_dt_61->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61->GetYaxis()->SetTitle("Counts");
    // h_dt_61->SetMaximum(18);
    // h_dt_61->SetLineColor(kGreen);
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
    h_dt_61_2->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_2->GetYaxis()->SetTitle("Counts");
    h_dt_61_2->Draw();
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
    /*
    TCanvas* c_bkgdint = new TCanvas("c_bkgdint", "Distribution of All Background Integral Values", 1200, 700);
    //c_bkgdint->SetLogy();
    c_bkgdint->cd();
    h_oldbkgdint->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_oldbkgdint->GetYaxis()->SetTitle("Counts");
    h_oldbkgdint->Draw();
    c_bkgdint->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_oldbkgdint_%i.png", run));
    // c_bkgdint->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_oldbkgdint.png");
    c_bkgdint->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_oldbkgdint.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_oldbkgdint.png"));
    h_oldbkgdint->Reset();
    */
    TCanvas* c_slice_vpintavg = new TCanvas("c_slice_vpintavg", "Distribution of Averaged Veto Panel Integral Values", 1200, 700);
    c_slice_vpintavg->SetLogy();
    c_slice_vpintavg->cd();
    h_slice_vpintavg->GetXaxis()->SetTitle("Integral (ADC)");
    h_slice_vpintavg->GetYaxis()->SetTitle("Counts");
    h_slice_vpintavg->Draw();
    c_slice_vpintavg->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_slice_vpintavg_%i.png", run));
    // c_slice_vpintavg->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_slice_vpintavg.png");
    c_slice_vpintavg->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_slice_vpintavg.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_slice_vpintavg.png"));
    h_slice_vpintavg->Reset();

    TCanvas* c_detint = new TCanvas("c_detint", "Distribution of All Ev61Det Integral Values", 1200, 700);
    c_detint->SetLogy();
    c_detint->cd();
    h_detint->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_detint->GetYaxis()->SetTitle("Counts");
    h_detint->SetMinimum(1);
    h_detint->SetMaximum(20000);
    h_detint->SetLineColor(kBlue);
    h_detint->Draw();
    TLatex *lat_detint = new TLatex(0.6, 0.7, Form("Average background energy = %.2f", avg_detint));
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
    h_detint_1->SetMaximum(5000);
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
    h_dt_v_detint->SetMarkerStyle(7);
    h_dt_v_detint->Draw();
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

    TCanvas* c_dt_v_tdet = new TCanvas("c_dt_v_tdet", "Event 61 Detector dt vs Detector dt", 1200, 700);
    // c_dt_v_tdet->SetLogy();
    c_dt_v_tdet->cd();
    h_dt_v_tdet->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_v_tdet->GetYaxis()->SetTitle("Time Since Last Detector Event (ns)");
    h_dt_v_tdet->SetMarkerStyle(7);
    h_dt_v_tdet->Draw();
    c_dt_v_tdet->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_dt_v_tdet_%i.png", run));
    // c_dt_v_tdet->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_dt_v_tdet.png");
    c_dt_v_tdet->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_tdet.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_tdet.png"));
    h_dt_v_tdet->Reset();

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

    /*TCanvas* c_sample_events = new TCanvas("c_sample_events", "Background in Sample Windows", 900, 700);
    // c_sample_events->SetLogy();
    c_sample_events->cd();
    h_sample_events->GetXaxis()->SetTitle("Sample Window Number");
    h_sample_events->GetYaxis()->SetTitle("Number of Events in Window");
    h_sample_events->Draw("HIST");
    c_sample_events->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sample_events_%i.png", run));
    // c_sample_events->SaveAs("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/plots/h_sample_events.png");
    c_sample_events->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sample_events_%i.png", run));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sample_events.png"));
    h_sample_events->Reset();*/

    auto p_sample_events = new TGraphErrors();
    p_sample_events->SetTitle("Background Statistics as a Function of Distance from Event 61");
    p_sample_events->GetXaxis()->SetTitle("Sample Window Number");
    p_sample_events->GetYaxis()->SetTitle("Number of Events in Window");
    for(int iSample = 0; iSample < bg_samples.size(); iSample++) {
        p_sample_events->AddPoint(iSample, bg_samples[iSample]);
        p_sample_events->SetPointError(iSample, 0, TMath::Sqrt(bg_samples[iSample]));}
    p_sample_events->SetMarkerStyle(20);
    p_sample_events->SetMarkerSize(0.8);
    p_sample_events->SetMarkerColor(kBlack);
    // p_sample_events->SetMinimum(0);
    p_sample_events->Draw("AP");
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/p_sample_events_%i.png", run));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/p_sample_events.png"));

    cout << "\n" << "End of code." << "\n" << endl;

    return 0;

}