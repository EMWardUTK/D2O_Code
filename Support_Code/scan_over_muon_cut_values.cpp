#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPad.h>
#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TSystem.h>
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"
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
const char* outputFileName = "PMTWaveformAnalysis";         // Why is output file going into ~ directory?
const char* outputStatsName = "PMTAnalysisStats";

std::vector<double> amplitudeToPE{ 26,28,30,29,27,12,30,28,30,30,24,30 };                // 12/3: run4176

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
    cout<<"run: "<< run <<" last run: " << last_run <<endl;

    int ADCSIZE = 45;
    TH1D *h_wf = new TH1D("h_wf", "Waveform", ADCSIZE, 0, ADCSIZE);

    /* Create vectors for plotting histograms */

    int muon_dt_min = 1 * pow(10, 4);

    int muon_dt_max = 2 * pow(10, 4);

    int muon_dt_step = 1 * pow(10, 3);

    int ev_61_dt_max = pow(10, 5);

    bool before_ev61_peak = true;

    bool bt_HL_n_LL = false;

    int num_events_low_dt = 0;

    /* Initialize histograms */                 // NUMBER OF BINS SHOULD BE A MULTIPLE OF THE SMALLEST AXIS UNIT!!! (OR LENGTH OF AXIS DIVIDED BY SMALLEST AXIS UNIT?)

    TH1D* h_dt_61 = new TH1D("h_dt_61", "Distribution of Time Separations Between Event 61 and Detector Peaks", 1875, - 15 * pow(10, 3), 15 * pow(10, 3));
    /*
    auto p_muon_dt_v_num_low_dt = new TGraph();
    TCanvas* c_muon_dt_v_num_low_dt = new TCanvas("c_muon_dt_v_num_low_dt", "Muon Cut Threshold vs Number of Low dt Events", 1200, 700);
    // c_muon_dt_v_num_low_dt->SetLogx();
    c_muon_dt_v_num_low_dt->cd();
    p_muon_dt_v_num_low_dt->SetTitle("Muon Cut Threshold vs Number of Low dt Events");
    p_muon_dt_v_num_low_dt->GetXaxis()->SetTitle("Delta-T (ns)");
    p_muon_dt_v_num_low_dt->GetYaxis()->SetTitle("Number of Events");
    */

    for(int iRun = run; iRun <= last_run ; iRun++){

        TFile *f; 
        //your root file location here
///        if(gSystem->AccessPathName(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v4.root", iRun))){
        if (gSystem->AccessPathName(Form("/data9/coherent/data/d2o/processedData/run%i_processed_v4.root", iRun))) {
///        if(gSystem->AccessPathName(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/run%i_processed_v4.root", iRun))){
            cout << "Could not open file" << endl;
            return -1; 
        } else{
///            f = new TFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v4.root", iRun));
            f = new TFile(Form("/data9/coherent/data/d2o/processedData/run%i_processed_v4.root", iRun));
///            f = new TFile(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/run%i_processed_v4.root", iRun));
        }

        // Create a text string, which is used to output the text file

        std::vector<double> integralToPE;

        integralToPE.clear();

        string int_val;

        // Read from the text file

        ifstream ReadSPIFile(Form("Single_Phe_Integral_Values_Run%i.txt", iRun));

        // Use a while loop together with the getline() function to read the file line by line

        while (getline(ReadSPIFile, int_val)) {

            // Output the text from the file

            integralToPE.push_back(stod(int_val));

        }

        // Close the file

        ReadSPIFile.close();

        // std::cout << "\n" << "integralToPE = { " << integralToPE[0] << " , " << integralToPE[1] << " , " << integralToPE[2] << " , " << integralToPE[3] << " , " << integralToPE[4] << " , " << integralToPE[5] << " , " << integralToPE[6] << " , " << integralToPE[7] << " , " << integralToPE[8] << " , " << integralToPE[9] << " , " << integralToPE[10] << " , " << integralToPE[11] << " }" << "\n";

///        TTree *t = (TTree *) f->Get("wf");
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

        // const char* outputFileName = Form("hdt61_val_start_run_%i", run);

        // Create output root file of processed waveform data
        TString outputFile;
        if (strcmp(outputFilePath, "") == 0) {
            outputFile.Form("%s.root", Form("hdt61_val_start_run_%i", iRun));
        } else {
            outputFile.Form("%s/%s.root", outputFilePath, Form("hdt61_val_start_run_%i", iRun));
        }
        TFile *fileOut = new TFile(outputFile, "RECREATE");
        TTree *hdt61Tree = new TTree("hdt61Tree", "hdt61Tree");
        hdt61Tree->SetDirectory(fileOut);

        // Data to keep track of for each hdt61Tree entry
        int *br_mudt;     /* Max muon dt cutoff value */
        int *br_numev;    /* Num of events in low dt region */

        hdt61Tree->Branch("mudt", &br_mudt, "mudt/I");
        hdt61Tree->Branch("numev", &br_numev, "numev/I");

        for (int iMuon_dt = muon_dt_min; iMuon_dt <= muon_dt_max ; iMuon_dt += muon_dt_step) {

            // Get statistics for up to 1 million entries from fileIn TTree T
            int numEntries = std::min((int)t->GetEntries(), MAX_NUM_ENTRIES);

            std::vector<double> chan_lengths;

            std::vector<double> peak_pos_RMS;

            std::vector<double> chan_start_no_outliers;

            double avg_peak_pos_RMS = 0.0;

            double var_peak_pos_RMS = 0.0;

            int ch15_on = 0;

            double last_muon_time = 0.0;

            double last_det_time = 0.0;

            double last_ev61_time = 0.0;

            double last_HL_time = 0.0;

            double last_LL_time = 0.0;
        
            // Can use either "t->GetEntries()" or "numEntries" or number
            for (int iEnt = 0; iEnt < t->GetEntries(); iEnt++) {

                // std::cout << "\n" << "Processing event " << iEnt + 1 << " of " << t->GetEntries() << "\n";

                Long64_t tentry = t->LoadTree(iEnt);

                b_eventID->GetEntry(tentry);
                b_nSamples->GetEntry(tentry);
                ///            b_adcTime->GetEntry(tentry);
                ///            b_adcSize->GetEntry(tentry);
                b_adcVal->GetEntry(tentry);
                ///            b_trigPattern->GetEntry(tentry);
                ///            b_timeStamp_extTrig->GetEntry(tentry);
                b_baselineMean->GetEntry(tentry);
                b_baselineRMS->GetEntry(tentry);
                ///            b_peakBin->GetEntry(tentry);
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
                double chan_15_start = 0.;

                std::vector<double> all_chan_start;
                std::vector<double> all_chan_end;
                std::vector<double> all_chan_peak;
                std::vector<double> all_chan_energy;
                std::vector<double> peak_energy;
                std::vector<double> ten_chan_peak;
                std::vector<double> all_chan_peakbin;
                std::vector<double> ten_chan_peakbin;
                std::vector<double> side_veto_panel_energy;
                std::vector<double> top_veto_panel_energy;
                std::vector<double> all_chan_start_vp;

                bool all_chan_beam = false;

                bool top_vp_event = false;

                int num_chan = 0;

                int num_chan_over_1phe = 0;

                int num_chan_simple = 0;

                double num_chan_vp = 0;

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

                chan_lengths.clear();

                for (int iChan = 0; iChan < 26; iChan++) {

                    if (iChan == 12 || iChan == 13 || iChan == 14) {

                        continue;

                    }

                    for (int i = 0; i < ADCSIZE; i++) {

                        h_wf->SetBinContent(i + 1, adcVal[iChan][i] - baselineMean[iChan]);

                    }

                    // h_wf->GetXaxis()->SetTitle("Time (1 bin = 16 ns)");
                    // h_wf->GetYaxis()->SetTitle("Counts");
                    // h_wf->Draw();
                    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/wf_image/wf_run_%i_event_%i_ch_%i.png", iRun, iEnt + 1, iChan));
                    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/wf_image/wf_run_%i_event_%i_ch_%i.png", iRun, iEnt + 1, iChan));

                    if (iChan == 15) {

                        for (int iBin = 1; iBin <= h_wf->GetNbinsX(); iBin++) {

                            double iBinContent = h_wf->GetBinContent(iBin);

                            Ev61Energy += iBinContent;

                        }

                        if (Ev61Energy > EV61_THRESHOLD) {

                            all_chan_beam = true;

                            ch15_on += 1;

                        }

                        Ev61Energy = 0.;

                    }

                    pulses_temp.clear();

                    onLastPulseTail = false;

                    bool peak_in_chan = false;

                    bool peak_in_chan_simple = false;

                    for (int iBin = 1; iBin <= h_wf->GetNbinsX(); iBin++) {

                        double iBinContent = h_wf->GetBinContent(iBin);
                        // std::cout << iBinContent << "\t";

                        if (iChan <= 11 && !peak_in_chan_simple && iBinContent >= PULSE_THRESHOLD) {

                            peak_in_chan_simple = true;

                            num_chan_simple += 1;

                        }
                        /*
                        if (iChan >= 16 && !peak_in_chan_veto && iBinContent >= PULSE_THRESHOLD) {

                            peak_in_chan_veto = true;

                            num_chan_veto += 1;

                        }
                        */
                        if (iBin > 15) {

                            AllPulseEnergy += iBinContent;

                        }

                        if (pulses_temp.size() > 0) {    // Could get rid of this check b/c event window is so small

                            onLastPulseTail = iBin * 16 - pulses_temp.back().start < tailWindow;

                        }

                        // Find pulse
                        if (!onPulse && !onLastPulseTail && iBinContent >= PULSE_THRESHOLD) {

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
                            if (iBinContent < BS_UNCERTAINTY) {

                                // std::cout << "\n" << "Found a pulse in Event " << iEnt + 1 << " Channel " << iChan << "\n";

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
                                    }
                                    // Record start of pulse (10% of peak)
                                    if (h_wf->GetBinContent(j) > peak * 0.1) {
                                        p.start = j * 16.;
                                    }
                                }

                                if (iChan <= 11) {

                                    // Record energy of pulse
                                    p.energy = pulseEnergy / integralToPE[iChan];

                                    all_chan_start.push_back(p.start);      // This is triggering for multiple bins above threshold? Even within single peak?
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

                                if (iChan == 15) {

                                    chan_15_start = p.start;

                                }

                                if (iChan >= 16) {

                                    all_chan_start_vp.push_back(p.start);

                                }

                                // Clear current pulse variables to look for new pulse
                                peak = 0.;
                                peakBin = 0;
                                pulseEnergy = 0.;
                                thresholdBin = 0;
                                onPulse = false;

                            }

                            if (iBin == ADCSIZE) {
                            
                                // Clear current pulse variables to look for new pulse
                                peak = 0.;
                                peakBin = 0;
                                pulseEnergy = 0.;
                                thresholdBin = 0;
                                onPulse = false;

                            }

                        }

                    }

                    if (iChan <= 11) {

                        all_chan_energy.push_back(AllPulseEnergy / integralToPE[iChan]);

                        if (AllPulseEnergy / integralToPE[iChan] > 1) {

                            num_chan_over_1phe += 1;

                        }

                        if (!peak_in_chan) {

                            ten_chan_peak.push_back(0);

                            ten_chan_peakbin.push_back(0);

                        }

                    }

                    if (iChan >= 16 && iChan <= 23) {

                        if (AllPulseEnergy > 200) {

                            num_chan_vp += 1;

                        }

                        if (iChan == 16) {

                            side_veto_panel_energy.push_back(AllPulseEnergy);

                            p_int_16 = AllPulseEnergy;

                        }

                        else if (iChan == 17) {

                            side_veto_panel_energy.push_back(AllPulseEnergy);

                            p_int_17 = AllPulseEnergy;

                        }

                        else if (iChan == 18) {

                            side_veto_panel_energy.push_back(AllPulseEnergy);

                            p_int_18 = AllPulseEnergy;

                        }

                        else if (iChan == 19) {

                            side_veto_panel_energy.push_back(AllPulseEnergy);

                            p_int_19 = AllPulseEnergy;

                        }

                        else if (iChan == 20) {

                            side_veto_panel_energy.push_back(AllPulseEnergy);

                            p_int_20 = AllPulseEnergy;

                        }

                        else if (iChan == 21) {

                            side_veto_panel_energy.push_back(AllPulseEnergy);

                            p_int_21 = AllPulseEnergy;

                        }

                        else if (iChan == 22) {

                            side_veto_panel_energy.push_back(AllPulseEnergy);

                            p_int_22 = AllPulseEnergy;

                        }

                        else if (iChan == 23) {

                            side_veto_panel_energy.push_back(AllPulseEnergy);

                            p_int_23 = AllPulseEnergy;

                        }

                    }

                    if (iChan >= 24) {

                        if (!top_vp_event && AllPulseEnergy > 200) {

                            num_chan_vp += 1;

                            top_vp_event = true;
                            
                        }

                        if (iChan == 24) {

                            top_veto_panel_energy.push_back(AllPulseEnergy);

                            p_int_24 = 1.07809 * AllPulseEnergy;

                        }

                        else if (iChan == 25) {

                            top_veto_panel_energy.push_back(AllPulseEnergy);

                            p_int_25 = AllPulseEnergy;

                        }

                    }

                    AllPulseEnergy = 0.;

                    h_wf->Reset();

                    chan_lengths.push_back(nSamples[iChan]);

                }   // Channel loop

                double pulse_start_time = mostFrequent(all_chan_start);

                double vp_start_time = nsTime + mostFrequent(all_chan_start_vp);

                if (num_chan_over_1phe >= 10) {

                    int pulse_start_index = mostFrequentIndex(all_chan_start);

                }

                //if (pulse_start_time == 0) {

                    //continue;

                //}

                struct pulse avg_pulse;

                avg_pulse.start = nsTime + pulse_start_time;                                                       // Need a way to weed out secondary pulses - use weighted average?
                avg_pulse.end = nsTime + mostFrequent(all_chan_end);
                avg_pulse.energy = std::accumulate(all_chan_energy.begin(), all_chan_energy.end(), 0.0);           // Integrals still sum over secondary peaks discarded in .start and .end values 
                avg_pulse.peak = std::accumulate(all_chan_peak.begin(), all_chan_peak.end(), 0.0);
                avg_pulse.number = num_chan_over_1phe;
                avg_pulse.beam = all_chan_beam;
                avg_pulse.trigger = triggerBits;
                avg_pulse.length = getAverage(chan_lengths);
                avg_pulse.side_vp_energy = std::accumulate(side_veto_panel_energy.begin(), side_veto_panel_energy.end(), 0.0);
                avg_pulse.top_vp_energy = std::accumulate(top_veto_panel_energy.begin(), top_veto_panel_energy.end(), 0.0);

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

                if (var_val < 5 * 16) {            // Pulses largely not consistent, need to weed out outlier and secondary pulses

                    avg_pulse.single = true;

                }

                else {

                    avg_pulse.single = false;

                }

                num_chan = 0;
                num_chan_over_1phe = 0;
                num_chan_simple = 0;

                all_chan_beam = false;

                top_vp_event = false;

                all_chan_start.clear();
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

                chan_lengths.clear();

                /*
                
                New Ev61-Det dt plot code loop
                
                */

                // Only looking for events between high light LED pulse and low light LED pulse

                if (avg_pulse.trigger == 8) {bt_HL_n_LL = true; last_HL_time = avg_pulse.start;}

                else if (avg_pulse.trigger == 16) {bt_HL_n_LL = false; before_ev61_peak = true; last_muon_time = 0.0; last_det_time = 0.0; last_ev61_time = 0.0; last_LL_time = avg_pulse.start;}

                // Stop looking (only) for detector pulses after Event 61 pulses

                if (avg_pulse.start - last_ev61_time > ev_61_dt_max) {before_ev61_peak = true;}

                // Is this event a muon?

                if (avg_pulse.energy >= 100 && (p_int_16 >= 750 || p_int_17 >= 950 || p_int_18 >= 1200 || p_int_19 >= 1375 || p_int_20 >= 525 || p_int_21 >= 700 || p_int_22 >= 700 || p_int_23 >= 500 || avg_pulse.top_vp_energy >= 450)) {

                    last_muon_time = vp_start_time;

                    all_chan_start_vp.clear();

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

                    continue;

                }

                // Is this event an internally triggered detector pulse?

                if (bt_HL_n_LL && avg_pulse.number >= 10 && avg_pulse.energy >= 100) {         // Should not limit when/where I am searching for internally triggered detector peaks? No "bt_HL_n_LL" check here?

                    // What if event is an Event 61 pulse AND an internally triggered detector pulse?
                    
                    if (avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {

                        if (p_int_16 <= 750 && p_int_17 <= 950 && p_int_18 <= 1200 && p_int_19 <= 1375 && p_int_20 <= 525 && p_int_21 <= 700 && p_int_22 <= 700 && p_int_23 <= 500 && avg_pulse.top_vp_energy <= 450) {

                            last_det_time = avg_pulse.start;

                            // Is this event ALSO an Event 61 event?

                            if (avg_pulse.beam) {

                                last_ev61_time = nsTime + chan_15_start;

                                before_ev61_peak = false;
                                
                            }

                            if (before_ev61_peak) {

                                all_chan_start_vp.clear();

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

                                continue;

                            }

                            else if (!before_ev61_peak) {           // Need another cut for first detector check to see if there was recent Event 61 peak?

                                if (last_det_time - last_muon_time < iMuon_dt) {

                                    before_ev61_peak = true;           // Delete this?

                                    all_chan_start_vp.clear();

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

                                    continue;

                                }

                                else if (last_ev61_time != 0.0 && (last_det_time - last_ev61_time) <= ev_61_dt_max) {

                                    if ((last_det_time - last_ev61_time) > -10000 && (last_det_time - last_ev61_time) < 10000) {num_events_low_dt += 1;}

                                    h_dt_61->Fill(last_det_time - last_ev61_time);

                                    before_ev61_peak = true;

                                    all_chan_start_vp.clear();

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

                                    continue;

                                }

                                else {

                                    before_ev61_peak = true;           // Delete this?

                                    all_chan_start_vp.clear();

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

                                    continue;

                                }

                            }

                        }

                    }

                }

                // Is this event an Event 61 pulse?

                if (bt_HL_n_LL && before_ev61_peak && avg_pulse.beam) {            // Does before_ev61_peak check make me miss ev61 events?

                    last_ev61_time = nsTime + chan_15_start;

                    if (last_ev61_time - last_muon_time < iMuon_dt) {

                        before_ev61_peak = false;           // Delete this?

                        all_chan_start_vp.clear();

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

                        continue;

                    }

                    else if (last_det_time != 0.0 && (last_ev61_time - last_det_time) <= ev_61_dt_max) {

                        if ((last_det_time - last_ev61_time) > -10000 && (last_det_time - last_ev61_time) < 10000) {num_events_low_dt += 1;}

                        h_dt_61->Fill(last_det_time - last_ev61_time);

                        before_ev61_peak = false;

                        all_chan_start_vp.clear();

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

                        continue;

                    }
                    
                    else {

                        before_ev61_peak = false;           // Delete this?

                        all_chan_start_vp.clear();

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
                        
                        continue;
                        
                    }
                    
                }

                all_chan_start_vp.clear();

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

                // Ev61Energy = 0.;

            } // Event loop

            ch15_on = 0;

            integralToPE.clear();

            avg_peak_pos_RMS = getAverage(peak_pos_RMS);

            var_peak_pos_RMS = variance(peak_pos_RMS);

            peak_pos_RMS.clear();

            // p_muon_dt_v_num_low_dt->AddPoint(iMuon_dt, num_events_low_dt);

            int mudt = iMuon_dt;
            int numev = num_events_low_dt;

            hdt61Tree->SetBranchAddress("mudt", &mudt);
            hdt61Tree->SetBranchAddress("numev", &numev);

            hdt61Tree->Fill();

            num_events_low_dt = 0;

        } // Muon cut loop

        hdt61Tree->Write();

        fileOut->Close();

    } // Run loop

    /* Fill & plot histograms */
    /*
    TCanvas* c_dt_61 = new TCanvas("c_dt_61", "Delta-T Between Event 61 and Detector Events", 900, 700);
    c_dt_61->SetLogy();
    c_dt_61->cd();
    h_dt_61->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61->GetYaxis()->SetTitle("Counts");
    h_dt_61->Draw();
    c_dt_61->SaveAs("h_dt_61.png");
    c_dt_61->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61.png"));
    h_dt_61->Reset();

    p_muon_dt_v_num_low_dt->Draw("AL*");
    c_muon_dt_v_num_low_dt->SaveAs("p_muon_dt_v_num_low_dt.png");
    c_muon_dt_v_num_low_dt->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/p_muon_dt_v_num_low_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/p_muon_dt_v_num_low_dt.png"));
    */
    std::cout << "\n" << "End of code." << "\n";

    return 0;

}
