#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPad.h>
#include <iostream>
#include <fstream>
#include <cstdio>
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
};

/** @brief Constants for pulse and pulse-edge detection */
const int PULSE_THRESHOLD = 30;    /* Pulse detected if read above this value */
const int BS_UNCERTAINTY = 5;      /* Baseline uncertainty */
const int PULSE_PE_THRESHOLD = 40; /* Large pulse detected if read above this value (given in units of photoelectrons) */

/** @brief Maximum number of waveforms to process from input root file */
const int MAX_NUM_ENTRIES = 1900000;        // run4144: 2700000 ; run4176: 2500000 ; run4193: 1900000

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

    int run_counter = 0;

    int hour12_run = 0;

    Long64_t hour12_st = 0.0;

    // std::vector<double> integralToPE;

    // std::vector<double> integralToPE_25hrs;

    /* Initialize histograms */

    TH1D* h_sphe0 = new TH1D("h_sphe0", "Single Photoelectron Integral Distribution, Channel 0", 100, -50, 400);

    TH1D* h_sphe1 = new TH1D("h_sphe1", "Single Photoelectron Integral Distribution, Channel 1", 100, -50, 400);

    TH1D* h_sphe2 = new TH1D("h_sphe2", "Single Photoelectron Integral Distribution, Channel 2", 100, -50, 400);

    TH1D* h_sphe3 = new TH1D("h_sphe3", "Single Photoelectron Integral Distribution, Channel 3", 100, -50, 400);

    TH1D* h_sphe4 = new TH1D("h_sphe4", "Single Photoelectron Integral Distribution, Channel 4", 100, -50, 400);

    TH1D* h_sphe5 = new TH1D("h_sphe5", "Single Photoelectron Integral Distribution, Channel 5", 100, -50, 400);

    TH1D* h_sphe6 = new TH1D("h_sphe6", "Single Photoelectron Integral Distribution, Channel 6", 100, -50, 400);

    TH1D* h_sphe7 = new TH1D("h_sphe7", "Single Photoelectron Integral Distribution, Channel 7", 100, -50, 400);

    TH1D* h_sphe8 = new TH1D("h_sphe8", "Single Photoelectron Integral Distribution, Channel 8", 100, -50, 400);

    TH1D* h_sphe9 = new TH1D("h_sphe9", "Single Photoelectron Integral Distribution, Channel 9", 100, -50, 400);

    TH1D* h_sphe10 = new TH1D("h_sphe10", "Single Photoelectron Integral Distribution, Channel 10", 100, -50, 400);

    TH1D* h_sphe11 = new TH1D("h_sphe11", "Single Photoelectron Integral Distribution, Channel 11", 100, -50, 400);

    // TH1D* h_sphe2_top = new TH1D("h_sphe2_top", "Channel 2 Sphe Integrals, Top Bimodal Population", 100, -50, 400);

    // TH1D* h_sphe2_bot = new TH1D("h_sphe2_bot", "Channel 2 Sphe Integrals, Bottom Bimodal Population", 100, -50, 400);

    // Create a text string, which is used to output the text file
    
    std::vector<int> runlist;

    string intval;

    // Read from the text file

    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/spring2025_magnet_days/spring2025_magnet_days_%i.txt", run));         // runs 21625 - 21649(?) ; run 21629 missing ; runs 21745 - 21768
    // ifstream ReadRunListFile("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/top_pop_runs_all.txt");
    // ifstream ReadRunListFile("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/bot_pop_runs_all.txt");
    // ifstream ReadRunListFile("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/may_pop_runs_all.txt");
    // ifstream ReadRunListFile("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/event_rate_over_time_runlist.txt");
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/yearoneevents/yearoneevents_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/sum2024goldenrunlist/sum2024goldenrunlist_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/2024allgoldenruns/2024allgoldenruns_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/weekly2024goldenruns/weekly2024goldenruns_%i.txt", run));
    ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/allweeklyruns/allweeklyruns_%i.txt", run));            // runs 7580-13136 = 7/21/23-3/10/24 = files 2-31  // runs 14072-14976 = 4/24/24-5/31/24 = 1713954001-1717208401 = files 36-41
    // ifstream ReadRunListFile("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/2024allgoldenruns.txt");

    // Use a while loop together with the getline() function to read the file line by linecd

    int iLine = 0;

    while (getline(ReadRunListFile, intval)) {

        // Output the text from the file

        runlist.push_back(stod(intval));

        iLine += 1;

        if (iLine == 24) {break;}
        
    }
    
    // Close the file

    ReadRunListFile.close();

    ReadRunListFile.clear();

    // Create and open a text file
    // ofstream SPIFile("sphe_int_aggregated_runs_test.txt");
    // ofstream SPIFile("sphe_int_incomplete.txt");
    ofstream SPIFile("sphe_int_day_avg.txt");
    // ofstream TPFile(Form("top_pop_runs_%i.txt", run));
    // ofstream BPFile(Form("bot_pop_runs_%i.txt", run));
    // ofstream MPFile(Form("may_pop_runs_%i.txt", run));

    // for(int iRun = run; iRun <= last_run ; iRun++){

    for (size_t iRun = 0; iRun < runlist.size(); iRun++) {

        // int run_iterable = iRun;

        int run_iterable = runlist[iRun];

        // int data_num = 9; int vers_num = 4;

        int data_num = 41; int vers_num = 5;

        if (run_iterable <= 14976) {data_num = 9; vers_num = 4;}

        else if (run_iterable > 14976 && run_iterable < 15696) {data_num = 41; vers_num = 4;}

        else if (run_iterable >= 15696) {data_num = 41; vers_num = 5;}

        TFile *f; 
        //your root file location here
///        if(gSystem->AccessPathName(Form("/data%i/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v%i.root", data_num, run_iterable, vers_num))){
        if (gSystem->AccessPathName(Form("/data%i/coherent/data/d2o/processedData/run%i_processed_v%i.root", data_num, run_iterable, vers_num))) {
///        if(gSystem->AccessPathName(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/run%i_processed_v%i.root", run_iterable, vers_num))){
            cout << "Could not open file " << run_iterable << endl;
            continue;           // return -1; 
        } else{
///            f = new TFile(Form("/data%i/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v%i.root", data_num, run_iterable, vers_num));
            f = new TFile(Form("/data%i/coherent/data/d2o/processedData/run%i_processed_v%i.root", data_num, run_iterable, vers_num));
///            f = new TFile(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/run%i_processed_v%i.root", run_iterable, vers_num));
        }

        auto tsstart = (TParameter<Long64_t> *) f->Get("starttime");
        
        Long64_t run_starttime = tsstart->GetVal();

///        TTree *t = (TTree *) f->Get("wf");
        TTree* t = (TTree*)f->Get("tree");

        // TChain t("tree");

        // t.Add(Form("/data%i/coherent/data/d2o/processedData/run%i_processed_v%i.root", data_num, run_iterable, vers_num))

        // Declaration of leaf types
        Int_t           eventID;
        Int_t           nSamples[32];
        Short_t         adcVal[32][45];
        Double_t        baselineMean[32];
        Double_t        baselineRMS[32];
        Double_t        pulseH[32];
        Int_t           peakPosition[32];
        Double_t        area[32];
        Long64_t        nsTime;
        Int_t           triggerBits;

        // List of branches
        TBranch        *b_eventID;
        TBranch        *b_nSamples;
        TBranch        *b_adcVal;
        TBranch        *b_baselineMean;
        TBranch        *b_baselineRMS;
        TBranch        *b_pulseH;
        TBranch        *b_peakPosition;
        TBranch        *b_area;
        TBranch        *b_nsTime;
        TBranch        *b_triggerBits;

        t->SetBranchAddress("eventID", &eventID, &b_eventID);
        t->SetBranchAddress("nSamples", &nSamples, &b_nSamples);
        t->SetBranchAddress("adcVal", adcVal, &b_adcVal);
        t->SetBranchAddress("baselineMean", baselineMean, &b_baselineMean);
        t->SetBranchAddress("baselineRMS", baselineRMS, &b_baselineRMS);
        t->SetBranchAddress("pulseH", pulseH, &b_pulseH);
        t->SetBranchAddress("peakPosition", &peakPosition, &b_peakPosition);
        t->SetBranchAddress("area", area, &b_area);
        t->SetBranchAddress("nsTime", &nsTime, &b_nsTime);
        t->SetBranchAddress("triggerBits", &triggerBits, &b_triggerBits);

        // Get statistics for up to 1 million entries from fileIn TTree T
        int numEntries = std::min((int)t->GetEntries(), MAX_NUM_ENTRIES);
    
        // Replaced "t->GetEntries()" with "numEntries"
        for (int iEnt = 0; iEnt < t->GetEntries(); iEnt++) {

	        // std::cout << "\n" << "Processing event " << iEnt + 1 << " of " << t->GetEntries() << "\n";

            Long64_t tentry = t->LoadTree(iEnt);

            b_eventID->GetEntry(tentry);
            b_nSamples->GetEntry(tentry);
            b_adcVal->GetEntry(tentry);
            b_baselineMean->GetEntry(tentry);
            b_baselineRMS->GetEntry(tentry);
            b_pulseH->GetEntry(tentry);
            b_peakPosition->GetEntry(tentry);
            b_area->GetEntry(tentry);
            b_nsTime->GetEntry(tentry);
            b_triggerBits->GetEntry(tentry);

            // Create variables to hold info of each current pulse
            double AllPulseEnergy = 0.;

            /* Record single Ph.e. integral values by observing "low light trigger" events */

            if (triggerBits == 16) {

                for (int iChan = 0; iChan < 12; iChan++) {

                    for (int i = 0; i < ADCSIZE; i++) {

                        h_wf->SetBinContent(i + 1, adcVal[iChan][i] - baselineMean[iChan]);

                    }

                    AllPulseEnergy = 0.;

                    for (int iBin = 1; iBin < h_wf->GetNbinsX(); iBin++) {

                        double iBinContent = h_wf->GetBinContent(iBin);

                        if (iBin >= 20) {

                            AllPulseEnergy += iBinContent;

                        }

                    }

                    if (iChan == 0) {

                        h_sphe0->Fill(AllPulseEnergy);

                    }

                    if (iChan == 1) {

                        h_sphe1->Fill(AllPulseEnergy);

                    }

                    if (iChan == 2) {

                        h_sphe2->Fill(AllPulseEnergy);

                    }

                    if (iChan == 3) {

                        h_sphe3->Fill(AllPulseEnergy);

                    }

                    if (iChan == 4) {

                        h_sphe4->Fill(AllPulseEnergy);

                    }

                    if (iChan == 5) {

                        h_sphe5->Fill(AllPulseEnergy);

                    }

                    if (iChan == 6) {

                        h_sphe6->Fill(AllPulseEnergy);

                    }

                    if (iChan == 7) {

                        h_sphe7->Fill(AllPulseEnergy);

                    }

                    if (iChan == 8) {

                        h_sphe8->Fill(AllPulseEnergy);

                    }

                    if (iChan == 9) {

                        h_sphe9->Fill(AllPulseEnergy);

                    }

                    if (iChan == 10) {

                        h_sphe10->Fill(AllPulseEnergy);

                    }

                    if (iChan == 11) {

                        h_sphe11->Fill(AllPulseEnergy);

                    }

                    AllPulseEnergy = 0.;

                    h_wf->Reset();

                }   // Channel loop

            }   // Trigger Type loop

        } // Event loop

        run_counter += 1;

        if (run_counter == 12) {hour12_st = run_starttime; hour12_run = run_iterable;}

        else if (run_counter == 24) {

            /* Fill & plot histograms */

            // std::vector<double> *integralToPE = new std::vector<double>();

            TCanvas* c_sphe0 = new TCanvas("c_sphe0", "Single Photoelectron Integral Conversion", 900, 700);
            c_sphe0->cd();

            TF1* func0 = new TF1("func0", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
            double Par0[] = {100, 350000, 0, 20, 90000, 70, 20, 20000, 150, 20, 1000, 200, 20};
            func0->SetParameters(Par0);
            func0->SetParName(0, "Noise");
            func0->SetParName(1, "ZPhe_Amplitude");
            func0->SetParName(2, "ZPhe_PeakCenter");
            func0->SetParName(3, "ZPhe_StandardDeviation");
            func0->SetParName(4, "SPhe_Amplitude");
            func0->SetParName(5, "SPhe_PeakCenter");
            func0->SetParName(6, "SPhe_StandardDeviation");
            func0->SetParName(7, "2Phe_Amplitude");
            func0->SetParName(8, "2Phe_PeakCenter");
            func0->SetParName(9, "2Phe_StandardDeviation");
            func0->SetParName(10, "3Phe_Amplitude");
            func0->SetParName(11, "3Phe_PeakCenter");
            func0->SetParName(12, "3Phe_StandardDeviation");
            h_sphe0->Fit("func0", "R");

            double_t param0 = func0->GetParameter(5);

            // integralToPE.push_back(param0);

            TLatex *lat_sphe0 = new TLatex(0.6, 0.75, Form("Avg Sphe Int = %.2f", param0));
            lat_sphe0->SetNDC();
            lat_sphe0->SetTextColor(1);
            lat_sphe0->SetTextSize(0.02);
            lat_sphe0->Draw();

            h_sphe0->GetXaxis()->SetTitle("Integral (ADC)");
            h_sphe0->GetYaxis()->SetTitle("Counts");
            h_sphe0->Draw();
            func0->Draw("same");
            c_sphe0->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe0/h_sphe0_run%i.png", hour12_run));
            c_sphe0->Close();
            // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe0.png"));
            // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe0.png"));
            h_sphe0->Reset();

            TCanvas* c_sphe1 = new TCanvas("c_sphe1", "Single Photoelectron Integral Conversion", 900, 700);
            c_sphe1->cd();

            TF1* func1 = new TF1("func1", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
            double Par1[] = {100, 300000, 0, 20, 80000, 100, 30, 20000, 180, 40, 1000, 300, 50};
            func1->SetParameters(Par1);
            func1->SetParName(0, "Noise");
            func1->SetParName(1, "ZPhe_Amplitude");
            func1->SetParName(2, "ZPhe_PeakCenter");
            func1->SetParName(3, "ZPhe_StandardDeviation");
            func1->SetParName(4, "SPhe_Amplitude");
            func1->SetParName(5, "SPhe_PeakCenter");
            func1->SetParName(6, "SPhe_StandardDeviation");
            func1->SetParName(7, "2Phe_Amplitude");
            func1->SetParName(8, "2Phe_PeakCenter");
            func1->SetParName(9, "2Phe_StandardDeviation");
            func1->SetParName(10, "3Phe_Amplitude");
            func1->SetParName(11, "3Phe_PeakCenter");
            func1->SetParName(12, "3Phe_StandardDeviation");
            h_sphe1->Fit("func1", "R");

            double_t param1 = func1->GetParameter(5);

            // integralToPE.push_back(param1);

            TLatex *lat_sphe1 = new TLatex(0.6, 0.75, Form("Avg Sphe Int = %.2f", param1));
            lat_sphe1->SetNDC();
            lat_sphe1->SetTextColor(1);
            lat_sphe1->SetTextSize(0.02);
            lat_sphe1->Draw();

            h_sphe1->GetXaxis()->SetTitle("Integral (ADC)");
            h_sphe1->GetYaxis()->SetTitle("Counts");
            h_sphe1->Draw();
            func1->Draw("same");
            c_sphe1->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe1/h_sphe1_run%i.png", hour12_run));
            c_sphe1->Close();
            // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe1.png"));
            // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe1.png"));
            h_sphe1->Reset();

            TCanvas* c_sphe2 = new TCanvas("c_sphe2", "Single Photoelectron Integral Conversion", 900, 700);
            c_sphe2->cd();

            TF1* func2 = new TF1("func2", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
            double Par2[] = {100, 300000, 0, 20, 90000, 100, 30, 20000, 160, 40, 1000, 250, 50};
            func2->SetParameters(Par2);
            func2->SetParName(0, "Noise");
            func2->SetParName(1, "ZPhe_Amplitude");
            func2->SetParName(2, "ZPhe_PeakCenter");
            func2->SetParName(3, "ZPhe_StandardDeviation");
            func2->SetParName(4, "SPhe_Amplitude");
            func2->SetParName(5, "SPhe_PeakCenter");
            func2->SetParName(6, "SPhe_StandardDeviation");
            func2->SetParName(7, "2Phe_Amplitude");
            func2->SetParName(8, "2Phe_PeakCenter");
            func2->SetParName(9, "2Phe_StandardDeviation");
            func2->SetParName(10, "3Phe_Amplitude");
            func2->SetParName(11, "3Phe_PeakCenter");
            func2->SetParName(12, "3Phe_StandardDeviation");
            h_sphe2->Fit("func2", "R");

            double_t param2 = func2->GetParameter(5);
            /*
            if (run_starttime <= 1710102002) {           // runs 7580-13136 = 7/21/23-3/10/24 = files 2-31   // runs 14072-14976 = 4/24/24-5/31/24 = 1713954001-1717208401 = files 36-41

                if (param2 >= 98) {TPFile << run_iterable << "\n";}

                else if (param2 <= 97) {BPFile << run_iterable << "\n";}

                flush(TPFile);

                flush(BPFile);
            
            }

            else if (run_starttime >= 1713954001 && run_starttime <= 1717208401) {

                if (param2 >= 96.5 && param2 <= 97.5) {MPFile << run_iterable << "\n";}

                flush(MPFile);

            }
            */
            // integralToPE.push_back(param2);

            TLatex *lat_sphe2 = new TLatex(0.6, 0.75, Form("Avg Sphe Int = %.2f", param2));
            lat_sphe2->SetNDC();
            lat_sphe2->SetTextColor(1);
            lat_sphe2->SetTextSize(0.02);
            lat_sphe2->Draw();

            h_sphe2->GetXaxis()->SetTitle("Integral (ADC)");
            h_sphe2->GetYaxis()->SetTitle("Counts");
            h_sphe2->Draw();
            func2->Draw("same");
            c_sphe2->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe2/h_sphe2_run%i.png", hour12_run));
            c_sphe2->Close();
            // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe2.png"));
            // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe2.png"));
            h_sphe2->Reset();

            TCanvas* c_sphe3 = new TCanvas("c_sphe3", "Single Photoelectron Integral Conversion", 900, 700);
            c_sphe3->cd();

            TF1* func3 = new TF1("func3", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
            double Par3[] = {100, 350000, 0, 20, 70000, 90, 30, 20000, 160, 40, 1000, 250, 50};
            func3->SetParameters(Par3);
            func3->SetParName(0, "Noise");
            func3->SetParName(1, "ZPhe_Amplitude");
            func3->SetParName(2, "ZPhe_PeakCenter");
            func3->SetParName(3, "ZPhe_StandardDeviation");
            func3->SetParName(4, "SPhe_Amplitude");
            func3->SetParName(5, "SPhe_PeakCenter");
            func3->SetParName(6, "SPhe_StandardDeviation");
            func3->SetParName(7, "2Phe_Amplitude");
            func3->SetParName(8, "2Phe_PeakCenter");
            func3->SetParName(9, "2Phe_StandardDeviation");
            func3->SetParName(10, "3Phe_Amplitude");
            func3->SetParName(11, "3Phe_PeakCenter");
            func3->SetParName(12, "3Phe_StandardDeviation");
            h_sphe3->Fit("func3", "R");

            double_t param3 = func3->GetParameter(5);

            // integralToPE.push_back(param3);

            TLatex *lat_sphe3 = new TLatex(0.6, 0.75, Form("Avg Sphe Int = %.2f", param3));
            lat_sphe3->SetNDC();
            lat_sphe3->SetTextColor(1);
            lat_sphe3->SetTextSize(0.02);
            lat_sphe3->Draw();

            h_sphe3->GetXaxis()->SetTitle("Integral (ADC)");
            h_sphe3->GetYaxis()->SetTitle("Counts");
            h_sphe3->Draw();
            func3->Draw("same");
            c_sphe3->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe3/h_sphe3_run%i.png", hour12_run));
            c_sphe3->Close();
            // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe3.png"));
            // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe3.png"));
            h_sphe3->Reset();

            TCanvas* c_sphe4 = new TCanvas("c_sphe4", "Single Photoelectron Integral Conversion", 900, 700);
            c_sphe4->cd();

            TF1* func4 = new TF1("func4", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
            double Par4[] = {100, 300000, 0, 20, 60000, 130, 50, 20000, 230, 50, 1000, 350, 50};
            func4->SetParameters(Par4);
            func4->SetParName(0, "Noise");
            func4->SetParName(1, "ZPhe_Amplitude");
            func4->SetParName(2, "ZPhe_PeakCenter");
            func4->SetParName(3, "ZPhe_StandardDeviation");
            func4->SetParName(4, "SPhe_Amplitude");
            func4->SetParName(5, "SPhe_PeakCenter");
            func4->SetParName(6, "SPhe_StandardDeviation");
            func4->SetParName(7, "2Phe_Amplitude");
            func4->SetParName(8, "2Phe_PeakCenter");
            func4->SetParName(9, "2Phe_StandardDeviation");
            func4->SetParName(10, "3Phe_Amplitude");
            func4->SetParName(11, "3Phe_PeakCenter");
            func4->SetParName(12, "3Phe_StandardDeviation");
            h_sphe4->Fit("func4", "R");

            double_t param4 = func4->GetParameter(5);

            // integralToPE.push_back(param4);

            TLatex *lat_sphe4 = new TLatex(0.6, 0.75, Form("Avg Sphe Int = %.2f", param4));
            lat_sphe4->SetNDC();
            lat_sphe4->SetTextColor(1);
            lat_sphe4->SetTextSize(0.02);
            lat_sphe4->Draw();

            h_sphe4->GetXaxis()->SetTitle("Integral (ADC)");
            h_sphe4->GetYaxis()->SetTitle("Counts");
            h_sphe4->Draw();
            func4->Draw("same");
            c_sphe4->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe4/h_sphe4_run%i.png", hour12_run));
            c_sphe4->Close();
            // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe4.png"));
            // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe4.png"));
            h_sphe4->Reset();

            TCanvas* c_sphe5 = new TCanvas("c_sphe5", "Single Photoelectron Integral Conversion", 900, 700);
            c_sphe5->cd();

            TF1* func5 = new TF1("func5", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
            double Par5[] = {100, 350000, 0, 20, 80000, 100, 30, 20000, 160, 40, 1000, 250, 50};
            func5->SetParameters(Par5);
            func5->SetParName(0, "Noise");
            func5->SetParName(1, "ZPhe_Amplitude");
            func5->SetParName(2, "ZPhe_PeakCenter");
            func5->SetParName(3, "ZPhe_StandardDeviation");
            func5->SetParName(4, "SPhe_Amplitude");
            func5->SetParName(5, "SPhe_PeakCenter");
            func5->SetParName(6, "SPhe_StandardDeviation");
            func5->SetParName(7, "2Phe_Amplitude");
            func5->SetParName(8, "2Phe_PeakCenter");
            func5->SetParName(9, "2Phe_StandardDeviation");
            func5->SetParName(10, "3Phe_Amplitude");
            func5->SetParName(11, "3Phe_PeakCenter");
            func5->SetParName(12, "3Phe_StandardDeviation");
            h_sphe5->Fit("func5", "R");

            double_t param5 = func5->GetParameter(5);

            // integralToPE.push_back(param5);

            TLatex *lat_sphe5 = new TLatex(0.6, 0.75, Form("Avg Sphe Int = %.2f", param5));
            lat_sphe5->SetNDC();
            lat_sphe5->SetTextColor(1);
            lat_sphe5->SetTextSize(0.02);
            lat_sphe5->Draw();

            h_sphe5->GetXaxis()->SetTitle("Integral (ADC)");
            h_sphe5->GetYaxis()->SetTitle("Counts");
            h_sphe5->Draw();
            func5->Draw("same");
            c_sphe5->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe5/h_sphe5_run%i.png", hour12_run));
            c_sphe5->Close();
            // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe5.png"));
            // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe5.png"));
            h_sphe5->Reset();

            TCanvas* c_sphe6 = new TCanvas("c_sphe6", "Single Photoelectron Integral Conversion", 900, 700);
            c_sphe6->cd();

            TF1* func6 = new TF1("func6", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
            double Par6[] = {100, 300000, 0, 20, 90000, 110, 30, 20000, 230, 40, 1000, 300, 50};
            func6->SetParameters(Par6);
            func6->SetParName(0, "Noise");
            func6->SetParName(1, "ZPhe_Amplitude");
            func6->SetParName(2, "ZPhe_PeakCenter");
            func6->SetParName(3, "ZPhe_StandardDeviation");
            func6->SetParName(4, "SPhe_Amplitude");
            func6->SetParName(5, "SPhe_PeakCenter");
            func6->SetParName(6, "SPhe_StandardDeviation");
            func6->SetParName(7, "2Phe_Amplitude");
            func6->SetParName(8, "2Phe_PeakCenter");
            func6->SetParName(9, "2Phe_StandardDeviation");
            func6->SetParName(10, "3Phe_Amplitude");
            func6->SetParName(11, "3Phe_PeakCenter");
            func6->SetParName(12, "3Phe_StandardDeviation");
            h_sphe6->Fit("func6", "R");

            double_t param6 = func6->GetParameter(5);

            // integralToPE.push_back(param6);

            TLatex *lat_sphe6 = new TLatex(0.6, 0.75, Form("Avg Sphe Int = %.2f", param6));
            lat_sphe6->SetNDC();
            lat_sphe6->SetTextColor(1);
            lat_sphe6->SetTextSize(0.02);
            lat_sphe6->Draw();

            h_sphe6->GetXaxis()->SetTitle("Integral (ADC)");
            h_sphe6->GetYaxis()->SetTitle("Counts");
            h_sphe6->Draw();
            func6->Draw("same");
            c_sphe6->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe6/h_sphe6_run%i.png", hour12_run));
            c_sphe6->Close();
            // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe6.png"));
            // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe6.png"));
            h_sphe6->Reset();

            TCanvas* c_sphe7 = new TCanvas("c_sphe7", "Single Photoelectron Integral Conversion", 900, 700);
            c_sphe7->cd();

            TF1* func7 = new TF1("func7", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
            double Par7[] = {100, 300000, 0, 20, 80000, 110, 30, 20000, 200, 40, 1000, 300, 50};
            func7->SetParameters(Par7);
            func7->SetParName(0, "Noise");
            func7->SetParName(1, "ZPhe_Amplitude");
            func7->SetParName(2, "ZPhe_PeakCenter");
            func7->SetParName(3, "ZPhe_StandardDeviation");
            func7->SetParName(4, "SPhe_Amplitude");
            func7->SetParName(5, "SPhe_PeakCenter");
            func7->SetParName(6, "SPhe_StandardDeviation");
            func7->SetParName(7, "2Phe_Amplitude");
            func7->SetParName(8, "2Phe_PeakCenter");
            func7->SetParName(9, "2Phe_StandardDeviation");
            func7->SetParName(10, "3Phe_Amplitude");
            func7->SetParName(11, "3Phe_PeakCenter");
            func7->SetParName(12, "3Phe_StandardDeviation");
            h_sphe7->Fit("func7", "R");

            double_t param7 = func7->GetParameter(5);

            // integralToPE.push_back(param7);

            TLatex *lat_sphe7 = new TLatex(0.6, 0.75, Form("Avg Sphe Int = %.2f", param7));
            lat_sphe7->SetNDC();
            lat_sphe7->SetTextColor(1);
            lat_sphe7->SetTextSize(0.02);
            lat_sphe7->Draw();

            h_sphe7->GetXaxis()->SetTitle("Integral (ADC)");
            h_sphe7->GetYaxis()->SetTitle("Counts");
            h_sphe7->Draw();
            func7->Draw("same");
            c_sphe7->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe7/h_sphe7_run%i.png", hour12_run));
            c_sphe7->Close();
            // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe7.png"));
            // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe7.png"));
            h_sphe7->Reset();

            TCanvas* c_sphe8 = new TCanvas("c_sphe8", "Single Photoelectron Integral Conversion", 900, 700);
            c_sphe8->cd();

            TF1* func8 = new TF1("func8", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
            double Par8[] = {100, 300000, 0, 20, 70000, 100, 30, 20000, 200, 40, 1000, 300, 50};
            func8->SetParameters(Par8);
            func8->SetParName(0, "Noise");
            func8->SetParName(1, "ZPhe_Amplitude");
            func8->SetParName(2, "ZPhe_PeakCenter");
            func8->SetParName(3, "ZPhe_StandardDeviation");
            func8->SetParName(4, "SPhe_Amplitude");
            func8->SetParName(5, "SPhe_PeakCenter");
            func8->SetParName(6, "SPhe_StandardDeviation");
            func8->SetParName(7, "2Phe_Amplitude");
            func8->SetParName(8, "2Phe_PeakCenter");
            func8->SetParName(9, "2Phe_StandardDeviation");
            func8->SetParName(10, "3Phe_Amplitude");
            func8->SetParName(11, "3Phe_PeakCenter");
            func8->SetParName(12, "3Phe_StandardDeviation");
            h_sphe8->Fit("func8", "R");

            double_t param8 = func8->GetParameter(5);

            // integralToPE.push_back(param8);

            TLatex *lat_sphe8 = new TLatex(0.6, 0.75, Form("Avg Sphe Int = %.2f", param8));
            lat_sphe8->SetNDC();
            lat_sphe8->SetTextColor(1);
            lat_sphe8->SetTextSize(0.02);
            lat_sphe8->Draw();

            h_sphe8->GetXaxis()->SetTitle("Integral (ADC)");
            h_sphe8->GetYaxis()->SetTitle("Counts");
            h_sphe8->Draw();
            func8->Draw("same");
            c_sphe8->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe8/h_sphe8_run%i.png", hour12_run));
            c_sphe8->Close();
            // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe8.png"));
            // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe8.png"));
            h_sphe8->Reset();

            TCanvas* c_sphe9 = new TCanvas("c_sphe9", "Single Photoelectron Integral Conversion", 900, 700);
            c_sphe9->cd();

            TF1* func9 = new TF1("func9", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
            double Par9[] = {100, 400000, 0, 20, 40000, 120, 50, 4000, 200, 50, 400, 300, 50};
            func9->SetParameters(Par9);
            func9->SetParName(0, "Noise");
            func9->SetParName(1, "ZPhe_Amplitude");
            func9->SetParName(2, "ZPhe_PeakCenter");
            func9->SetParName(3, "ZPhe_StandardDeviation");
            func9->SetParName(4, "SPhe_Amplitude");
            func9->SetParName(5, "SPhe_PeakCenter");
            func9->SetParName(6, "SPhe_StandardDeviation");
            func9->SetParName(7, "2Phe_Amplitude");
            func9->SetParName(8, "2Phe_PeakCenter");
            func9->SetParName(9, "2Phe_StandardDeviation");
            func9->SetParName(10, "3Phe_Amplitude");
            func9->SetParName(11, "3Phe_PeakCenter");
            func9->SetParName(12, "3Phe_StandardDeviation");
            h_sphe9->Fit("func9", "R");

            double_t param9 = func9->GetParameter(5);

            // integralToPE.push_back(param9);

            TLatex *lat_sphe9 = new TLatex(0.6, 0.75, Form("Avg Sphe Int = %.2f", param9));
            lat_sphe9->SetNDC();
            lat_sphe9->SetTextColor(1);
            lat_sphe9->SetTextSize(0.02);
            lat_sphe9->Draw();

            h_sphe9->GetXaxis()->SetTitle("Integral (ADC)");
            h_sphe9->GetYaxis()->SetTitle("Counts");
            h_sphe9->Draw();
            func9->Draw("same");
            c_sphe9->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe9/h_sphe9_run%i.png", hour12_run));
            c_sphe9->Close();
            // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe9.png"));
            // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe9.png"));
            h_sphe9->Reset();

            TCanvas* c_sphe10 = new TCanvas("c_sphe10", "Single Photoelectron Integral Conversion", 900, 700);
            c_sphe10->cd();

            TF1* func10 = new TF1("func10", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
            double Par10[] = {100, 300000, 0, 20, 70000, 120, 40, 20000, 200, 50, 1000, 300, 50};
            func10->SetParameters(Par10);
            func10->SetParName(0, "Noise");
            func10->SetParName(1, "ZPhe_Amplitude");
            func10->SetParName(2, "ZPhe_PeakCenter");
            func10->SetParName(3, "ZPhe_StandardDeviation");
            func10->SetParName(4, "SPhe_Amplitude");
            func10->SetParName(5, "SPhe_PeakCenter");
            func10->SetParName(6, "SPhe_StandardDeviation");
            func10->SetParName(7, "2Phe_Amplitude");
            func10->SetParName(8, "2Phe_PeakCenter");
            func10->SetParName(9, "2Phe_StandardDeviation");
            func10->SetParName(10, "3Phe_Amplitude");
            func10->SetParName(11, "3Phe_PeakCenter");
            func10->SetParName(12, "3Phe_StandardDeviation");
            h_sphe10->Fit("func10", "R");

            double_t param10 = func10->GetParameter(5);

            // integralToPE.push_back(param10);

            TLatex *lat_sphe10 = new TLatex(0.6, 0.75, Form("Avg Sphe Int = %.2f", param10));
            lat_sphe10->SetNDC();
            lat_sphe10->SetTextColor(1);
            lat_sphe10->SetTextSize(0.02);
            lat_sphe10->Draw();

            h_sphe10->GetXaxis()->SetTitle("Integral (ADC)");
            h_sphe10->GetYaxis()->SetTitle("Counts");
            h_sphe10->Draw();
            func10->Draw("same");
            c_sphe10->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe10/h_sphe10_run%i.png", hour12_run));
            c_sphe10->Close();
            // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe10.png"));
            // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe10.png"));
            h_sphe10->Reset();

            TCanvas* c_sphe11 = new TCanvas("c_sphe11", "Single Photoelectron Integral Conversion", 900, 700);
            c_sphe11->cd();

            TF1* func11 = new TF1("func11", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
            double Par11[] = {100, 300000, 0, 20, 60000, 110, 40, 20000, 200, 50, 1000, 300, 50};
            func11->SetParameters(Par11);
            func11->SetParName(0, "Noise");
            func11->SetParName(1, "ZPhe_Amplitude");
            func11->SetParName(2, "ZPhe_PeakCenter");
            func11->SetParName(3, "ZPhe_StandardDeviation");
            func11->SetParName(4, "SPhe_Amplitude");
            func11->SetParName(5, "SPhe_PeakCenter");
            func11->SetParName(6, "SPhe_StandardDeviation");
            func11->SetParName(7, "2Phe_Amplitude");
            func11->SetParName(8, "2Phe_PeakCenter");
            func11->SetParName(9, "2Phe_StandardDeviation");
            func11->SetParName(10, "3Phe_Amplitude");
            func11->SetParName(11, "3Phe_PeakCenter");
            func11->SetParName(12, "3Phe_StandardDeviation");
            h_sphe11->Fit("func11", "R");

            double_t param11 = func11->GetParameter(5);

            // integralToPE.push_back(param11);

            TLatex *lat_sphe11 = new TLatex(0.6, 0.75, Form("Avg Sphe Int = %.2f", param11));
            lat_sphe11->SetNDC();
            lat_sphe11->SetTextColor(1);
            lat_sphe11->SetTextSize(0.02);
            lat_sphe11->Draw();

            h_sphe11->GetXaxis()->SetTitle("Integral (ADC)");
            h_sphe11->GetYaxis()->SetTitle("Counts");
            h_sphe11->Draw();
            func11->Draw("same");
            c_sphe11->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe11/h_sphe11_run%i.png", hour12_run));
            c_sphe11->Close();
            // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe11.png"));
            // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe11.png"));
            h_sphe11->Reset();

            std::cout << "\n" << "integralToPE, Run " << run << " = { " << param0 << "," << param1 << "," << param2 << "," << param3 << "," << param4 << "," << param5 << "," << param6 << "," << param7 << "," << param8 << "," << param9 << "," << param10 << "," << param11 << " }" << "\n";

            // for (int iVec = 0; iVec < integralToPE.size(); iVec++) {integralToPE_25hrs.push_back(integralToPE[iVec]);}

            // Write to the file
            SPIFile << hour12_st << "\t";
            
            SPIFile << param0 << "\t";

            SPIFile << param1 << "\t";

            SPIFile << param2 << "\t";

            SPIFile << param3 << "\t";

            SPIFile << param4 << "\t";

            SPIFile << param5 << "\t";

            SPIFile << param6 << "\t";

            SPIFile << param7 << "\t";

            SPIFile << param8 << "\t";

            SPIFile << param9 << "\t";

            SPIFile << param10 << "\t";

            SPIFile << param11 << "\n";

            flush(SPIFile);

            cout << "\n" << "Wrote to text file, Run " << run_iterable << endl;

            // delete integralToPE;

            // integralToPE.erase(integralToPE.begin(), integralToPE.end());

            run_counter = 0;

            // integralToPE.clear();

            // integralToPE.shrink_to_fit();

        }

        f->Close();

    } // Run loop
    /*
    TCanvas* c_sphe2_top = new TCanvas("c_sphe2_top", "Single Photoelectron Integral Conversion", 900, 700);
    c_sphe2_top->cd();

    TF1* func2t = new TF1("func2t", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
    double Par2t[] = {100, 300000, 0, 20, 90000, 100, 30, 20000, 160, 40, 1000, 250, 50};
    func2t->SetParameters(Par2t);
    func2t->SetParName(0, "Noise");
    func2t->SetParName(1, "ZPhe_Amplitude");
    func2t->SetParName(2, "ZPhe_PeakCenter");
    func2t->SetParName(3, "ZPhe_StandardDeviation");
    func2t->SetParName(4, "SPhe_Amplitude");
    func2t->SetParName(5, "SPhe_PeakCenter");
    func2t->SetParName(6, "SPhe_StandardDeviation");
    func2t->SetParName(7, "2Phe_Amplitude");
    func2t->SetParName(8, "2Phe_PeakCenter");
    func2t->SetParName(9, "2Phe_StandardDeviation");
    func2t->SetParName(10, "3Phe_Amplitude");
    func2t->SetParName(11, "3Phe_PeakCenter");
    func2t->SetParName(12, "3Phe_StandardDeviation");
    h_sphe2_top->Fit("func2t", "R");

    h_sphe2_top->GetXaxis()->SetTitle("Integral (ADC)");
    h_sphe2_top->GetYaxis()->SetTitle("Counts");
    h_sphe2_top->Draw();
    func2t->Draw("same");
    c_sphe2_top->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe2/h_sphe2_top_run%i.png", hour12_run));
    c_sphe2_top->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe2_top.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe2_top.png"));
    h_sphe2_top->Reset();

    TCanvas* c_sphe2_bot = new TCanvas("c_sphe2_bot", "Single Photoelectron Integral Conversion", 900, 700);
    c_sphe2_bot->cd();

    TF1* func2b = new TF1("func2b", "[0]+[1]*exp(-pow((x-[2]),2)/(2*pow([3],2)))+[4]*exp(-pow((x-[5]),2)/(2*pow([6],2)))+[7]*exp(-pow((x-[8]),2)/(2*pow([9],2)))+[10]*exp(-pow((x-[11]),2)/(2*pow([12],2)))", -50, 400);
    double Par2b[] = {100, 300000, 0, 20, 90000, 100, 30, 20000, 160, 40, 1000, 250, 50};
    func2b->SetParameters(Par2b);
    func2b->SetParName(0, "Noise");
    func2b->SetParName(1, "ZPhe_Amplitude");
    func2b->SetParName(2, "ZPhe_PeakCenter");
    func2b->SetParName(3, "ZPhe_StandardDeviation");
    func2b->SetParName(4, "SPhe_Amplitude");
    func2b->SetParName(5, "SPhe_PeakCenter");
    func2b->SetParName(6, "SPhe_StandardDeviation");
    func2b->SetParName(7, "2Phe_Amplitude");
    func2b->SetParName(8, "2Phe_PeakCenter");
    func2b->SetParName(9, "2Phe_StandardDeviation");
    func2b->SetParName(10, "3Phe_Amplitude");
    func2b->SetParName(11, "3Phe_PeakCenter");
    func2b->SetParName(12, "3Phe_StandardDeviation");
    h_sphe2_bot->Fit("func2b", "R");

    h_sphe2_bot->GetXaxis()->SetTitle("Integral (ADC)");
    h_sphe2_bot->GetYaxis()->SetTitle("Counts");
    h_sphe2_bot->Draw();
    func2b->Draw("same");
    c_sphe2_bot->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe2/h_sphe2_bot_run%i.png", hour12_run));
    c_sphe2_bot->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_sphe2_bot.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_sphe2_bot.png"));
    h_sphe2_bot->Reset();
    */
    // Close the file

    SPIFile.close();

    // TPFile.close();

    // BPFile.close();

    // MPFile.close();
    /*
    int PMT_count = 0;

    std::cout << "\n" << "integralToPE_25hrs = { ";

    for (int iVec = 0; iVec < integralToPE_25hrs.size(); iVec++) {
        
        std::cout << integralToPE_25hrs[iVec]; PMT_count +=1;
    
        if (PMT_count == 12) {std::cout << " }" << "\n" << "\n" << "integralToPE_25hrs = { "; PMT_count = 0;}

        else {std::cout << " , ";}
        
    }

    std::cout << "}" << "\n" << "\n";
    */
    return 0;

}
