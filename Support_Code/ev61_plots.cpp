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
#include "TTimeStamp.h"
#include <stdio.h>
#include <string.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <bits/stdc++.h>
#include <chrono>
#include <ctime>
#include <iomanip>

using std::cout; 
using std::endl;
using namespace std;

/** @brief Constants for pulse and pulse-edge detection */
const int PULSE_THRESHOLD = 30;    /* Pulse detected if read above this value */
const int BS_UNCERTAINTY = 5;      /* Baseline uncertainty */
const int PULSE_PE_THRESHOLD = 40; /* Large pulse detected if read above this value (given in units of photoelectrons) */
const int EV61_THRESHOLD = 1000;    /* SNS Beam assumed to be on if Event 61 read above this value (given in units of ADC) */

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

    // std::vector<double> integralToPE;

    // std::vector<double> integralToPE_25hrs;

    /* Initialize Plots */

    TH1D* h_intev61 = new TH1D("h_intev61", "Event 61 Integral Values", 2200, -200, 2000);

    auto p_ev61_avg_int = new TGraph();

    // TH1D* h_sphe2_top = new TH1D("h_sphe2_top", "Channel 2 Sphe Integrals, Top Bimodal Population", 100, -50, 400);

    // TH1D* h_sphe2_bot = new TH1D("h_sphe2_bot", "Channel 2 Sphe Integrals, Bottom Bimodal Population", 100, -50, 400);

    // Create a text string, which is used to output the text file
    
    // std::vector<int> runlist;

    // string intval;

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
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/allweeklyruns/allweeklyruns_%i.txt", run));            // runs 7580-13136 = 7/21/23-3/10/24 = files 2-31  // runs 14072-14976 = 4/24/24-5/31/24 = 1713954001-1717208401 = files 36-41
    // ifstream ReadRunListFile("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/2024allgoldenruns.txt");
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/F25-S26_GR_Runs_4084GR/F25-S26_GR_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/F25-S26_GR_Weekly_4084GR/F25-S26_GR_Weekly_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/S25-S26_GR_Runs_6235GR/S25-S26_GR_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/S25-S26_GR_Weekly_6149GR/S25-S26_GR_Weekly_%i.txt", run));
    // ifstream ReadRunListFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/S25F25S26_GR_3_rand_days_%i.txt", run));

    // Use a while loop together with the getline() function to read the file line by linecd

    // int iLine = 0;

    // while (getline(ReadRunListFile, intval)) {

        // Output the text from the file

        // runlist.push_back(stod(intval));

        // iLine += 1;

        // if (iLine == 24) {break;}
        
    // }
    
    // Close the file

    // ReadRunListFile.close();

    // ReadRunListFile.clear();

    for(int iRun = run; iRun <= last_run ; iRun++){

        // if (iRun % 50 != 0) {continue;}

        // for (size_t iRun = 0; iRun < runlist.size(); iRun++) {

        int run_iterable = iRun;

        // int run_iterable = runlist[iRun];

        // int data_num = 9; int vers_num = 4;

        int data_num = 41; int vers_num = 5;

        if (run_iterable <= 14976) {data_num = 9; vers_num = 4;}

        else if (run_iterable > 14976 && run_iterable < 15696) {data_num = 41; vers_num = 4;}

        else if (run_iterable >= 15696) {data_num = 41; vers_num = 5;}

        int ev61_min = 1200; int ev61_max = 2000;

        if (run_iterable < 24813) {ev61_min = 1200; ev61_max = 2000;}

        else if (run_iterable >= 24813) {ev61_min = 700; ev61_max = 1000;}

        bool old_channel_map = false; int Event61Chan = 15; int SiPM1Chan = 16; int SiPM2Chan = 17; int SiPM3Chan = 18; int SiPM4Chan = 19; int SiPM5Chan = 20; int SiPM6Chan = 21; int SiPM7Chan = 22; int SiPM8Chan = 23; int SiPM9Chan = 24; int SiPM10Chan = 25;

        if (vers_num < 5) {old_channel_map = true; Event61Chan = 15; SiPM1Chan = 16; SiPM2Chan = 17; SiPM3Chan = 18; SiPM4Chan = 19; SiPM5Chan = 20; SiPM6Chan = 21; SiPM7Chan = 22; SiPM8Chan = 23; SiPM9Chan = 24; SiPM10Chan = 25;}

        else if (vers_num == 5) {old_channel_map = false; SiPM1Chan = 12; SiPM2Chan = 13; SiPM3Chan = 14; SiPM4Chan = 15; SiPM5Chan = 16; SiPM6Chan = 17; SiPM7Chan = 18; SiPM8Chan = 19; SiPM9Chan = 20; SiPM10Chan = 21; Event61Chan = 22;}

        TFile *f; 
        //your root file location here
        // if(gSystem->AccessPathName(Form("/data%i/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v%i.root", data_num, run_iterable, vers_num))){
        if (gSystem->AccessPathName(Form("/data%i/coherent/data/d2o/processedData/run%i_processed_v%i.root", data_num, run_iterable, vers_num))) {
            // if(gSystem->AccessPathName(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/run%i_processed_v%i.root", run_iterable, vers_num))){
            cout << "Could not open file " << run_iterable << endl;
            continue;           // return -1; 
        } else{
            // f = new TFile(Form("/data%i/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v%i.root", data_num, run_iterable, vers_num));
            f = new TFile(Form("/data%i/coherent/data/d2o/processedData/run%i_processed_v%i.root", data_num, run_iterable, vers_num));
            // f = new TFile(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/run%i_processed_v%i.root", run_iterable, vers_num));
        }

        auto tsstart = (TParameter<Long64_t> *) f->Get("starttime");
        
        Long64_t run_starttime = tsstart->GetVal();

        // TTree *t = (TTree *) f->Get("wf");
        TTree* t = (TTree*)f->Get("tree");

        // TChain t("tree");

        // t.Add(Form("/data%i/coherent/data/d2o/processedData/run%i_processed_v%i.root", data_num, run_iterable, vers_num))

        // Declaration of leaf types
        // Int_t           eventID;
        // Int_t           nSamples[32];
        Short_t         adcVal[32][45];
        Double_t        baselineMean[32];
        // Double_t        baselineRMS[32];
        // Double_t        pulseH[32];
        // Int_t           peakPosition[32];
        // Double_t        area[32];
        // Long64_t        nsTime;
        // Int_t           triggerBits;

        // List of branches
        // TBranch        *b_eventID;
        // TBranch        *b_nSamples;
        TBranch        *b_adcVal;
        TBranch        *b_baselineMean;
        // TBranch        *b_baselineRMS;
        // TBranch        *b_pulseH;
        // TBranch        *b_peakPosition;
        // TBranch        *b_area;
        // TBranch        *b_nsTime;
        // TBranch        *b_triggerBits;

        // t->SetBranchAddress("eventID", &eventID, &b_eventID);
        // t->SetBranchAddress("nSamples", &nSamples, &b_nSamples);
        t->SetBranchAddress("adcVal", adcVal, &b_adcVal);
        t->SetBranchAddress("baselineMean", baselineMean, &b_baselineMean);
        // t->SetBranchAddress("baselineRMS", baselineRMS, &b_baselineRMS);
        // t->SetBranchAddress("pulseH", pulseH, &b_pulseH);
        // t->SetBranchAddress("peakPosition", &peakPosition, &b_peakPosition);
        // t->SetBranchAddress("area", area, &b_area);
        // t->SetBranchAddress("nsTime", &nsTime, &b_nsTime);
        // t->SetBranchAddress("triggerBits", &triggerBits, &b_triggerBits);

        // Get statistics for up to 1 million entries from fileIn TTree T
        // int numEntries = std::min((int)t->GetEntries(), MAX_NUM_ENTRIES);

        // Event 61 average integral for each run
        double Ev61IntTot = 0.0;

        int Ev61Count = 0;
    
        // Replaced "t->GetEntries()" with "numEntries"
        for (int iEnt = 0; iEnt < t->GetEntries(); iEnt++) {

	        // std::cout << "\n" << "Processing event " << iEnt + 1 << " of " << t->GetEntries() << "\n";

            Long64_t tentry = t->LoadTree(iEnt);

            // b_eventID->GetEntry(tentry);
            // b_nSamples->GetEntry(tentry);
            b_adcVal->GetEntry(tentry);
            b_baselineMean->GetEntry(tentry);
            // b_baselineRMS->GetEntry(tentry);
            // b_pulseH->GetEntry(tentry);
            // b_peakPosition->GetEntry(tentry);
            // b_area->GetEntry(tentry);
            // b_nsTime->GetEntry(tentry);
            // b_triggerBits->GetEntry(tentry);

            // Create variables to hold info of each current pulse
            double AllPulseEnergy = 0.;

            for (int iChan = 15; iChan < 23; iChan++) {

                if (iChan != Event61Chan) {continue;}

                for (int i = 0; i < ADCSIZE; i++) {

                    h_wf->SetBinContent(i + 1, adcVal[iChan][i] - baselineMean[iChan]);

                }

                AllPulseEnergy = 0.;

                for (int iBin = 1; iBin < h_wf->GetNbinsX(); iBin++) {
                    // for (int iBin = 23; iBin <= 39; iBin++) {

                    double iBinContent = h_wf->GetBinContent(iBin);

                    // if (iBin >= 23 && iBin <= 39) {

                    AllPulseEnergy += iBinContent;

                    // }

                }

                if (iChan == Event61Chan && AllPulseEnergy > 1200 && AllPulseEnergy < 1400) {h_intev61->Fill(AllPulseEnergy);}

                if (iChan == Event61Chan && AllPulseEnergy > 750) {

                    Ev61IntTot += AllPulseEnergy;

                    Ev61Count += 1;

                }

                AllPulseEnergy = 0.;

                h_wf->Reset();

            }   // Channel loop

        } // Event loop

        if (Ev61Count > 1000) {
            
            double Ev61IntAvg = Ev61IntTot / Ev61Count;
            
            p_ev61_avg_int->AddPoint(run_starttime, Ev61IntAvg);

            if (Ev61IntAvg < 1200 && run_starttime > 1735603261 && run_starttime < 1740787261) {cout << "\n" << "Low Event 61 integral run in January/February 2025: " << run_iterable << "\n" << endl;}
        
        }

        f->Close();

        // run_counter += 1;

        // if (run_counter % 10 == 0) {cout << "\n" << "Finished run " << run_counter << endl;}

        TTimeStamp ts((unsigned int)run_starttime);

        TCanvas* c_intev61 = new TCanvas("c_intev61", "Distribution of Event 61 Integral Values", 1200, 700);
        c_intev61->SetLogy();
        c_intev61->cd();
        h_intev61->SetTitle(Form("Event 61 Integral Values for %s", ts.AsString("s")));
        h_intev61->GetXaxis()->SetTitle("Integral (ADC)");
        h_intev61->GetYaxis()->SetTitle("Counts");
        h_intev61->Draw();
        c_intev61->SaveAs(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/h_intev61_%i.png", run_iterable));
        c_intev61->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_intev61.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_intev61.png"));
        h_intev61->Reset();

    } // Run loop

    p_ev61_avg_int->SetTitle("Average Event 61 Integral Over Time");
    p_ev61_avg_int->GetXaxis()->SetTitle("Date");
    p_ev61_avg_int->GetYaxis()->SetTitle("Integral (ADC)");
    p_ev61_avg_int->Draw("A*");
    // Tell the axis its values are Unix timestamps, and set the offset to epoch (0)
    p_ev61_avg_int->GetXaxis()->SetTimeDisplay(1);
    p_ev61_avg_int->GetXaxis()->SetTimeOffset(0, "gmt");   // "gmt" avoids local timezone shifting
    p_ev61_avg_int->GetXaxis()->SetTimeFormat("%m-%d-%Y");
    // p_ev61_avg_int->GetXaxis()->SetLabelSize(0.03);
    // p_ev61_avg_int->GetXaxis()->LabelsOption("v");  // vertical labels, helps avoid crowding
    // p_ev61_avg_int->GetXaxis()->SetLimits(2600, 3600);
    // p_ev61_avg_int->SetMinimum(50);
    // p_ev61_avg_int->SetMaximum(75);
    gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots/p_ev61_avg_int.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/p_ev61_avg_int.png"));

    return 0;

}