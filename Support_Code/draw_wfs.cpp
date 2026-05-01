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

int main(int argc, char *argv[])
{
    int run {0};
    int event {0};
    if (argc == 3) {
        sscanf(argv[1], "%i", &run);
        sscanf(argv[2], "%i", &event);
    } else if(argc == 2){ 
        sscanf(argv[1], "%i", &run);
        event = run;
    } else {
        cout << "Usase: "<< argv[0] <<" [run]" << " [last run] (optional)" << endl;
        return -1; 
    }   
    cout << "\n" << "Run: " << run << " Event: " << event << endl;

    int ADCSIZE = 45;
    TH1D *h_wf = new TH1D("h_wf", "Waveform", ADCSIZE, 0, ADCSIZE);

    // int data_num = 9; int vers_num = 4;

    int data_num = 41; int vers_num = 5;

    bool old_channel_map = false;

    if (vers_num < 5) {old_channel_map = true;}

    else if (vers_num == 5) {old_channel_map = false;}

    TFile *f; 
    //your root file location here
///        if(gSystem->AccessPathName(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v4.root", iRun))){
    if (gSystem->AccessPathName(Form("/data%i/coherent/data/d2o/processedData/run%i_processed_v%i.root", data_num, run, vers_num))) {
///        if(gSystem->AccessPathName(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/run%i_processed_v4.root", iRun))){
        cout << "Could not open file" << endl;
        return -1; 
    } else{
///            f = new TFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v4.root", iRun));
        f = new TFile(Form("/data%i/coherent/data/d2o/processedData/run%i_processed_v%i.root", data_num, run, vers_num));
///            f = new TFile(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/run%i_processed_v4.root", iRun));
    }

    TTree* t = (TTree*)f->Get("tree");

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

    Long64_t tentry = t->LoadTree(event);           // event - 1

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

    for (int iChan = 0; iChan < 26; iChan++) {

        if (old_channel_map) {if (iChan == 12 || iChan == 13 || iChan == 14) {continue;}}

        else if (!old_channel_map) {if (iChan == 23 || iChan == 24 || iChan == 25) {continue;}}

        for (int i = 0; i < ADCSIZE; i++) {

            h_wf->SetBinContent(i + 1, adcVal[iChan][i] - baselineMean[iChan]);

        }

        h_wf->GetXaxis()->SetTitle("Time (1 bin = 16 ns)");
        h_wf->GetYaxis()->SetTitle("Amplitude (ADC)");
        // h_wf->SetMaximum(320);
        // h_wf->SetMinimum(-20);
        h_wf->SetLineColor(kBlack);
        h_wf->SetLineWidth(3);
        h_wf->Draw();
        // if (run == 17006) {h_wf->SetLineColor(kGreen);}
        // else if (run == 17032) {h_wf->SetLineColor(kRed);}
        // else if (run == 17085) {h_wf->SetLineColor(kBlue);}
        // else if (run == 17146) {h_wf->SetLineColor(kBlack);}
        // else if (run == 17172) {h_wf->SetLineColor(kAzure);}
        // else if (run == 17197) {h_wf->SetLineColor(kViolet);}
        gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/wf_image/wf_run_%i_event_%i_ch_%i.png", run, event, iChan));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/wf_image/wf_run_%i_event_%i_ch_%i.png", run, event, iChan));
        h_wf->Reset();

    }

    f->Close();

    cout << "\n" << "Run, Event Number & triggerBits Value = " << run << " " << event << " " << triggerBits << endl;

    cout << "\n" << "End of code." << endl;
    
}
