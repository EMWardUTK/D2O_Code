#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TF1.h>
#include "TH1.h"
// RooFit includes (replaces Fit/Fitter.h and Math/Functor.h)
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooMsgService.h"
#include "RooGaussian.h"
#include "RooProfileLL.h"
#include "RooPlot.h"
#include "RooAbsReal.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <filesystem>

using std::cout; 
using std::endl;
using std::string;
using std::ifstream;

// Global pointers to histograms
TH2D* g_hist_vD = nullptr;   // MC: Neutrino-Deuterium
TH2D* g_hist_vO = nullptr;   // MC: Neutrino-Oxygen
TH2D* g_hist_BRN = nullptr;  // MC: Beam-Related Neutrons
TH2D* g_hist_vFe = nullptr;  // MC: Neutrino-Iron
TH2D* g_hist_vPb = nullptr;  // MC: Neutrino-Lead
TH2D* g_hist_SSB = nullptr;  // RD: Steady-State Background
TH2D* g_hist_UB = nullptr;   // RD: Unblinded Beam Spill Data
TH2D* g_target = nullptr;    // Pseudo-Real Data: MC & SSB Events

//Added separate function for loading simulation files. Returns an efficiency of applying
// TODO: Check this is very similar to what Justin was doing
// - nPMT cut
// - Min/max energy cut
// - Integration time
// - Whether events show up in the veto panel
double create2DHist(TH2D* hist,std::string simFilename, std::string timingFilename,
                    double integration_time = 416,int minPMTs_req=10,int minHits_req=30,int maxHits_req=500) {
    /////////////////////////////////
    //Get timing distribution first//
    /////////////////////////////////
    string line;
    Long64_t nTimingEvents = 0;
    double minTime_ns = hist->GetXaxis()->GetXmin();                        // THIS LINE IS NEW
    double maxTime_ns = hist->GetXaxis()->GetXmax();
    double binTime_ns = static_cast<int>(maxTime_ns - minTime_ns);          // THIS LINE IS NEW
    TH1D timingHist("timingHist", "",binTime_ns,minTime_ns,maxTime_ns);     // THIS LINE IS CHANGED
    ifstream txt_file(static_cast<TString>(timingFilename));
    std::vector<double> timing_events = {};
    while (getline(txt_file, line)) {
        timingHist.Fill(std::stoi(line));
        nTimingEvents++;
    }
    txt_file.close();

    ///////////////////////////////
    //Now get energy and fill TH2//
    ///////////////////////////////
    TFile* inpFile = new TFile(static_cast<TString>(simFilename),"READ");
    TTree* tree = (TTree*)inpFile->Get("Sim_Tree");
    long nEntries = tree->GetEntries();

    TTreeReader reader(tree);
    TTreeReaderValue<Int_t> eventNumber(reader, "eventNumber");
    TTreeReaderArray<Int_t> pmtNum(reader, "pmtHits.pmtNum");
    TTreeReaderArray<Double_t> eventTime(reader, "pmtHits.eventTime");
    TTreeReaderValue<Bool_t> veto_tag(reader, "veto_tag");

    Long64_t nValid=0;
    while (reader.Next()) {
        if ((*eventNumber) % 50000 == 0) {
            std::cout<<"On event "<<(*eventNumber)<<" of "<<nEntries<<std::endl;
        }

        int nHits=0; //Stores the number of PEs (total) for this event
        vector<int> hit_pmts = {}; //Stores a unique list of the hit PMTs within this event
        vector<double> hit_times = {}; //Stores the individual hit times of the event, only used for finding the first event time

        if (*veto_tag) continue;

        //Step through hits
        for (int ihit=0;ihit<eventTime.GetSize();ihit++) {
            //Check the hit occurs within integration time of the simulation start
            if (eventTime[ihit]<integration_time) {
                nHits++;
                hit_times.push_back(eventTime[ihit]);
                //Check if we have already registered this PMT in our hit_pmts vector
                if (std::find(hit_pmts.begin(), hit_pmts.end(), pmtNum[ihit]) == hit_pmts.end()){
                    hit_pmts.push_back(pmtNum[ihit]);
                }
            }
        }
        
        //Apply our cuts
        if ((hit_pmts.size() >= minPMTs_req) && (nHits >= minHits_req) && (nHits < maxHits_req)) {
            //If this passes the cuts, get our first hit time to add to neutrino timing distribution
            sort(hit_times.begin(),hit_times.end());
            double startTime = hit_times.at(0);

            int hitBin = hist->GetYaxis()->FindBin(nHits);

            for (int timingBin = 1; timingBin <= timingHist.GetNbinsX(); timingBin++) {
                double count = timingHist.GetBinContent(timingBin);
                if (count == 0) continue;

                double shiftedTime = timingHist.GetXaxis()->GetBinCenter(timingBin) + startTime;

                if (shiftedTime >= maxTime_ns) continue;

                int timeBin = hist->GetXaxis()->FindBin(shiftedTime);

                hist->AddBinContent(hist->GetBin(timeBin, hitBin),count);

                nValid += count;
            }
        }
    }
    double efficiency = static_cast<float>(nValid)/static_cast<float>(nEntries*nTimingEvents);
    inpFile->Close();
    return efficiency;
}

int fillSSBHist(TH2D* hist,std::string fname,int minHits_req=30,int maxHits_req=500) {
    double time;
    double energy;

    std::ifstream file(static_cast<TString>(fname));

    int nEntries = 0;
    while (file >> time >> energy) {
        if (time >= hist->GetXaxis()->GetXmin() && time < hist->GetXaxis()->GetXmax() && energy >= minHits_req && energy < maxHits_req) {
            hist->Fill(time, energy);  // x = energy, y = time
            nEntries++;
        }
    }
    return nEntries;
}

void main_fitting_asimov() {

    bool regenerate_pdfs = false; // Set to true to generate PDFs. Set to false to load from file.
    bool regenerate_plls = false; // Set to true to generate PLLs. Set to false to skip generation.

    // RooFit always uses NLL (extended likelihood) - fit_method kept for
    // compatibility with plotting/labeling code inherited from the original
    // chi2-vs-NLL version of this script.
    int fit_method = 2;

    // Define dimensions of all 2D histograms
    int tbins = 40;
    int tmin = 0;
    int tmax = 10000;
    int ebins = 40;
    int emin = 30;
    int emax = 500;

    int minPMTs_req = 10;
    double integration_time_ns = 416;

    // Define dimensions of event rate histograms
    int vDbins = 60;
    int vDmin = 300;
    int vDmax = 900;
    int vObins = 50;
    int vOmin = -100;
    int vOmax = 400;
    int BRNbins = 40;
    int BRNmin = 0;
    int BRNmax = 400;
    int vFebins = 60;
    int vFemin = -100;
    int vFemax = 500;
    int vPbbins = 40;
    int vPbmin = 0;
    int vPbmax = 400;
    int SSBbins = 60;
    int SSBmin = 1400;
    int SSBmax = 2000;
    int NLLbins = 100;
    int NLLmin = -20;
    int NLLmax = -10;

    // Expected event ranges
    vDbins = 60; vDmin = 300; vDmax = 900; vObins = 70; vOmin = -200; vOmax = 500; BRNbins = 40; BRNmin = 200; BRNmax = 600; vFebins = 70; vFemin = -200; vFemax = 500; vPbbins = 60; vPbmin = 0; vPbmax = 600; SSBbins = 60; SSBmin = 3300; SSBmax = 3900; NLLbins = 100; NLLmin = 30000; NLLmax = 40000;

    // Unblinded event ranges
    vDbins= 60; vDmin = 400; vDmax = 1000; BRNbins = 40; BRNmin = 200; BRNmax = 2200; NLLbins = 100; NLLmin = 35000; NLLmax = 45000;

    double numscale_F24 = 0.47403734;       // F24: 3250 MWhrs of “golden” SNS beam time; S25: 3606 MWhrs of “golden” SNS beam time
    double numscale_S25 = 0.52596266;
    
    // Define number of events in target, without efficiencies. From Yuri's calculations
    double vDnum = 581.4; 
    double vOnum = 94.6;
    double BRNnum = 654; //This is an estimate, scaled to match the number reported in the original code with similar effs applied.
    double vFenum = 501;
    double vPbnum = 17244.5;

    double peak_spread_rms_eff = 0.998;     //Afterglow cut, so equivalent to a dead time i.e. uniform efficiency
    double time_integration_eff = 0.987;    //Only impacts neutrino PDFs
    double dead_time_eff = 0.996;           //Uniform

    double vD_sim_eff;
    double vO_sim_eff;
    double vFe_sim_eff;
    double vPb_sim_eff;
    double BRN_sim_eff;
    int nSSB_x100;
    int nUB;

    if (regenerate_pdfs) {
        g_hist_vD = new TH2D("h_vD", "Neutrino-Deuterium Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);
        g_hist_vO = new TH2D("h_vO", "Neutrino-Oxygen Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);
        g_hist_BRN = new TH2D("h_BRN", "Beam-Related Neutrons Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);
        g_hist_vFe = new TH2D("h_vFe", "Neutrino-Iron Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);
        g_hist_vPb = new TH2D("h_vPb", "Neutrino-Lead Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);
        g_hist_SSB = new TH2D("h_SSB", "Out-of-Beam-Window Steady State Background Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);
        g_hist_UB = new TH2D("h_UB", "Unblinded Beam Spill Data Time vs Energy Distribution", tbins, tmin, tmax, ebins, emin, emax);
        
        vD_sim_eff = create2DHist(g_hist_vD,"sims/vD_Energy.root","sims/vD_Time.txt",integration_time_ns,minPMTs_req,emin,emax);
        vO_sim_eff = create2DHist(g_hist_vO,"sims/vO_Energy.root","sims/vO_Time.txt",integration_time_ns,minPMTs_req,emin,emax);
        vFe_sim_eff = create2DHist(g_hist_vFe,"sims/vFe_Energy.root","sims/vFe_Time.txt",integration_time_ns,minPMTs_req,emin,emax);
        vPb_sim_eff = create2DHist(g_hist_vPb,"sims/vPb_Energy.root","sims/vPb_Time.txt",integration_time_ns,minPMTs_req,emin,emax);
        BRN_sim_eff = create2DHist(g_hist_BRN,"sims/BRN_Energy.root","sims/BRN_Time.txt",integration_time_ns,minPMTs_req,emin,emax);
        nSSB_x100 = fillSSBHist(g_hist_SSB,"data/SSB_Time_And_Energy.txt",emin,emax);
        nUB = fillSSBHist(g_hist_UB,"data/Unblinded_BeamSpill_Time_And_Energy.txt",emin,emax);

        TFile* outFile = new TFile("pdfs.root","RECREATE");
        //We are storing the simulation efficiency as the integral to preserve this number in the outout
        g_hist_vD->Scale(vD_sim_eff/g_hist_vD->Integral());
        g_hist_vO->Scale(vO_sim_eff/g_hist_vO->Integral());
        g_hist_vFe->Scale(vFe_sim_eff/g_hist_vFe->Integral());
        g_hist_vPb->Scale(vPb_sim_eff/g_hist_vPb->Integral());
        g_hist_BRN->Scale(BRN_sim_eff/g_hist_BRN->Integral());

        g_hist_vD->Write();
        g_hist_vO->Write();
        g_hist_vFe->Write();
        g_hist_vPb->Write();
        g_hist_BRN->Write();
        g_hist_SSB->Write();
        g_hist_UB->Write();
        outFile->Close();
    }
    else {
        TFile* inFile = new TFile("pdfs.root","READ");
        g_hist_vD = (TH2D*)inFile->Get("h_vD");
        g_hist_vO = (TH2D*)inFile->Get("h_vO");
        g_hist_vFe = (TH2D*)inFile->Get("h_vFe");
        g_hist_vPb = (TH2D*)inFile->Get("h_vPb");
        g_hist_BRN = (TH2D*)inFile->Get("h_BRN");
        g_hist_SSB = (TH2D*)inFile->Get("h_SSB");
        g_hist_UB = (TH2D*)inFile->Get("h_UB");

        vD_sim_eff = g_hist_vD->Integral();
        vO_sim_eff = g_hist_vO->Integral();
        vFe_sim_eff = g_hist_vFe->Integral();
        vPb_sim_eff = g_hist_vPb->Integral();
        BRN_sim_eff = g_hist_BRN->Integral();
        nSSB_x100 = g_hist_SSB->Integral();
        nUB = g_hist_UB->Integral();
    }

    double SSBnum = 0.01*nSSB_x100;
    std::cout<<"vD simulation eff (not including fixed effs) is "<<vD_sim_eff<<std::endl;
    std::cout<<"vO simulation eff (not including fixed effs) is "<<vO_sim_eff<<std::endl;
    std::cout<<"vPb simulation eff (not including fixed effs) is "<<vPb_sim_eff<<std::endl;
    std::cout<<"vFe simulation eff (not including fixed effs) is "<<vFe_sim_eff<<std::endl;
    std::cout<<"BRN simulation ef  (not including fixed effs) is "<<BRN_sim_eff<<std::endl;
    std::cout<<"Number of SSB events in time window (x100) is "<<nSSB_x100<<std::endl;
    
    //Include fixed efficiencies
    vDnum *= vD_sim_eff*peak_spread_rms_eff*time_integration_eff*dead_time_eff;
    vOnum *= vO_sim_eff*peak_spread_rms_eff*time_integration_eff*dead_time_eff;
    vFenum *= vFe_sim_eff*peak_spread_rms_eff*time_integration_eff*dead_time_eff;
    vPbnum *= vPb_sim_eff*peak_spread_rms_eff*time_integration_eff*dead_time_eff;
    BRNnum *= BRN_sim_eff*peak_spread_rms_eff*dead_time_eff; 
    double total_num = vDnum + vOnum + vFenum + vPbnum + BRNnum + SSBnum;

    std::cout<<"Nominal number of vD events is "<<vDnum<<std::endl;
    std::cout<<"Nominal number of vO events is "<<vOnum<<std::endl;
    std::cout<<"Nominal number of vFe events is "<<vFenum<<std::endl;
    std::cout<<"Nominal number of vPb events is "<<vPbnum<<std::endl;
    std::cout<<"Nominal number of BRN events is "<<BRNnum<<std::endl;
    std::cout<<"Nominal number of SSB events is "<<SSBnum<<std::endl;
    std::cout<<"Nominal number of UB events is "<<nUB<<std::endl;

    ///////////////////////////////////
    //For plotting output of toy fits//
    ///////////////////////////////////
    TH2D* h_fitted_2D = new TH2D("h_fitted_2D", "Fitted MC (weighted sum)", tbins, tmin, tmax, ebins, emin, emax);

    // g_target holds the generated toy dataset for the current loop iteration,
    // as a TH2D, so that the (restored) result-extraction and plotting code
    // below -- written against a TH2D "target" histogram -- continues to work.
    g_target = new TH2D("h_target", "Beam Spill Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);

    TH1D* h_vD_num = new TH1D("h_vD_num", "Number of vD Events Found", vDbins, vDmin, vDmax);
    TH1D* h_vO_num = new TH1D("h_vO_num", "Number of vO Events Found", vObins, vOmin, vOmax);
    TH1D* h_BRN_num = new TH1D("h_BRN_num", "Number of BRN Events Found", BRNbins, BRNmin, BRNmax);
    TH1D* h_vFe_num = new TH1D("h_vFe_num", "Number of vFe Events Found", vFebins, vFemin, vFemax);
    TH1D* h_vPb_num = new TH1D("h_vPb_num", "Number of vPb Events Found", vPbbins, vPbmin, vPbmax);
    TH1D* h_SSB_num = new TH1D("h_SSB_num", "Number of SSB Events Found", SSBbins, SSBmin, SSBmax);

    TH2D* h_vD_vO = new TH2D("h_vD_vO", "vD vs vO Recovered Event Rates", vDbins, vDmin, vDmax, vObins, vOmin, vOmax);
    TH2D* h_vD_vFe = new TH2D("h_vD_vFe", "vD vs vFe Recovered Event Rates", vDbins, vDmin, vDmax, vFebins, vFemin, vFemax);
    TH2D* h_vD_vPb = new TH2D("h_vD_vPb", "vD vs vPb Recovered Event Rates", vDbins, vDmin, vDmax, vPbbins, vPbmin, vPbmax);
    TH2D* h_vO_vFe = new TH2D("h_vO_vFe", "vO vs vFe Recovered Event Rates", vObins, vOmin, vOmax, vFebins, vFemin, vFemax);
    TH2D* h_vO_vPb = new TH2D("h_vO_vPb", "vO vs vPb Recovered Event Rates", vObins, vOmin, vOmax, vPbbins, vPbmin, vPbmax);
    TH2D* h_vFe_vPb = new TH2D("h_vFe_vPb", "vFe vs vPb Recovered Event Rates", vFebins, vFemin, vFemax, vPbbins, vPbmin, vPbmax);

    TH1D* h_chi2_per_ndf = new TH1D("h_chi2_per_ndf", "Negative Log-Likelihood per Bins Used", NLLbins, NLLmin, NLLmax);

    TH2D* h_background_true = new TH2D("h_background_true", "Out-of-Beam-Window Steady State Background Data", tbins, tmin, tmax, ebins, emin, emax);
    TH1D* h_bg_stats = new TH1D("h_bg_stats", "Number of Background Events Put into Target", SSBbins, SSBmin, SSBmax);
    TH2D* h_bg_vO_trend = new TH2D("h_bg_vO_trend", "Number of Background Events in Target vs Recovered Number of vO Events", SSBbins, SSBmin, SSBmax, vObins, vOmin, vOmax);
    TH1D* h_bg_energy = new TH1D("h_bg_energy", "Out-of-Beam-Window Steady State Background Energy", ebins, emin, emax);
    TH1D* h_sg_energy = new TH1D("h_sg_energy", "In-Beam Window Signal Energy", ebins, emin, emax);

    std::filesystem::create_directory("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/RooFit_Code/Fitting_Results/Plots");

    /////////////////
    //RooFit set-up//
    /////////////////
    // Suppress RooFit output - only show errors
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    
    //RooRealVars defined outside the loop
    RooRealVar roo_time("roo_time", "Time", tmin, tmax);
    RooRealVar roo_energy("roo_energy", "Energy", emin, emax);
    roo_time.setBins(tbins); 
    roo_energy.setBins(ebins);
    RooArgSet obs_set(roo_time, roo_energy);

    // Create RooDataHist for each MC template from TH2D
    RooDataHist rdh_vD ("rdh_vD",  "vD template",  RooArgList(roo_time, roo_energy), g_hist_vD);
    RooDataHist rdh_vO ("rdh_vO",  "vO template",  RooArgList(roo_time, roo_energy), g_hist_vO);
    RooDataHist rdh_BRN("rdh_BRN", "BRN template", RooArgList(roo_time, roo_energy), g_hist_BRN);
    RooDataHist rdh_vFe("rdh_vFe", "vFe template", RooArgList(roo_time, roo_energy), g_hist_vFe);
    RooDataHist rdh_vPb("rdh_vPb", "vPb template", RooArgList(roo_time, roo_energy), g_hist_vPb);
    RooDataHist rdh_SSB("rdh_SSB", "SSB template", RooArgList(roo_time, roo_energy), g_hist_SSB);
    RooDataHist rdh_UB ("rdh_UB", "Unblinded Data", RooArgList(roo_time, roo_energy), g_hist_UB);

    // Create a shape PDF for each template (normalized to unit area internally)
    RooHistPdf pdf_vD ("pdf_vD",  "vD PDF",  obs_set, rdh_vD,0);
    RooHistPdf pdf_vO ("pdf_vO",  "vO PDF",  obs_set, rdh_vO,0);
    RooHistPdf pdf_BRN("pdf_BRN", "BRN PDF", obs_set, rdh_BRN,0);
    RooHistPdf pdf_vFe("pdf_vFe", "vFe PDF", obs_set, rdh_vFe,0);
    RooHistPdf pdf_vPb("pdf_vPb", "vPb PDF", obs_set, rdh_vPb,0);
    RooHistPdf pdf_SSB("pdf_SSB", "SSB PDF", obs_set, rdh_SSB,0);

    RooRealVar N_vD ("N_vD",  "vD yield",  vDnum,  0.0, 1e8);
    RooRealVar N_vO ("N_vO",  "vO yield",  vOnum,  0.0, 1e8);
    RooRealVar N_BRN("N_BRN", "BRN yield", BRNnum, 0.0, 1e8);
    RooRealVar N_vFe("N_vFe", "vFe yield", vFenum, 0.0, 1e8); 
    RooRealVar N_vPb("N_vPb", "vPb yield", vPbnum, 0.0, 1e8);
    RooRealVar N_SSB("N_SSB", "SSB yield", SSBnum, 0.0, 1e8);

    // Build extended sum PDF: model = N_vD*pdf_vD + N_vO*pdf_vO + ...
    // Passing N yields (not N-1 fractions) makes this an extended PDF
    RooAddPdf model("model", "Signal + Background Model",
        RooArgList(pdf_vD, pdf_vO, pdf_vFe, pdf_vPb, pdf_BRN, pdf_SSB),
        RooArgList(N_vD, N_vO, N_vFe, N_vPb, N_BRN, N_SSB));

    //Random seed set outside the loop
    RooRandom::randomGenerator()->SetSeed(12345);

    ///////////////
    //Constraints//
    ///////////////
    double ssb_mean  = nSSB_x100 / 100.0;
    double ssb_sigma = std::sqrt(nSSB_x100) / 100.0;
    std::cout<<ssb_mean<<","<<ssb_sigma<<std::endl;
    RooConstVar N_SSB_mean("N_SSB_mean","",ssb_mean);
    RooConstVar N_SSB_sigma("N_SSB_sigma","",ssb_sigma);
    RooGaussian SSB_constraint("SSB_constraint","",N_SSB,N_SSB_mean,N_SSB_sigma);
    
    RooConstVar N_vFe_mean("N_vFe_mean","",vFenum);
    RooConstVar N_vFe_sigma("N_vFe_sigma","",vFenum*0.37);
    RooGaussian vFe_constraint("vFe_constraint","",N_vFe,N_vFe_mean,N_vFe_sigma);

    // RooConstVar N_vO_mean("N_vO_mean","",vOnum);
    // RooConstVar N_vO_sigma("N_vO_sigma","",vOnum*1.0);
    // RooGaussian vO_constraint("vO_constraint","",N_vO,N_vO_mean,N_vO_sigma);

    RooProdPdf* model_c = new RooProdPdf("model_c","model_c",RooArgList(model,SSB_constraint,vFe_constraint));

    ///////////////
    // MAIN LOOP //
    ///////////////
    // Best fit values as starting points for NLL
    int max_loop = 1;
    int min_loop = 0;
    int sample_loop = (max_loop - 1);
    int nFailed=0;
    for (int iLoop = min_loop; iLoop < max_loop; iLoop++) {

        cout << "\n========================================" << endl;
        cout << "================ " << "LOOP " << iLoop + 1 << " ================" << endl;
        cout << "========================================" << endl;

        //Reset amplitudes to initial guesses each time
        N_vD.setVal(vDnum);
        N_vO.setVal(vOnum);
        N_vFe.setVal(vFenum);
        N_vPb.setVal(vPbnum);
        N_BRN.setVal(BRNnum);
        N_SSB.setVal(SSBnum);
            
        //Generate toy data set with poisson statistics
        //auto* data = model.generate(obs_set,total_num); //Causing an error?
        auto* data = model.generateBinned(obs_set,
                                          RooFit::NumEvents(static_cast<int>(std::round(total_num))),
                                          RooFit::Extended(true));

        // Populate g_target (a TH2D) from the generated toy dataset so that
        // downstream result-extraction and plotting code (written against a
        // TH2D "target" histogram) continues to work unmodified.
        g_target->Reset();
        for (int ibin = 0; ibin < data->numEntries(); ibin++) {
            data->get(ibin);
            g_target->Fill(roo_time.getVal(), roo_energy.getVal(), data->weight());
        }

        // Perform extended NLL fit with pull terms and MINOS errors
        // Extended(kTRUE):           Poisson constraint on total yield
        // ExternalConstraints(...):  pull terms on vFe and SSB components
        // Minos(RooArgSet(...)):     profile likelihood errors on vD, vO, vPb only
        // Save(kTRUE):               returns RooFitResult with full covariance
        // PrintLevel(-1):            suppresses MINUIT printout
        RooFitResult* fitResult = model_c->fitTo(
            rdh_UB,         // *data
            RooFit::Extended(kTRUE),
            RooFit::Save(kTRUE),
            RooFit::PrintLevel(-1),
            RooFit::Strategy(2),
            RooFit::Minos(RooArgSet(N_vD)),
            RooFit::Constrain(RooArgSet(N_SSB,N_vFe)),
            RooFit::MaxCalls(6*100000)
        );

        fitResult->Print("v");
        fitResult->correlationMatrix().Print();

        bool fit_success = (fitResult->status() == 0);

        if (!fit_success) {
            cout << "\nERROR: Fit failed! Status: " << fitResult->status() << endl;
            for (int i = 0; i < fitResult->numStatusHistory(); ++i) {
                std::cout << fitResult->statusLabelHistory(i) << " = " << fitResult->statusCodeHistory(i) << '\n';   
            };
            delete fitResult;
            nFailed++;
            continue;
        }
        else {
        }

        // ===================================================================
        // EXTRACT RESULTS
        // Convert RooFit yields back to weights for compatibility with
        // all downstream code: weight_i = N_i / g_hist_i->Integral()
        // This is the inverse of: N_i = weight_i * g_hist_i->Integral()
        // ===================================================================
        double param[6], param_err[6];
        param[0]     = N_vD.getVal()   / g_hist_vD->Integral();
        param[1]     = N_vO.getVal()   / g_hist_vO->Integral();
        param[2]     = N_BRN.getVal()  / g_hist_BRN->Integral();
        param[3]     = N_vFe.getVal()  / g_hist_vFe->Integral();
        param[4]     = N_vPb.getVal()  / g_hist_vPb->Integral();
        param[5]     = N_SSB.getVal()  / g_hist_SSB->Integral();
        // HESSE errors: symmetric, from second derivative of NLL at minimum
        param_err[0] = N_vD.getError()  / g_hist_vD->Integral();
        param_err[1] = N_vO.getError()  / g_hist_vO->Integral();
        param_err[2] = N_BRN.getError() / g_hist_BRN->Integral();
        param_err[3] = N_vFe.getError() / g_hist_vFe->Integral();
        param_err[4] = N_vPb.getError() / g_hist_vPb->Integral();
        param_err[5] = N_SSB.getError() / g_hist_SSB->Integral();

        // MINOS asymmetric errors for vD, vO, vPb (profile likelihood intervals)
        // getErrorHi() = upper (positive) deviation from best-fit value
        // getErrorLo() = lower (negative) deviation from best-fit value
        // These are more statistically rigorous than HESSE for low-statistics fits
        // double minos_vD_hi  = N_vD.getErrorHi();  // / g_hist_vD->Integral();
        // double minos_vD_lo  = N_vD.getErrorLo();  // / g_hist_vD->Integral();
        // double minos_vO_hi  = N_vO.getErrorHi();  // / g_hist_vO->Integral();
        // double minos_vO_lo  = N_vO.getErrorLo();  // / g_hist_vO->Integral();
        // double minos_vPb_hi = N_vPb.getErrorHi(); // / g_hist_vPb->Integral();
        // double minos_vPb_lo = N_vPb.getErrorLo(); // / g_hist_vPb->Integral();

        // minNll() returns the minimum negative log-likelihood found
        double minFcnValue = fitResult->minNll();
        // covQual() >= 2 means covariance matrix is at least approximately correct
        bool fitIsValid = (fitResult->status() == 0 && fitResult->covQual() >= 2);

        delete fitResult;

        cout << "\nFitted event rates (HESSE symmetric errors):" << endl;
        cout << "  events_vD  = " << N_vD.getVal() << " +/- " << N_vD.getError() << endl;
        cout << "  events_vO  = " << N_vO.getVal() << " +/- " << N_vO.getError() << endl;
        cout << "  events_BRN = " << N_BRN.getVal() << " +/- " << N_BRN.getError() << endl;
        cout << "  events_vFe = " << N_vFe.getVal() << " +/- " << N_vFe.getError() << endl;
        cout << "  events_vPb = " << N_vPb.getVal() << " +/- " << N_vPb.getError() << endl;
        cout << "  events_SSB = " << N_SSB.getVal() << " +/- " << N_SSB.getError() << endl;

        // cout << "\nFitted event rates (MINOS asymmetric errors for vD, vO, vPb):" << endl;
        // cout << "  events_vD  = " << N_vD.getVal() << " +" << minos_vD_hi  << " / " << minos_vD_lo  << endl;
        // cout << "  events_vO  = " << N_vO.getVal() << " +" << minos_vO_hi  << " / " << minos_vO_lo  << endl;
        // cout << "  events_vPb = " << N_vPb.getVal() << " +" << minos_vPb_hi << " / " << minos_vPb_lo << endl;

        h_vD_num->Fill(param[0] * g_hist_vD->Integral());
        h_vO_num->Fill(param[1] * g_hist_vO->Integral());
        h_BRN_num->Fill(param[2] * g_hist_BRN->Integral());
        h_vFe_num->Fill(param[3] * g_hist_vFe->Integral());
        h_vPb_num->Fill(param[4] * g_hist_vPb->Integral());
        h_SSB_num->Fill(param[5] * g_hist_SSB->Integral());

        h_vD_vO->Fill(param[0] * g_hist_vD->Integral(), param[1] * g_hist_vO->Integral());
        h_vD_vFe->Fill(param[0] * g_hist_vD->Integral(), param[3] * g_hist_vFe->Integral());
        h_vD_vPb->Fill(param[0] * g_hist_vD->Integral(), param[4] * g_hist_vPb->Integral());
        h_vO_vFe->Fill(param[1] * g_hist_vO->Integral(), param[3] * g_hist_vFe->Integral());
        h_vO_vPb->Fill(param[1] * g_hist_vO->Integral(), param[4] * g_hist_vPb->Integral());
        h_vFe_vPb->Fill(param[3] * g_hist_vFe->Integral(), param[4] * g_hist_vPb->Integral());

        // h_bg_vO_trend depended on "count", a per-event tally of how many
        // true background events landed in the toy target under the old
        // manual event-splitting generation scheme. generateBinned() draws
        // Poisson/multinomial bin counts directly from the total model and
        // does not track which component each generated count came from, so
        // this quantity no longer exists and this fill cannot be restored
        // as-is. Left disabled; see also canvas c7 below.
        // h_bg_vO_trend->Fill(count, param[1] * g_hist_vO->Integral());
        
        cout << "\n========================================" << endl;
        cout << "Fit Results" << endl;
        cout << "========================================" << endl;
        cout << "Fit status: " << (fitIsValid ? "VALID" : "INVALID") << endl;
        cout << "Chi-squared/Negative Log-Likelihood: " << minFcnValue << endl;
        
        int ndf = 0;
        int nbins_used = 0;
        if (fit_method == 1) {
            // Calculate number of degrees of freedom (count 2D bins with data)
            for (int ix = 1; ix <= g_target->GetNbinsX(); ++ix) {
                for (int iy = 1; iy <= g_target->GetNbinsY(); ++iy) {
                    if (g_target->GetBinError(ix, iy) > 0) ndf++;
                }
            }
            ndf -= 6;  // Subtract number of fit parameters
            
            cout << "NDF: " << ndf << endl;
            cout << "Chi2/NDF: " << minFcnValue / ndf << endl;

            h_chi2_per_ndf->Fill(minFcnValue / ndf);
        }
        else if (fit_method == 2) {
            // For likelihood fits, we can't directly compute chi2/NDF the same way
            // But we can report the number of bins used
            for (int ix = 1; ix <= g_target->GetNbinsX(); ++ix) {
                for (int iy = 1; iy <= g_target->GetNbinsY(); ++iy) {
                    if (g_target->GetBinContent(ix, iy) > 0) nbins_used++;
                }
            }
            
            cout << "Number of bins with data: " << nbins_used << endl;
            cout << "Number of parameters: 6" << endl;
            cout << "Effective NDF: " << (nbins_used - 6) << endl;
            cout << "NLL per bin: " << minFcnValue / nbins_used << endl;

            h_chi2_per_ndf->Fill(minFcnValue / nbins_used);
        }
        
        // Create the fitted MC histogram (weighted sum) - 2D
        for (int ix = 1; ix <= h_fitted_2D->GetNbinsX(); ++ix) {
            for (int iy = 1; iy <= h_fitted_2D->GetNbinsY(); ++iy) {
                double val = param[0] * g_hist_vD->GetBinContent(ix, iy) +
                             param[1] * g_hist_vO->GetBinContent(ix, iy) +
                             param[2] * g_hist_BRN->GetBinContent(ix, iy) +
                             param[3] * g_hist_vFe->GetBinContent(ix, iy) +
                             param[4] * g_hist_vPb->GetBinContent(ix, iy) +
                             param[5] * g_hist_SSB->GetBinContent(ix, iy);
                h_fitted_2D->SetBinContent(ix, iy, val);
            }
        }

        if (iLoop == sample_loop) {
            // ===================================================================
            // PLOTTING
            // ===================================================================
            cout << "\n========================================" << endl;
            cout << "Creating plots..." << endl;
            cout << "========================================" << endl;

            TCanvas* c2 = new TCanvas("c2", "2D Event Rates", 1600, 1200);
            c2->Divide(4, 2);

            // ===================================================================
            // Plot 1: Target
            // ===================================================================
            c2->cd(1);
            TH2D *h_target_copy = (TH2D*)g_target->Clone("h_target_copy");
            h_target_copy->SetTitle("Target");
            h_target_copy->GetXaxis()->SetTitle("Time");
            h_target_copy->GetYaxis()->SetTitle("Energy");
            h_target_copy->Draw("COLZ");
            
            // ===================================================================
            // Plot 2: vD
            // ===================================================================
            c2->cd(2);
            TH2D *h_vD_copy = (TH2D*)g_hist_vD->Clone("h_vD_copy");
            h_vD_copy->SetTitle("vD");
            h_vD_copy->GetXaxis()->SetTitle("Time");
            h_vD_copy->GetYaxis()->SetTitle("Energy");
            h_vD_copy->Draw("COLZ");
            
            // ===================================================================
            // Plot 3: vO
            // ===================================================================
            c2->cd(3);
            TH2D *h_vO_copy = (TH2D*)g_hist_vO->Clone("h_vO_copy");
            h_vO_copy->SetTitle("vO");
            h_vO_copy->GetXaxis()->SetTitle("Time");
            h_vO_copy->GetYaxis()->SetTitle("Energy");
            h_vO_copy->Draw("COLZ");
            
            // ===================================================================
            // Plot 4: BRN
            // ===================================================================
            c2->cd(4);
            TH2D *h_BRN_copy = (TH2D*)g_hist_BRN->Clone("h_BRN_copy");
            h_BRN_copy->SetTitle("BRN");
            h_BRN_copy->GetXaxis()->SetTitle("Time");
            h_BRN_copy->GetYaxis()->SetTitle("Energy");
            h_BRN_copy->Draw("COLZ");

            // ===================================================================
            // Plot 5: vFe
            // ===================================================================
            c2->cd(5);
            TH2D *h_vFe_copy = (TH2D*)g_hist_vFe->Clone("h_vFe_copy");
            h_vFe_copy->SetTitle("vFe");
            h_vFe_copy->GetXaxis()->SetTitle("Time");
            h_vFe_copy->GetYaxis()->SetTitle("Energy");
            h_vFe_copy->Draw("COLZ");

            // ===================================================================
            // Plot 6: vPb
            // ===================================================================
            c2->cd(6);
            TH2D *h_vPb_copy = (TH2D*)g_hist_vPb->Clone("h_vPb_copy");
            h_vPb_copy->SetTitle("vPb");
            h_vPb_copy->GetXaxis()->SetTitle("Time");
            h_vPb_copy->GetYaxis()->SetTitle("Energy");
            h_vPb_copy->Draw("COLZ");

            // ===================================================================
            // Plot 7: SSB
            // ===================================================================
            c2->cd(7);
            TH2D *h_SSB_copy = (TH2D*)g_hist_SSB->Clone("h_SSB_copy");
            h_SSB_copy->SetTitle("SSB");
            h_SSB_copy->GetXaxis()->SetTitle("Time");
            h_SSB_copy->GetYaxis()->SetTitle("Energy");
            h_SSB_copy->Draw("COLZ");

            c2->Update();

            c2->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/RooFit_Code/Fitting_Results/Plots/2D_Event_Rates.png");
            
            TCanvas* c3 = new TCanvas("c3", "MC to Data Fit Results - 2D", 1600, 1200);
            c3->Divide(3, 2);
            
            // ===================================================================
            // Plot 1: Data 2D histogram
            // ===================================================================
            c3->cd(1);
            h_target_copy->SetTitle("Fake Data (Time vs Energy)");
            h_target_copy->GetXaxis()->SetTitle("Time");
            h_target_copy->GetYaxis()->SetTitle("Energy");
            h_target_copy->Draw("COLZ");
            
            // ===================================================================
            // Plot 2: Fitted MC 2D histogram
            // ===================================================================
            c3->cd(2);
            h_fitted_2D->SetTitle("Fitted MC (Time vs Energy)");
            h_fitted_2D->GetXaxis()->SetTitle("Time");
            h_fitted_2D->GetYaxis()->SetTitle("Energy");
            h_fitted_2D->Draw("COLZ");
            
            // ===================================================================
            // Plot 3: Residuals 2D
            // ===================================================================
            c3->cd(3);
            TH2D* h_residuals_2D = (TH2D*)g_target->Clone("h_residuals_2D");
            h_residuals_2D->Add(h_fitted_2D, -1);
            h_residuals_2D->SetTitle("Residuals (Data - MC)");
            h_residuals_2D->GetXaxis()->SetTitle("Time");
            h_residuals_2D->GetYaxis()->SetTitle("Energy");
            h_residuals_2D->Draw("COLZ");
            
            // ===================================================================
            // Plot 4: Energy projection comparison
            // ===================================================================
            c3->cd(4);
            gPad->SetLogy();
            
            // Project onto Y axis (Energy)
            TH1D* h_data_energy = g_target->ProjectionY("h_data_energy");
            TH1D* h_fitted_energy = h_fitted_2D->ProjectionY("h_fitted_energy");
            
            h_data_energy->SetMarkerStyle(20);
            h_data_energy->SetMarkerSize(0.8);
            h_data_energy->SetMarkerColor(kBlack);
            h_data_energy->SetLineColor(kBlack);
            h_data_energy->SetTitle("Energy Projection");
            h_data_energy->GetXaxis()->SetTitle("Energy");
            h_data_energy->GetYaxis()->SetTitle("Events");
            h_data_energy->Draw("E");
            
            h_fitted_energy->SetLineColor(kRed);
            h_fitted_energy->SetLineWidth(2);
            h_fitted_energy->Draw("HIST SAME");
            
            TLegend* leg1 = new TLegend(0.6, 0.75, 0.88, 0.88);
            leg1->AddEntry(h_data_energy, "Real Data", "lep");
            leg1->AddEntry(h_fitted_energy, "Fitted MC", "l");
            leg1->Draw();
            
            if (fit_method == 1) {
                // Add chi2/ndf text
                TLatex* latex1 = new TLatex();
                latex1->SetNDC();
                latex1->SetTextSize(0.04);
                latex1->DrawLatex(0.6, 0.70, Form("#chi^{2}/NDF = %.2f/%d", minFcnValue, ndf));
                latex1->DrawLatex(0.6, 0.65, Form("= %.2f", minFcnValue/ndf));
            }
            else if (fit_method == 2) {
                // Add NLL info
                TLatex* latex1 = new TLatex();
                latex1->SetNDC();
                latex1->SetTextSize(0.04);
                latex1->DrawLatex(0.6, 0.70, Form("-log(L) = %.1f", minFcnValue));
                latex1->DrawLatex(0.6, 0.65, Form("bins = %d", nbins_used));
            }
            
            // ===================================================================
            // Plot 5: Time projection comparison
            // ===================================================================
            c3->cd(5);
            gPad->SetLogy();
            
            // Project onto X axis (Time)
            TH1D* h_data_time = g_target->ProjectionX("h_data_time");
            TH1D* h_fitted_time = h_fitted_2D->ProjectionX("h_fitted_time");
            
            h_data_time->SetMarkerStyle(20);
            h_data_time->SetMarkerSize(0.8);
            h_data_time->SetMarkerColor(kBlack);
            h_data_time->SetLineColor(kBlack);
            h_data_time->SetTitle("Time Projection");
            h_data_time->GetXaxis()->SetTitle("Time");
            h_data_time->GetYaxis()->SetTitle("Events");
            h_data_time->Draw("E");
            
            h_fitted_time->SetLineColor(kRed);
            h_fitted_time->SetLineWidth(2);
            h_fitted_time->Draw("HIST SAME");
            
            TLegend* leg2 = new TLegend(0.6, 0.75, 0.88, 0.88);
            leg2->AddEntry(h_data_time, "Real Data", "lep");
            leg2->AddEntry(h_fitted_time, "Fitted MC", "l");
            leg2->Draw();
            
            // ===================================================================
            // Plot 6: Pull distribution
            // ===================================================================
            c3->cd(6);
            
            TH1D* h_pulls = new TH1D("h_pulls", "Pull Distribution", 50, -5, 5);
            if (fit_method == 1) {
                // Chi-squared pull calculation
                for (int ix = 1; ix <= g_target->GetNbinsX(); ++ix) {
                    for (int iy = 1; iy <= g_target->GetNbinsY(); ++iy) {
                        double data_val = g_target->GetBinContent(ix, iy);
                        double mc_val = h_fitted_2D->GetBinContent(ix, iy);
                        double error = g_target->GetBinError(ix, iy);
                        if (error > 0) {
                            double pull = (data_val - mc_val) / error;
                            h_pulls->Fill(pull);
                        }
                    }
                }
            }
            else if (fit_method == 2) {
                // Likelihood pull calculation
                for (int ix = 1; ix <= g_target->GetNbinsX(); ++ix) {
                    for (int iy = 1; iy <= g_target->GetNbinsY(); ++iy) {
                        double data_val = g_target->GetBinContent(ix, iy);
                        double mc_val = h_fitted_2D->GetBinContent(ix, iy);
                        
                        // For Poisson, error = sqrt(observed count)
                        if (data_val > 0) {
                            double error = TMath::Sqrt(data_val);
                            double pull = (data_val - mc_val) / error;
                            h_pulls->Fill(pull);
                        }
                    }
                }
            }
            
            h_pulls->GetXaxis()->SetTitle("Pull");
            h_pulls->GetYaxis()->SetTitle("Entries");
            h_pulls->SetFillColor(kAzure-9);
            h_pulls->Draw();
            
            // Fit pull distribution with Gaussian
            TF1* gauss = new TF1("gauss", "gaus", -5, 5);
            h_pulls->Fit(gauss, "Q");
            
            TLatex* latex2 = new TLatex();
            latex2->SetNDC();
            latex2->SetTextSize(0.04);
            latex2->DrawLatex(0.15, 0.85, Form("Mean = %.2f", gauss->GetParameter(1)));
            latex2->DrawLatex(0.15, 0.80, Form("Sigma = %.2f", gauss->GetParameter(2)));
            
            c3->Update();

            c3->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/RooFit_Code/Fitting_Results/Plots/MC_to_Data_Fit_Results-2D.png");
            
            // ===================================================================
            // Create a second canvas with MC components
            // ===================================================================
            TCanvas* c4 = new TCanvas("c4", "MC Components", 2400, 800);
            c4->Divide(7, 2);
            
            // Plot weighted MC components - Energy projections
            c4->cd(1);
            gPad->SetLogy();
            TH2D* h_vD_scaled = (TH2D*)g_hist_vD->Clone("h_vD_scaled");
            h_vD_scaled->Scale(param[0]);
            TH1D* h_vD_energy = h_vD_scaled->ProjectionY("h_vD_energy");
            h_vD_energy->SetLineColor(kRed);
            h_vD_energy->SetLineWidth(2);
            h_vD_energy->SetTitle(Form("#nu-D Energy (w=%.5f)", param[0]));
            h_vD_energy->GetXaxis()->SetTitle("Energy");
            h_vD_energy->Draw("HIST");
            
            c4->cd(2);
            gPad->SetLogy();
            TH2D* h_vO_scaled = (TH2D*)g_hist_vO->Clone("h_vO_scaled");
            h_vO_scaled->Scale(param[1]);
            TH1D* h_vO_energy = h_vO_scaled->ProjectionY("h_vO_energy");
            h_vO_energy->SetLineColor(kBlue);
            h_vO_energy->SetLineWidth(2);
            h_vO_energy->SetTitle(Form("#nu-O Energy (w=%.5f)", param[1]));
            h_vO_energy->GetXaxis()->SetTitle("Energy");
            h_vO_energy->Draw("HIST");
            
            c4->cd(3);
            gPad->SetLogy();
            TH2D* h_BRN_scaled = (TH2D*)g_hist_BRN->Clone("h_BRN_scaled");
            h_BRN_scaled->Scale(param[2]);
            TH1D* h_BRN_energy = h_BRN_scaled->ProjectionY("h_BRN_energy");
            h_BRN_energy->SetLineColor(kGreen+2);
            h_BRN_energy->SetLineWidth(2);
            h_BRN_energy->SetTitle(Form("BRN Energy (w=%.5f)", param[2]));
            h_BRN_energy->GetXaxis()->SetTitle("Energy");
            h_BRN_energy->Draw("HIST");

            c4->cd(4);
            gPad->SetLogy();
            TH2D* h_vFe_scaled = (TH2D*)g_hist_vFe->Clone("h_vFe_scaled");
            h_vFe_scaled->Scale(param[3]);
            TH1D* h_vFe_energy = h_vFe_scaled->ProjectionY("h_vFe_energy");
            h_vFe_energy->SetLineColor(kCyan+2);
            h_vFe_energy->SetLineWidth(2);
            h_vFe_energy->SetTitle(Form("vFe Energy (w=%.5f)", param[3]));
            h_vFe_energy->GetXaxis()->SetTitle("Energy");
            h_vFe_energy->Draw("HIST");

            c4->cd(5);
            gPad->SetLogy();
            TH2D* h_vPb_scaled = (TH2D*)g_hist_vPb->Clone("h_vPb_scaled");
            h_vPb_scaled->Scale(param[4]);
            TH1D* h_vPb_energy = h_vPb_scaled->ProjectionY("h_vPb_energy");
            h_vPb_energy->SetLineColor(kYellow+3);
            h_vPb_energy->SetLineWidth(2);
            h_vPb_energy->SetTitle(Form("vPb Energy (w=%.5f)", param[4]));
            h_vPb_energy->GetXaxis()->SetTitle("Energy");
            h_vPb_energy->Draw("HIST");

            c4->cd(6);
            gPad->SetLogy();
            TH2D* h_SSB_scaled = (TH2D*)g_hist_SSB->Clone("h_SSB_scaled");
            h_SSB_scaled->Scale(param[5]);
            TH1D* h_SSB_energy = h_SSB_scaled->ProjectionY("h_SSB_energy");
            h_SSB_energy->SetLineColor(kOrange);
            h_SSB_energy->SetLineWidth(2);
            h_SSB_energy->SetTitle(Form("SSB Energy (w=%.5f)", param[5]));
            h_SSB_energy->GetXaxis()->SetTitle("Energy");
            h_SSB_energy->Draw("HIST");

            c4->cd(7);
            gPad->SetLogy();
            TH2D* h_Me_scaled = (TH2D*)g_target->Clone("h_Me_scaled");
            TH1D* h_Me_energy = h_Me_scaled->ProjectionY("h_Me_energy");
            h_Me_energy->SetLineColor(kMagenta+2);
            h_Me_energy->SetLineWidth(2);
            h_Me_energy->SetTitle(Form("Target Energy"));
            h_Me_energy->GetXaxis()->SetTitle("Energy");
            h_Me_energy->Draw("HIST");
            
            // Plot weighted MC components - Time projections
            c4->cd(8);
            gPad->SetLogy();
            TH1D* h_vD_time = h_vD_scaled->ProjectionX("h_vD_time");
            h_vD_time->SetLineColor(kRed);
            h_vD_time->SetLineWidth(2);
            h_vD_time->SetTitle(Form("#nu-D Time (w=%.5f)", param[0]));
            h_vD_time->GetXaxis()->SetTitle("Time");
            h_vD_time->Draw("HIST");
            
            c4->cd(9);
            gPad->SetLogy();
            TH1D* h_vO_time = h_vO_scaled->ProjectionX("h_vO_time");
            h_vO_time->SetLineColor(kBlue);
            h_vO_time->SetLineWidth(2);
            h_vO_time->SetTitle(Form("#nu-O Time (w=%.5f)", param[1]));
            h_vO_time->GetXaxis()->SetTitle("Time");
            h_vO_time->Draw("HIST");
            
            c4->cd(10);
            gPad->SetLogy();
            TH1D* h_BRN_time = h_BRN_scaled->ProjectionX("h_BRN_time");
            h_BRN_time->SetLineColor(kGreen+2);
            h_BRN_time->SetLineWidth(2);
            h_BRN_time->SetTitle(Form("BRN Time (w=%.5f)", param[2]));
            h_BRN_time->GetXaxis()->SetTitle("Time");
            h_BRN_time->Draw("HIST");

            c4->cd(11);
            gPad->SetLogy();
            TH1D* h_vFe_time = h_vFe_scaled->ProjectionX("h_vFe_time");
            h_vFe_time->SetLineColor(kCyan+2);
            h_vFe_time->SetLineWidth(2);
            h_vFe_time->SetTitle(Form("vFe Time (w=%.5f)", param[3]));
            h_vFe_time->GetXaxis()->SetTitle("Time");
            h_vFe_time->Draw("HIST");

            c4->cd(12);
            gPad->SetLogy();
            TH1D* h_vPb_time = h_vPb_scaled->ProjectionX("h_vPb_time");
            h_vPb_time->SetLineColor(kYellow+3);
            h_vPb_time->SetLineWidth(2);
            h_vPb_time->SetTitle(Form("vPb Time (w=%.5f)", param[4]));
            h_vPb_time->GetXaxis()->SetTitle("Time");
            h_vPb_time->Draw("HIST");

            c4->cd(13);
            gPad->SetLogy();
            TH1D* h_SSB_time = h_SSB_scaled->ProjectionX("h_SSB_time");
            h_SSB_time->SetLineColor(kOrange);
            h_SSB_time->SetLineWidth(2);
            h_SSB_time->SetTitle(Form("SSB Time (w=%.5f)", param[5]));
            h_SSB_time->GetXaxis()->SetTitle("Time");
            h_SSB_time->Draw("HIST");

            c4->cd(14);
            gPad->SetLogy();
            TH1D* h_Me_time = h_Me_scaled->ProjectionX("h_Me_time");
            h_Me_time->SetLineColor(kMagenta+2);
            h_Me_time->SetLineWidth(2);
            h_Me_time->SetTitle(Form("Target Time"));
            h_Me_time->GetXaxis()->SetTitle("Time");
            h_Me_time->Draw("HIST");
            
            c4->Update();

            c4->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/RooFit_Code/Fitting_Results/Plots/MC_Components.png");

            // ===================================================================
            // Create a fourth canvas showing different combinations of signals
            // ===================================================================
            TCanvas* c6 = new TCanvas("c6", "1D Plots", 1600, 1200);
            c6->Divide(2, 1);

            // c6->cd(1);
            gPad->SetLogy();
            TH1D* h_BRN_sim_energy = g_hist_BRN->ProjectionY("h_BRN_sim_energy");
            h_BRN_sim_energy->SetLineColor(kRed);
            h_BRN_sim_energy->SetLineWidth(2);
            // h_BRN_sim_energy->SetTitle("Original BRN Energy");
            h_BRN_sim_energy->GetXaxis()->SetTitle("Energy");
            // h_BRN_sim_energy->Draw("HIST");

            // c6->cd(2);
            gPad->SetLogy();
            TH1D* h_vO_sim_energy = g_hist_vO->ProjectionY("h_vO_sim_energy");
            TH1D* h_BRN_vO_sim_energy = (TH1D*)h_BRN_sim_energy->Clone("h_BRN_vO_sim_energy");
            h_BRN_vO_sim_energy->Add(h_vO_sim_energy);
            h_BRN_vO_sim_energy->SetLineColor(kBlue);
            h_BRN_vO_sim_energy->SetLineWidth(2);
            // h_BRN_vO_sim_energy->SetTitle("Original BRN + #nu-O Energy");
            h_BRN_vO_sim_energy->GetXaxis()->SetTitle("Energy");
            // h_BRN_vO_sim_energy->Draw("HIST");

            c6->cd(1);
            gPad->SetLogy();
            TH1D* h_vD_sim_energy = g_hist_vD->ProjectionY("h_vD_sim_energy");
            TH1D* h_BRN_vO_vD_sim_energy = (TH1D*)h_BRN_vO_sim_energy->Clone("h_BRN_vO_vD_sim_energy");
            h_BRN_vO_vD_sim_energy->Add(h_vD_sim_energy);
            h_BRN_vO_vD_sim_energy->SetLineColor(kRed);
            h_BRN_vO_vD_sim_energy->SetLineWidth(2);
            h_BRN_vO_vD_sim_energy->SetTitle(Form("Original BRN + #nu-O + #nu-D Energy"));
            h_BRN_vO_vD_sim_energy->GetXaxis()->SetTitle("Energy");
            h_BRN_vO_vD_sim_energy->Draw("HIST");
            h_BRN_vO_sim_energy->SetLineColor(kBlue);
            h_BRN_vO_sim_energy->SetLineWidth(2);
            h_BRN_vO_sim_energy->Draw("HIST same");
            h_BRN_sim_energy->SetLineColor(kGreen+2);
            h_BRN_sim_energy->SetLineWidth(2);
            h_BRN_sim_energy->Draw("HIST same");
            TLegend* leg3 = new TLegend(0.6, 0.75, 0.88, 0.88);
            leg3->AddEntry(h_BRN_sim_energy, "BRN", "l");
            leg3->AddEntry(h_BRN_vO_sim_energy, "BRN + #nu-O", "l");
            leg3->AddEntry(h_BRN_vO_vD_sim_energy, "BRN + #nu-O + #nu-D", "l");
            leg3->Draw();

            // c6->cd(4);
            gPad->SetLogy();
            h_BRN_energy->SetLineColor(kRed);
            h_BRN_energy->SetLineWidth(2);
            // h_BRN_energy->SetTitle("Weighted BRN Energy");
            h_BRN_energy->GetXaxis()->SetTitle("Energy");
            // h_BRN_energy->Draw("HIST");

            // c6->cd(5);
            gPad->SetLogy();
            TH1D* h_BRN_vO_energy = (TH1D*)h_BRN_energy->Clone("h_BRN_vO_energy");
            h_BRN_vO_energy->Add(h_vO_energy);
            h_BRN_vO_energy->SetLineColor(kBlue);
            h_BRN_vO_energy->SetLineWidth(2);
            // h_BRN_vO_energy->SetTitle("Weighted BRN + #nu-O Energy");
            h_BRN_vO_energy->GetXaxis()->SetTitle("Energy");
            // h_BRN_vO_energy->Draw("HIST");

            c6->cd(2);
            gPad->SetLogy();
            TH1D* h_BRN_vO_vD_energy = (TH1D*)h_BRN_vO_energy->Clone("h_BRN_vO_vD_energy");
            h_BRN_vO_vD_energy->Add(h_vD_energy);
            h_BRN_vO_vD_energy->SetLineColor(kRed);
            h_BRN_vO_vD_energy->SetLineWidth(2);
            h_BRN_vO_vD_energy->SetTitle(Form("Weighted BRN + #nu-O + #nu-D Energy"));
            h_BRN_vO_vD_energy->GetXaxis()->SetTitle("Energy");
            h_BRN_vO_vD_energy->Draw("HIST");
            h_BRN_vO_energy->SetLineColor(kBlue);
            h_BRN_vO_energy->SetLineWidth(2);
            h_BRN_vO_energy->Draw("HIST same");
            h_BRN_energy->SetLineColor(kGreen+2);
            h_BRN_energy->SetLineWidth(2);
            h_BRN_energy->Draw("HIST same");
            TLegend* leg4 = new TLegend(0.6, 0.75, 0.88, 0.88);
            leg4->AddEntry(h_BRN_energy, "BRN", "l");
            leg4->AddEntry(h_BRN_vO_energy, "BRN + #nu-O", "l");
            leg4->AddEntry(h_BRN_vO_vD_energy, "BRN + #nu-O + #nu-D", "l");
            leg4->Draw();

            c6->Update();

            c6->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/RooFit_Code/Fitting_Results/Plots/1D_Plots.png");

            // ===================================================================
            // Canvases c7 (Steady State Background) and c9 (Residual Errors)
            // are DISABLED. Both depend on histograms (h_bg_stats,
            // h_bg_vO_trend, h_sg_energy, h_bg_energy) that were filled
            // event-by-event in the old manual toy-generation code (tracking
            // which individual events were true background vs. true signal).
            // Under the new scheme, create2DHist()/fillSSBHist() build
            // aggregate MC/background templates once outside the loop, and
            // model.generateBinned() draws toy data directly from the total
            // model's bin counts -- there is no per-event species label to
            // tally, so these four histograms are never filled and these
            // canvases would only show empty axes. Left disabled; a
            // meaningful replacement would need a different diagnostic
            // (e.g. inspecting fitResult->correlationMatrix() for N_SSB).
            // ===================================================================
            /*
            TCanvas* c7 = new TCanvas("c7", "Steady State Background", 1600, 1200);
            c7->Divide(2, 1);

            c7->cd(1);
            // gPad->SetLogy();
            h_bg_stats->GetXaxis()->SetTitle("Number of Events");
            h_bg_stats->GetYaxis()->SetTitle("Counts");
            h_bg_stats->Draw();

            c7->cd(2);
            // gPad->SetLogy();
            h_bg_vO_trend->GetXaxis()->SetTitle("Background in Target");
            h_bg_vO_trend->GetYaxis()->SetTitle("vO Events");
            h_bg_vO_trend->Draw("COLZ");

            c7->Update();

            c7->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/RooFit_Code/Fitting_Results/Plots/Steady_State_Background.png");

            // ===================================================================
            // Create a canvas showing true and scaled background plots
            // ===================================================================
            TCanvas* c9 = new TCanvas("c9", "Residual Errors", 2000, 1000);
            c9->Divide(2, 1);

            c9->cd(1);
            h_sg_energy->GetXaxis()->SetTitle("Number of PEs");
            h_sg_energy->GetYaxis()->SetTitle("Counts per 25 PEs");
            h_sg_energy->SetMarkerStyle(20);
            h_sg_energy->SetMarkerSize(1.0);
            h_sg_energy->Draw("E1 X0 P");

            c9->cd(2);
            TH1D* h_res_energy = (TH1D*)h_sg_energy->Clone("h_res_energy");
            h_bg_energy->Scale(0.01);
            h_res_energy->Add(h_bg_energy, -1);
            int numBins = h_sg_energy->GetNbinsX();
            for (int i = 1; i <= numBins; i++) {
                double res_counts = h_sg_energy->GetBinContent(i);
                double res_error  = (res_counts > 0) ? std::sqrt(res_counts) : 0.0;
                h_res_energy->SetBinError(i, res_error);
                // std::cout << "h_sg_energy, Bin " << i << ": Counts = " << h_sg_energy->GetBinContent(i) << ", Error = " << h_sg_energy->GetBinError(i) << std::endl;
                // std::cout << "h_res_energy, Bin " << i << ": Counts = " << h_res_energy->GetBinContent(i) << ", Error = " << h_res_energy->GetBinError(i) << std::endl;
            }
            h_res_energy->SetTitle("In-Beam Window Residual Energy");
            h_res_energy->GetXaxis()->SetTitle("Number of PEs");
            h_res_energy->GetYaxis()->SetTitle("Counts per 25 PEs");
            h_res_energy->SetMarkerStyle(20);
            h_res_energy->SetMarkerSize(1.0);
            h_res_energy->Draw("E1 X0 P");

            c9->Update();

            c9->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/RooFit_Code/Fitting_Results/Plots/Residual_Errors.png");
            */
        }

    }

    // ===================================================================
    // PROFILE LIKELIHOOD SCAN (on the Asimov dataset)
    // Reuses model_c and the N_vD/N_vO/N_BRN/N_vFe/N_vPb/N_SSB objects
    // already built and still in scope from the main loop above, rather
    // than rebuilding a parallel "pl_"-prefixed copy of everything.
    // Runs one clean fit, builds the NLL, then plots the profile likelihood
    // curve for vD, vO, and vPb - the unconstrained parameters of interest.
    // The curves cross PLL=0.5 at the 1-sigma CI and PLL=2.0 at 2-sigma.
    //
    // NOTE on Asimov data: Extended() and ExpectedData() are options for
    // GENERATING data (generate()/generateBinned()), not for fitTo() or
    // createNLL() -- those two methods just fit whatever RooAbsData object
    // they're given, so passing ExpectedData() to them does not change
    // which dataset is used and is simply not a recognized option there.
    // The correct way to get the (non-fluctuated) Asimov dataset is to
    // generate it explicitly with both Extended(true) and ExpectedData(true):
    // Extended(true) tells the generator to use the model's own extended
    // yield sum as the expected total; ExpectedData(true) tells it to fill
    // every bin to that exact expectation instead of drawing a
    // Poisson/multinomial-fluctuated toy. (kTRUE and true are the same
    // value here -- kTRUE is just ROOT's old Bool_t alias for true.)
    // ===================================================================
    cout << "\n===================================================" << endl;
    cout << "Running Profile Likelihood Scan (Asimov dataset)..." << endl;
    cout << "===================================================" << endl;

    // Reset yields to nominal before generating/fitting the Asimov dataset,
    // for a clean, reproducible starting point (matches how the main loop
    // resets these at the start of every iteration).
    N_vD.setVal(vDnum);
    N_vO.setVal(vOnum);
    N_BRN.setVal(BRNnum);
    N_vFe.setVal(vFenum);
    N_vPb.setVal(vPbnum);
    N_SSB.setVal(SSBnum);

    // N_vO's fit range (0, 1e8) doesn't allow the scan below zero. Widen it
    // temporarily just for this diagnostic scan so the profile likelihood
    // curve can still show whether zero is excluded -- nothing later in the
    // script depends on N_vO's original range, since the main loop is done.
    // N_vO.setRange(-200.0, 1e8);

    RooDataHist* asimov_data = model.generateBinned(obs_set,
                                                     RooFit::Extended(true),
                                                     RooFit::ExpectedData(true));

    // Run the fit once to find the best-fit point on the Asimov dataset.
    // With no statistical fluctuation in the data, this should recover
    // ~the nominal yields (up to any pull from the vFe/SSB constraints) --
    // this is the "median expected" fit, useful for judging expected
    // sensitivity rather than one particular toy's outcome.
    model_c->fitTo(rdh_UB,          // *asimov_data
        RooFit::Extended(kTRUE),
        RooFit::PrintLevel(-1),
        RooFit::Strategy(2),
        RooFit::Constrain(RooArgSet(N_SSB, N_vFe)),
        RooFit::MaxCalls(6*100000));

    cout << "  Profile scan best-fit point (Asimov):" << endl;
    cout << "    N_vD  = " << N_vD.getVal()  << endl;
    cout << "    N_vO  = " << N_vO.getVal()  << endl;
    cout << "    N_vPb = " << N_vPb.getVal() << endl;
    cout << "    N_vFe = " << N_vFe.getVal() << endl;
    cout << "    N_BRN = " << N_BRN.getVal() << endl;
    cout << "    N_SSB = " << N_SSB.getVal() << endl;

    // Build the NLL function including the same pull constraints as the
    // main loop's fit. This is the function RooProfileLL minimizes over
    // nuisance parameters at each fixed value of the parameter of interest.
    RooAbsReal* pl_nll = model_c->createNLL(rdh_UB,         // *asimov_data
        RooFit::Extended(kTRUE),
        RooFit::Constrain(RooArgSet(N_SSB, N_vFe)));

    // Build profile likelihood objects for vD, vO, and vPb, directly on the
    // same RooRealVars used throughout the rest of the script.
    RooProfileLL pll_vD ("pll_vD",  "Profile NLL: N_{vD}",  *pl_nll, RooArgSet(N_vD));
    RooProfileLL pll_vO ("pll_vO",  "Profile NLL: N_{vO}",  *pl_nll, RooArgSet(N_vO));
    RooProfileLL pll_vPb("pll_vPb", "Profile NLL: N_{vPb}", *pl_nll, RooArgSet(N_vPb));
    RooProfileLL pll_vFe("pll_vFe", "Profile NLL: N_{vFe}", *pl_nll, RooArgSet(N_vFe));
    RooProfileLL pll_BRN("pll_BRN", "Profile NLL: N_{BRN}", *pl_nll, RooArgSet(N_BRN));
    RooProfileLL pll_SSB("pll_SSB", "Profile NLL: N_{SSB}", *pl_nll, RooArgSet(N_SSB));

    // ===================================================================
    // Plot the profile likelihood curves
    // ===================================================================
    TCanvas* c_profile = new TCanvas("c_profile", "Profile Likelihood", 1800, 1200);
    c_profile->Divide(3, 2);

    // Scan range: ±50% (vD, vPb) or ±200 events (vO) around the best-fit value
    double vD_val  = N_vD.getVal();
    double vO_val  = N_vO.getVal();
    double vPb_val = N_vPb.getVal();
    double vFe_val = N_vFe.getVal();
    double BRN_val = N_BRN.getVal();
    double SSB_val = N_SSB.getVal();

    // ===================================================================
    // Profile: N_vD
    // ===================================================================
    c_profile->cd(1);

    N_vD.setRange("scan_vD", vD_val * 0.5, vD_val * 1.5);
    RooPlot* frame_vD = N_vD.frame(
        RooFit::Range("scan_vD"),
        RooFit::Title("Profile Likelihood: N_{#nu-D}"),
        RooFit::Bins(50));
    frame_vD->GetXaxis()->SetTitle("N_{#nu-D} (events)");
    frame_vD->GetYaxis()->SetTitle("-#Delta log(L)");

    pll_vD.plotOn(frame_vD, RooFit::LineColor(kBlue), RooFit::LineWidth(2));

    frame_vD->SetMinimum(0.0);
    frame_vD->SetMaximum(5.0);
    frame_vD->Draw();

    // 1-sigma line at 0.5, 2-sigma line at 2.0
    TLine* line1_vD = new TLine(vD_val * 0.5, 0.5, vD_val * 1.5, 0.5);
    TLine* line2_vD = new TLine(vD_val * 0.5, 2.0, vD_val * 1.5, 2.0);
    line1_vD->SetLineColor(kRed);   line1_vD->SetLineStyle(2); line1_vD->SetLineWidth(2);
    line2_vD->SetLineColor(kOrange); line2_vD->SetLineStyle(2); line2_vD->SetLineWidth(2);
    line1_vD->Draw("same");
    line2_vD->Draw("same");

    TLatex* pl_latex_vD = new TLatex();
    pl_latex_vD->SetNDC();
    pl_latex_vD->SetTextSize(0.04);
    pl_latex_vD->SetTextColor(kBlue);
    pl_latex_vD->DrawLatex(0.15, 0.85, Form("Best fit: %.1f", vD_val));
    pl_latex_vD->SetTextColor(kRed);
    pl_latex_vD->DrawLatex(0.15, 0.78, "-- 1#sigma (PLL = 0.5)");
    pl_latex_vD->SetTextColor(kOrange+1);
    pl_latex_vD->DrawLatex(0.15, 0.71, "-- 2#sigma (PLL = 2.0)");

    // ===================================================================
    // Profile: N_vO
    // ===================================================================
    c_profile->cd(2);

    N_vO.setRange("scan_vO", vO_val - 200, vO_val + 200);
    RooPlot* frame_vO = N_vO.frame(
        RooFit::Range("scan_vO"),
        RooFit::Title("Profile Likelihood: N_{#nu-O}"),
        RooFit::Bins(50));
    frame_vO->GetXaxis()->SetTitle("N_{#nu-O} (events)");
    frame_vO->GetYaxis()->SetTitle("-#Delta log(L)");

    pll_vO.plotOn(frame_vO, RooFit::LineColor(kBlue), RooFit::LineWidth(2));

    frame_vO->SetMinimum(0.0);
    frame_vO->SetMaximum(5.0);
    frame_vO->Draw();

    TLine* line1_vO = new TLine(vO_val - 200, 0.5, vO_val + 200, 0.5);
    TLine* line2_vO = new TLine(vO_val - 200, 2.0, vO_val + 200, 2.0);
    line1_vO->SetLineColor(kRed);    line1_vO->SetLineStyle(2); line1_vO->SetLineWidth(2);
    line2_vO->SetLineColor(kOrange); line2_vO->SetLineStyle(2); line2_vO->SetLineWidth(2);
    line1_vO->Draw("same");
    line2_vO->Draw("same");

    TLatex* pl_latex_vO = new TLatex();
    pl_latex_vO->SetNDC();
    pl_latex_vO->SetTextSize(0.04);
    pl_latex_vO->SetTextColor(kBlue);
    pl_latex_vO->DrawLatex(0.15, 0.85, Form("Best fit: %.1f", vO_val));
    pl_latex_vO->SetTextColor(kRed);
    pl_latex_vO->DrawLatex(0.15, 0.78, "-- 1#sigma (PLL = 0.5)");
    pl_latex_vO->SetTextColor(kOrange+1);
    pl_latex_vO->DrawLatex(0.15, 0.71, "-- 2#sigma (PLL = 2.0)");

    // ===================================================================
    // Profile: N_vPb
    // ===================================================================
    c_profile->cd(3);

    N_vPb.setRange("scan_vPb", vPb_val * 0.5, vPb_val * 1.5);
    RooPlot* frame_vPb = N_vPb.frame(
        RooFit::Range("scan_vPb"),
        RooFit::Title("Profile Likelihood: N_{#nu-Pb}"),
        RooFit::Bins(50));
    frame_vPb->GetXaxis()->SetTitle("N_{#nu-Pb} (events)");
    frame_vPb->GetYaxis()->SetTitle("-#Delta log(L)");

    pll_vPb.plotOn(frame_vPb, RooFit::LineColor(kBlue), RooFit::LineWidth(2));

    frame_vPb->SetMinimum(0.0);
    frame_vPb->SetMaximum(5.0);
    frame_vPb->Draw();

    TLine* line1_vPb = new TLine(vPb_val * 0.5, 0.5, vPb_val * 1.5, 0.5);
    TLine* line2_vPb = new TLine(vPb_val * 0.5, 2.0, vPb_val * 1.5, 2.0);
    line1_vPb->SetLineColor(kRed);    line1_vPb->SetLineStyle(2); line1_vPb->SetLineWidth(2);
    line2_vPb->SetLineColor(kOrange); line2_vPb->SetLineStyle(2); line2_vPb->SetLineWidth(2);
    line1_vPb->Draw("same");
    line2_vPb->Draw("same");

    TLatex* pl_latex_vPb = new TLatex();
    pl_latex_vPb->SetNDC();
    pl_latex_vPb->SetTextSize(0.04);
    pl_latex_vPb->SetTextColor(kBlue);
    pl_latex_vPb->DrawLatex(0.15, 0.85, Form("Best fit: %.1f", vPb_val));
    pl_latex_vPb->SetTextColor(kRed);
    pl_latex_vPb->DrawLatex(0.15, 0.78, "-- 1#sigma (PLL = 0.5)");
    pl_latex_vPb->SetTextColor(kOrange+1);
    pl_latex_vPb->DrawLatex(0.15, 0.71, "-- 2#sigma (PLL = 2.0)");

    // ===================================================================
    // Profile: N_vFe
    // ===================================================================
    c_profile->cd(4);

    N_vFe.setRange("scan_vFe", vFe_val * 0.5, vFe_val * 1.5);
    RooPlot* frame_vFe = N_vFe.frame(
        RooFit::Range("scan_vFe"),
        RooFit::Title("Profile Likelihood: N_{#nu-Fe}"),
        RooFit::Bins(50));
    frame_vFe->GetXaxis()->SetTitle("N_{#nu-Fe} (events)");
    frame_vFe->GetYaxis()->SetTitle("-#Delta log(L)");

    pll_vFe.plotOn(frame_vFe, RooFit::LineColor(kBlue), RooFit::LineWidth(2));

    frame_vFe->SetMinimum(0.0);
    frame_vFe->SetMaximum(5.0);
    frame_vFe->Draw();

    TLine* line1_vFe = new TLine(vFe_val * 0.5, 0.5, vFe_val * 1.5, 0.5);
    TLine* line2_vFe = new TLine(vFe_val * 0.5, 2.0, vFe_val * 1.5, 2.0);
    line1_vFe->SetLineColor(kRed);    line1_vFe->SetLineStyle(2); line1_vFe->SetLineWidth(2);
    line2_vFe->SetLineColor(kOrange); line2_vFe->SetLineStyle(2); line2_vFe->SetLineWidth(2);
    line1_vFe->Draw("same");
    line2_vFe->Draw("same");

    TLatex* pl_latex_vFe = new TLatex();
    pl_latex_vFe->SetNDC();
    pl_latex_vFe->SetTextSize(0.04);
    pl_latex_vFe->SetTextColor(kBlue);
    pl_latex_vFe->DrawLatex(0.15, 0.85, Form("Best fit: %.1f", vFe_val));
    pl_latex_vFe->SetTextColor(kRed);
    pl_latex_vFe->DrawLatex(0.15, 0.78, "-- 1#sigma (PLL = 0.5)");
    pl_latex_vFe->SetTextColor(kOrange+1);
    pl_latex_vFe->DrawLatex(0.15, 0.71, "-- 2#sigma (PLL = 2.0)");

    // ===================================================================
    // Profile: N_BRN
    // ===================================================================
    c_profile->cd(5);

    N_BRN.setRange("scan_BRN", BRN_val * 0.5, BRN_val * 1.5);
    RooPlot* frame_BRN = N_BRN.frame(
        RooFit::Range("scan_BRN"),
        RooFit::Title("Profile Likelihood: N_BRN"),
        RooFit::Bins(50));
    frame_BRN->GetXaxis()->SetTitle("N_BRN (events)");
    frame_BRN->GetYaxis()->SetTitle("-#Delta log(L)");

    pll_BRN.plotOn(frame_BRN, RooFit::LineColor(kBlue), RooFit::LineWidth(2));

    frame_BRN->SetMinimum(0.0);
    frame_BRN->SetMaximum(5.0);
    frame_BRN->Draw();

    TLine* line1_BRN = new TLine(BRN_val * 0.5, 0.5, BRN_val * 1.5, 0.5);
    TLine* line2_BRN = new TLine(BRN_val * 0.5, 2.0, BRN_val * 1.5, 2.0);
    line1_BRN->SetLineColor(kRed);    line1_BRN->SetLineStyle(2); line1_BRN->SetLineWidth(2);
    line2_BRN->SetLineColor(kOrange); line2_BRN->SetLineStyle(2); line2_BRN->SetLineWidth(2);
    line1_BRN->Draw("same");
    line2_BRN->Draw("same");

    TLatex* pl_latex_BRN = new TLatex();
    pl_latex_BRN->SetNDC();
    pl_latex_BRN->SetTextSize(0.04);
    pl_latex_BRN->SetTextColor(kBlue);
    pl_latex_BRN->DrawLatex(0.15, 0.85, Form("Best fit: %.1f", BRN_val));
    pl_latex_BRN->SetTextColor(kRed);
    pl_latex_BRN->DrawLatex(0.15, 0.78, "-- 1#sigma (PLL = 0.5)");
    pl_latex_BRN->SetTextColor(kOrange+1);
    pl_latex_BRN->DrawLatex(0.15, 0.71, "-- 2#sigma (PLL = 2.0)");

    // ===================================================================
    // Profile: N_SSB
    // ===================================================================
    c_profile->cd(6);

    N_SSB.setRange("scan_SSB", SSB_val * 0.5, SSB_val * 1.5);
    RooPlot* frame_SSB = N_SSB.frame(
        RooFit::Range("scan_SSB"),
        RooFit::Title("Profile Likelihood: N_SSB"),
        RooFit::Bins(50));
    frame_SSB->GetXaxis()->SetTitle("N_SSB (events)");
    frame_SSB->GetYaxis()->SetTitle("-#Delta log(L)");

    pll_SSB.plotOn(frame_SSB, RooFit::LineColor(kBlue), RooFit::LineWidth(2));

    frame_SSB->SetMinimum(0.0);
    frame_SSB->SetMaximum(5.0);
    frame_SSB->Draw();

    TLine* line1_SSB = new TLine(SSB_val * 0.5, 0.5, SSB_val * 1.5, 0.5);
    TLine* line2_SSB = new TLine(SSB_val * 0.5, 2.0, SSB_val * 1.5, 2.0);
    line1_SSB->SetLineColor(kRed);    line1_SSB->SetLineStyle(2); line1_SSB->SetLineWidth(2);
    line2_SSB->SetLineColor(kOrange); line2_SSB->SetLineStyle(2); line2_SSB->SetLineWidth(2);
    line1_SSB->Draw("same");
    line2_SSB->Draw("same");

    TLatex* pl_latex_SSB = new TLatex();
    pl_latex_SSB->SetNDC();
    pl_latex_SSB->SetTextSize(0.04);
    pl_latex_SSB->SetTextColor(kBlue);
    pl_latex_SSB->DrawLatex(0.15, 0.85, Form("Best fit: %.1f", SSB_val));
    pl_latex_SSB->SetTextColor(kRed);
    pl_latex_SSB->DrawLatex(0.15, 0.78, "-- 1#sigma (PLL = 0.5)");
    pl_latex_SSB->SetTextColor(kOrange+1);
    pl_latex_SSB->DrawLatex(0.15, 0.71, "-- 2#sigma (PLL = 2.0)");

    c_profile->Update();
    c_profile->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/RooFit_Code/Fitting_Results/Plots/Profile_Likelihood.png");

    cout << "  Profile likelihood scan complete." << endl;
    cout << "  Canvas saved to Profile_Likelihood.png" << endl;

    // Clean up heap-allocated objects
    delete pl_nll;
    delete asimov_data;

    // ===================================================================
    // PLOTTING
    // ===================================================================
    cout << "\n========================================" << endl;
    cout << "Creating plots..." << endl;
    cout << "========================================" << endl;

    // ===================================================================
    // Create a canvas showing fitting results over many loops
    // ===================================================================
    TCanvas* c1 = new TCanvas("c1", "Interaction Event Rates", 1800, 1200);
    c1->Divide(4, 2);
    
    // ===================================================================
    // Plot 1: vD
    // ===================================================================
    c1->cd(1);
    h_vD_num->SetTitle("vD");
    h_vD_num->GetXaxis()->SetTitle("Number of Events");
    h_vD_num->GetYaxis()->SetTitle("Number of Fit Attempts");
    h_vD_num->Draw();

    // Draw line at predicted event rate
    TLine *line1 = new TLine(vDnum, 0, vDnum, h_vD_num->GetMaximum());
    line1->SetLineColor(kBlack);
    line1->SetLineWidth(3);
    line1->Draw("same");

    // Fit distribution with Gaussian
    TF1* gauss1 = new TF1("gauss", "gaus", 560, 690);
    h_vD_num->Fit(gauss1, "Q");
    
    TLatex* latex1 = new TLatex();
    latex1->SetNDC();
    latex1->SetTextSize(0.04);
    latex1->DrawLatex(0.15, 0.85, Form("Mean = %.2f", gauss1->GetParameter(1)));
    latex1->DrawLatex(0.15, 0.80, Form("Sigma = %.2f", gauss1->GetParameter(2)));
    
    // ===================================================================
    // Plot 2: vO
    // ===================================================================
    c1->cd(2);
    h_vO_num->SetTitle("vO");
    h_vO_num->GetXaxis()->SetTitle("Number of Events");
    h_vO_num->GetYaxis()->SetTitle("Number of Fit Attempts");
    h_vO_num->Draw();

    // Draw line at predicted event rate
    TLine *line2 = new TLine(vOnum, 0, vOnum, h_vO_num->GetMaximum());
    line2->SetLineColor(kBlack);
    line2->SetLineWidth(3);
    line2->Draw("same");

    // Fit distribution with Gaussian
    TF1* gauss2 = new TF1("gauss", "gaus", 60, 190);
    h_vO_num->Fit(gauss2, "Q");
    
    TLatex* latex2 = new TLatex();
    latex2->SetNDC();
    latex2->SetTextSize(0.04);
    latex2->DrawLatex(0.15, 0.85, Form("Mean = %.2f", gauss2->GetParameter(1)));
    latex2->DrawLatex(0.15, 0.80, Form("Sigma = %.2f", gauss2->GetParameter(2)));
    
    // ===================================================================
    // Plot 3: BRN
    // ===================================================================
    c1->cd(3);
    h_BRN_num->SetTitle("BRN");
    h_BRN_num->GetXaxis()->SetTitle("Number of Events");
    h_BRN_num->GetYaxis()->SetTitle("Number of Fit Attempts");
    h_BRN_num->Draw();

    // Draw line at predicted event rate
    TLine *line3 = new TLine(BRNnum, 0, BRNnum, h_BRN_num->GetMaximum());
    line3->SetLineColor(kBlack);
    line3->SetLineWidth(3);
    line3->Draw("same");

    // Fit distribution with Gaussian
    TF1* gauss3 = new TF1("gauss", "gaus", 110, 150);
    h_BRN_num->Fit(gauss3, "Q");
    
    TLatex* latex3 = new TLatex();
    latex3->SetNDC();
    latex3->SetTextSize(0.04);
    latex3->DrawLatex(0.15, 0.85, Form("Mean = %.2f", gauss3->GetParameter(1)));
    latex3->DrawLatex(0.15, 0.80, Form("Sigma = %.2f", gauss3->GetParameter(2)));

    // ===================================================================
    // Plot 4: vFe
    // ===================================================================
    c1->cd(4);
    h_vFe_num->SetTitle("vFe");
    h_vFe_num->GetXaxis()->SetTitle("Number of Events");
    h_vFe_num->GetYaxis()->SetTitle("Number of Fit Attempts");
    h_vFe_num->Draw();

    // Draw line at predicted event rate
    TLine *line4 = new TLine(vFenum, 0, vFenum, h_vFe_num->GetMaximum());
    line4->SetLineColor(kBlack);
    line4->SetLineWidth(3);
    line4->Draw("same");

    // Fit distribution with Gaussian
    TF1* gauss4 = new TF1("gauss", "gaus", 110, 150);
    h_vFe_num->Fit(gauss4, "Q");
    
    TLatex* latex4 = new TLatex();
    latex4->SetNDC();
    latex4->SetTextSize(0.04);
    latex4->DrawLatex(0.15, 0.85, Form("Mean = %.2f", gauss4->GetParameter(1)));
    latex4->DrawLatex(0.15, 0.80, Form("Sigma = %.2f", gauss4->GetParameter(2)));

    // ===================================================================
    // Plot 5: vPb
    // ===================================================================
    c1->cd(5);
    h_vPb_num->SetTitle("vPb");
    h_vPb_num->GetXaxis()->SetTitle("Number of Events");
    h_vPb_num->GetYaxis()->SetTitle("Number of Fit Attempts");
    h_vPb_num->Draw();

    // Draw line at predicted event rate
    TLine *line5 = new TLine(vPbnum, 0, vPbnum, h_vPb_num->GetMaximum());
    line5->SetLineColor(kBlack);
    line5->SetLineWidth(3);
    line5->Draw("same");

    // Fit distribution with Gaussian
    TF1* gauss5 = new TF1("gauss", "gaus", 110, 150);
    h_vPb_num->Fit(gauss5, "Q");
    
    TLatex* latex5 = new TLatex();
    latex5->SetNDC();
    latex5->SetTextSize(0.04);
    latex5->DrawLatex(0.15, 0.85, Form("Mean = %.2f", gauss5->GetParameter(1)));
    latex5->DrawLatex(0.15, 0.80, Form("Sigma = %.2f", gauss5->GetParameter(2)));

    // ===================================================================
    // Plot 6: SSB
    // ===================================================================
    c1->cd(6);
    h_SSB_num->SetTitle("SSB");
    h_SSB_num->GetXaxis()->SetTitle("Number of Events");
    h_SSB_num->GetYaxis()->SetTitle("Number of Fit Attempts");
    h_SSB_num->Draw();

    // Draw line at predicted event rate
    TLine *line6 = new TLine(SSBnum, 0, SSBnum, h_SSB_num->GetMaximum());
    line6->SetLineColor(kBlack);
    line6->SetLineWidth(3);
    line6->Draw("same");

    // Fit distribution with Gaussian
    TF1* gauss6 = new TF1("gauss", "gaus", 110, 150);
    h_SSB_num->Fit(gauss6, "Q");
    
    TLatex* latex6 = new TLatex();
    latex6->SetNDC();
    latex6->SetTextSize(0.04);
    latex6->DrawLatex(0.15, 0.85, Form("Mean = %.2f", gauss6->GetParameter(1)));
    latex6->DrawLatex(0.15, 0.80, Form("Sigma = %.2f", gauss6->GetParameter(2)));

    // ===================================================================
    // Plot 7: NLL per bin
    // ===================================================================

    c1->cd(7);
    // h_chi2_per_ndf->SetTitle("NLL per Number of Bins Used");
    if (fit_method == 1) {h_chi2_per_ndf->GetXaxis()->SetTitle("Chi2/NDF");}
    if (fit_method == 2) {h_chi2_per_ndf->GetXaxis()->SetTitle("NLL/Bins");}
    h_chi2_per_ndf->GetYaxis()->SetTitle("Number of Fit Attempts");
    h_chi2_per_ndf->Draw();

    // Fit distribution with Gaussian
    TF1* gauss7 = new TF1("gauss", "gaus", -29, -27);
    h_chi2_per_ndf->Fit(gauss7, "Q");

    TLatex* latex7 = new TLatex();
    latex7->SetNDC();
    latex7->SetTextSize(0.04);
    latex7->DrawLatex(0.15, 0.85, Form("Mean = %.2f", gauss7->GetParameter(1)));
    latex7->DrawLatex(0.15, 0.80, Form("Sigma = %.2f", gauss7->GetParameter(2)));

    c1->Update();

    c1->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/RooFit_Code/Fitting_Results/Plots/Interaction_Event_Rates.png");

    // ===================================================================
    // Create a canvas showing correlations between different weights
    // ===================================================================
    TCanvas* c8 = new TCanvas("c8", "Recovered Weight Correlations", 1600, 1200);
    c8->Divide(3, 2);
    
    // ===================================================================
    // Plot 1: vD-vO
    // ===================================================================
    
    c8->cd(1);
    h_vD_vO->SetTitle("vD-vO");
    h_vD_vO->GetXaxis()->SetTitle("vD Event Rate");
    h_vD_vO->GetYaxis()->SetTitle("vO Event Rate");
    h_vD_vO->Draw("COLZ");

    // ===================================================================
    // Plot 2: vD-vFe
    // ===================================================================
    
    c8->cd(2);
    h_vD_vFe->SetTitle("vD-vFe");
    h_vD_vFe->GetXaxis()->SetTitle("vD Event Rate");
    h_vD_vFe->GetYaxis()->SetTitle("vFe Event Rate");
    h_vD_vFe->Draw("COLZ");

    // ===================================================================
    // Plot 3: vD-vPb
    // ===================================================================
    
    c8->cd(3);
    h_vD_vPb->SetTitle("vD-vPb");
    h_vD_vPb->GetXaxis()->SetTitle("vD Event Rate");
    h_vD_vPb->GetYaxis()->SetTitle("vPb Event Rate");
    h_vD_vPb->Draw("COLZ");

    // ===================================================================
    // Plot 4: vO-vFe
    // ===================================================================
    
    c8->cd(4);
    h_vO_vFe->SetTitle("vO-vFe");
    h_vO_vFe->GetXaxis()->SetTitle("vO Event Rate");
    h_vO_vFe->GetYaxis()->SetTitle("vFe Event Rate");
    h_vO_vFe->Draw("COLZ");

    // ===================================================================
    // Plot 5: vO-vPb
    // ===================================================================
    
    c8->cd(5);
    h_vO_vPb->SetTitle("vO-vPb");
    h_vO_vPb->GetXaxis()->SetTitle("vO Event Rate");
    h_vO_vPb->GetYaxis()->SetTitle("vPb Event Rate");
    h_vO_vPb->Draw("COLZ");

    // ===================================================================
    // Plot 6: vFe-vPb
    // ===================================================================
    
    c8->cd(6);
    h_vFe_vPb->SetTitle("vFe-vPb");
    h_vFe_vPb->GetXaxis()->SetTitle("vFe Event Rate");
    h_vFe_vPb->GetYaxis()->SetTitle("vPb Event Rate");
    h_vFe_vPb->Draw("COLZ");

    c8->Update();

    c8->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/RooFit_Code/Fitting_Results/Plots/Recovered_Weight_Correlations.png");
    
    cout << "\nPlots created successfully!" << endl;
    cout << "Canvas 1: Number of Fitted Events of Each Interaction Type" << endl;
    cout << "Canvas 2: Data vs MC comparison with 2D histograms and projections" << endl;
    cout << "Canvas 3: Individual MC components (energy and time projections)" << endl;
    cout << "========================================\n" << endl;

    std::cout<<nFailed<<" fits failed"<<std::endl;
}
    
