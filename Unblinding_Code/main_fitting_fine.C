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
using namespace RooFit;

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
    double maxTime_ns = hist->GetXaxis()->GetXmax(); 
    TH1D timingHist("timingHist", "",static_cast<int>(maxTime_ns),0,maxTime_ns);
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

//TODO THIS IGNORES 5/6 OF THE STEADY STATE DATA
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

void main_fitting_fine() {

    bool regenerate_pdfs = false; //Set to true to generate PDFs. Set to false to load from file

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

    // Define number of events in target, without efficiencies. From Yuri's calculations
    double vDnum = 581.4; 
    double vOnum = 94.6;
    double BRNnum = 654; //This is an estimate, scaled to match the number reported in the original code with similar effs applied.
    double vFenum = 501;
    double vPbnum = 17244.5;

    double peak_spread_rms_eff = 0.99;  //Afterglow cut, so equivalent to a dead time i.e. uniform efficiency
    double time_integration_eff = 0.987; //Only impacts neutrino PDFs
    double dead_time_eff = 0.996; //Uniform

    double vD_sim_eff;
    double vO_sim_eff;
    double vFe_sim_eff;
    double vPb_sim_eff;
    double brn_sim_eff;
    int nSSB_x100;
    if (regenerate_pdfs) {
        g_hist_vD = new TH2D("h_vD", "Neutrino-Deuterium Time vs Energy Distribution", tbins, tmin, tmax, ebins, emin, emax);
        g_hist_vO = new TH2D("h_vO", "Neutrino-Oxygen Time vs Energy Distribution", tbins, tmin, tmax, ebins, emin, emax);
        g_hist_BRN = new TH2D("h_BRN", "Beam-Related Neutrons Time vs Energy Distribution", tbins, tmin, tmax, ebins, emin, emax);
        g_hist_vFe = new TH2D("h_vFe", "Neutrino-Iron Time vs Energy Distribution", tbins, tmin, tmax, ebins, emin, emax);
        g_hist_vPb = new TH2D("h_vPb", "Neutrino-Lead Time vs Energy Distribution", tbins, tmin, tmax, ebins, emin, emax);
        g_hist_SSB = new TH2D("h_SSB", "Out-of-Beam-Window Steady State Background Time vs Energy Distribution", tbins, tmin, tmax, ebins, emin, emax);
        g_hist_UB = new TH2D("h_UB", "Unblinded Beam Spill Data Time vs Energy Distribution", tbins, tmin, tmax, ebins, emin, emax);
        
        vD_sim_eff = create2DHist(g_hist_vD,"sims/vD_Energy.root","sims/neutrino_timing.txt",integration_time_ns,minPMTs_req,emin,emax);
        vO_sim_eff = create2DHist(g_hist_vO,"sims/vO_Energy.root","sims/neutrino_timing.txt",integration_time_ns,minPMTs_req,emin,emax);
        vFe_sim_eff = create2DHist(g_hist_vFe,"sims/vFe_Energy.root","sims/neutrino_timing.txt",integration_time_ns,minPMTs_req,emin,emax);
        vPb_sim_eff = create2DHist(g_hist_vPb,"sims/vPb_Energy.root","sims/neutrino_timing.txt",integration_time_ns,minPMTs_req,emin,emax);
        brn_sim_eff = create2DHist(g_hist_BRN,"sims/BRN_Energy.root","sims/BRN_timing.txt",integration_time_ns,minPMTs_req,emin,emax);
        nSSB_x100 = fillSSBHist(g_hist_SSB,"sims/SSB_Time_And_Energy.txt",emin,emax);
        nUB = fillSSBHist(g_hist_UB,"data/UB_Time_And_Energy.txt",emin,emax);

        TFile* outFile = new TFile("pdfs.root","RECREATE");
        //We are storing the simulation efficiency as the integral to preserve this number in the outout
        g_hist_vD->Scale(vD_sim_eff/g_hist_vD->Integral());
        g_hist_vO->Scale(vO_sim_eff/g_hist_vO->Integral());
        g_hist_vFe->Scale(vFe_sim_eff/g_hist_vFe->Integral());
        g_hist_vPb->Scale(vPb_sim_eff/g_hist_vPb->Integral());
        g_hist_BRN->Scale(brn_sim_eff/g_hist_BRN->Integral());

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
        brn_sim_eff = g_hist_BRN->Integral();
        nSSB_x100 = g_hist_SSB->Integral();
        nUB = g_hist_UB->Integral();
    }

    double SSBnum = 0.01*nSSB_x100;
    std::cout<<"vD simulation eff (not including fixed effs) is "<<vD_sim_eff<<std::endl;
    std::cout<<"vO simulation eff (not including fixed effs) is "<<vO_sim_eff<<std::endl;
    std::cout<<"vPb simulation eff (not including fixed effs) is "<<vPb_sim_eff<<std::endl;
    std::cout<<"vFe simulation eff (not including fixed effs) is "<<vFe_sim_eff<<std::endl;
    std::cout<<"BRN simulation ef  (not including fixed effs) is "<<brn_sim_eff<<std::endl;
    std::cout<<"Number of SSB events in time window (x100) is "<<nSSB_x100<<std::endl;
    
    //Include fixed efficiencies
    vDnum *= vD_sim_eff*peak_spread_rms_eff*time_integration_eff*dead_time_eff;
    vOnum *= vO_sim_eff*peak_spread_rms_eff*time_integration_eff*dead_time_eff;
    vFenum *= vFe_sim_eff*peak_spread_rms_eff*time_integration_eff*dead_time_eff;
    vPbnum *= vPb_sim_eff*peak_spread_rms_eff*time_integration_eff*dead_time_eff;
    BRNnum *= brn_sim_eff*peak_spread_rms_eff*dead_time_eff; 
    double total_num = vDnum + vOnum + vFenum + vPbnum + BRNnum + SSBnum;

    std::cout<<"Nominal number of vD events is "<<vDnum<<std::endl;
    std::cout<<"Nominal number of vO events is "<<vOnum<<std::endl;
    std::cout<<"Nominal number of vFe events is "<<vFenum<<std::endl;
    std::cout<<"Nominal number of vPb events is "<<vPbnum<<std::endl;
    std::cout<<"Nominal number of BRN events is "<<BRNnum<<std::endl;
    std::cout<<"Nominal number of SSB events is "<<SSBnum<<std::endl;

    ///////////////////////////////////
    //For plotting output of toy fits//
    ///////////////////////////////////
    TH2D* h_fitted_2D = new TH2D("h_fitted_2D", "Fitted MC (weighted sum)", tbins, tmin, tmax, ebins, emin, emax);

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

    std::filesystem::create_directory("plots");

    /////////////////
    //RooFit set-up//
    /////////////////
    // Suppress RooFit output - only show errors
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    
    //RooRealVars defined outside the loop
    RooRealVar roo_time("roo_time", "Time", tmin, tmax);
    RooRealVar roo_energy("roo_energy", "Energy", emin, emax);
    RooArgSet obs_set(roo_time, roo_energy);

    // Create RooDataHist for each MC template from TH2D
    RooDataHist rdh_vD ("rdh_vD",  "vD template",  RooArgList(roo_time, roo_energy), g_hist_vD);
    RooDataHist rdh_vO ("rdh_vO",  "vO template",  RooArgList(roo_time, roo_energy), g_hist_vO);
    RooDataHist rdh_BRN("rdh_BRN", "BRN template", RooArgList(roo_time, roo_energy), g_hist_BRN);
    RooDataHist rdh_vFe("rdh_vFe", "vFe template", RooArgList(roo_time, roo_energy), g_hist_vFe);
    RooDataHist rdh_vPb("rdh_vPb", "vPb template", RooArgList(roo_time, roo_energy), g_hist_vPb);
    RooDataHist rdh_SSB("rdh_SSB", "SSB template", RooArgList(roo_time, roo_energy), g_hist_SSB);
    RooDataHist rdh_UB ("rdh_UB",  "Unblinded Data",  RooArgList(roo_time, roo_energy), g_hist_UB);

    // Create a shape PDF for each template (normalized to unit area internally)
    RooHistPdf pdf_vD ("pdf_vD",  "vD PDF",  obs_set, rdh_vD,0);
    RooHistPdf pdf_vO ("pdf_vO",  "vO PDF",  obs_set, rdh_vO,0);
    RooHistPdf pdf_BRN("pdf_BRN", "BRN PDF", obs_set, rdh_BRN,0);
    RooHistPdf pdf_vFe("pdf_vFe", "vFe PDF", obs_set, rdh_vFe,0);
    RooHistPdf pdf_vPb("pdf_vPb", "vPb PDF", obs_set, rdh_vPb,0);
    RooHistPdf pdf_SSB("pdf_SSB", "SSB PDF", obs_set, rdh_SSB,0);

    RooRealVar N_vD ("N_vD",  "vD yield",  vDnum,  0, 1e6);
    RooRealVar N_vO ("N_vO",  "vO yield",  vOnum,  0, 1e6);
    RooRealVar N_BRN("N_BRN", "BRN yield", BRNnum, 0, 1e6);
    RooRealVar N_vFe("N_vFe", "vFe yield", vFenum, 0, 1e6); 
    RooRealVar N_vPb("N_vPb", "vPb yield", vPbnum, 0, 1e6);
    RooRealVar N_SSB("N_SSB", "SSB yield", SSBnum, 0, 1e6);

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
    RooGaussian Fe_constraint("Fe_constraint","",N_vFe,N_vFe_mean,N_vFe_sigma);

    RooProdPdf* model_c = new RooProdPdf("model_c","model_c",RooArgList(model,SSB_constraint,Fe_constraint));

    ///////////////
    // MAIN LOOP //
    ///////////////
    //best fit values as starting points for NLL
    int max_loop = 1;
    int min_loop = 0;
    int sample_loop = (max_loop - 1);
    int nFailed=0;
    for (int iLoop = min_loop; iLoop < max_loop; iLoop++) {

        //cout << "\n========================================" << endl;
        //cout << "================ " << "LOOP " << iLoop + 1 << " ================" << endl;
        //cout << "========================================" << endl;

        //Reset amplitudes to initial guesses each time
        N_vD.setVal(vDnum);
        N_vO.setVal(vOnum);
        N_vFe.setVal(vFenum);
        N_vPb.setVal(vPbnum);
        N_BRN.setVal(BRNnum);
        N_SSB.setVal(SSBnum);
            
        auto* data = model.generateBinned(obs_set,
                                          RooFit::NumEvents(static_cast<int>(std::round(total_num))),
                                          RooFit::Extended(true),
                                          RooFit::ExpectedData(true));

        RooFitResult* fitResult = model_c->fitTo(
            *rdh_UB,            // *data,
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

        // ============================================================
        // This is a code block to do 1D and 2D profiles for every pair of variables.
        // Plotting in ROOT is (in my opinion) fugly, so the philosophy is just to dump text files and plot them in Python.
        // There's a bit of effort to simplify the code here to not repeat many for loops. 
        // The gist of it is to make a vector of the RooRealVars, a vector of their names, and a vector of their scan ranges. Then loop over those vectors to do the profiling.
        // ============================================================
        
        std::filesystem::create_directory("pll_output");

        // Make vectors of variables & names
        std::vector<RooRealVar*> params = { &N_vD, &N_vO, &N_vFe, &N_vPb, &N_BRN, &N_SSB };
        std::vector<std::string> pnames = { "vD",  "vO",  "vFe",  "vPb",  "BRN",  "SSB"  };

        /*// Scan range: ±50% (vD, vPb) or ±200 events (vO) around the best-fit value
        double vD_val  = vDnum;            // N_vD.getVal();
        double vO_val  = vOnum;            // N_vO.getVal();
        double vFe_val = vFenum;           // N_vFe.getVal();
        double vPb_val = vPbnum;           // N_vPb.getVal();
        double BRN_val = BRNnum;           // N_BRN.getVal();
        double SSB_val = SSBnum;           // N_SSB.getVal();

        // TODO: if you need to change the ranges in the text/plots you can do it here.
        std::vector<std::pair<double,double>> pranges = {
            { vD_val * 0.5, vD_val * 1.5 },     // vD
            { 0, vO_val + 200 },                // vO       // vO_val - 200
            { vFe_val * 0.5, vFe_val * 1.5 },   // vFe
            { vPb_val * 0.5, vPb_val * 1.5 },   // vPb
            { BRN_val * 0.5, BRN_val * 1.5 },   // BRN
            { SSB_val * 0.5, SSB_val * 1.5 }    // SSB
        };*/

        const size_t nPar = params.size();

        // Do initial fit, ensure all are floating
        for (auto* p : params) { p->setConstant(false); }
        std::unique_ptr<RooAbsReal> nll{model_c->createNLL(*rdh_UB, NumCPU(2))};      // Can just delete ", NumCPU(2)"            // *data,
        for (auto* p : params) { p->setConstant(false); }
            
        RooMinimizer m(*nll);
        m.setPrintLevel(-1);
        m.migrad();

        const double nll_min = nll->getVal();

        //Set profile ranges to +-3 sigma (hesse errors), restricting to positive values
        std::vector<std::pair<double,double>> pranges;
        for (auto* p : params) {
            double lower = p->getVal() - 3.0*p->getError();
            double upper = p->getVal() + 3.0*p->getError();
            if (lower < 0) lower = 0;
            pranges.push_back({lower, upper});
        }

        // 1D
        const int n1d = 20;
        for (size_t i = 0; i < nPar; ++i) {
            RooRealVar* par = params[i];
            const double lo = pranges[i].first, hi = pranges[i].second;

            std::ofstream f1d("pll_output/pll_1D_" + pnames[i] + ".txt");
            f1d << "# N_" << pnames[i] << "  twoDeltaNLL\n";

            for (int ix = 0; ix < n1d; ++ix) {
                double val = lo + (hi - lo) * (ix + 0.5) / n1d;

                //Set all to initial values
                N_vD.setVal(vDnum);
                N_vO.setVal(vOnum);
                N_vFe.setVal(vFenum);
                N_vPb.setVal(vPbnum);
                N_BRN.setVal(BRNnum);
                N_SSB.setVal(SSBnum);

                par->setVal(val);
                par->setConstant(true);

                RooMinimizer m(*nll);
                m.setPrintLevel(-1);
                m.migrad();

                double twoDeltaNLL = 2.0 * (nll->getVal() - nll_min);       // x2 so we can compare to chi2 values
                if (twoDeltaNLL < 0 && twoDeltaNLL > -1e-8) twoDeltaNLL = 0.0;
                f1d << val << " " << twoDeltaNLL << "\n";

                par->setConstant(false);
            }
            f1d.close();
            std::cout << "Wrote pll_output/pll_1D_" << pnames[i] << ".txt" << std::endl;
        }

        // 2D
        const int n2d = 20;
        for (size_t i = 0; i < nPar; ++i) {
            for (size_t j = i + 1; j < nPar; ++j) {
                RooRealVar* pX = params[i];
                RooRealVar* pY = params[j];
                const double xlo = pranges[i].first, xhi = pranges[i].second;
                const double ylo = pranges[j].first, yhi = pranges[j].second;

                std::ofstream f2d("pll_output/pll_2D_" + pnames[i] + "_" + pnames[j] + ".txt");
                f2d << "# N_" << pnames[i] << "  N_" << pnames[j] << "  twoDeltaNLL\n";

                for (int ix = 0; ix < n2d; ++ix) {
                    double xVal = xlo + (xhi - xlo) * (ix + 0.5) / n2d;
                    for (int iy = 0; iy < n2d; ++iy) {
                        double yVal = ylo + (yhi - ylo) * (iy + 0.5) / n2d;

                        //Set all to initial values
                        N_vD.setVal(vDnum);
                        N_vO.setVal(vOnum);
                        N_vFe.setVal(vFenum);
                        N_vPb.setVal(vPbnum);
                        N_BRN.setVal(BRNnum);
                        N_SSB.setVal(SSBnum);

                        pX->setVal(xVal); pY->setVal(yVal);
                        pX->setConstant(true); pY->setConstant(true);

                        RooMinimizer m(*nll);
                        m.setPrintLevel(-1);
                        m.migrad();

                        double twoDeltaNLL = 2.0 * (nll->getVal() - nll_min);
                        if (twoDeltaNLL < 0 && twoDeltaNLL > -1e-8) twoDeltaNLL = 0.0;
                        f2d << xVal << " " << yVal << " " << twoDeltaNLL << "\n";

                        pX->setConstant(false); pY->setConstant(false);
                    }
                }
                f2d.close();
                std::cout << "Wrote pll_output/pll_2D_" << pnames[i]
                          << "_" << pnames[j] << ".txt" << std::endl;
            }
        }

    }
    std::cout<<nFailed<<" fits failed"<<std::endl;
}