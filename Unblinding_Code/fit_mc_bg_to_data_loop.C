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
#include "Fit/Fitter.h"
#include "Math/Functor.h"
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

// Global pointers to histograms (needed for chi-squared function)
TH2D* g_hist_vD = nullptr;   // MC: Neutrino-Deuterium
TH2D* g_hist_vO = nullptr;   // MC: Neutrino-Oxygen
TH2D* g_hist_BRN = nullptr;  // MC: Beam-Related Neutrons
TH2D* g_hist_vFe = nullptr;  // MC: Neutrino-Iron
TH2D* g_hist_vPb = nullptr;  // MC: Neutrino-Lead
TH2D* g_hist_SSB = nullptr;  // RD: Steady-State Background
TH2D* g_target = nullptr;    // Pseudo-Real Data: MC & SSB Events

// Chi-squared function that will be minimized
// Parameters: p[0] = weight for vD, p[1] = weight for vO, p[2] = weight for BRN
double ChiSquaredFCN(const double* p) {
    double chi2 = 0.0;
    
    // Loop over all bins in the 2D histogram (time and energy)
    int nbinsX = g_target->GetNbinsX();  // Time bins
    int nbinsY = g_target->GetNbinsY();  // Energy bins
    
    for (int ix = 1; ix <= nbinsX; ++ix) {
        for (int iy = 1; iy <= nbinsY; ++iy) {
            // Calculate the weighted sum of MC histograms for this bin
            double mc_sum = p[0] * g_hist_vD->GetBinContent(ix, iy) + 
                            p[1] * g_hist_vO->GetBinContent(ix, iy) + 
                            p[2] * g_hist_BRN->GetBinContent(ix, iy) +
                            p[3] * g_hist_vFe->GetBinContent(ix, iy) +
                            p[4] * g_hist_vPb->GetBinContent(ix, iy) +
                            p[5] * g_hist_SSB->GetBinContent(ix, iy);
            
            // Real data value and error
            double data_val = g_target->GetBinContent(ix, iy);
            double data_err = 1 + TMath::Sqrt(data_val + 0.75); // g_target->GetBinError(ix, iy);
            
            // Calculate chi-squared value
            if (data_val > 0 && data_err > 0) {
                // Standard Pearson chi-squared for non-empty bins
                double diff = mc_sum - data_val;
                chi2 += (diff * diff) / (data_err * data_err);
            } /*else if (data_val == 0 && mc_sum > 0) {
                // Poisson chi-squared for empty bins
                // For n=0, chi² = 2×μ (from Poisson likelihood)
                chi2 += 2.0 * mc_sum;
            } else if (data_val > 0 && data_err == 0) {
                // Edge case: data present but no error assigned
                // Use Poisson error: σ = √n
                double poisson_err = TMath::Sqrt(data_val);
                double diff = mc_sum - data_val;
                chi2 += (diff * diff) / (poisson_err * poisson_err);
            }
            if (mc_sum > 0) {  // Avoid division by small numbers
                // χ² = (observed - expected)² / expected
                double diff = mc_sum - data_val;
                chi2 += (diff * diff) / mc_sum;
            }*/
            // If both are zero, no contribution
        }
    }
    
    return chi2;  // Return chi-squared value
}

// Negative Log-Likelihood function that will be minimized
// Parameters: p[0] = weight for vD, p[1] = weight for vO, p[2] = weight for BRN
// For Poisson statistics: -log(L) = Σ[μ - n×log(μ)]
// where μ = expected (MC model), n = observed (data)
double NegLogLikelihoodFCN(const double* p) {
    double nll = 0.0;
    
    // Loop over all bins in the 2D histogram (time and energy)
    int nbinsX = g_target->GetNbinsX();  // Time bins
    int nbinsY = g_target->GetNbinsY();  // Energy bins
    
    for (int ix = 1; ix <= nbinsX; ++ix) {
        for (int iy = 1; iy <= nbinsY; ++iy) {
            // Calculate the weighted sum of MC histograms for this bin (expected μ)
            double mu = p[0] * g_hist_vD->GetBinContent(ix, iy) + 
                        p[1] * g_hist_vO->GetBinContent(ix, iy) + 
                        p[2] * g_hist_BRN->GetBinContent(ix, iy) +
                        p[3] * g_hist_vFe->GetBinContent(ix, iy) +
                        p[4] * g_hist_vPb->GetBinContent(ix, iy) +
                        p[5] * g_hist_SSB->GetBinContent(ix, iy);
            
            // Observed count (data)
            double n = g_target->GetBinContent(ix, iy);

            // Poisson negative log-likelihood: -log(L) = μ - n×log(μ)
            if (n == 0 || mu == 0) {
                nll += mu;
            } else {
                nll += mu - n*TMath::Log(mu);
            }
            
        }
    }
    
    return nll;  // Return negative log-likelihood value
}

void fit_mc_bg_to_data_loop() {

    // Define dimensions of all 2D histograms
    int tbins = 20;
    int tmin = 0;
    int tmax = 10000;
    int ebins = 20;
    int emin = 0;
    int emax = 500;

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

    // Define number of events in target
    int vDnum = 615;
    int vOnum = 95;
    int BRNnum = 206;
    int vFenum = 132;
    int vPbnum = 215;
    int SSBnum = 1679;

    // Define selection cuts
    int minEnergy = 50;
    int maxEnergy = 500;

    int sel_cut_case = 3;
    // 1: >=50phe, <= 500phe, mult >= 10
    if (sel_cut_case == 1) {vDnum = 615; vOnum = 95; BRNnum = 206; vFenum = 132; vPbnum = 215; minEnergy = 50; maxEnergy = 500; emax = 500; BRNbins = 40; BRNmin = 0; BRNmax = 400; vPbbins = 40; vPbmin = 0; vPbmax = 400; SSBbins = 60; SSBmin = 1400; SSBmax = 2000; NLLbins = 100; NLLmin = -20; NLLmax = -10;}
    // 2: >=30phe, <= 500phe, mult >= 10
    if (sel_cut_case == 2) {vDnum = 620; vOnum = 104; BRNnum = 717; vFenum = 161; vPbnum = 530; minEnergy = 30; maxEnergy = 500; emax = 500; BRNbins = 40; BRNmin = 500; BRNmax = 900; vPbbins = 40; vPbmin = 350; vPbmax = 750; SSBbins = 60; SSBmin = 3200; SSBmax = 3800; NLLbins = 100; NLLmin = -55; NLLmax = -45;}
    // 3: >=50phe, <= 600phe, mult >= 10
    if (sel_cut_case == 3) {vDnum = 618; vOnum = 95; BRNnum = 206; vFenum = 133; vPbnum = 215; minEnergy = 50; maxEnergy = 600; emax = 600; BRNbins = 40; BRNmin = 0; BRNmax = 400; vPbbins = 40; vPbmin = 0; vPbmax = 400; SSBbins = 60; SSBmin = 1400; SSBmax = 2000; NLLbins = 100; NLLmin = -20; NLLmax = -10;}
    // 4: >=50phe, <= 500phe, mult == 12
    if (sel_cut_case == 4) {vDnum = 611; vOnum = 92; BRNnum = 157; vFenum = 126; vPbnum = 175; minEnergy = 50; maxEnergy = 500; emax = 500; BRNbins = 40; BRNmin = 0; BRNmax = 400; vPbbins = 40; vPbmin = 0; vPbmax = 400; SSBbins = 60; SSBmin = 900; SSBmax = 1500; NLLbins = 100; NLLmin = -15; NLLmax = -5;}

    // Define function to be minimized
    int fit_method = 2;

    g_target = new TH2D("h_target", "Beam Spill Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);
    g_hist_vD = new TH2D("h_vD", "Neutrino-Deuterium Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);
    g_hist_vO = new TH2D("h_vO", "Neutrino-Oxygen Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);
    g_hist_BRN = new TH2D("h_BRN", "Beam-Related Neutrons Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);
    g_hist_vFe = new TH2D("h_vFe", "Neutrino-Iron Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);
    g_hist_vPb = new TH2D("h_vPb", "Neutrino-Lead Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);
    g_hist_SSB = new TH2D("h_SSB", "Out-of-Beam-Window Steady State Background Time vs Energy Distributions", tbins, tmin, tmax, ebins, emin, emax);

    TH2D* h_fitted_2D = new TH2D("h_fitted_2D", "Fitted MC (weighted sum)", tbins, tmin, tmax, ebins, emin, emax);

    TH1D* h_vD_num = new TH1D("h_vD_num", "Number of vD Events Found", vDbins, vDmin, vDmax);
    TH1D* h_vO_num = new TH1D("h_vO_num", "Number of vO Events Found", vObins, vOmin, vOmax);
    TH1D* h_BRN_num = new TH1D("h_BRN_num", "Number of BRN Events Found", BRNbins, BRNmin, BRNmax);
    TH1D* h_vFe_num = new TH1D("h_vFe_num", "Number of vFe Events Found", vFebins, vFemin, vFemax);
    TH1D* h_vPb_num = new TH1D("h_vPb_num", "Number of vPb Events Found", vPbbins, vPbmin, vPbmax);
    TH1D* h_SSB_num = new TH1D("h_SSB_num", "Number of SSB Events Found", SSBbins, SSBmin, SSBmax);

    TH2D* h_vD_vO = new TH2D("h_vD_vO", "vD vs vO Recovered Event Rates", vDbins, vDmin, vDmax, vObins, vOmin, vOmax);
    TH2D* h_vD_BRN = new TH2D("h_vD_BRN", "vD vs BRN Recovered Event Rates", vDbins, vDmin, vDmax, BRNbins, BRNmin, BRNmax);
    TH2D* h_vD_SSB = new TH2D("h_vD_SSB", "vD vs SSB Recovered Event Rates", vDbins, vDmin, vDmax, SSBbins, SSBmin, SSBmax);
    TH2D* h_vO_BRN = new TH2D("h_vO_BRN", "vO vs BRN Recovered Event Rates", vObins, vOmin, vOmax, BRNbins, BRNmin, BRNmax);
    TH2D* h_vO_SSB = new TH2D("h_vO_SSB", "vO vs SSB Recovered Event Rates", vObins, vOmin, vOmax, SSBbins, SSBmin, SSBmax);
    TH2D* h_BRN_SSB = new TH2D("h_BRN_SSB", "BRN vs SSB Recovered Event Rates", BRNbins, BRNmin, BRNmax, SSBbins, SSBmin, SSBmax);

    TH1D* h_chi2_per_ndf;
    if (fit_method == 1) {h_chi2_per_ndf = new TH1D("h_chi2_per_ndf", "Chi-Squared per Degrees of Freedom", 20, 0, 2);}
    if (fit_method == 2) {h_chi2_per_ndf = new TH1D("h_chi2_per_ndf", "Negative Log-Likelihood per Bins Used", NLLbins, NLLmin, NLLmax);}

    TH2D* h_background_true = new TH2D("h_background_true", "Out-of-Beam-Window Steady State Background Data", tbins, tmin, tmax, ebins, emin, emax);
    TH1D* h_bg_stats = new TH1D("h_bg_stats", "Number of Background Events Put into Target", SSBbins, SSBmin, SSBmax);
    TH2D* h_bg_vO_trend = new TH2D("h_bg_vO_trend", "Number of Background Events in Target vs Recovered Number of vO Events", SSBbins, SSBmin, SSBmax, vObins, vOmin, vOmax);
    TH1D* h_bg_energy = new TH1D("h_bg_energy", "Out-of-Beam-Window Steady State Background Energy", ebins, emin, emax);
    TH1D* h_sg_energy = new TH1D("h_sg_energy", "In-Beam Window Signal Energy", ebins, emin, emax);

    std::filesystem::create_directory("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/fitting_results/four_test_cases/plots");

    int max_loop = 1000;
    int min_loop = 0;
    int sample_loop = (max_loop - 1);

    for (int iLoop = min_loop; iLoop < max_loop; iLoop++) {

        cout << "\n========================================" << endl;
        cout << "================ " << "LOOP " << iLoop + 1 << " ================" << endl;
        cout << "========================================" << endl;

        // Generate random number seed
        TRandom3 rng(0);

        // Create a text string, which is used to output the text file
        std::vector<double> bg_time;
        std::vector<double> bg_energy;
        std::vector<double> num_vec;
        string intval;
        string num_str;

        // Clear histograms and vectors
        g_target->Reset();
        g_hist_vD->Reset();
        g_hist_vO->Reset();
        g_hist_BRN->Reset();
        g_hist_vFe->Reset();
        g_hist_vPb->Reset();
        g_hist_SSB->Reset();
        h_background_true->Reset();
        h_bg_energy->Reset();
        h_sg_energy->Reset();
        bg_time.clear();
        bg_energy.clear();
    
        // cout << "\n========================================" << endl;
        // cout << "Loading Data from ROOT Files" << endl;
        // cout << "========================================" << endl;
        
        // ===================================================================
        // LOAD MC HISTOGRAM 1: Neutrino-Deuterium (vD)
        // ===================================================================
        // cout << "\nLoading MC data: Neutrino-Deuterium..." << endl;
        
        TFile *f1 = new TFile("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/nud_1M.root");
        if (!f1 || f1->IsZombie()) {
            cout << "Error: Could not open nud_1M.root" << endl;
            return;
        }
        
        TTree* t1 = (TTree*)f1->Get("Sim_Tree");
        if (!t1) {
            cout << "Error: Could not find tree Sim_Tree" << endl;
            return;
        }
        
        Int_t numHits;
        TBranch *b_numHits;
        t1->SetBranchAddress("numHits", &numHits, &b_numHits);
        
        // Create a text string, which is used to output the text file
        std::vector<int> vD_time_dist;
        string intval1;

        // Read from the text file
        ifstream ReadTextFile1("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/SNS_v_Time_Dist_1M_10us_1.txt");

        // Use a while loop together with the getline() function to read the file line by line
        while (getline(ReadTextFile1, intval1)) {
            // Output the text from the file
            vD_time_dist.push_back(stod(intval1));
        }

        // Close the file
        ReadTextFile1.close();
        ReadTextFile1.clear();

        // Generate a random number from a Gaussian distribution
        double randMean  = 0.0;
        double randSigma = 3.33;
        int randVal1 = std::round(rng.Gaus(randMean, randSigma));

        std::vector<double> rand_vD_val;
        for(int iRand = 0; iRand < vDnum + randVal1; iRand++) {
            rand_vD_val.push_back(std::round(rng.Uniform(0.0, (t1->GetEntries()) * 0.98)));
        }
        
        int count = 0;
        for(int iEvent = 0; iEvent < t1->GetEntries(); iEvent++){
            Long64_t tentry1 = t1->LoadTree(iEvent);
            b_numHits->GetEntry(tentry1);
            if (numHits > minEnergy && numHits < maxEnergy) {
                if (std::find(rand_vD_val.begin(), rand_vD_val.end(), count) != rand_vD_val.end()) {
                    g_target->Fill(vD_time_dist[iEvent], numHits);
                    h_sg_energy->Fill(numHits);
                }
                else {
                    g_hist_vD->Fill(vD_time_dist[iEvent], numHits);
                }
                /*if (count >= iLoop * 650 && count < iLoop * 650 + vDnum + randVal1) {
                    g_target->Fill(vD_time_dist[iEvent], numHits);
                }
                else if (count < iLoop * 650 || count >= iLoop * 650 + vDnum + randVal1) {
                    g_hist_vD->Fill(vD_time_dist[iEvent], numHits);
                }*/
                count += 1;
            }
        }
        
        // cout << "  Loaded " << t1->GetEntries() << " events into h_vD" << endl;
        // cout << "  Histogram integral: " << g_hist_vD->Integral() << endl;
        // cout << "  h_sg_energy integral: " << h_sg_energy->Integral() << endl;
        
        f1->Close();
        delete f1;
        
        // ===================================================================
        // LOAD MC HISTOGRAM 2: Neutrino-Oxygen (vO)
        // ===================================================================
        // cout << "\nLoading MC data: Neutrino-Oxygen..." << endl;
        
        TFile *f2 = new TFile("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/nuO_1M_pyNewton.root");
        if (!f2 || f2->IsZombie()) {
            cout << "Error: Could not open nuO_1M_pyNewton.root" << endl;
            return;
        }
        
        TTree* t2 = (TTree*)f2->Get("Sim_Tree");
        if (!t2) {
            cout << "Error: Could not find tree Sim_Tree" << endl;
            return;
        }
        
        // Int_t numHits;
        // TBranch *b_numHits;
        t2->SetBranchAddress("numHits", &numHits, &b_numHits);

        // Create a text string, which is used to output the text file
        std::vector<int> vO_time_dist;
        string intval2;

        // Read from the text file
        ifstream ReadTextFile2("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/SNS_v_Time_Dist_1M_10us_2.txt");

        // Use a while loop together with the getline() function to read the file line by line
        while (getline(ReadTextFile2, intval2)) {
            // Output the text from the file
            vO_time_dist.push_back(stod(intval2));
        }

        // Close the file
        ReadTextFile2.close();
        ReadTextFile2.clear();

        int randVal2 = std::round(rng.Gaus(randMean, randSigma));

        std::vector<double> rand_vO_val;
        for(int iRand = 0; iRand < vOnum + randVal2; iRand++) {
            rand_vO_val.push_back(std::round(rng.Uniform(0.0, (t2->GetEntries()) * 0.82)));
        }

        count = 0;
        for(int iEvent = 0; iEvent < t2->GetEntries(); iEvent++){
            Long64_t tentry2 = t2->LoadTree(iEvent);
            b_numHits->GetEntry(tentry2);
            if (numHits > minEnergy && numHits < maxEnergy) {
                if (std::find(rand_vO_val.begin(), rand_vO_val.end(), count) != rand_vO_val.end()) {
                    g_target->Fill(vO_time_dist[iEvent], numHits);
                    h_sg_energy->Fill(numHits);
                }
                else {
                    g_hist_vO->Fill(vO_time_dist[iEvent], numHits);
                }
                /*if (count >= iLoop * 650 && count < iLoop * 650 + vOnum + randVal2) {
                    g_target->Fill(vO_time_dist[iEvent], numHits);
                }
                else if (count < iLoop * 650 || count >= iLoop * 650 + vOnum + randVal2) {
                    g_hist_vO->Fill(vO_time_dist[iEvent], numHits);
                }*/
                count += 1;
            }
        }
        
        // cout << "  Loaded " << t2->GetEntries() << " events into h_vO" << endl;
        // cout << "  Histogram integral: " << g_hist_vO->Integral() << endl;
        // cout << "  h_sg_energy integral: " << h_sg_energy->Integral() << endl;
        
        f2->Close();
        delete f2;
        
        // ===================================================================
        // LOAD MC HISTOGRAM 3: Beam-Related Neutrons (BRN)
        // ===================================================================
        // cout << "\nLoading MC data: Beam-Related Neutrons..." << endl;

        TFile *f3;
        if (sel_cut_case == 2) {f3 = new TFile("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/BRN_40MeV_150k_30phe.root");}
        else {f3 = new TFile("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/BRN_40MeV_151k.root");}
        if (!f3 || f3->IsZombie()) {
            cout << "Error: Could not open BRN ROOT file" << endl;
            return;
        }
        
        TTree* t3 = (TTree*)f3->Get("Sim_Tree");
        if (!t3) {
            cout << "Error: Could not find tree Sim_Tree" << endl;
            return;
        }
        
        // Int_t numHits;
        // TBranch *b_numHits;
        t3->SetBranchAddress("numHits", &numHits, &b_numHits);

        // Create a text string, which is used to output the text file
        std::vector<int> BRN_time_dist;
        string intval3;

        // Read from the text file
        ifstream ReadTextFile3("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/BRN_Time_Dist_1M.txt");

        // Use a while loop together with the getline() function to read the file line by line
        while (getline(ReadTextFile3, intval3)) {
            // Output the text from the file
            BRN_time_dist.push_back(stod(intval3));
        }

        // Close the file
        ReadTextFile3.close();
        ReadTextFile3.clear();

        int randVal3 = std::round(rng.Gaus(randMean, randSigma));

        std::vector<double> rand_BRN_val;
        for(int iRand = 0; iRand < BRNnum + randVal3; iRand++) {
            rand_BRN_val.push_back(std::round(rng.Uniform(0.0, (t3->GetEntries()) * 0.93)));
        }

        count = 0;
        for(int iEvent = 0; iEvent < t3->GetEntries(); iEvent++){
            Long64_t tentry3 = t3->LoadTree(iEvent);
            b_numHits->GetEntry(tentry3);
            if (numHits > minEnergy && numHits < maxEnergy) {
                if (std::find(rand_BRN_val.begin(), rand_BRN_val.end(), count) != rand_BRN_val.end()) {
                    g_target->Fill(BRN_time_dist[iEvent], numHits);
                }
                else {
                    g_hist_BRN->Fill(BRN_time_dist[iEvent], numHits);
                }
                /*if (count >= iLoop * 150 && count < iLoop * 150 + BRNnum + randVal3) {
                    g_target->Fill(BRN_time_dist[iEvent], numHits);
                }
                else if (count < iLoop * 150 || count >= iLoop * 150 + BRNnum + randVal3) {
                    g_hist_BRN->Fill(BRN_time_dist[iEvent], numHits);
                }*/
                count += 1;
            }
        }
        
        // cout << "  Loaded " << t3->GetEntries() << " events into h_BRN" << endl;
        // cout << "  Histogram integral: " << g_hist_BRN->Integral() << endl;
        
        f3->Close();
        delete f3;

        // ===================================================================
        // LOAD MC HISTOGRAM 4: Neutrino-Iron (vFe)
        // ===================================================================
        // cout << "\nLoading MC data: Neutrino-Iron..." << endl;
        
        TFile *f4 = new TFile("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/nuFe_100k.root");
        if (!f4 || f4->IsZombie()) {
            cout << "Error: Could not open nuFe_100k.root" << endl;
            return;
        }
        
        TTree* t4 = (TTree*)f4->Get("Sim_Tree");
        if (!t4) {
            cout << "Error: Could not find tree Sim_Tree" << endl;
            return;
        }
        
        // Int_t numHits;
        // TBranch *b_numHits;
        t4->SetBranchAddress("numHits", &numHits, &b_numHits);
        
        // Create a text string, which is used to output the text file
        std::vector<int> vFe_time_dist;
        string intval4;

        // Read from the text file
        ifstream ReadTextFile4("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/SNS_v_Time_Dist_1M_10us_3.txt");

        // Use a while loop together with the getline() function to read the file line by line
        while (getline(ReadTextFile4, intval4)) {
            // Output the text from the file
            vFe_time_dist.push_back(stod(intval4));
        }

        // Close the file
        ReadTextFile4.close();
        ReadTextFile4.clear();

        int randVal4 = std::round(rng.Gaus(randMean, randSigma));

        std::vector<double> rand_vFe_val;
        for(int iRand = 0; iRand < vFenum + randVal4; iRand++) {
            rand_vFe_val.push_back(std::round(rng.Uniform(0.0, (t4->GetEntries()) * 0.23)));
        }
        
        count = 0;
        for(int iEvent = 0; iEvent < t4->GetEntries(); iEvent++){
            Long64_t tentry4 = t4->LoadTree(iEvent);
            b_numHits->GetEntry(tentry4);
            if (numHits > minEnergy && numHits < maxEnergy) {
                if (std::find(rand_vFe_val.begin(), rand_vFe_val.end(), count) != rand_vFe_val.end()) {
                    g_target->Fill(vFe_time_dist[iEvent], numHits);
                    h_sg_energy->Fill(numHits);
                }
                else {
                    g_hist_vFe->Fill(vFe_time_dist[iEvent], numHits);
                }
                /*if (count >= iLoop * 150 && count < iLoop * 150 + vFenum + randVal4) {
                    g_target->Fill(vFe_time_dist[iEvent], numHits);
                }
                else if (count < iLoop * 150 || count >= iLoop * 150 + vFenum + randVal4) {
                    g_hist_vFe->Fill(vFe_time_dist[iEvent], numHits);
                }*/
                count += 1;
            }
        }
        
        // cout << "  Loaded " << t4->GetEntries() << " events into h_vFe" << endl;
        // cout << "  Histogram integral: " << g_hist_vFe->Integral() << endl;
        // cout << "  h_sg_energy integral: " << h_sg_energy->Integral() << endl;
        
        f4->Close();
        delete f4;

        // ===================================================================
        // LOAD MC HISTOGRAM 5: Neutrino-Lead (vPb)
        // ===================================================================
        // cout << "\nLoading MC data: Neutrino-Lead..." << endl;
        
        TFile *f5 = new TFile("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/nuPb_1M.root");
        if (!f5 || f5->IsZombie()) {
            cout << "Error: Could not open nuPb_1M.root" << endl;
            return;
        }
        
        TTree* t5 = (TTree*)f5->Get("Sim_Tree");
        if (!t5) {
            cout << "Error: Could not find tree Sim_Tree" << endl;
            return;
        }
        
        // Int_t numHits;
        // TBranch *b_numHits;
        t5->SetBranchAddress("numHits", &numHits, &b_numHits);
        
        // Create a text string, which is used to output the text file
        std::vector<int> vPb_time_dist;
        string intval5;

        // Read from the text file
        ifstream ReadTextFile5("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/SNS_v_Time_Dist_1M_10us_4.txt");

        // Use a while loop together with the getline() function to read the file line by line
        while (getline(ReadTextFile5, intval5)) {
            // Output the text from the file
            vPb_time_dist.push_back(stod(intval5));
        }

        // Close the file
        ReadTextFile5.close();
        ReadTextFile5.clear();

        int randVal5 = std::round(rng.Gaus(randMean, randSigma));

        std::vector<double> rand_vPb_val;
        for(int iRand = 0; iRand < vPbnum + randVal5; iRand++) {
            rand_vPb_val.push_back(std::round(rng.Uniform(0.0, (t5->GetEntries()) * 0.010)));
        }
        
        count = 0;
        for(int iEvent = 0; iEvent < t5->GetEntries(); iEvent++){
            Long64_t tentry5 = t5->LoadTree(iEvent);
            b_numHits->GetEntry(tentry5);
            if (numHits > minEnergy && numHits < maxEnergy) {
                if (std::find(rand_vPb_val.begin(), rand_vPb_val.end(), count) != rand_vPb_val.end()) {
                    g_target->Fill(vPb_time_dist[iEvent], numHits);
                    h_sg_energy->Fill(numHits);
                }
                else {
                    g_hist_vPb->Fill(vPb_time_dist[iEvent], numHits);
                }
                /*if (count >= iLoop * 250 && count < iLoop * 250 + vPbnum + randVal5) {
                    g_target->Fill(vPb_time_dist[iEvent], numHits);
                }
                else if (count < iLoop * 250 || count >= iLoop * 250 + vPbnum + randVal5) {
                    g_hist_vPb->Fill(vPb_time_dist[iEvent], numHits);
                }*/
                count += 1;
            }
        }
        
        // cout << "  Loaded " << t5->GetEntries() << " events into h_vPb" << endl;
        // cout << "  Histogram integral: " << g_hist_vPb->Integral() << endl;
        // cout << "  h_sg_energy integral: " << h_sg_energy->Integral() << endl;
        
        f5->Close();
        delete f5;
        
        // ===================================================================
        // LOAD REAL DATA HISTOGRAM: Steady State Background (SSB)
        // ===================================================================
        // cout << "\nLoading real data: Background Events..." << endl;

        // Read from the background data text file
        ifstream ReadBackgroundFile;
        if (sel_cut_case == 1) {ReadBackgroundFile.open("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/2026/mar_26/week_mar_23/background_int_dt_txt_files/Background_Int_dt_4216_standard.txt");}
        else if (sel_cut_case == 2) {ReadBackgroundFile.open("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/2026/mar_26/week_mar_23/background_int_dt_txt_files/Background_Int_dt_4216_30phe.txt");}
        else if (sel_cut_case == 3) {ReadBackgroundFile.open("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/2026/mar_26/week_mar_23/background_int_dt_txt_files/Background_Int_dt_4216_600phe.txt");}
        else if (sel_cut_case == 4) {ReadBackgroundFile.open("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/2026/mar_26/week_mar_23/background_int_dt_txt_files/Background_Int_dt_4216_mult12.txt");}
        if (!ReadBackgroundFile.is_open()) {
            cout << "Error: Could not open background file" << endl;
            return;
        }

        count = 0;
        // Use a while loop together with the getline() function to read the file line by line
        while (getline(ReadBackgroundFile, intval)) {
            num_vec.clear();
            num_str.clear();
            intval.push_back('\t');
            for (char ch : intval) {           
                if (ch != '\t') {num_str.push_back(ch);}            
                else if (ch == '\t') {num_vec.push_back(stod(num_str)); num_str.clear();}        
            }
            if (num_vec[0] >= 0.0 && num_vec[0] <= 10000.0 && num_vec[1] > minEnergy && num_vec[1] <= maxEnergy) {
                double rand_num = std::floor(rng.Uniform(1.00001, 100.99999));
                if (rand_num == 50.0) {
                    // if (count < SSBnum) {
                    g_target->Fill(num_vec[0], num_vec[1]);
                    h_sg_energy->Fill(num_vec[1]);
                    count += 1;
                }
                g_hist_SSB->Fill(num_vec[0], num_vec[1]);
                h_background_true->Fill(num_vec[0], num_vec[1]);
                h_bg_energy->Fill(num_vec[1]);
                bg_time.push_back(num_vec[0]);
                bg_energy.push_back(num_vec[1]);
            }
        }

        h_bg_stats->Fill(count);
        SSBnum = 0.01 * g_hist_SSB->Integral();

        // Close the file
        ReadBackgroundFile.close();
        ReadBackgroundFile.clear();
        
        // cout << "  Loaded " << bg_time.size() << " events into h_SSB" << endl;
        // cout << "  Histogram integral: " << g_hist_SSB->Integral() << endl;
        // cout << "  h_sg_energy integral: " << h_sg_energy->Integral() << endl;

        cout << "\nFilled histograms:" << endl;
        cout << "  Loaded " << g_target->Integral() << " events into h_target" << endl;
        cout << "  Loaded " << g_hist_vD->Integral() << " events into h_vD" << endl;
        cout << "  Loaded " << g_hist_vO->Integral() << " events into h_vO" << endl;
        cout << "  Loaded " << g_hist_BRN->Integral() << " events into h_BRN" << endl;
        cout << "  Loaded " << g_hist_vFe->Integral() << " events into h_vFe" << endl;
        cout << "  Loaded " << g_hist_vPb->Integral() << " events into h_vPb" << endl;
        cout << "  Loaded " << g_hist_SSB->Integral() << " events into h_SSB" << endl;

        cout << "\nPredicted event rates:" << endl;
        cout << "  vDnum = " << vDnum << endl;
        cout << "  vOnum = " << vOnum << endl;
        cout << "  BRNnum = " << BRNnum << endl;
        cout << "  vFenum = " << vFenum << endl;
        cout << "  vPbnum = " << vPbnum << endl;
        cout << "  SSBnum = " << SSBnum << endl;
        
        // cout << "\n========================================" << endl;
        // cout << "Setting up Fitting Function" << endl;
        // cout << "========================================" << endl;

        // Normalize MC histograms to match target statistics
        double target_total = g_target->Integral();
        double vD_total = g_hist_vD->Integral();
        double vO_total = g_hist_vO->Integral();
        double BRN_total = g_hist_BRN->Integral();
        double vFe_total = g_hist_vFe->Integral();
        double vPb_total = g_hist_vPb->Integral();
        double SSB_total = g_hist_SSB->Integral();
        double scale_vD = ((vDnum + randVal1) / vD_total);
        double scale_vO = ((vOnum + randVal2) / vO_total);
        double scale_BRN = ((BRNnum + randVal3) / BRN_total);
        double scale_vFe = ((vFenum + randVal4) / vFe_total);
        double scale_vPb = ((vPbnum + randVal5) / vPb_total);
        double scale_SSB = (count / SSB_total);

        // g_hist_vD->Scale(scale_vD);
        // g_hist_vO->Scale(scale_vO);
        // g_hist_BRN->Scale(scale_BRN);
        // g_hist_vFe->Scale(scale_vFe);
        // g_hist_vPb->Scale(scale_vPb);
        // g_hist_SSB->Scale(scale_SSB);

        // Choose chi2 or NLL fitting method (1 = ChiSquaredFCN, 2 = NegLogLikelihoodFCN)
        // fit_method = 1;
        
        // Set up the fitter
        ROOT::Fit::Fitter fitter;
        
        // Set the function to minimize (our chi-squared or NLL function)
        ROOT::Math::Functor fcn(NegLogLikelihoodFCN, 6);
        fitter.SetFCN(fcn);

        // Set initial parameter values - start with random weights
        // TRandom3 rng(0);  // Use 0 for random seed based on system time, or use a fixed number for reproducibility
        double initial_params[6];
        // initial_params[0] = rng.Uniform(0.05, 10.5);  // Random weight for vD between 0.05 and 10.5
        // initial_params[1] = rng.Uniform(0.05, 10.5);  // Random weight for vO between 0.05 and 10.5
        // initial_params[2] = rng.Uniform(0.05, 10.5);  // Random weight for BRN between 0.05 and 10.5
        // initial_params[3] = rng.Uniform(0.05, 10.5);  // Random weight for vFe between 0.05 and 10.5
        // initial_params[4] = rng.Uniform(0.05, 10.5);  // Random weight for vPb between 0.05 and 10.5
        // initial_params[5] = rng.Uniform(0.05, 10.5);  // Random weight for SSB between 0.05 and 10.5

        // Estimate initial weights from your target composition
        initial_params[0] = scale_vD;
        initial_params[1] = scale_vO;
        initial_params[2] = scale_BRN;
        initial_params[3] = scale_vFe;
        initial_params[4] = scale_vPb;
        initial_params[5] = scale_SSB;

        fitter.Config().SetParamsSettings(6, initial_params);
        
        // Set parameter names
        fitter.Config().ParSettings(0).SetName("weight_vD");
        fitter.Config().ParSettings(1).SetName("weight_vO");
        fitter.Config().ParSettings(2).SetName("weight_BRN");
        fitter.Config().ParSettings(3).SetName("weight_vFe");
        fitter.Config().ParSettings(4).SetName("weight_vPb");
        fitter.Config().ParSettings(5).SetName("weight_SSB");
        
        // Set parameter limits (weights should be non-negative)
        fitter.Config().ParSettings(0).SetLimits(-1e6, 1e6);
        fitter.Config().ParSettings(1).SetLimits(-1e6, 1e6);
        fitter.Config().ParSettings(2).SetLimits(-1e6, 1e6);
        fitter.Config().ParSettings(3).SetLimits(scale_vFe * 0.7, scale_vFe * 1.3);
        fitter.Config().ParSettings(4).SetLimits(-1e6, 1e6);
        fitter.Config().ParSettings(5).SetLimits(-1e6, 1e6);
        
        cout << "\nInitial parameters:" << endl;
        cout << "  weight_vD  = " << initial_params[0] << endl;
        cout << "  weight_vO  = " << initial_params[1] << endl;
        cout << "  weight_BRN = " << initial_params[2] << endl;
        cout << "  weight_vFe = " << initial_params[3] << endl;
        cout << "  weight_vPb = " << initial_params[4] << endl;
        cout << "  weight_SSB = " << initial_params[5] << endl;

        // fitter.SetFitType(1); // chi2
        // fitter.SetFitType(3); // NLL
        
        // Perform the fit
        // cout << "\nPerforming function minimization..." << endl;
        bool fit_success = fitter.FitFCN();
        
        if (!fit_success) {
            cout << "\nERROR: Fit failed!" << endl;
            continue;
        }
        
        // Get results
        const ROOT::Fit::FitResult& result = fitter.Result();

        cout << "\nFitted weights:" << endl;
        cout << "  weight_vD  = " << result.Parameter(0) << " +/- " << result.ParError(0) << endl;
        cout << "  weight_vO  = " << result.Parameter(1) << " +/- " << result.ParError(1) << endl;
        cout << "  weight_BRN = " << result.Parameter(2) << " +/- " << result.ParError(2) << endl;
        cout << "  weight_vFe = " << result.Parameter(3) << " +/- " << result.ParError(3) << endl;
        cout << "  weight_vPb = " << result.Parameter(4) << " +/- " << result.ParError(4) << endl;
        cout << "  weight_SSB = " << result.Parameter(5) << " +/- " << result.ParError(5) << endl;

        h_vD_num->Fill(result.Parameter(0) * g_hist_vD->Integral());
        h_vO_num->Fill(result.Parameter(1) * g_hist_vO->Integral());
        h_BRN_num->Fill(result.Parameter(2) * g_hist_BRN->Integral());
        h_vFe_num->Fill(result.Parameter(3) * g_hist_vFe->Integral());
        h_vPb_num->Fill(result.Parameter(4) * g_hist_vPb->Integral());
        h_SSB_num->Fill(result.Parameter(5) * g_hist_SSB->Integral());

        h_vD_vO->Fill(result.Parameter(0) * g_hist_vD->Integral(), result.Parameter(1) * g_hist_vO->Integral());
        h_vD_BRN->Fill(result.Parameter(0) * g_hist_vD->Integral(), result.Parameter(2) * g_hist_BRN->Integral());
        h_vD_SSB->Fill(result.Parameter(0) * g_hist_vD->Integral(), result.Parameter(5) * g_hist_SSB->Integral());
        h_vO_BRN->Fill(result.Parameter(1) * g_hist_vO->Integral(), result.Parameter(2) * g_hist_BRN->Integral());
        h_vO_SSB->Fill(result.Parameter(1) * g_hist_vO->Integral(), result.Parameter(5) * g_hist_SSB->Integral());
        h_BRN_SSB->Fill(result.Parameter(2) * g_hist_BRN->Integral(), result.Parameter(5) * g_hist_SSB->Integral());
        
        h_bg_vO_trend->Fill(count, result.Parameter(1) * g_hist_vO->Integral());

        cout << "\nRecovered event rates:" << endl;
        cout << "  events_vD  = " << result.Parameter(0) * g_hist_vD->Integral() << " +/- " << result.ParError(0) * g_hist_vD->Integral() << endl;
        cout << "  events_vO  = " << result.Parameter(1) * g_hist_vO->Integral() << " +/- " << result.ParError(1) * g_hist_vO->Integral() << endl;
        cout << "  events_BRN = " << result.Parameter(2) * g_hist_BRN->Integral() << " +/- " << result.ParError(2) * g_hist_BRN->Integral() << endl;
        cout << "  events_vFe = " << result.Parameter(3) * g_hist_vFe->Integral() << " +/- " << result.ParError(3) * g_hist_vFe->Integral() << endl;
        cout << "  events_vPb = " << result.Parameter(4) * g_hist_vPb->Integral() << " +/- " << result.ParError(4) * g_hist_vPb->Integral() << endl;
        cout << "  events_SSB = " << result.Parameter(5) * g_hist_SSB->Integral() << " +/- " << result.ParError(5) * g_hist_SSB->Integral() << endl;
        
        cout << "\n========================================" << endl;
        cout << "Fit Results" << endl;
        cout << "========================================" << endl;
        cout << "Fit status: " << (result.IsValid() ? "VALID" : "INVALID") << endl;
        cout << "Chi-squared/Negative Log-Likelihood: " << result.MinFcnValue() << endl;
        
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
            cout << "Chi2/NDF: " << result.MinFcnValue() / ndf << endl;

            h_chi2_per_ndf->Fill(result.MinFcnValue() / ndf);
        }
        else if (fit_method == 2) {
            // For likelihood fits, we can't directly compute χ²/NDF the same way
            // But we can report the number of bins used
            for (int ix = 1; ix <= g_target->GetNbinsX(); ++ix) {
                for (int iy = 1; iy <= g_target->GetNbinsY(); ++iy) {
                    if (g_target->GetBinContent(ix, iy) > 0) nbins_used++;
                }
            }
            
            cout << "Number of bins with data: " << nbins_used << endl;
            cout << "Number of parameters: 6" << endl;
            cout << "Effective NDF: " << (nbins_used - 6) << endl;
            cout << "NLL per bin: " << result.MinFcnValue() / nbins_used << endl;

            h_chi2_per_ndf->Fill(result.MinFcnValue() / nbins_used);
        }
        
        // Create the fitted MC histogram (weighted sum) - 2D
        for (int ix = 1; ix <= h_fitted_2D->GetNbinsX(); ++ix) {
            for (int iy = 1; iy <= h_fitted_2D->GetNbinsY(); ++iy) {
                double val = result.Parameter(0) * g_hist_vD->GetBinContent(ix, iy) +
                             result.Parameter(1) * g_hist_vO->GetBinContent(ix, iy) +
                             result.Parameter(2) * g_hist_BRN->GetBinContent(ix, iy) +
                             result.Parameter(3) * g_hist_vFe->GetBinContent(ix, iy) +
                             result.Parameter(4) * g_hist_vPb->GetBinContent(ix, iy) +
                             result.Parameter(5) * g_hist_SSB->GetBinContent(ix, iy);
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

            c2->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/fitting_results/four_test_cases/plots/2D_Event_Rates.png");
            
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
                latex1->DrawLatex(0.6, 0.70, Form("#chi^{2}/NDF = %.2f/%d", result.MinFcnValue(), ndf));
                latex1->DrawLatex(0.6, 0.65, Form("= %.2f", result.MinFcnValue()/ndf));
            }
            else if (fit_method == 2) {
                // Add NLL info
                TLatex* latex1 = new TLatex();
                latex1->SetNDC();
                latex1->SetTextSize(0.04);
                latex1->DrawLatex(0.6, 0.70, Form("-log(L) = %.1f", result.MinFcnValue()));
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

            c3->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/fitting_results/four_test_cases/plots/MC_to_Data_Fit_Results-2D.png");
            
            // ===================================================================
            // Create a second canvas with MC components
            // ===================================================================
            TCanvas* c4 = new TCanvas("c4", "MC Components", 2400, 800);
            c4->Divide(7, 2);
            
            // Plot weighted MC components - Energy projections
            c4->cd(1);
            gPad->SetLogy();
            TH2D* h_vD_scaled = (TH2D*)g_hist_vD->Clone("h_vD_scaled");
            h_vD_scaled->Scale(result.Parameter(0));
            TH1D* h_vD_energy = h_vD_scaled->ProjectionY("h_vD_energy");
            h_vD_energy->SetLineColor(kRed);
            h_vD_energy->SetLineWidth(2);
            h_vD_energy->SetTitle(Form("#nu-D Energy (w=%.5f)", result.Parameter(0)));
            h_vD_energy->GetXaxis()->SetTitle("Energy");
            h_vD_energy->Draw("HIST");
            
            c4->cd(2);
            gPad->SetLogy();
            TH2D* h_vO_scaled = (TH2D*)g_hist_vO->Clone("h_vO_scaled");
            h_vO_scaled->Scale(result.Parameter(1));
            TH1D* h_vO_energy = h_vO_scaled->ProjectionY("h_vO_energy");
            h_vO_energy->SetLineColor(kBlue);
            h_vO_energy->SetLineWidth(2);
            h_vO_energy->SetTitle(Form("#nu-O Energy (w=%.5f)", result.Parameter(1)));
            h_vO_energy->GetXaxis()->SetTitle("Energy");
            h_vO_energy->Draw("HIST");
            
            c4->cd(3);
            gPad->SetLogy();
            TH2D* h_BRN_scaled = (TH2D*)g_hist_BRN->Clone("h_BRN_scaled");
            h_BRN_scaled->Scale(result.Parameter(2));
            TH1D* h_BRN_energy = h_BRN_scaled->ProjectionY("h_BRN_energy");
            h_BRN_energy->SetLineColor(kGreen+2);
            h_BRN_energy->SetLineWidth(2);
            h_BRN_energy->SetTitle(Form("BRN Energy (w=%.5f)", result.Parameter(2)));
            h_BRN_energy->GetXaxis()->SetTitle("Energy");
            h_BRN_energy->Draw("HIST");

            c4->cd(4);
            gPad->SetLogy();
            TH2D* h_vFe_scaled = (TH2D*)g_hist_vFe->Clone("h_vFe_scaled");
            h_vFe_scaled->Scale(result.Parameter(3));
            TH1D* h_vFe_energy = h_vFe_scaled->ProjectionY("h_vFe_energy");
            h_vFe_energy->SetLineColor(kCyan+2);
            h_vFe_energy->SetLineWidth(2);
            h_vFe_energy->SetTitle(Form("vFe Energy (w=%.5f)", result.Parameter(3)));
            h_vFe_energy->GetXaxis()->SetTitle("Energy");
            h_vFe_energy->Draw("HIST");

            c4->cd(5);
            gPad->SetLogy();
            TH2D* h_vPb_scaled = (TH2D*)g_hist_vPb->Clone("h_vPb_scaled");
            h_vPb_scaled->Scale(result.Parameter(4));
            TH1D* h_vPb_energy = h_vPb_scaled->ProjectionY("h_vPb_energy");
            h_vPb_energy->SetLineColor(kYellow+3);
            h_vPb_energy->SetLineWidth(2);
            h_vPb_energy->SetTitle(Form("vPb Energy (w=%.5f)", result.Parameter(4)));
            h_vPb_energy->GetXaxis()->SetTitle("Energy");
            h_vPb_energy->Draw("HIST");

            c4->cd(6);
            gPad->SetLogy();
            TH2D* h_SSB_scaled = (TH2D*)g_hist_SSB->Clone("h_SSB_scaled");
            h_SSB_scaled->Scale(result.Parameter(5));
            TH1D* h_SSB_energy = h_SSB_scaled->ProjectionY("h_SSB_energy");
            h_SSB_energy->SetLineColor(kOrange);
            h_SSB_energy->SetLineWidth(2);
            h_SSB_energy->SetTitle(Form("SSB Energy (w=%.5f)", result.Parameter(5)));
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
            h_vD_time->SetTitle(Form("#nu-D Time (w=%.5f)", result.Parameter(0)));
            h_vD_time->GetXaxis()->SetTitle("Time");
            h_vD_time->Draw("HIST");
            
            c4->cd(9);
            gPad->SetLogy();
            TH1D* h_vO_time = h_vO_scaled->ProjectionX("h_vO_time");
            h_vO_time->SetLineColor(kBlue);
            h_vO_time->SetLineWidth(2);
            h_vO_time->SetTitle(Form("#nu-O Time (w=%.5f)", result.Parameter(1)));
            h_vO_time->GetXaxis()->SetTitle("Time");
            h_vO_time->Draw("HIST");
            
            c4->cd(10);
            gPad->SetLogy();
            TH1D* h_BRN_time = h_BRN_scaled->ProjectionX("h_BRN_time");
            h_BRN_time->SetLineColor(kGreen+2);
            h_BRN_time->SetLineWidth(2);
            h_BRN_time->SetTitle(Form("BRN Time (w=%.5f)", result.Parameter(2)));
            h_BRN_time->GetXaxis()->SetTitle("Time");
            h_BRN_time->Draw("HIST");

            c4->cd(11);
            gPad->SetLogy();
            TH1D* h_vFe_time = h_vFe_scaled->ProjectionX("h_vFe_time");
            h_vFe_time->SetLineColor(kCyan+2);
            h_vFe_time->SetLineWidth(2);
            h_vFe_time->SetTitle(Form("vFe Time (w=%.5f)", result.Parameter(3)));
            h_vFe_time->GetXaxis()->SetTitle("Time");
            h_vFe_time->Draw("HIST");

            c4->cd(12);
            gPad->SetLogy();
            TH1D* h_vPb_time = h_vPb_scaled->ProjectionX("h_vPb_time");
            h_vPb_time->SetLineColor(kYellow+3);
            h_vPb_time->SetLineWidth(2);
            h_vPb_time->SetTitle(Form("vPb Time (w=%.5f)", result.Parameter(4)));
            h_vPb_time->GetXaxis()->SetTitle("Time");
            h_vPb_time->Draw("HIST");

            c4->cd(13);
            gPad->SetLogy();
            TH1D* h_SSB_time = h_SSB_scaled->ProjectionX("h_SSB_time");
            h_SSB_time->SetLineColor(kOrange);
            h_SSB_time->SetLineWidth(2);
            h_SSB_time->SetTitle(Form("SSB Time (w=%.5f)", result.Parameter(5)));
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

            c4->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/fitting_results/four_test_cases/plots/MC_Components.png");

            // ===================================================================
            // Create a third canvas with original 2D histograms
            // ===================================================================
            /*TCanvas* c5 = new TCanvas("c5", "2D Plots", 1600, 1200);
            c5->Divide(2, 2);

            // vD 2D histogram
            c5->cd(1);
            g_hist_vD->SetTitle("Original vD Data");
            g_hist_vD->GetXaxis()->SetTitle("Time");
            g_hist_vD->GetYaxis()->SetTitle("Energy");
            g_hist_vD->Draw("COLZ");
            
            // vO 2D histogram
            c5->cd(2);
            g_hist_vO->SetTitle("Original vO Data");
            g_hist_vO->GetXaxis()->SetTitle("Time");
            g_hist_vO->GetYaxis()->SetTitle("Energy");
            g_hist_vO->Draw("COLZ");

            // BRN 2D histogram
            c5->cd(3);
            g_hist_BRN->SetTitle("Original BRN Data");
            g_hist_BRN->GetXaxis()->SetTitle("Time");
            g_hist_BRN->GetYaxis()->SetTitle("Energy");
            g_hist_BRN->Draw("COLZ");

            // Michels 2D histogram
            c5->cd(4);
            // g_target->SetTitle("Original Michel Data");
            g_target->GetXaxis()->SetTitle("Time");
            g_target->GetYaxis()->SetTitle("Energy");
            g_target->Draw("COLZ");

            c5->Update();*/

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

            c6->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/fitting_results/four_test_cases/plots/1D_Plots.png");

            // ===================================================================
            // Create a canvas showing true and scaled background plots
            // ===================================================================
            TCanvas* c7 = new TCanvas("c7", "Steady State Background", 1600, 1200);
            c7->Divide(2, 1);

            /*c7->cd(1);
            // gPad->SetLogy();
            h_background_true->SetTitle("Original Background");
            h_background_true->GetXaxis()->SetTitle("Time");
            h_background_true->GetYaxis()->SetTitle("Energy");
            h_background_true->Draw("COLZ");

            c7->cd(2);
            // gPad->SetLogy();
            g_hist_SSB->SetTitle("100x Scaled Background");
            g_hist_SSB->GetXaxis()->SetTitle("Time");
            g_hist_SSB->GetYaxis()->SetTitle("Energy");
            g_hist_SSB->Draw("COLZ");*/

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

            c7->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/fitting_results/four_test_cases/plots/Steady_State_Background.png");

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

            c9->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/fitting_results/four_test_cases/plots/Residual_Errors.png");

        }

    }

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
    // Plot 7: Chi2/NDF
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

    c1->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/fitting_results/four_test_cases/plots/Interaction_Event_Rates.png");

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
    // Plot 1: vD-BRN
    // ===================================================================
    
    c8->cd(2);
    h_vD_BRN->SetTitle("vD-BRN");
    h_vD_BRN->GetXaxis()->SetTitle("vD Event Rate");
    h_vD_BRN->GetYaxis()->SetTitle("BRN Event Rate");
    h_vD_BRN->Draw("COLZ");

    // ===================================================================
    // Plot 1: vD-SSB
    // ===================================================================
    
    c8->cd(3);
    h_vD_SSB->SetTitle("vD-SSB");
    h_vD_SSB->GetXaxis()->SetTitle("vD Event Rate");
    h_vD_SSB->GetYaxis()->SetTitle("SSB Event Rate");
    h_vD_SSB->Draw("COLZ");

    // ===================================================================
    // Plot 1: vO-BRN
    // ===================================================================
    
    c8->cd(4);
    h_vO_BRN->SetTitle("vO-BRN");
    h_vO_BRN->GetXaxis()->SetTitle("vO Event Rate");
    h_vO_BRN->GetYaxis()->SetTitle("BRN Event Rate");
    h_vO_BRN->Draw("COLZ");

    // ===================================================================
    // Plot 1: vO-SSB
    // ===================================================================
    
    c8->cd(5);
    h_vO_SSB->SetTitle("vO-SSB");
    h_vO_SSB->GetXaxis()->SetTitle("vO Event Rate");
    h_vO_SSB->GetYaxis()->SetTitle("SSB Event Rate");
    h_vO_SSB->Draw("COLZ");

    // ===================================================================
    // Plot 1: BRN-SSB
    // ===================================================================
    
    c8->cd(6);
    h_BRN_SSB->SetTitle("BRN-SSB");
    h_BRN_SSB->GetXaxis()->SetTitle("BRN Event Rate");
    h_BRN_SSB->GetYaxis()->SetTitle("SSB Event Rate");
    h_BRN_SSB->Draw("COLZ");

    c8->Update();

    c8->SaveAs("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/fitting_results/four_test_cases/plots/Recovered_Weight_Correlations.png");

    /*TCanvas* c2 = new TCanvas("c2", "Interaction Event Rates", 1600, 1200);
    c2->Divide(4, 1);

    // ===================================================================
    // Plot 1: Target
    // ===================================================================
    c2->cd(1);
    g_target->SetTitle("Target");
    g_target->GetXaxis()->SetTitle("Time");
    g_target->GetYaxis()->SetTitle("Energy");
    g_target->Draw();
    
    // ===================================================================
    // Plot 1: vD
    // ===================================================================
    c2->cd(2);
    g_hist_vD->SetTitle("vD");
    g_hist_vD->GetXaxis()->SetTitle("Time");
    g_hist_vD->GetYaxis()->SetTitle("Energy");
    g_hist_vD->Draw();
    
    // ===================================================================
    // Plot 2: vO
    // ===================================================================
    c2->cd(3);
    g_hist_vO->SetTitle("vO");
    g_hist_vO->GetXaxis()->SetTitle("Time");
    g_hist_vO->GetYaxis()->SetTitle("Energy");
    g_hist_vO->Draw();
    
    // ===================================================================
    // Plot 3: BRN
    // ===================================================================
    c2->cd(4);
    g_hist_BRN->SetTitle("BRN");
    g_hist_BRN->GetXaxis()->SetTitle("Time");
    g_hist_BRN->GetYaxis()->SetTitle("Energy");
    g_hist_BRN->Draw();
    
    TCanvas* c3 = new TCanvas("c3", "MC to Data Fit Results - 2D", 1600, 1200);
    c3->Divide(3, 2);
    
    // ===================================================================
    // Plot 1: Data 2D histogram
    // ===================================================================
    c3->cd(1);
    g_target->SetTitle("Fake Data (Time vs Energy)");
    g_target->GetXaxis()->SetTitle("Time");
    g_target->GetYaxis()->SetTitle("Energy");
    g_target->Draw("COLZ");
    
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
        latex1->DrawLatex(0.6, 0.70, Form("#chi^{2}/NDF = %.2f/%d", result.MinFcnValue(), ndf));
        latex1->DrawLatex(0.6, 0.65, Form("= %.2f", result.MinFcnValue()/ndf));
    }
    else if (fit_method == 2) {
        // Add NLL info
        TLatex* latex1 = new TLatex();
        latex1->SetNDC();
        latex1->SetTextSize(0.04);
        latex1->DrawLatex(0.6, 0.70, Form("-log(L) = %.1f", result.MinFcnValue()));
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
    
    // ===================================================================
    // Create a second canvas with MC components
    // ===================================================================
    TCanvas* c4 = new TCanvas("c4", "MC Components", 1600, 800);
    c4->Divide(4, 2);
    
    // Plot weighted MC components - Energy projections
    c4->cd(1);
    gPad->SetLogy();
    TH2D* h_vD_scaled = (TH2D*)g_hist_vD->Clone("h_vD_scaled");
    h_vD_scaled->Scale(result.Parameter(0));
    TH1D* h_vD_energy = h_vD_scaled->ProjectionY("h_vD_energy");
    h_vD_energy->SetLineColor(kRed);
    h_vD_energy->SetLineWidth(2);
    h_vD_energy->SetTitle(Form("#nu-D Energy (w=%.3f)", result.Parameter(0)));
    h_vD_energy->GetXaxis()->SetTitle("Energy");
    h_vD_energy->Draw("HIST");
    
    c4->cd(2);
    gPad->SetLogy();
    TH2D* h_vO_scaled = (TH2D*)g_hist_vO->Clone("h_vO_scaled");
    h_vO_scaled->Scale(result.Parameter(1));
    TH1D* h_vO_energy = h_vO_scaled->ProjectionY("h_vO_energy");
    h_vO_energy->SetLineColor(kBlue);
    h_vO_energy->SetLineWidth(2);
    h_vO_energy->SetTitle(Form("#nu-O Energy (w=%.3f)", result.Parameter(1)));
    h_vO_energy->GetXaxis()->SetTitle("Energy");
    h_vO_energy->Draw("HIST");
    
    c4->cd(3);
    gPad->SetLogy();
    TH2D* h_BRN_scaled = (TH2D*)g_hist_BRN->Clone("h_BRN_scaled");
    h_BRN_scaled->Scale(result.Parameter(2));
    TH1D* h_BRN_energy = h_BRN_scaled->ProjectionY("h_BRN_energy");
    h_BRN_energy->SetLineColor(kGreen+2);
    h_BRN_energy->SetLineWidth(2);
    h_BRN_energy->SetTitle(Form("BRN Energy (w=%.3f)", result.Parameter(2)));
    h_BRN_energy->GetXaxis()->SetTitle("Energy");
    h_BRN_energy->Draw("HIST");

    c4->cd(4);
    gPad->SetLogy();
    TH2D* h_Me_scaled = (TH2D*)g_target->Clone("h_Me_scaled");
    TH1D* h_Me_energy = h_Me_scaled->ProjectionY("h_Me_energy");
    h_Me_energy->SetLineColor(kMagenta+2);
    h_Me_energy->SetLineWidth(2);
    h_Me_energy->SetTitle(Form("Target Energy"));
    h_Me_energy->GetXaxis()->SetTitle("Energy");
    h_Me_energy->Draw("HIST");
    
    // Plot weighted MC components - Time projections
    c4->cd(5);
    gPad->SetLogy();
    TH1D* h_vD_time = h_vD_scaled->ProjectionX("h_vD_time");
    h_vD_time->SetLineColor(kRed);
    h_vD_time->SetLineWidth(2);
    h_vD_time->SetTitle(Form("#nu-D Time (w=%.3f)", result.Parameter(0)));
    h_vD_time->GetXaxis()->SetTitle("Time");
    h_vD_time->Draw("HIST");
    
    c4->cd(6);
    gPad->SetLogy();
    TH1D* h_vO_time = h_vO_scaled->ProjectionX("h_vO_time");
    h_vO_time->SetLineColor(kBlue);
    h_vO_time->SetLineWidth(2);
    h_vO_time->SetTitle(Form("#nu-O Time (w=%.3f)", result.Parameter(1)));
    h_vO_time->GetXaxis()->SetTitle("Time");
    h_vO_time->Draw("HIST");
    
    c4->cd(7);
    gPad->SetLogy();
    TH1D* h_BRN_time = h_BRN_scaled->ProjectionX("h_BRN_time");
    h_BRN_time->SetLineColor(kGreen+2);
    h_BRN_time->SetLineWidth(2);
    h_BRN_time->SetTitle(Form("BRN Time (w=%.3f)", result.Parameter(2)));
    h_BRN_time->GetXaxis()->SetTitle("Time");
    h_BRN_time->Draw("HIST");

    c4->cd(8);
    gPad->SetLogy();
    TH1D* h_Me_time = h_Me_scaled->ProjectionX("h_Me_time");
    h_Me_time->SetLineColor(kMagenta+2);
    h_Me_time->SetLineWidth(2);
    h_Me_time->SetTitle(Form("Target Time"));
    h_Me_time->GetXaxis()->SetTitle("Time");
    h_Me_time->Draw("HIST");
    
    c4->Update();

    // ===================================================================
    // Create a third canvas with original 2D histograms
    // ===================================================================
    TCanvas* c5 = new TCanvas("c5", "2D Plots", 1600, 1200);
    c5->Divide(2, 2);

    // vD 2D histogram
    c5->cd(1);
    g_hist_vD->SetTitle("Original vD Data");
    g_hist_vD->GetXaxis()->SetTitle("Time");
    g_hist_vD->GetYaxis()->SetTitle("Energy");
    g_hist_vD->Draw("COLZ");
    
    // vO 2D histogram
    c5->cd(2);
    g_hist_vO->SetTitle("Original vO Data");
    g_hist_vO->GetXaxis()->SetTitle("Time");
    g_hist_vO->GetYaxis()->SetTitle("Energy");
    g_hist_vO->Draw("COLZ");

    // BRN 2D histogram
    c5->cd(3);
    g_hist_BRN->SetTitle("Original BRN Data");
    g_hist_BRN->GetXaxis()->SetTitle("Time");
    g_hist_BRN->GetYaxis()->SetTitle("Energy");
    g_hist_BRN->Draw("COLZ");

    // Michels 2D histogram
    c5->cd(4);
    // g_target->SetTitle("Original Michel Data");
    g_target->GetXaxis()->SetTitle("Time");
    g_target->GetYaxis()->SetTitle("Energy");
    g_target->Draw("COLZ");

    c5->Update();

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
    h_BRN_sim_energy->SetTitle("Original BRN Energy");
    h_BRN_sim_energy->GetXaxis()->SetTitle("Energy");
    // h_BRN_sim_energy->Draw("HIST");

    // c6->cd(2);
    gPad->SetLogy();
    TH1D* h_vO_sim_energy = g_hist_vO->ProjectionY("h_vO_sim_energy");
    TH1D* h_BRN_vO_sim_energy = (TH1D*)h_BRN_sim_energy->Clone("h_BRN_vO_sim_energy");
    h_BRN_vO_sim_energy->Add(h_vO_sim_energy);
    h_BRN_vO_sim_energy->SetLineColor(kBlue);
    h_BRN_vO_sim_energy->SetLineWidth(2);
    h_BRN_vO_sim_energy->SetTitle("Original BRN + #nu-O Energy");
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
    h_BRN_energy->SetTitle("Weighted BRN Energy");
    h_BRN_energy->GetXaxis()->SetTitle("Energy");
    // h_BRN_energy->Draw("HIST");

    // c6->cd(5);
    gPad->SetLogy();
    TH1D* h_BRN_vO_energy = (TH1D*)h_BRN_energy->Clone("h_BRN_vO_energy");
    h_BRN_vO_energy->Add(h_vO_energy);
    h_BRN_vO_energy->SetLineColor(kBlue);
    h_BRN_vO_energy->SetLineWidth(2);
    h_BRN_vO_energy->SetTitle("Weighted BRN + #nu-O Energy");
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

    c6->Update();*/
    
    cout << "\nPlots created successfully!" << endl;
    cout << "Canvas 1: Number of Fitted Events of Each Interaction Type" << endl;
    cout << "Canvas 2: Data vs MC comparison with 2D histograms and projections" << endl;
    cout << "Canvas 3: Individual MC components (energy and time projections)" << endl;
    cout << "========================================\n" << endl;
}