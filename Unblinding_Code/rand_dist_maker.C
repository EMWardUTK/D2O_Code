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
#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

using std::cout; 
using std::endl;
using std::string;
using std::ifstream;
using namespace std;

void rand_dist_maker() {

    // ===================================================================
    // CREATE RANDOM LIST OF SNS TIME VALUES
    // ===================================================================
    cout << "\n========================================" << endl;
    cout << "Creating First List of Random SNS v Times" << endl;
    cout << "========================================" << endl;

    TH1D* h_v_time_1 = new TH1D("h_v_time_1", "Random SNS v Time Distribution 1", 100, 0, 5000);
    TH1D* h_v_time_2 = new TH1D("h_v_time_2", "Random SNS v Time Distribution 2", 100, 0, 5000);
    TH1D* h_brn_time = new TH1D("h_brn_time", "Random BRN Time Distribution", 100, 0, 1200);

    int num_data_pts = 1000000;

    std::vector<double> sns_v_time_spec;
    string intval1;
    ifstream ReadTextFile1("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/SNS_v_Time_Plot_w_Fluctuations.txt");
    while (getline(ReadTextFile1, intval1)) {
        sns_v_time_spec.push_back(stod(intval1));
    }
    ReadTextFile1.close();
    ReadTextFile1.clear();

    TRandom3 rng(0);
    std::vector<double> sns_v_time_dist;
    double rand_t_val;
    double rand_c_val;
    int i = 0;
    while (i < num_data_pts) {
        rand_t_val = rand() % 5001;
        rand_c_val = rng.Uniform(0.001, *max_element(sns_v_time_spec.begin(), sns_v_time_spec.end()));
        if (rand_t_val >= 0 && rand_c_val < sns_v_time_spec[rand_t_val]) {
            sns_v_time_dist.push_back(rand_t_val);
            h_v_time_1->Fill(rand_t_val);
            i += 1;
        }
    }

    ofstream OutputTextFile1("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/SNS_v_Time_Dist_1M_5us_1.txt");

    for(int k = 0; k < sns_v_time_dist.size(); k++){
        OutputTextFile1 << sns_v_time_dist[k] << "\n";
    }

    OutputTextFile1.close();

    cout << "sns_v_time_dist.size() = " << sns_v_time_dist.size() << endl;

    cout << "\n========================================" << endl;
    cout << "Creating Second List of Random SNS v Times" << endl;
    cout << "========================================" << endl;

    sns_v_time_dist.clear();
    i = 0;
    while (i < num_data_pts) {
        rand_t_val = rand() % 5001;
        rand_c_val = rng.Uniform(0.001, *max_element(sns_v_time_spec.begin(), sns_v_time_spec.end()));
        if (rand_t_val >= 0 && rand_c_val < sns_v_time_spec[rand_t_val]) {
            sns_v_time_dist.push_back(rand_t_val);
            h_v_time_2->Fill(rand_t_val);
            i += 1;
        }
    }

    ofstream OutputTextFile2("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/SNS_v_Time_Dist_1M_5us_2.txt");

    for(int k = 0; k < sns_v_time_dist.size(); k++){
        OutputTextFile2 << sns_v_time_dist[k] << "\n";
    }

    OutputTextFile2.close();

    cout << "sns_v_time_dist.size() = " << sns_v_time_dist.size() << endl;

    cout << "\n========================================" << endl;
    cout << "Creating One List of Random SNS BRN Times" << endl;
    cout << "========================================" << endl;

    std::vector<double> brn_time_spec;
    string intval2;
    ifstream ReadTextFile2("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/BRN_Time_Plot_w_Fluctuations.txt");
    while (getline(ReadTextFile2, intval2)) {
        brn_time_spec.push_back(stod(intval2));
    }
    ReadTextFile2.close();
    ReadTextFile2.clear();

    std::vector<double> brn_time_dist;
    i = 0;
    while (i < num_data_pts) {
        rand_t_val = rand() % 1201;
        rand_c_val = rng.Uniform(0.001, *max_element(brn_time_spec.begin(), brn_time_spec.end()));
        if (rand_c_val < brn_time_spec[rand_t_val]) {
            brn_time_dist.push_back(rand_t_val);
            h_brn_time->Fill(rand_t_val);
            i += 1;
        }
    }

    ofstream OutputTextFile3("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/Data_MC_Comp/BRN_Time_Dist_1M.txt");

    for(int k = 0; k < brn_time_dist.size(); k++){
        OutputTextFile3 << brn_time_dist[k] << "\n";
    }

    OutputTextFile3.close();

    cout << "brn_time_dist.size() = " << brn_time_dist.size() << endl;

    // ===================================================================
    // Create a canvas showing random SNS time distributions
    // ===================================================================
    TCanvas* c1 = new TCanvas("c1", "SNS Time Distributions", 1600, 800);
    c1->Divide(3, 1);
    
    c1->cd(1);
    gPad->SetLogy();
    h_v_time_1->SetLineColor(kRed);
    h_v_time_1->SetLineWidth(2);
    h_v_time_1->SetTitle("SNS v Time Distribution");
    h_v_time_1->GetXaxis()->SetTitle("Time (ns)");
    h_v_time_1->GetYaxis()->SetTitle("Counts");
    h_v_time_1->Draw();
    
    c1->cd(2);
    gPad->SetLogy();
    h_v_time_2->SetLineColor(kBlue);
    h_v_time_2->SetLineWidth(2);
    h_v_time_2->SetTitle("SNS v Time Distribution");
    h_v_time_2->GetXaxis()->SetTitle("Time (ns)");
    h_v_time_2->GetYaxis()->SetTitle("Counts");
    h_v_time_2->Draw();

    c1->cd(3);
    // gPad->SetLogy();
    h_brn_time->SetLineColor(kGreen+2);
    h_brn_time->SetLineWidth(2);
    h_brn_time->SetTitle("SNS BRN Time Distribution");
    h_brn_time->GetXaxis()->SetTitle("Time (ns)");
    h_brn_time->GetYaxis()->SetTitle("Counts");
    h_brn_time->Draw();

    c1->Update();

}
