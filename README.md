# D2O_Code
All source codes for analysis and unblinding of COHERENT's heavy water Cherenkov detector.

Before running the core analysis code on newly obtained data, energy calibration and corrections need to be performed on the data. The calibration code is currently set up such that it needs to be run in a specific order, given below.

Energy Calibration Procedure:
  1.) Single phe int conversion code: single_phe_int.cpp
  2.) SiPM noise pedestal cuts code: main_michels_vp_noise_cuts.cpp
  3.) Minimum muon cuts code: main_min_muons.cpp
  4.) Vertical muon integral correction code: main_vert_muons.cpp
