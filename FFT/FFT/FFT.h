// FFT.h - Contains declarations of functions for computing FFT of spectral data
// This header file declares some functions to perform FFT calculations

#pragma once

double convert_dBm_to_mW(double dBm_val);

double convert_mW_to_dBm(double mW_val);

void compute_FFT(unsigned long& N_spctr_data, double& spctr_hor_spacing, std::vector<double>& spctr_data, int& N_fft_data, std::vector<double>& fft_data, std::vector<double>& fft_abcissae);

extern "C" _declspec(dllexport) void Spectrum_FFT(unsigned long& N_spctr_data, double spctr_hor_spacing, double spctr_data[], int& N_fft_data, double fft_data[], double fft_abcissae[]);
